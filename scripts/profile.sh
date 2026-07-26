#!/usr/bin/env bash
#
# Profile Vortex-Transport with linux perf.
#
# Usage:
#   scripts/profile.sh                          # runs ./Vortex-Transport-prof input_files/input
#   scripts/profile.sh ./Vortex-Transport-prof my_input
#
# Run from the repository root. Build the profiling binary first:
#   make profile ENABLE_OPENMP=TRUE     # profile the OpenMP build (what production runs use)
#   make profile                        # profile the serial build
# Run `make clean` first when switching between the two, and note that
# ENABLE_OPENMP is a separate knob: plain `make profile` silently builds a
# SERIAL binary. This script warns if the binary it is given has no OpenMP.
#
# Outputs (written to the current directory):
#   perf.data        raw samples          (interactive: perf report -i perf.data)
#   perf_report.txt  text version of the perf report
#   flamegraph.svg   flamegraph, only if a flamegraph tool is installed
#
# ---------------------------------------------------------------------------
# Installing perf (Debian/Ubuntu):
#   sudo apt-get install linux-tools-common linux-tools-generic linux-tools-$(uname -r)
#
# On WSL2 there is usually no linux-tools package matching the kernel
# ($(uname -r) ends in "-microsoft-standard-WSL2"). Options:
#   1. Install linux-tools-generic and call the versioned binary directly:
#        sudo apt-get install linux-tools-generic
#        ls /usr/lib/linux-tools/*/perf     # this script finds it automatically
#   2. Or build perf from the WSL2 kernel sources:
#        git clone --depth 1 https://github.com/microsoft/WSL2-Linux-Kernel
#        make -C WSL2-Linux-Kernel/tools/perf && sudo cp WSL2-Linux-Kernel/tools/perf/perf /usr/local/bin/
#
# If perf record fails with a permissions error, allow profiling with:
#   sudo sysctl kernel.perf_event_paranoid=1
#
# Optional flamegraph tools (either one works):
#   cargo install inferno                                    # inferno-collapse-perf / inferno-flamegraph
#   git clone https://github.com/brendangregg/FlameGraph     # stackcollapse-perf.pl / flamegraph.pl (add to PATH)
#
# Complementary counters (not covered by this script):
#   perf stat ./Vortex-Transport-prof input_files/input
# The "CPUs utilized" line is the quickest OpenMP sanity check: it should be
# close to OMP_NUM_THREADS if threads are busy. On WSL2 hardware counters
# (cycles, instructions, cache-misses) are usually <not supported>; perf then
# samples on task-clock, which is fine for time-based profiles.
#
# Reading OpenMP profiles:
#   - Outlined parallel regions show up as "func(...) [clone ._omp_fn.N]".
#     Worker-thread samples root at gomp_thread_start/start_thread, NOT at
#     main: only the master thread's chain goes through GOMP_parallel back to
#     the caller. That is inherent to OpenMP thread pools, not a perf bug.
#   - Time in gomp_team_barrier_wait_end / do_wait / do_spin is threads
#     SPINNING at barriers (idle CPU burn, load imbalance or too little work
#     per region) - it is not real work in the function that contains the
#     parallel loop. Resolving these names needs libgomp debug symbols
#     (libgomp1-dbgsym); this script warns if they are missing.
#   - To measure compute without the spin noise, rerun with
#     OMP_WAIT_POLICY=passive (threads sleep instead of spinning).
# ---------------------------------------------------------------------------

set -euo pipefail

BINARY="${1:-./Vortex-Transport-prof}"
shift || true
ARGS=("${@:-input_files/input}")

# ---------------------------------------------------------------------------
# locate perf: PATH first, then the versioned linux-tools binaries (WSL2 case,
# where /usr/bin/perf refuses to run because no tool matches the kernel)
# ---------------------------------------------------------------------------
PERF=""
if command -v perf >/dev/null 2>&1 && perf --version >/dev/null 2>&1; then
    PERF="perf"
else
    for candidate in /usr/lib/linux-tools/*/perf; do
        if [ -x "$candidate" ]; then
            PERF="$candidate"
            break
        fi
    done
fi

if [ -z "$PERF" ]; then
    echo "ERROR: perf not found. See the install instructions in the comments of this script." >&2
    exit 1
fi
echo "using perf: $PERF ($($PERF --version))"

if [ ! -x "$BINARY" ]; then
    echo "ERROR: '$BINARY' not found or not executable." >&2
    echo "Build it first with:  make profile ENABLE_OPENMP=TRUE" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# OpenMP sanity checks: is this the parallel binary, and can libgomp's
# internal barrier/spin functions be resolved to names in the report?
# ---------------------------------------------------------------------------
GOMP_PATH="$(ldd "$BINARY" 2>/dev/null | awk '/libgomp/ {print $3}')"
if [ -n "$GOMP_PATH" ]; then
    echo "binary is an OpenMP build (links $GOMP_PATH); perf records all OpenMP threads"
    BUILD_ID="$(readelf -n "$GOMP_PATH" 2>/dev/null | awk '/Build ID/ {print $NF}')"
    DEBUG_FILE="/usr/lib/debug/.build-id/${BUILD_ID:0:2}/${BUILD_ID:2}.debug"
    if [ -n "$BUILD_ID" ] && [ ! -f "$DEBUG_FILE" ]; then
        echo "WARNING: no debug symbols for libgomp (missing $DEBUG_FILE)." >&2
        echo "         Barrier/spin-wait time will show as raw hex addresses in libgomp.so.1" >&2
        echo "         instead of names like gomp_team_barrier_wait_end. Install the" >&2
        echo "         libgomp1-dbgsym package matching 'dpkg -l libgomp1' to fix this." >&2
    fi
else
    echo "WARNING: '$BINARY' is a SERIAL build (does not link libgomp)." >&2
    echo "         If you meant to profile the OpenMP code, rebuild first:" >&2
    echo "             make clean && make profile ENABLE_OPENMP=TRUE" >&2
fi

# ---------------------------------------------------------------------------
# record: -g captures call graphs with frame-pointer unwinding. This works
# here because the profiling build uses -fno-omit-frame-pointer and Ubuntu
# 24.04+ system libraries (libgomp, libc) are built with frame pointers too;
# --call-graph dwarf would only cost ~10x more perf.data for the same chains.
# -F 499 samples ~499 times/second PER THREAD (perf follows all threads the
# child spawns, so no -a needed).
# ---------------------------------------------------------------------------
echo "recording: $PERF record -g -F 499 -o perf.data -- $BINARY ${ARGS[*]}"
RECORD_START=$(date +%s.%N)
"$PERF" record -g -F 499 -o perf.data -- "$BINARY" "${ARGS[@]}"
RECORD_END=$(date +%s.%N)
echo "program wall time under perf: $(echo "$RECORD_END - $RECORD_START" | bc) s (compare against a plain run to gauge perf overhead)"

# ---------------------------------------------------------------------------
# text report
# ---------------------------------------------------------------------------
"$PERF" report -i perf.data --stdio > perf_report.txt
echo ""
echo "wrote perf_report.txt (top of the report below)"
echo "for the interactive view run: $PERF report -i perf.data"
echo "----------------------------------------------------------------------"
grep -v '^#' perf_report.txt | head -n 25 || true
echo "----------------------------------------------------------------------"

# ---------------------------------------------------------------------------
# flamegraph, if a tool is available (inferno or Brendan Gregg's FlameGraph);
# falls back gracefully with a hint if neither is installed
# ---------------------------------------------------------------------------
if command -v inferno-collapse-perf >/dev/null 2>&1 && command -v inferno-flamegraph >/dev/null 2>&1; then
    "$PERF" script -i perf.data | inferno-collapse-perf | inferno-flamegraph > flamegraph.svg
    echo "wrote flamegraph.svg (inferno)"
elif command -v stackcollapse-perf.pl >/dev/null 2>&1 && command -v flamegraph.pl >/dev/null 2>&1; then
    "$PERF" script -i perf.data | stackcollapse-perf.pl | flamegraph.pl > flamegraph.svg
    echo "wrote flamegraph.svg (FlameGraph)"
else
    echo "no flamegraph tool found, skipping flamegraph.svg"
    echo "  install one with:  cargo install inferno"
    echo "  or clone https://github.com/brendangregg/FlameGraph and add it to PATH"
fi

echo ""
echo "done. artifacts: perf.data, perf_report.txt$( [ -f flamegraph.svg ] && echo ', flamegraph.svg' )"
