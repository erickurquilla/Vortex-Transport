#!/usr/bin/env bash
#
# Profile Vortex-Transport with linux perf.
#
# Usage:
#   scripts/profile.sh                          # runs ./Vortex-Transport-prof input
#   scripts/profile.sh ./Vortex-Transport-prof my_input
#
# Run from the repository root. Build the profiling binary first:
#   make profile        (or: make PROFILE=1)
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
# Complementary hardware counters (not covered by this script):
#   perf stat -e cycles,instructions,cache-misses,page-faults ./Vortex-Transport-prof input
# ---------------------------------------------------------------------------

set -euo pipefail

BINARY="${1:-./Vortex-Transport-prof}"
shift || true
ARGS=("${@:-input}")

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
    echo "Build it first with:  make profile" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# record: -g captures call graphs (works well because the profiling build
# uses -fno-omit-frame-pointer), -F 499 samples ~499 times/second
# ---------------------------------------------------------------------------
echo "recording: $PERF record -g -F 499 -o perf.data -- $BINARY ${ARGS[*]}"
"$PERF" record -g -F 499 -o perf.data -- "$BINARY" "${ARGS[@]}"

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
