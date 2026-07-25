# Vortex-Transport

2D Discontinuous Galerkin solver for the compressible Euler equations on triangular elements, set up for the isentropic vortex transport benchmark.

## Building and running

```bash
make                       # builds ./Vortex-Transport
./Vortex-Transport input   # runs the simulation, writes grid/ and output/
```

Requires `g++` (C++11) and the Eigen headers (`/usr/include/eigen3`, e.g. `sudo apt-get install libeigen3-dev`).

## Profiling

Profiling is opt-in and does not affect the default build. It has three layers: `perf` sampling (where is the time spent, down to the line), an instrumented timing table (explicit per-function totals), and a peak-memory report.

### 1. Build with profiling enabled

```bash
make profile               # or: make PROFILE=1
```

This builds a separate executable, `Vortex-Transport-prof`, from separate object files, compiled with `-O2 -g -fno-omit-frame-pointer` (realistic optimization + accurate perf call stacks) and `-DENABLE_TIMERS` (enables the timing macros in `Profiling.H`). The regular `make` target is unchanged, and the two builds can coexist.

### 2. Run under perf

```bash
scripts/profile.sh                              # runs ./Vortex-Transport-prof input
scripts/profile.sh ./Vortex-Transport-prof input
```

The script produces:

- `perf.data` — raw samples; explore interactively with `perf report -i perf.data`
- `perf_report.txt` — the same report as text
- `flamegraph.svg` — only if `inferno-flamegraph` or `flamegraph.pl` is on your PATH (the script skips it otherwise and tells you how to install one)

If `perf` is missing, see the install instructions in the comments at the top of `scripts/profile.sh` (including WSL2-specific notes). If recording fails with a permission error, run `sudo sysctl kernel.perf_event_paranoid=1`.

For hardware counters (cache behavior, page faults) use `perf stat` instead of `perf record`:

```bash
perf stat -e cycles,instructions,cache-misses,page-faults ./Vortex-Transport-prof input
```

For detailed heap analysis beyond peak RSS, `valgrind --tool=massif` is the right tool (much slower, but attributes memory to allocation sites).

### 3. Interpreting the output

- **perf report / flamegraph** — statistical samples of where CPU time is spent, including code you did not instrument (e.g. `numerical_flux`, `libm`, allocator). Use the flamegraph for the big picture, `perf report` (and `perf annotate`) to drill into a specific function. This is the primary tool for deciding what to parallelize.
- **Timing tables** — the profiling binary prints two tables to stderr at exit. Every function in the codebase is instrumented; functions that were never called still appear with 0 calls and 0%.
  - **By self time** (exclusive: elapsed time minus time spent in nested instrumented calls) — answers "where is CPU time actually spent". Self percentages sum to ~100% in a single-threaded run; with multiple threads the sum can exceed 100% of wall time. Sort key for deciding what to optimize or parallelize.
  - **By inclusive time** (children included) — answers "which high-level code path is expensive". Inclusive percentages do not sum to 100% because nested calls overlap (`rk4` contains `evolve_elem`, which contains the `compute_*` kernels); this is expected. Also shows min/max time per call.
  - The timers are thread-safe by construction: each thread accumulates into thread-local storage (no locks on the hot path) and results are merged at program exit, so the numbers stay valid after parallelizing with OpenMP. For very hot micro-functions (e.g. `numerical_flux`, millions of calls) the timer itself adds some overhead — trust perf for those.
- **Peak memory** — the same summary reports `VmHWM` (peak resident set size) and `VmPeak` (peak virtual memory) read from `/proc/self/status`.

To instrument a new function, add `PROFILE_SCOPE("name");` at the top of the function body and `PROFILE_DECLARE("name");` at file scope (the declare makes the name show up in the table even when the function is never called). Both macros compile to nothing in the default build; see `Profiling.H`.

## Post-processing

```bash
python3 plottools/gridplot.py                  # plottools/grid.pdf
python3 plottools/hydrodynamicstateplot.py     # plottools/step_*.pdf + animation.mp4
python3 plottools/l2stateerror.py              # prints L2 error initial vs final state
```

Run these from the repository root; they need `numpy`, `pandas`, `matplotlib` and `ffmpeg` (for the animation).
