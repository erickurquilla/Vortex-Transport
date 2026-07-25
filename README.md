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
- **Timing table** — the profiling binary prints a table to stderr at exit: function name, call count, total seconds, and % of wall time, for the instrumented functions (the DG kernels in `Evolve.cpp`, the time steppers in `Timestepping.cpp`, and `Element::write_data`). These are exact inclusive wall times, useful for before/after comparisons when optimizing. Nested timers overlap (`rk4` contains `evolve_elem`, which contains the `compute_*` kernels), so percentages do not sum to 100.
- **Peak memory** — the same summary reports `VmHWM` (peak resident set size) and `VmPeak` (peak virtual memory) read from `/proc/self/status`.

To instrument additional functions, add `PROFILE_SCOPE("name");` at the top of the function (include `Profiling.H`). The macro compiles to nothing in the default build.

## Post-processing

```bash
python3 plottools/gridplot.py                  # plottools/grid.pdf
python3 plottools/hydrodynamicstateplot.py     # plottools/step_*.pdf + animation.mp4
python3 plottools/l2stateerror.py              # prints L2 error initial vs final state
```

Run these from the repository root; they need `numpy`, `pandas`, `matplotlib` and `ffmpeg` (for the animation).
