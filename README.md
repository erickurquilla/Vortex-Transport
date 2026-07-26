# Vortex-Transport

2D nodal discontinuous Galerkin (DG) solver for the compressible Euler equations on a structured triangulation, configured for the classic isentropic vortex transport problem.

The conserved state is \(\mathbf{U} = (\rho, \rho u, \rho v, \rho E)\). Each triangular element stores nodal degrees of freedom for a Lagrange basis of order `p` on the reference triangle with vertices \((0,0)\), \((1,0)\), \((0,1)\). Interface coupling uses a Roe numerical flux (`numerical_flux` in `Numericalflux.cpp`). Time integration is explicit forward Euler or classical RK4.

---

## Key features

**What the code does today**

- Solves the 2D compressible Euler equations with \(\gamma = 1.4\) (hardcoded).
- Nodal DG on triangles: mass and stiffness operators built once in reference space (`Preevolve.cpp`), then mapped to physical space per element via the affine Jacobian (`Element::build_jacobians`, `build_mass_matrix_inverse`, `build_stiffness_matrix`).
- Structured mesh: each Cartesian cell is split into two triangles (types 0 / 1: right-angle vertex down / up). Domain is centered at the origin, \([-\texttt{domain\_x}/2, +\texttt{domain\_x}/2] \times [-\texttt{domain\_y}/2, +\texttt{domain\_y}/2]\).
- Neighbor connectivity wraps at the outer edges (periodic-style boundaries); see `Meshgeneration.cpp`.
- Optional interior-node grid perturbation via `perturb_grid` (boundary nodes stay fixed).
- Isentropic vortex initial condition, hardcoded in `Element::initialize_hydrodinamics` (parameters such as \(\varepsilon=0.3\), \(M_\infty=0.5\), \(U_\infty=V_\infty=1/\sqrt{2}\), center \((0,0)\)).
- Two ways to place the IC on the DG nodes: direct nodal interpolation (`U_initialization_type=0`) or \(L^2\) / least-squares projection (`U_initialization_type=1`).
- Explicit time stepping: forward Euler (`stepping_method=0`) or RK4 (`stepping_method=1`).
- Opt-in OpenMP over element-level loops (`ENABLE_OPENMP=TRUE`), with `schedule(runtime)` so `OMP_SCHEDULE` is honored.
- Opt-in `perf` + instrumented timing build (`make profile`); archived profile snapshots live under `performance_history/`.
- Python post-processing: mesh plot, hydrodynamic field plots / animation, and an L2 state-error comparison of first vs last dump.

**What it does not do (current limits)**

- Only polynomial orders `p = 0,1,2,3` (hardcoded `switch` in `Lagrangebasis.cpp`).
- Only the isentropic vortex IC; vortex constants and \(\gamma\) are not input-file parameters.
- No adaptive / CFL-based \(\Delta t\): `time_step = simulation_time / number_time_steps`.
- No MPI / distributed memory.
- Default `make` binary is compiled without `-O2` (optimization is enabled on the profiling binary only).
- Quadrature `integration_order` must be in `0..19` (lookup tables in `Quadraturerule.cpp`).

---

## Requirements

| Dependency | Role |
|---|---|
| `g++` with C++11 (`-std=c++11`) | Compiler used by the `makefile` |
| Eigen 3 headers at `/usr/include/eigen3` | Mass-matrix inverse in `Preevolve.cpp` only (`#include <Eigen/Dense>`). Install e.g. `sudo apt-get install libeigen3-dev` |
| OpenMP (optional) | Enabled with `make ENABLE_OPENMP=TRUE` (`-fopenmp -DVORTEX_USE_OPENMP`) |
| Linux `perf` (optional) | Used by `scripts/profile.sh` |
| Python 3 + `numpy`, `pandas`, `matplotlib` (optional) | `plottools/` |
| `ffmpeg` (optional) | Animation in `hydrodynamicstateplot.py` |

The makefile hardcodes `-I/usr/include/eigen3`. There is no CMake build.

---

## Build

```bash
# Serial build (default) → ./Vortex-Transport
make

# OpenMP build → ./Vortex-Transport
make clean          # required when switching OpenMP on/off (same .o names)
make ENABLE_OPENMP=TRUE

# Profiling binary → ./Vortex-Transport-prof  (-O2 -g -fno-omit-frame-pointer -DENABLE_TIMERS)
make profile
make profile ENABLE_OPENMP=TRUE   # profile the OpenMP configuration

make clean          # also removes grid/, output/, perf artifacts, plottools/*.pdf|*.mp4
```

Notes from the makefile:

- `ENABLE_OPENMP` defaults to `FALSE`. Plain `make profile` produces a **serial** profiling binary; use `ENABLE_OPENMP=TRUE` if you want to profile the parallel code.
- The short git hash is baked into `main` at compile time (`GIT_COMMIT_HASH`) and printed at startup.
- Regular objects use `CFLAGS = -std=c++11 -Wall` only; `-O2` is part of `PROFILE_FLAGS` for `*.prof.o` / `Vortex-Transport-prof`.

---

## Usage

```bash
./Vortex-Transport <input_file>
# examples:
./Vortex-Transport input_files/input
./Vortex-Transport input_files/input_omp_performance
```

Exactly one CLI argument is required (`Parameters.cpp`). On start, the program prints the git commit, OpenMP/CPU diagnostics, then runs.

### Input file format

`key=value` lines. Comments are lines without `=` (or text after a value on the same line is kept as part of the value string for some keys — prefer comments on their own lines as in the shipped files). Keys read by `read_input_files`:

| Key | Meaning |
|---|---|
| `p` | Lagrange polynomial order (supported: 0–3) |
| `domain_x`, `domain_y` | Full domain width/height (cm in the sample comments); mesh spans \(\pm\) half of each |
| `num_element_in_x`, `num_element_in_y` | Cartesian cells; total triangles = `2 * num_element_in_x * num_element_in_y` |
| `perturb_grid` | Interior-node mesh perturbation strength (0 = Cartesian) |
| `integration_order` | Index into Gauss quadrature tables (0–19) for line and area integrals |
| `simulation_time` | Final physical time \(T\) |
| `number_time_steps` | Number of steps \(N\); \(\Delta t = T/N\) (computed, not read as its own key) |
| `write_every_steps` | Write `output/step_<k>/` when `k % write_every_steps == 0` (plus always step 0) |
| `U_initialization_type` | `0` direct interpolation; `1` least-squares projection |
| `stepping_method` | `0` forward Euler; `1` RK4 |

Shipped inputs:

- `input_files/input` — `p=3`, \(8\times8\) cells (128 triangles), \(T \approx 14.14213562373\), 10000 RK4 steps, write every 250. That \(T\) equals \(10\sqrt{2}\), i.e. one diagonal transit of a vortex with speed \(1/\sqrt{2}\) across a side-10 domain (the usual “return to start” check under periodic wrap).
- `input_files/input_omp_performance` — large mesh (`256×256` cells), tiny \(T\) and 10 steps; intended for OpenMP / profiling throughput tests, not a full vortex period.

### Program flow (`main.cpp`)

1. **Initialize** — read parameters; `generate_mesh` (writes `grid/element_*.txt`); build line/area Gauss rules; reference nodes; reference inverse mass matrix and stiffness matrix; construct each `Element` (Jacobians, physical mass/stiffness, IC); write step 0; construct each `Evolve_element` and pre-evaluate basis values at face quadrature points.
2. **Time loop** — for `a = 1 .. number_time_steps`: call `forward_euler` or `rk4`. Each stage of those methods calls `evolve_elem`, which for every element computes face states \(U^\pm\), Roe flux, face flux integral, stiffness vector, residual \(\mathbf{R} = \mathbf{S} - \int \hat{\mathbf{F}}\), then \(\dot{\mathbf{U}} = M^{-1}\mathbf{R}\). After the step, update nodal \(\mathbf{U}\) and recompute the Euler flux \(\mathbf{F}(\mathbf{U})\).
3. **Output** — when due, `clean_create_directory("output/step_<a>")` then `write_output` → one `element_<id>.txt` per triangle. Step 0 also writes `JMS_element_*.txt` (Jacobians, inverse mass, stiffness, outward unit normals).
4. **Timing** — prints greppable `[PHASE]` / `[TIMING]` lines for init, timestepping, teardown, and total wall time.

### Output layout

```
grid/element_<i>.txt          # mesh: vertices, type, right/left/vertical neighbors
output/step_<k>/element_<i>.txt
# header: node_number time x y u0 u1 u2 u3 fx0..fx3 fy0..fy3
# u0..u3 = ρ, ρu, ρv, ρE ; fx*/fy* = Euler flux components
output/step_0/JMS_element_<i>.txt   # operators / normals (step 0 only)
```

Each run recreates `output/` and `grid/` from scratch (`clean_create_directory` uses `rm -rf`).

---

## Directory / source map

| Path | Role |
|---|---|
| `main.cpp` | Driver: diagnostics, init, time loop, phase timings |
| `Parameters.{H,cpp}` | `parameters` struct + `read_input_files` |
| `Meshgeneration.{H,cpp}` | Structured triangulation, neighbor wrap, reference nodes, `grid/` I/O |
| `Element.{H,cpp}` | Per-element geometry, operators, IC, state \(\mathbf{U}\)/\(\mathbf{F}\), write |
| `Lagrangebasis.{H,cpp}` | \(\phi_i(\xi,\eta)\), \(\nabla\phi_i\), reference→physical map (`p≤3`) |
| `Quadraturerule.{H,cpp}` | Gauss line/area rules vs `integration_order` |
| `Preevolve.{H,cpp}` | Reference \(M^{-1}\) (via Eigen) and reference stiffness |
| `Evolve.{H,cpp}` | Per-element DG residual / \(\dot{\mathbf{U}}\) pipeline |
| `Numericalflux.{H,cpp}` | Roe flux (in-place, no heap alloc on the hot path) |
| `Timestepping.{H,cpp}` | `evolve_elem`, `forward_euler`, `rk4` |
| `Utilities.{H,cpp}` | File/dir helpers, parallel `write_output` |
| `Profiling.H` | Opt-in `PROFILE_SCOPE` / `PROFILE_DECLARE` timers + peak RSS |
| `makefile` | Serial / OpenMP / profile targets |
| `input_files/` | Example parameter files |
| `scripts/profile.sh` | `perf record` + report (+ optional flamegraph) |
| `performance_history/` | Manually archived profile run logs |
| `plottools/` | Mesh / state plots, animation, L2 error |
| `grid/`, `output/` | Runtime artifacts (gitignored; deleted by `make clean`) |

---

## Parallelization (OpenMP)

Enabled only when built with `ENABLE_OPENMP=TRUE` (defines `VORTEX_USE_OPENMP` and links `-fopenmp`).

**Parallel loops** (all `#pragma omp parallel for default(shared) schedule(runtime)`):

| Location | What is parallelized |
|---|---|
| `main.cpp` | Element construction / Jacobians / operators / IC |
| `main.cpp` | `Evolve_element` construction + `evaluate_basis_in_quadrature_poits` |
| `Timestepping.cpp` → `evolve_elem` | Per-element DG residual / \(\dot U\) (hot path each RK stage) |
| `Utilities.cpp` → `write_output` | Per-element `write_data` |

**Still serial:** mesh generation, reference-space operator assembly, and the RK4/Euler loops that advance nodal \(\mathbf{U}\) and refresh \(\mathbf{F}\) after `evolve_elem`.

At startup, OpenMP builds print `_OPENMP`, `omp_get_max_threads()`, a probed active thread count, and the environment variables below.

### Controlling OpenMP at runtime

Because loops use `schedule(runtime)`, the schedule is taken from the environment:

```bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DYNAMIC=false
export OMP_SCHEDULE=static          # or static,1 / dynamic / guided, ...
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_WAIT_POLICY=passive      # optional: sleep at barriers instead of spin (cleaner perf profiles)

make clean && make ENABLE_OPENMP=TRUE
./Vortex-Transport input_files/input_omp_performance
```

`main` reports `OMP_NUM_THREADS`, `OMP_PROC_BIND`, `OMP_PLACES`, `OMP_DYNAMIC`, `OMP_SCHEDULE`, and `OMP_MAX_ACTIVE_LEVELS` (or `"not set"`).

---

## Performance and profiling

Profiling is opt-in and does not change the default `make` binary. There are three layers: `perf` sampling, instrumented timing tables (`Profiling.H`), and a peak-memory report (`VmHWM` / `VmPeak` from `/proc/self/status`).

### 1. Build with profiling enabled

```bash
make clean
make profile ENABLE_OPENMP=TRUE     # or omit ENABLE_OPENMP for a serial profile build
```

Produces `Vortex-Transport-prof` from `*.prof.o` with `-O2 -g -fno-omit-frame-pointer -DENABLE_TIMERS -pthread`.

### 2. Run under perf

```bash
scripts/profile.sh                                          # defaults: ./Vortex-Transport-prof input_files/input
scripts/profile.sh ./Vortex-Transport-prof input_files/input_omp_performance

OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores \
  OMP_DYNAMIC=false OMP_SCHEDULE=static OMP_MAX_ACTIVE_LEVELS=1 \
  scripts/profile.sh
```

Artifacts in the repo root: `perf.data`, `perf_report.txt`, and optionally `flamegraph.svg` (if `inferno-flamegraph` or `flamegraph.pl` is on `PATH`). Install / WSL2 notes for `perf` are in the header comments of `scripts/profile.sh`. If recording fails on permissions: `sudo sysctl kernel.perf_event_paranoid=1`.

Complementary:

```bash
perf stat -e cycles,instructions,cache-misses,page-faults ./Vortex-Transport-prof input_files/input
# heap attribution (slow): valgrind --tool=massif ./Vortex-Transport-prof ...
```

### 3. Interpreting the output

- **perf / flamegraph** — statistical CPU samples, including uninstrumented code. Primary guide for what to optimize or parallelize. OpenMP outlined regions appear as `func(...) [clone ._omp_fn.N]`; barrier spin shows up in libgomp (install `libgomp1-dbgsym` for symbol names). See comments in `scripts/profile.sh`.
- **Timing tables** (stderr at exit of the profiling binary) — every `PROFILE_DECLARE`d function appears, including never-called ones at 0. **Self** time ≈ where CPU is spent; **inclusive** time ≈ which high-level path is expensive (nested scopes overlap, so inclusive % does not sum to 100%). Timers are thread-local and merged at exit (safe with OpenMP). For micro-functions with millions of calls (`numerical_flux`), prefer perf — timer overhead is measurable.
- **Peak memory** — `VmHWM` / `VmPeak` in the same summary.

To instrument a new function: `PROFILE_DECLARE("name");` at file scope and `PROFILE_SCOPE("name");` at the top of the body. Both macros are no-ops without `-DENABLE_TIMERS`.

### `performance_history/`

Hand-archived copies of profile runs (not written automatically by `profile.sh`), named

```text
performance_history/profile_<git_short_hash>_<MM_DD_YYYY>_<HH_MM>.txt
```

They typically contain the program stdout (steps, `[TIMING]`), the instrumented timer tables, and a truncated `perf report`. Useful for comparing commits after hot-path changes.

---

## Post-processing

Run from the repository root after a simulation has produced `grid/` and `output/`:

```bash
python3 plottools/gridplot.py                 # → plottools/grid.pdf
python3 plottools/hydrodynamicstateplot.py    # → plottools/step_*.pdf + animation.mp4
python3 plottools/l2stateerror.py             # prints L2 state error, first dump vs last
```

Needs `numpy`, `pandas`, `matplotlib`; the animation path also needs `ffmpeg`. `l2stateerror.py` currently hardcodes `p = 3` and a domain-normalization factor consistent with the default \(10\times10\) setup.

---

## Known limitations / caveats

- **Hardcoded physics constants** — \(\gamma\), vortex amplitude/Mach/freestream, and vortex center live in `Element::initialize_hydrodinamics` and are duplicated in the Euler flux updates inside `Timestepping.cpp` / `Numericalflux.cpp`.
- **Supported `p` and quadrature** — see “Key features”; unsupported values `exit` with an error message.
- **Default build is unoptimized** — use `Vortex-Transport-prof` or add your own `-O2`/`-O3` if comparing wall times for production-style runs.
- **OpenMP coverage is partial** — residual evaluation is parallel; RK stage state updates are not.
- **I/O wipes previous results** — every run deletes and recreates `output/` (and `generate_mesh` recreates `grid/`).
- **No CONTRIBUTING / multi-package install** — single Makefile, Eigen path fixed to `/usr/include/eigen3`.
