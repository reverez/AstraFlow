# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

AstraFlow is a 2D axisymmetric compressible Euler/Navier–Stokes finite-volume solver in C++20 and CUDA (target `sm_120`, RTX 5070 Laptop), with a deterministic CPU reference backend, a Catch2 verification suite, a CLI, an ImGui GUI and a Python analysis layer. It is a scientific-computing portfolio project run as an **architect-directed, AI-assisted workflow**: the user owns the mathematics, numerical method, acceptance criteria and review; agents implement under those specifications. Read `README.md` for the physics and the authorship model before touching solver code.

## Build

Development happens in WSL2 Ubuntu 24.04. CMake ≥ 3.24, Ninja, GCC 13. CMake FetchContent downloads pinned nlohmann/json, Catch2, GLFW, ImGui, ImPlot at first configure.

```sh
# CPU-only (what CI runs)
cmake -S . -B build-cpu -G Ninja -DCMAKE_BUILD_TYPE=Release -DASTRAFLOW_ENABLE_CUDA=OFF
cmake --build build-cpu -j6

# CUDA (project-local toolchain; no system CUDA install)
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CUDA_COMPILER="$PWD/.toolchains/cuda-12.8.1/bin/nvcc" \
  -DASTRAFLOW_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=120
cmake --build build -j6
./build/astraflow_cli --device-info     # runs a kernel and reports compiled arch

# Optional targets
-DASTRAFLOW_BUILD_GUI=ON          # astraflow_gui (run scripts/bootstrap_gui.sh first if X11/GL headers are missing)
-DASTRAFLOW_BUILD_BENCHMARKS=ON   # astraflow_benchmark, astraflow_residual_probe, astraflow_diagnose_closure
```

If `.toolchains/cuda-12.8.1` is missing: `python3 -m venv .venv && .venv/bin/python scripts/bootstrap_cuda.py` (SHA256-checked NVIDIA redist packages, no sudo). The CUDA build only works on the target GPU machine; the GPU is not reachable from some agent sandboxes.

Sanitizer configuration used for acceptance: CPU build with `-fsanitize=address,undefined` and leak detection (see `docs/verification/verification_report.md`).

Existing build trees: `build` (CUDA), `build-cpu`, `build-sanitize`. The `*-before-move` / `*-before-work-root-move` directories are stale copies from repository relocations and can be ignored.

## Tests

```sh
ctest --test-dir build --output-on-failure        # CUDA build: 10 entries
ctest --test-dir build-cpu --output-on-failure    # CPU build: 8 entries

# Single CTest entry (names: unit, grid_refinement, viscous_verification, nozzle_verification,
# verification, cuda_verification, cuda_device, cli_help, cli_version, closure_analysis)
ctest --test-dir build -R nozzle_verification --output-on-failure

# Single Catch2 case or tag directly
./build/astraflow_tests "[steady]"
./build/astraflow_tests "HLLC contact and expansion fallback"
./build/astraflow_tests --list-tests

# Python analysis tests (run in .venv; system Python skips the 3 NumPy-dependent cases)
.venv/bin/python -m unittest discover -s tests/analysis -v
```

Catch2 tags: `[unit]`, `[verification]`, `[nozzle]`, `[viscous]`, `[refinement]`, `[cuda]`, `[balance]`, `[convergence]`, `[steady]`, `[parallel]`, `[wall]`, `[io]`. Anything tagged `[cuda]` is only compiled into the CUDA build.

OpenMP: libgomp's default spin-wait starves this hybrid-core WSL2 laptop (16 threads were 20× slower than 1). Run anything multi-threaded with `OMP_WAIT_POLICY=passive OMP_PROC_BIND=true OMP_PLACES=cores` and at most 8 threads; ctest sets this for the `unit` entry. The CLI defaults to one thread; `--threads N` opts in.

## Running

```sh
./build/astraflow_cli --config examples/rocket_nozzle/config.json --backend cuda \
    [--precision float|double] [--output runs/name] [--max-iterations N] [--initial-state FILE]
./build/astraflow_gui --config examples/rocket_nozzle/config.json [--smoke-test] [--capture out.ppm]
.venv/bin/python scripts/check_run.py runs/<dir>     # stdlib-only JSON/CSV/VTK validation
```

Examples: `sod`, `isentropic_nozzle`, `viscous_channel`, `rocket_nozzle/config.json` (3,000-step FP32 demo), `rocket_nozzle/scientific_closure.json` (V1.1 MC study, FP64), `rocket_nozzle/v12_closure.json` (V1.2 van Albada study with the numerical mass gate) and `rocket_nozzle/cooled_wall.json` (isothermal 600 K wall, Sutherland viscosity, pseudo-time).

Study drivers: `scripts/run_scientific_closure.py --config examples/rocket_nozzle/v12_closure.json --root runs/v12-closure --grids 128x32 256x64 512x128 --iteration-limit ...`, `scripts/run_pseudo_root_study.py run|analyze`, `scripts/sweep_nozzle.py`, `scripts/plot_readme_figures.py {sod,field,scaling,wall,sweep} SOURCE`, `scripts/make_banner.py`. An explicit `--output` directory must be new or empty; the writer refuses to overwrite.

Python: `.venv` (WSL Python 3.12). `requirements-analysis.txt` is matplotlib only, needed for figures; orchestration and statistics scripts use the standard library.

## Architecture

Single library target `astraflow` (`include/astraflow/`, `src/`) linked by `astraflow_cli`, `astraflow_gui`, `astraflow_tests` and the benchmark/diagnostic executables. Dependency direction is strict and should stay that way:

```
physics/state.hpp  →  numerics/{flux,euler1d_ops,finite_volume}.hpp  →  cpu/ , cuda/  →  core/Simulation  →  analysis/ , io/  →  apps/
```

**Shared header-only numerics, instantiated for float and double on both backends.** `physics/state.hpp` (ideal gas, `State<T>`, SoA `View<T>`), `numerics/flux.hpp` (limiters, rotated HLLC with HLLE fallback) and `numerics/finite_volume.hpp` (boundary ghost states, MUSCL slopes, corrected transport gradients, viscous flux, residual with axisymmetric source, CFL timestep) are plain `AF_HD` (`__host__ __device__`) templates. `src/cpu/solver.cpp` and `src/cuda/solver.cu` call *the same functions*; they differ only in scheduling and memory. This is what makes CPU/CUDA parity at ~1e-14 (FP64) a meaningful test. Any change to a numerics header changes both backends at once, and parity tests must still pass. `src/*/euler1d.*` is the older 1D Euler path kept for Sod verification.

**Backends.** `Solver` is the abstract interface (`initialize`, `step(max_dt)`, `state`, `stats`). `StepStats` carries the four RMS residuals and the numerical inlet/outlet mass fluxes from the same spatial evaluation. The CUDA step makes exactly one device-to-host transfer (one `Diagnostics` record: sums, timestep, boundary fluxes, invalid index, fallback counts); keep it that way when adding diagnostics. The CPU backend has OpenMP loops gated by `parallel_threshold` (2,048 cells) and the opt-in thread count; residual and boundary sums are serial so results are thread-count independent. `make_cpu_solver` / `make_cuda_solver` choose precision and `SteadySettings`. Each backend has two step paths: physical SSP-RK2 (default) and optional local pseudo-time RK2 (`steady_solver.enabled`), which reuses the identical spatial operator and reports `R(U)` at the completed iterate via a third residual evaluation. The CUDA backend uses RAII `cuda::Buffer`/`cuda::Event` (`src/cuda/memory/resources.cuh`), allocates nothing inside the step loop, reduces timestep and the four residual sums with CUB, and surfaces invalid-state cell indices and fallback counts as diagnostics. Optional CUDA-event stage profiling is behind `numerics.profile`.

**`core/Simulation` owns nondimensionalisation.** It builds the physical mesh, derives `Scales` (length = throat radius or channel height; p0, T0, sqrt(R T0)), rescales mesh, state, gas and boundary settings into solver units, selects the backend, and restores physical units in `state()`/`stats()`. Consequence: **residuals in CSV/summary are RMS time derivatives in nondimensional solver units; fields, times and engineering results are in configuration units.** Do not mix them. `Simulation::step` also samples engineering integrals, the numerical boundary-flux balance and the `ConvergenceMonitor` every `convergence.sampling_interval` iterations, and `termination()` distinguishes `steady_converged`, `iteration_limit`, `time_limit`, `user_stop`, `invalid_state`.

**Config** (`io/config.cpp`) parses JSON strictly: unknown keys in any section throw; every numeric is validated. `Limiter::FirstOrderDiagnostic` is deliberately unreachable from config (diagnostic executables only). Limiters: `minmod`, `van_leer`, `mc`, `van_albada`. Transport: `gas.viscosity_model` = `constant` | `sutherland` | `power_law` with `reference_temperature`, `sutherland_temperature`, `viscosity_exponent`. Convergence: `mass_gate` = `interior_estimate` (V1 legacy default) | `numerical_boundary`, and `balance_interval` (0 = host balance only at start/finish).

**Output** (`io/output.cpp`, `RunOutput`): per run `config.json`, `metadata.json`, `mesh.vts`, `convergence.csv`, `engineering.csv`, `balance.jsonl`, `initial_balance.json`, `final_balance.json`, `final_state.vts`, `summary.json`, `performance.json`, `wall.csv` for viscous no-slip cases (+ `initialization.json` with `--initial-state`). `export_state`/`import_state` implement the schema-1 state files used by `--write-state`/`--initial-state`. VTK is written by hand, no VTK library. Default location is a timestamped directory under git-ignored `runs/`.

**Analysis** (`analysis/`): `engineering.cpp` (mass flow, thrust, Isp, exit state from physical snapshots), `convergence.cpp` (rolling-window steady monitor: residual thresholds, observable relative range, mass-mismatch gates), `balance.cpp` (true four-equation boundary/source budgets and conservative 2:1 prolongation used by V1.1 diagnostics), `wall.cpp` (wall traction and heat flux from the production viscous flux; integrated loads), `fields.cpp` (GUI scalar fields).

**GUI** (`apps/gui/main.cpp` + `core/Controller`): a `std::jthread` worker owns the `Simulation`; it publishes immutable `Snapshot` objects (mesh, physical state, gas, stats, history) every ten iterations and at Run/Pause/Step/Reset boundaries. The renderer never touches live CUDA memory, and no physics lives in GUI code. `--smoke-test --capture` is the scripted verification path.

**Scripts** (`scripts/`): `run_scientific_closure.py` (sequential, lock-protected, resumable grid study that reuses verified attempts via config/binary/output SHA256), `analyze_scientific_closure.py` (Richardson/GCI, window statistics, CPU/CUDA comparison), `run_residual_closure.py` / `analyze_residual_closure.py` (C2b), `run_steady_acceleration.py` / `analyze_steady_acceleration.py` (C2c), `plot_scientific_closure.py`, `residual_probe.cpp` and `diagnose_closure.cpp` (C++ diagnostics built with benchmarks ON).

## Floating-point determinism

`-ffp-contract=off` (CXX) and `--fmad=false` (CUDA) are set on purpose so CPU and GPU round identically; FP32 parity bounds were met only after this. Keep them. Diagnostic executables set the same flag.

## Scientific-integrity rules that constrain code changes

These come from the README, `docs/PHASE1_STATUS.md` and the `docs/verification/*.md` reports and are treated as hard requirements, not style:

- Acceptance bounds and convergence gates are preset before a run. **Never loosen a tolerance, change a gate, or redefine an estimator after seeing a result** to make it pass. If a criterion is believed wrong, document the evidence and propose a pre-registered revision as a new milestone.
- Altered schemes (Van Leer, first-order, pseudo-time, warm starts) are **diagnostic** unless they pass the unchanged gates and any extra equivalence criteria. Do not report them as accepted solutions.
- Never rerun or overwrite a recorded scientific case to regenerate a report; the analysis scripts read preserved outputs. Runners refuse to overwrite and use new attempt/prefix directories.
- Reports distinguish production residual (second RK stage) from probe residual `R(U)` at the completed state, interior-cell mass estimate from numerical boundary-flux mismatch, and physical time from pseudo-time. Keep those distinctions explicit in any new output or doc.
- Do not claim Nsight/compute-sanitizer results, experimental validation, or grid independence; none exist.
- Numerical results are accepted only against analytical solutions, conservation, convergence order or CPU/GPU parity. Compiling and producing plausible fields is not acceptance.

## Repository state (as of 2026-10-21, V2.0)

- `main` = tag `v1.0.0`. All later work is on `v1.1-scientific-closure`; `docs/STATUS.md` is the live status record and `docs/PHASE1_STATUS.md` the V1.0/V1.1 history.
- V1.1 (MC limiter) closed only the 128×32 rocket grid. V1.2/V1.3 (van Albada, numerical mass gate, pre-registered protocols in `docs/verification/v12_closure_protocol.md`, `v13_continuation_protocol.md`) closed 128×32 (49,480 it.) and 256×64 (760,760 it.) and supported the pseudo-time same-root hypothesis; 512×128 reached its 2.5M cap still decaying and is continued under `v14_fine_continuation_protocol.md`. Results: `docs/verification/rocket_closure_v12.md`. Preserved runs live under git-ignored `runs/v12-closure`, `runs/v13-closure`, `runs/v14-closure`, `runs/v12-pseudo-root`.
- Never rerun a recorded study case to regenerate a report; the analysis and plot scripts read the preserved directories. Budgets/caps for a study are fixed in its protocol before launch; a new budget means a new pre-registered protocol (that is how V1.3 and V1.4 came about).
- Not done, deliberately: plume, chemistry, turbulence, real gas, multi-GPU, CUDA/OpenGL interop (no verification path yet).

## Conventions

- Formatting: `.clang-format` (LLVM base, 4-space indent, 100 columns). Code is dense; functions are short and comments are reserved for non-obvious numerical semantics.
- Commit messages follow `feat:`/`test:`/`docs:`/`perf:`/`chore:` prefixes with one milestone-sized commit per gate.
- Docs under `docs/` are the record of acceptance; update `docs/PHASE1_STATUS.md` and the relevant verification report when a milestone or diagnostic result changes.
- If running tools from the Windows side against this WSL path (`\\wsl.localhost\...`), git refuses the repository as "dubious ownership"; run git inside the distro, e.g. `wsl.exe -d Ubuntu-24.04 -e bash -lc 'cd "/home/arnav/Arnavs Folder/dev/projects/portfolio/AstraFlow" && git status'`.
