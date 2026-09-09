# AstraFlow Phase 1 Status

## Environment

See [environment.md](environment.md). Project-local CUDA 12.8.93 installed; SM120 kernel execution passed.

## Milestones

- [x] M0 Repository / Toolchain
- [x] M1 CPU Euler
- [x] M2 CUDA Euler
- [x] M3 2D finite volumes
- [x] M4 Axisymmetric nozzle
- [x] M5 Viscous Navier–Stokes
- [x] M6 Engineering analysis / output
- [x] M7 Interactive GUI
- [x] M8 Performance
- [x] M9 Final verification

## Current blocker

None. The user configured origin as git@github.com:reverez/AstraFlow.git; SSH fetch succeeds. Remote README commits were merged without rewriting history, retaining the user's structure and authorship text.

## Latest validation

Final Release CPU 6/6 and CUDA 8/8 CTest entries pass. CPU ASan/UBSan with leak detection passes 6/6 outside the tracing sandbox. Entropy-wave refinement orders are 1.81/2.00; Couette temperature orders 2.00/2.00. CUDA long-time Couette passes in both precisions. The perturbed nozzle recovers throat Mach 0.803 -> 0.997.

Both 3,000-step 128x32 rocket CLI runs and independent JSON/CSV/VTK validation pass. CUDA FP32 estimates thrust 171.05877 N with 0.04655% inlet/exit mass-flow mismatch; termination is the iteration limit, not steady convergence. Actual WSLg GUI control/rendering smoke and inspected capture pass.

Equal-FP32 benchmarks cover 128x32 through 1024x256. Largest GPU 1.98721 ms/iteration versus single-thread CPU 121.539 ms, measured 61.16x; device buffers 101,831,975 bytes. Event profiling and raw before/after measurements are in docs/performance/. Nsight and compute-sanitizer were unavailable.

## Completion

M0–M9 gates passed. Development continued at the user-relocated `/home/arnav/dev/projects/portfolio/AstraFlow`; the user commit and old build trees were preserved. Final results are consolidated in the verification and benchmark reports.

## Deferred

Only specified Phase 2 features: plume, chemistry, turbulence, advanced thermodynamics, CUDA/OpenGL interop and advanced nozzle optimisation.

## V1.1 scientific closure

- [x] C0: annotated v1.0.0 preserves ecd0337; development on v1.1-scientific-closure; main unchanged; baseline Release suite passes.
- [x] C1: convergence monitor, full engineering window, explicit termination/output and dedicated FP64 configuration; CPU 6/6, CUDA 8/8, legacy 3,000-step iteration-limit check and analysis unit tests pass.
- [ ] C2: 128x32 rocket steady closure.
- [ ] C3: independently converged rocket grid study and CPU reference cross-check.
- [ ] C4: Richardson/GCI analysis where valid.
- [ ] C5: scientific figures and closure report.
- [ ] C6: final acceptance pass and branch publication (no merge).

Next: measure baseline FP64 convergence, then refine only after all steady criteria pass. Details and exact criteria are in [rocket_scientific_closure.md](verification/rocket_scientific_closure.md).
