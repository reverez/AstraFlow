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

Phase 1 has no blocker. V1.1 refined-grid steady closure remains unresolved; see the V1.1 section below. The user configured origin as git@github.com:reverez/AstraFlow.git; SSH fetch succeeds. Remote README commits were merged without rewriting history, retaining the user's structure and authorship text.

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
- [x] C2: 128x32 CUDA FP64 steady closure at 82,480 iterations, 64.7 s wall; all residual/window/conservation gates pass, zero fallbacks.
- [ ] C3: blocked scientifically. CPU FP64 cross-check passes at 82,480 iterations (max field difference 1.20e-14). 256x64 reaches 250,000 steps without residual/mass-estimate closure; finer grids were not started.
- [ ] C4: analysis mathematics tests pass; numerical Richardson/GCI is not valid with only one accepted grid.
- [ ] C5: diagnostic closure report and four eligible baseline figures plus rejected-grid residual figure complete; grid-convergence figure withheld because the prerequisite failed.
- [ ] C6: final software pass is green (CPU 7/7, CUDA 9/9, analysis tests and output integrity); scientific acceptance remains false because C3 did not pass. Review branch published without merging main.

Scientific acceptance is unresolved. CFL sensitivity and an additional 50,000-step diagnostic continuation did not resolve the medium-grid residual. Existing cell-based mass estimate differs from the numerical boundary-flux balance; neither criterion was changed. Details and exact criteria are in [rocket_scientific_closure.md](verification/rocket_scientific_closure.md).

### V1.1-C2b residual investigation

- [x] True four-equation numerical boundary/global balances, signed local convective/transport/source decomposition, fixed-cell primitive/gradient/limiter traces and targeted regressions.
- [x] Conservative 128×32-to-256×64 initialization: maximum parent integral error 5.00e-16, no positivity limiting. Fresh CUDA FP64 production run ends at 250,000 steps with energy residual 4.57e-5 and cell mass mismatch 0.00117821; not converged.
- [x] Equal-time 100.397418 µs MC, Van Leer and first-order diagnostics from the same saved medium field; all diagnostic only. Numerical balance assembly error <=6.88e-15; signed near-axis residual cycles and neighboring limiter switches measured.
- [x] Evidence, three figures and full-cadence focus traces recorded in [residual_closure_c2b.md](verification/residual_closure_c2b.md). An implicit pseudo-transient proposal preserves R(U)=0; no acceleration implemented.
- [ ] Scientific acceptance remains unresolved: no accepted medium/fine pair or grid uncertainty. Accepted 128×32 and original 250k cold-medium cases were not rerun; 512×128 was not started. Production mathematics and all convergence gates remain unchanged.
