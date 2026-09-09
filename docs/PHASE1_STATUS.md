# AstraFlow Phase 1 Status

## Environment

See [environment.md](environment.md). Project-local CUDA 12.8.93 installed; SM120 kernel execution passed.

## Milestones

- [x] M0 Repository / Toolchain
- [x] M1 CPU Euler
- [x] M2 CUDA Euler
- [ ] M3 2D finite volumes
- [ ] M4 Axisymmetric nozzle
- [ ] M5 Viscous Navier–Stokes
- [ ] M6 Engineering analysis / output
- [ ] M7 Interactive GUI
- [ ] M8 Performance
- [ ] M9 Final verification

## Current blocker

GitHub CLI has no authenticated host; remote creation deferred. Local development is unblocked.

## Latest validation

M2 Release CUDA: all five CTest entries pass. GPU Sod/uniform and CPU/GPU parity pass in FP32/FP64 (4,008 assertions). FP32 max/L1/L2 differences: 4.61936e-6 / 2.49782e-7 / 5.32646e-7. FP64: 4.21885e-14 / 7.72314e-16 / 3.19405e-15. Matching contraction settings fixed the initially failed FP32 parity without changing tolerances. Native SM120 device probe passes.

## Execution plan

Follow M0–M9 in order; commit only at passed milestone gates, use targeted development tests and one comprehensive final pass. Keep physics/numerics independent of backend and UI. Measure actual errors and performance.

## Deferred

Only specified Phase 2 features: plume, chemistry, turbulence, advanced thermodynamics, CUDA/OpenGL interop and advanced nozzle optimisation.
