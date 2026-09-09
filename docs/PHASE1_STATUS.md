# AstraFlow Phase 1 Status

## Environment

See [environment.md](environment.md). Project-local CUDA 12.8.93 installed; SM120 kernel execution passed.

## Milestones

- [x] M0 Repository / Toolchain
- [ ] M1 CPU Euler
- [ ] M2 CUDA Euler
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

Release CPU and CUDA configure/build passed; both CTest suites passed (2 smoke tests each). RTX 5070 Laptop GPU, CC 12.0, 8,518,041,600 bytes, runtime 12080; device kernel reported sm_120. No scientific claims yet.

## Execution plan

Follow M0–M9 in order; commit only at passed milestone gates, use targeted development tests and one comprehensive final pass. Keep physics/numerics independent of backend and UI. Measure actual errors and performance.

## Deferred

Only specified Phase 2 features: plume, chemistry, turbulence, advanced thermodynamics, CUDA/OpenGL interop and advanced nozzle optimisation.
