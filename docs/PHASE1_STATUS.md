# AstraFlow Phase 1 Status

## Environment

See [environment.md](environment.md). Project-local CUDA 12.8.93 installed; SM120 kernel execution passed.

## Milestones

- [x] M0 Repository / Toolchain
- [x] M1 CPU Euler
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

M1 Release CPU: 3 CTest entries passed; 296 unit assertions and 407 Sod assertions. Sod nx=400, t=0.2: L1 density 0.00144607, axial velocity 0.00254609, pressure 0.000890535; mass error 3.33067e-15. No reconstruction corrections. Uniform state preserved. M0 native SM120 runtime gate passed previously.

## Execution plan

Follow M0–M9 in order; commit only at passed milestone gates, use targeted development tests and one comprehensive final pass. Keep physics/numerics independent of backend and UI. Measure actual errors and performance.

## Deferred

Only specified Phase 2 features: plume, chemistry, turbulence, advanced thermodynamics, CUDA/OpenGL interop and advanced nozzle optimisation.
