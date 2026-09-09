# AstraFlow Phase 1 Status

## Environment

See [environment.md](environment.md). Project-local CUDA 12.8.93 installed; SM120 kernel execution passed.

## Milestones

- [x] M0 Repository / Toolchain
- [x] M1 CPU Euler
- [x] M2 CUDA Euler
- [x] M3 2D finite volumes
- [ ] M4 Axisymmetric nozzle
- [ ] M5 Viscous Navier–Stokes
- [ ] M6 Engineering analysis / output
- [ ] M7 Interactive GUI
- [ ] M8 Performance
- [ ] M9 Final verification

## Current blocker

GitHub CLI has no authenticated host; remote creation deferred. Local development is unblocked.

## Latest validation

M3 Release CPU/CUDA builds and complete milestone CTest suites pass. Rectangular mesh validity/face closure, 2D periodic uniform preservation, conservative density-wave evolution and FP64 parity pass. All four integrated conserved quantities drift by less than 1e-12. Earlier Euler gates remain green.

## Execution plan

Follow M0–M9 in order; commit only at passed milestone gates, use targeted development tests and one comprehensive final pass. Keep physics/numerics independent of backend and UI. Measure actual errors and performance.

## Deferred

Only specified Phase 2 features: plume, chemistry, turbulence, advanced thermodynamics, CUDA/OpenGL interop and advanced nozzle optimisation.
