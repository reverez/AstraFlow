# AstraFlow Phase 1 Status

## Environment

See [environment.md](environment.md). Project-local CUDA 12.8.93 installed; SM120 kernel execution passed.

## Milestones

- [x] M0 Repository / Toolchain
- [x] M1 CPU Euler
- [x] M2 CUDA Euler
- [x] M3 2D finite volumes
- [x] M4 Axisymmetric nozzle
- [ ] M5 Viscous Navier–Stokes
- [ ] M6 Engineering analysis / output
- [ ] M7 Interactive GUI
- [ ] M8 Performance
- [ ] M9 Final verification

## Current blocker

GitHub CLI has no authenticated host; remote creation deferred. Local development is unblocked.

## Latest validation

M4 complete CUDA-build suite: 6/6 CTest entries pass, including CPU and GPU checks. Axis rest drift 2.22045e-15 after fixing inlet root precision. Nozzle 96x12, 4,000 steps: throat Mach 1.00957, mean Mach error 0.000319245, mass-flow spread 0.000774078. CPU/GPU nozzle parity passes 1e-10. No invalid states or relaxed tolerances.

## Execution plan

Follow M0–M9 in order; commit only at passed milestone gates, use targeted development tests and one comprehensive final pass. Keep physics/numerics independent of backend and UI. Measure actual errors and performance.

## Deferred

Only specified Phase 2 features: plume, chemistry, turbulence, advanced thermodynamics, CUDA/OpenGL interop and advanced nozzle optimisation.
