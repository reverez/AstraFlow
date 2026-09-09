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
- [ ] M7 Interactive GUI
- [ ] M8 Performance
- [ ] M9 Final verification

## Current blocker

GitHub CLI has no authenticated host; remote creation deferred. Local development is unblocked.

## Latest validation

M6 Release CPU/CUDA full suites pass. Real 128x32 viscous rocket CLI runs completed 300 steps on CPU FP64 and GPU FP32. Strict configuration parsing/roundtrip, dimensional scaling and annular engineering integrals pass. Both runs pass independent standard-library JSON/CSV/VTK structural and physical-field checks. GPU example: mass flow 0.09681747 kg/s, exit Mach 2.30520, estimated thrust 172.67063 N. These are iteration-limited transient runs, not convergence claims.

## Execution plan

Follow M0–M9 in order; commit only at passed milestone gates, use targeted development tests and one comprehensive final pass. Keep physics/numerics independent of backend and UI. Measure actual errors and performance.

## Deferred

Only specified Phase 2 features: plume, chemistry, turbulence, advanced thermodynamics, CUDA/OpenGL interop and advanced nozzle optimisation.
