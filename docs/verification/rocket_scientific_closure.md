# Rocket scientific closure — V1.1

## Configuration

The dedicated `examples/rocket_nozzle/scientific_closure.json` retains the Phase 1 physical case: chamber/throat/exit radii 0.006/0.0045/0.00675 m, chamber/contraction/expansion lengths 0.03/0.06/0.12 m, P0=2 MPa, T0=2800 K, ambient pressure 101325 Pa, gamma=1.4, R=287 J/(kg K), constant viscosity 4e-5 Pa s and Pr=0.72. Walls are stationary, no-slip and adiabatic. MC/MUSCL, HLLC/HLLE, SSP-RK2 and CFL=0.4 are unchanged. Scientific runs use FP64 throughout; no solver-kernel mathematics was changed.

The first grid is 128x32. An initial safety allowance of 150,000 iterations permits observation of convergence; it is not a claimed convergence time. Each run stops automatically only when the criteria below pass or an explicit limit/failure occurs.

## Convergence definition

The monitor samples the first completed step and then every 20 iterations. Initial residuals mean the first completed step's nondimensional RMS time derivatives, not the unadvanced zero diagnostics. After at least 5,000 iterations, the final 2,000-iteration window must simultaneously satisfy:

- Every sampled mass, axial-momentum, radial-momentum and energy residual <= 1e-6.
- Relative range `(max-min)/abs(mean)` <= 5e-4 for mass flow, thrust, Isp, throat Mach, exit Mach and exit pressure.
- Inlet/exit relative mass-flow mismatch <= 1e-3 and station relative mass-flow spread <= 2e-3 at every sampled point.

The population relative standard deviation and least-squares relative drift (slope times window span divided by absolute mean) are also recorded. Zero means use a 1e-30 denominator floor; undefined Isp cannot pass. An optional additional relative residual target is supported; zero disables that extra condition, never the absolute target. Reduction factors are initial/final (larger means more reduction); undefined zero-denominator factors are null.

Summary fields include `termination_reason`, `converged`, `initial_residuals`, `final_residuals`, `residual_reduction_factors`, `convergence_window` and `observable_stability`. Existing `termination` and `residuals_nondimensional` keys remain. A legacy residual-only stop is no longer labelled steady convergence. SIGINT/SIGTERM request `user_stop`; numerical-step exceptions produce `invalid_state` metadata without exporting an invalid final field. CLI configuration/startup failures remain errors, not completed simulations.

## Reproduction and run preservation

```sh
.venv/bin/python scripts/run_scientific_closure.py --grids 128x32
# Only after the first grid is independently verified:
.venv/bin/python scripts/run_scientific_closure.py --grids 128x32 256x64 512x128
.venv/bin/python scripts/analyze_scientific_closure.py runs/scientific-closure
.venv/bin/python -m unittest discover -s tests/analysis -v
```

The runner invokes the CLI sequentially and keeps each attempt under its grid directory. Completed runs are reused only after independent convergence checks and matching configuration, executable SHA256 and output checksums. Interrupted, failed or incompatible attempts are preserved. An exclusive lock prevents overlapping orchestrators. Resume restarts an incomplete attempt from its original initial condition in a new directory; it does not claim checkpoint restart. No unconverged case may advance to the next grid.

## Analysis definitions

Consecutive-grid relative differences use the finer value as denominator. Goals are 0.5% for mass flow, thrust, Isp, throat Mach and exit Mach, and 1% for exit pressure. Profiles use the first radial cell row as the centreline approximation and last row as the wall-pressure approximation. Axial profiles are linearly interpolated onto 2,049 common x/L points within their shared cell-centre extent; no endpoint extrapolation hides differences near the throat.

For r=2 and monotonic, decreasing, non-roundoff differences, observed order is `log(abs((coarse-medium)/(medium-fine)))/log(2)`. Richardson extrapolation is `fine+(fine-medium)/(2^p-1)` and fine-grid GCI is `100*1.25*abs((fine-medium)/fine)/(2^p-1)` percent. Oscillatory/non-monotonic, non-decreasing or roundoff-level sequences are rejected. These estimates assume a leading single power of h; three samples alone do not independently prove the asymptotic regime. Synthetic first-, second- and third-order sequences with both signs test the formulas.

## Current results

Convergence infrastructure passed all six CPU and eight CUDA Release CTest entries. New rolling-window, drift/range, conservation-limit and termination tests passed. The unchanged 3,000-step demonstration remains `iteration_limit`, `converged=false`, with identical final FP32 residuals to Phase 1. Five independent Python analysis tests pass. Long-run closure and grid acceptance have not yet been established.

## Limitations

Grid independence and convergence of this ideal-gas axisymmetric nozzle calculation do not constitute experimental validation of a real rocket engine. This study does not add physics or establish accuracy for arbitrary geometries/backpressures. Wall shear is not exported because V1 has no independently verified wall-shear output subsystem.
