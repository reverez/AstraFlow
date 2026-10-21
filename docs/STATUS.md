# AstraFlow status — V2.0 (2026-10-21)

Supersedes [PHASE1_STATUS.md](PHASE1_STATUS.md) as the live status record; that file is kept as
the V1.0/V1.1 acceptance history.

## Release state

| Release | State | Evidence |
|---|---|---|
| V1.0 (`v1.0.0`) | Accepted | [verification_report.md](verification/verification_report.md), [benchmark_report.md](performance/benchmark_report.md) |
| V1.1 | Closed as a diagnostic record; MC medium grid never converged | [rocket_scientific_closure.md](verification/rocket_scientific_closure.md), C2b, C2c |
| V1.2 | Partially accepted per its protocol: coarse closure, CPU parity, pseudo-time root study pass; medium cap under-budgeted | [rocket_closure_v12.md](verification/rocket_closure_v12.md) |
| V1.3 | Partially accepted: medium grid closed at 760,760 iterations; fine grid reached its 2.5M cap still decaying | same report |
| V1.4 | Pre-registered fine-grid continuation (4.5M budget); see the report for the outcome | [v14_fine_continuation_protocol.md](verification/v14_fine_continuation_protocol.md) |
| V2.0 | Complete | this file, `CHANGELOG` section below |

Release suites: CPU 8/8 CTest entries (37 Catch2 cases, 7,501 assertions), CUDA 10/10, 11 Python
analysis tests under the project venv. The V2 executable reproduces the V1.1 executable's
2,000-iteration FP64 rocket field bitwise.

## V2.0 changes

Optimisation (fields unchanged):
- One device-to-host transfer per CUDA step; the RK kernels read the reduced timestep on the
  device. 128×32: 0.36 → 0.22 ms/iteration; 512×128: 0.55 → 0.34 ms; 1024×256 unchanged
  (bandwidth-bound). Old and new executables benchmarked alternately in the same machine state.
- Solver-side numerical inlet/outlet mass fluxes (both backends, same sampling stage as the
  residual). The FP64 scientific configuration runs 2.3× faster end to end because the host
  finite-volume balance is no longer evaluated every 20 iterations.
- OpenMP CPU reference backend, opt-in (`--threads`), serial fixed-order reductions so results are
  thread-count independent; 3.6–3.9× at 1024×256 with 8 bound passive threads. libgomp's default
  spin-wait is pathological on this hybrid-core WSL2 host and is documented.

Numerics and physics:
- Van Albada limiter as a production option (Sod, entropy-wave order, parity gates pass).
- Sutherland and power-law viscosity with k = μ cp / Pr, verified on compressible Couette flow
  with the exact solution for arbitrary μ(T).
- Wall traction and conductive heat-flux analysis from the production viscous flux (`wall.csv`,
  integrated loads); verified on Couette (shear, dissipation balance) and on the four-equation
  energy balance of an isothermal nozzle wall.
- Convergence monitor gate on the numerical boundary mass flux (`mass_gate`), `balance_interval`.
- Isothermal walls colder than half the stagnation temperature are rejected at configuration
  (the mirrored ghost temperature would be negative). The cooled-wall example uses 1500 K.

Tooling:
- `--version`, `--threads`, `--write-state`; schema-1 state export/import in the library.
- Study drivers: `run_pseudo_root_study.py`, `run_v13_continuation.py`,
  `run_v14_fine_continuation.py`, `sweep_nozzle.py`, `plot_readme_figures.py`, `make_banner.py`.
- CI: clang-format check, analysis tests under a venv, end-to-end Sod run validation.

## Scientific results of this release

- 128×32 van Albada: `steady_converged` at 49,480 iterations (MC needed 82,480); CPU and CUDA FP64
  agree to 1.3e-14 at the same iteration.
- Pseudo-time root hypothesis supported: physical/pseudo endpoint difference 1.6e-6 → 5.1e-7 →
  1.6e-7 for tolerances 1e-6 → 3e-7 → 1e-7 (ratio 10.2).
- 256×64 van Albada: `steady_converged` at 760,760 total iterations, residuals < 1e-6, numerical
  mass mismatch 1.1e-9, zero fallbacks. First closure of the medium rocket grid in the project.
- Coarse-to-medium changes: mass flow 0.50%, thrust 0.44%, Isp 0.06%, exit pressure 0.04%, exit
  Mach 0.8%.
- 512×128: 2.5M iterations, energy residual 2.2e-5 and decaying exponentially (e-folding
  1.1–1.3M iterations); not accepted; continuation pre-registered as V1.4.

## Known limitations and open items

- The inviscid slip-wall nozzle plateaus near a 1e-5 residual under local pseudo-time; not
  investigated. The expansion sweep therefore uses the viscous no-slip configuration.
- Isothermal walls need 2 T_wall ≥ T0 with the mirrored ghost formulation; a positivity-preserving
  wall ghost with a diagnostic counter would remove that restriction.
- The "throat Mach" observable samples the minimum-area station whose cell-centre position moves
  under refinement; it is reported but not used as a convergence measure.
- Fine-grid physical-time convergence is slow (decay length grew ~8× from medium to fine, against
  the 2× budgeted); pseudo-time or an implicit steady solver would be the next acceleration step.
- Wall-time comparisons between runs of the campaign are affected by laptop load variance and are
  reported as such; iteration counts are the reliable measure.
- Not attempted in V2 because they cannot be verified to the project's standard here: external
  plume, reacting chemistry, turbulence modelling, real-gas thermodynamics, multi-GPU, CUDA/OpenGL
  interoperability.
