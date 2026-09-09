# Verification report

## M1 CPU Euler (Release, FP64)

Sod: unit interval, discontinuity at x=0.5, gamma=1.4, (rho,u,p)L=(1,0,1), (rho,u,p)R=(0.125,0,0.1), t=0.2, 400 cells, MC reconstruction, CFL=0.4. Errors compare cell-centre numerical primitives to the independent exact self-similar solution.

| Quantity | Mean absolute error | Acceptance bound |
|---|---:|---:|
| Density | 0.00144607 | 0.006 |
| Axial velocity | 0.00254609 | 0.008 |
| Pressure | 0.000890535 | 0.005 |
| Total mass drift | 3.33067e-15 | 1e-12 |

Bounds allow expected shock/contact discretisation error at this resolution while remaining below one percent of reference scales. No reconstruction fallback was required. Unit tests cover equal-state flux, face rotation, stationary contact, expansion fallback, limiters, thermodynamics, invalid states, and uniform preservation.


## M2 CUDA Euler (Release)

Both precisions pass the same exact-Sod bounds and uniform-state checks. At 400 cells and t=0.2:

| Precision | CPU/GPU max absolute | L1 | L2 |
|---|---:|---:|---:|
| FP32 | 4.61936e-6 | 2.49782e-7 | 5.32646e-7 |
| FP64 | 4.21885e-14 | 7.72314e-16 | 3.19405e-15 |

Norms cover all four conservative fields, nondimensional. The initial FP32 maximum 6.10352e-5 exceeded the preset 5e-5 bound. A controlled build with CPU `-ffp-contract=off` and CUDA `--fmad=false` reduced that discrepancy; this matched-rounding baseline is retained. No tolerance was loosened. GPU checks contain 4,008 passing assertions. Multidimensional verification follows below.

## M3 2D finite volumes

Rectangular cell volumes sum to the domain volume; outward face vectors close within 1e-14. Invalid dimensions/volumes are rejected. CPU and GPU uniform periodic flow (16x8, 20 steps) preserves every conservative component within 1e-13. A nonuniform periodic density wave (32x16, 30 steps) preserves all four integrated conserved quantities within 1e-12 and passes the preset FP64 parity bound 2e-11. Both complete Release CTest configurations pass, including prior 1D checks.

## M4 axisymmetric nozzle

The default smooth contour has chamber/throat/exit radii 0.2/0.15/0.225 and chamber/contraction/expansion lengths 1/2/4, with nondimensional P0=T0=R=1, gamma=1.4. A quasi-1D initial guess includes streamline radial velocity; it is then evolved by the full 2D flux/source solver for 4,000 iterations on 96x12 cells, reaching t=13.7749. It is not held fixed or overwritten by the analytical solution.

| Check | Measured | Preset bound |
|---|---:|---:|
| Near-throat area-averaged axial Mach | 1.00957 | 0.85–1.15 |
| Axial mean absolute Mach difference from quasi-1D | 0.000319245 | 0.12 |
| Station mass-flow spread / maximum | 0.000774078 | 0.05 |
| Stationary pressure state max drift, 50 steps | 2.22045e-15 | 1e-11 |

The nozzle's small wall slopes make quasi-1D theory a useful approximation. These bounds accommodate legitimate multidimensional effects; they do not assert exact quasi-1D behavior for steep contours or shocked flow. Geometry identities, C1 joins, rejected invalid radii, inlet stagnation enthalpy and supersonic outlet extrapolation also pass. The initial stationary-state failure was traced to the 28-step FP64 inlet root solve; increasing the root solve to FP64 precision repaired it without changing the test threshold. GPU parity is checked after 100 nozzle steps on 32x8.

## M5 compressible Navier–Stokes

The canonical case is exact steady compressible Couette flow with constant viscosity, periodic x, stationary lower wall, upper-wall U=1, and equal isothermal wall temperatures Tw=1. With H=1, R=1, gamma=1.4, mu=0.05, Pr=0.72, the analytical solution is u=y, v=0, T=1+mu y(1-y)/(2k), constant pressure and rho=p/(R T). Heat conduction balances viscous heating. This tests the full compressible energy equation, not an incompressible substitute.

The 4x24 grid starts from u=0.8y, T=1 and evolves to t=15 (9,924 steps). L1 velocity error is 5.14813e-5 (bound 0.002); L1 temperature error is 4.46196e-5 (bound 0.001). A 48x12 adiabatic no-slip nozzle remains physical for 1,000 steps to t=3.553 with no reconstruction fallback. On a 32x8 viscous nozzle after 100 steps, maximum conservative CPU/GPU differences are 3.9968e-15 (FP64) and 4.26173e-6 (FP32), within unchanged 1e-10/5e-5 bounds.

## M6 production workflow

Strict JSON parsing and roundtrip, dimensional scaling/restoration and uniform analytical exit-area/mass-flow/thrust/Isp integration pass unit tests. Actual 128x32 viscous rocket runs completed 300 steps on CPU FP64 and CUDA FP32. The CUDA run produced mdot=0.0968174695 kg/s, exit Mach=2.30519878, exit p=154822.8336 Pa and idealised thrust=172.6706309 N. Its mass-flow inlet/exit mismatch was 0.00312177. These are transient iteration-limit outputs, not steady convergence or experimental validation.

`scripts/check_run.py` independently parsed both runs' VTK XML and JSON/CSV, checked grid/array lengths, finite positive fields, last iteration consistency, and the thrust/Isp identity. All checks passed. No Python packages were installed.

## M7 GUI verification

Actual WSLg startup and OpenGL rendering passed. `astraflow_gui --smoke-test` drives the same control handlers used by the buttons, verifies exactly one requested step and reset-to-zero behavior, regenerates the mesh, cycles all nine fields, and exercises residual/engineering updates. This is scripted handler testing, not a claim of manual mouse testing. Independent worker tests cover pause boundaries, immutable old snapshots, gas-property retention and rejected invalid geometry.

The captured framebuffer exposed a convergence-axis bug (x stayed at 0–100 despite later iterations); it was repaired and the GUI smoke/capture rerun. The verified screenshot at `docs/screenshots/astraflow.png` shows the complete current history and engineering dashboard. No display limitation was encountered. Both relocated Release build suites also passed.

## Final grid-refinement studies

CPU FP64, MC reconstruction, SSP-RK2. Bounds were set before running these studies. The smooth periodic entropy wave has rho=1+0.2 sin(2 pi (x+y-0.6 t)), u=0.4, v=0.2 and p=1, evolved to t=0.25. Errors compare density at cell centres with the advected analytical wave.

| Grid | L1 | L2 | Linf | Integrated mass drift |
|---|---:|---:|---:|---:|
| 20x10 | 0.00424178 | 0.00494423 | 0.0101353 | 4.44e-16 |
| 40x20 | 0.00120883 | 0.00151437 | 0.00367988 | 1.55e-15 |
| 80x40 | 0.000301796 | 0.000452457 | 0.00139742 | 2.78e-15 |

L1 observed orders are 1.81106 and 2.00196, exceeding the preset 1.5 bound; finest L1 is below 0.002. Limiting at extrema affects the coarser-grid rate.

The Couette configuration above is evolved to t=20 with 12, 24 and 48 radial cells. Axial variation is absent, so four periodic axial cells suffice. The nonzero steady temperature error measures diffusion/energy discretisation; velocity has a small remaining transient.

| Grid | Temperature L1 | Temperature L2 | Temperature Linf | Velocity L1 |
|---|---:|---:|---:|---:|
| 4x12 | 0.000178571 | 0.000178571 | 0.000178821 | 4.64866e-6 |
| 4x24 | 0.0000446426 | 0.0000446429 | 0.0000448759 | 4.42355e-6 |
| 4x48 | 0.0000111605 | 0.0000111617 | 0.0000113896 | 4.36888e-6 |

Temperature L1 orders are 2.00001 and 2.00002, exceeding the preset 1.7 bound; finest error is below 3e-5. This is analytical verification of the canonical cases, not a grid-independence claim for the rocket example.

## Additional final verification

The 64x8 inviscid nozzle starts with all initial velocities reduced by 20% and evolves to t=20. Throat Mach recovers from 0.802966 to 0.996907, with station mass-flow spread 0.00157557. Preset checks require Mach to rise by at least 0.1, finish within 0.15 of unity and mass spread below 0.05. This verifies recovery rather than merely preserving a quasi-1D initial guess.

Long-time CUDA Couette runs use the same 4x24 analytical case, end time 15 and bounds as CPU, each taking 9,924 steps:

| Backend | Velocity L1 | Temperature L1 |
|---|---:|---:|
| CPU FP64 | 5.14813e-5 | 4.46196e-5 |
| CUDA FP64 | 5.14813e-5 | 4.46196e-5 |
| CUDA FP32 | 5.94457e-5 | 3.43532e-5 |

Final robustness checks reject nonfinite mesh centres/normals, invalid connectivity, overflowing cell counts, invalid reservoir/wall settings and unsupported moving nozzle walls. Chamber pressure uses annular-area weighting, covered by an independent integration test.

## Final production runs and acceptance pass

Release CPU and CUDA builds completed in the relocated repository on 2026-09-09. All six CPU CTest entries and all eight CUDA CTest entries passed. The final CUDA suite includes both long-time Couette precisions, all parity checks and the SM120 device probe.

Both production runs use the checked-in 128x32 viscous rocket example for 3,000 iterations. They reach the iteration limit, not the configured 1e-6 residual target. Reported values are transient ideal-gas estimates.

| Quantity | CPU FP64 | CUDA FP32 |
|---|---:|---:|
| Simulated time (s) | 0.000147368906 | 0.000147368944 |
| Mass flow (kg/s) | 0.0966167097 | 0.0966168406 |
| Throat Mach | 1.0069978 | 1.0069983 |
| Exit Mach | 2.22319146 | 2.22319197 |
| Exit pressure (Pa) | 162185.19 | 162185.369 |
| Estimated thrust (N) | 171.058507 | 171.05877 |
| Specific impulse (s) | 180.539302 | 180.539334 |
| Inlet/exit mass-flow relative mismatch | 0.000465017777 | 0.000465525019 |
| Station mass-flow relative spread | 0.00170668327 | 0.00170655231 |

Both runs have zero flux/reconstruction fallbacks. `scripts/check_run.py` passed independently for both directories, including VTK structure, finite positive state fields, configuration/grid consistency, last CSV iteration and the thrust/Isp identity. Generated data remain under ignored `runs/final-cpu` and `runs/final-cuda`.

The final actual WSLg GUI smoke/capture passed all nine fields, Run/Pause/Step/Reset/Regenerate, live residuals and engineering/performance updates. The captured 1480x980 image was inspected. Smoke mode raises the runtime limit and disables automatic convergence so a fast GPU cannot finish before scripted interaction checks; normal user configurations retain their limits. This tests handlers and rendering, not physical mouse events.

All six CPU CTest entries passed with AddressSanitizer, UndefinedBehaviorSanitizer and leak detection enabled (102.87 seconds). No sanitizer findings were reported. The initial sandboxed attempt passed unit assertions but LeakSanitizer reported its incompatibility with ptrace. That attempt was stopped; the complete successful pass used the same instrumented binaries outside the tracing sandbox, retaining leak detection.

To reproduce the CPU memory/runtime check (with leak detection outside a tracing sandbox):

```sh
cmake -S . -B build-sanitize -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DASTRAFLOW_ENABLE_CUDA=OFF -DASTRAFLOW_BUILD_GUI=OFF \
  '-DCMAKE_CXX_FLAGS=-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all' \
  '-DCMAKE_CXX_FLAGS_RELWITHDEBINFO=-O1 -g -DNDEBUG' \
  '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined'
cmake --build build-sanitize -j 4
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
  ctest --test-dir build-sanitize --output-on-failure
```

The instrumented suite is substantially slower than Release because it checks the full refinement integrations. Nsight and compute-sanitizer were not installed; no GPU memory-sanitizer pass is claimed. CPU CI is configured; the local validation reported here does not claim a completed hosted GitHub Actions run.
