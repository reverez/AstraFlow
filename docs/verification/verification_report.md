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

CPU CTest: 3 entries passed; 296 unit and 407 verification assertions. CUDA, nozzle, viscous and grid-refinement verification are not yet claimed.

## M2 CUDA Euler (Release)

Both precisions pass the same exact-Sod bounds and uniform-state checks. At 400 cells and t=0.2:

| Precision | CPU/GPU max absolute | L1 | L2 |
|---|---:|---:|---:|
| FP32 | 4.61936e-6 | 2.49782e-7 | 5.32646e-7 |
| FP64 | 4.21885e-14 | 7.72314e-16 | 3.19405e-15 |

Norms cover all four conservative fields, nondimensional. The initial FP32 maximum 6.10352e-5 exceeded the preset 5e-5 bound. A controlled build with CPU `-ffp-contract=off` and CUDA `--fmad=false` reduced that discrepancy; this matched-rounding baseline is retained. No tolerance was loosened. GPU checks contain 4,008 passing assertions. This does not yet verify a multidimensional solver.

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

The 4x24 grid starts from u=0.8y, T=1 and evolves to t=15 (9,924 steps). L1 velocity error is 5.14813e-5 (bound 0.002); L1 temperature error is 4.46196e-5 (bound 0.001). A 48x12 adiabatic no-slip nozzle remains physical for 1,000 steps to t=3.553 with no reconstruction fallback. On a 32x8 viscous nozzle after 100 steps, maximum conservative CPU/GPU differences are 3.9968e-15 (FP64) and 4.26173e-6 (FP32), within unchanged 1e-10/5e-5 bounds. All seven CTest entries pass.
