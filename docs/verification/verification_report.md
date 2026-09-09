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
