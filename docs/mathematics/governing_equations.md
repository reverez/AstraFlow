# Governing equations

The implemented 1D reference solves dU/dt + dF/dx = 0, with U = (rho, rho u, rho v, rho E). F = (rho u, rho u²+p, rho u v, u(rho E+p)). Transverse momentum is retained for reuse in multidimensional face rotations.

The calorically perfect ideal gas uses p=(gamma-1)[rho E-rho(u²+v²)/2], T=p/(rho R), a=sqrt(gamma p/rho). Gas properties are configurable. Sod uses nondimensional reference density, pressure and length of one; its R=1. Floats and doubles share templated mathematical functions.

The production solver extends these equations to weighted axisymmetric finite volumes, including Newtonian stress and Fourier conduction. See the separate axisymmetric formulation and numerical method documents.

## Reference scales and engineering output

Dimensional nozzle cases use Lref=throat radius, Tref=T0, pref=P0, rhoref=P0/(R T0), vref=sqrt(R T0), tref=Lref/vref. Cartesian cases use domain height as Lref. The scaled gas has R*=1, P0*=T0*=1 and mu*=mu/(rhoref vref Lref). Primitive floors scale consistently. State components scale by (rhoref,rhoref vref,rhoref vref,pref). Mesh measures scale by Lref, Lref² or Lref³ as appropriate; the radial source coefficient scales inversely with length. Export reverses these transformations; effective scales are saved in performance.json.

Engineering exit integration uses actual exit annular-face areas times 2 pi, with the last interior state approximating the exit face state. Mass flow is sum rho u A; thrust is sum [rho u²+(p-pa)]A; Isp=F/(mdot g0), g0=9.80665 m/s². Isp is null at nonpositive/vanishing mass flow. Pressure, temperature and Mach are area averaged; axial velocity is mass-flow weighted. Throat axial Mach uses the minimum-area cell column. Station mass-flow diagnostics use the cell-averaged section area V/dx. These quadratures approximate boundary/section values at finite resolution; they are not experimental engine predictions. Cartesian integrals are per unit depth.

VTK total_energy is specific total energy E; velocities and geometry are physical, and the domain is the meridional plane at z=0. Display/rendering never changes the exported field data.
