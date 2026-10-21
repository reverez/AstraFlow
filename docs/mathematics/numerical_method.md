# Numerical method

Cell-centred conservative finite volumes store four separate contiguous state planes (SoA). Primitive MUSCL reconstruction uses selectable minmod, van Leer or MC slopes; MC is the tested baseline. End-cell slopes vanish at transmissive boundaries.

Face states are rotated into normal/tangential coordinates for HLLC. Davis minimum/maximum acoustic wave bounds bracket the contact estimate. Nonfinite or nonphysical intermediate states trigger conservative HLLE using the same wave bounds. Invalid reconstructed states revert that face to first order and increment a diagnostic; the conservative solution itself is never clipped.

SSP-RK2 computes U1=U+dt L(U), then Unew=(U+U1+dt L(U1))/2. The 1D timestep is CFL*dx/max(|u|+a), capped by requested end time; CFL must be <=0.5. Every stage is checked for finite, positive density, pressure and temperature. Invalid stages throw with cell and iteration context. Fluxes are stored once per face, ensuring cancellation between neighbouring cells.

The exact Sod reference independently solves the pressure wave-curve equation by 100 bisections and samples analytical shock, contact and rarefaction regions. The reference is test-only. At shocks, second-order global convergence is not expected.

## Viscous discretisation

Constant dynamic viscosity is the guaranteed model. k=mu cp/Pr with cp=gamma R/(gamma-1). Central logical-coordinate gradients of velocity and temperature are mapped to physical coordinates: d/dr=(logical radial difference)/dr and d/dx=(logical axial difference)/dx-skew*d/dr. Face gradients average adjacent gradients, then correct along the centre-to-centre vector to reproduce the actual primitive jump. This correction controls nonorthogonal-grid diffusion and gives a consistent wall-normal derivative from mirrored ghost states.

The explicit limit includes D=max(4mu/(3rho),mu gamma/(Pr rho)), accounting for momentum diffusion and energy/temperature diffusion. A conservative diffusion limit dt_v=0.5 CFL/[D(1/dx²+(2+skew²)/dr²)] combines harmonically with the convective limit. The additional radial factor bounds cylindrical curvature and keeps the baseline conservative. No implicit viscous solve is used.

## V2 additions

### Van Albada limiter

The symmetric van Albada slope `sigma = a b (a + b) / (a^2 + b^2)` for `a b > 0` (zero otherwise)
is available as a production limiter alongside minmod, van Leer and MC. It equals `a` when the
one-sided differences agree, lies between minmod and MC in diffusivity, and is differentiable for
`a b > 0`, so it does not switch branches discretely the way MC does. The V1.2 rocket study uses
it because the C2b investigation traced the medium-grid stall of the MC scheme to branch
switching at near-axis cells. On the smooth entropy wave its observed L1 orders are 1.61 and
1.94 (MC: 1.81 and 2.00); on Sod at 400 cells its L1 errors are 1.83e-3 (rho), 3.33e-3 (u) and
1.22e-3 (p), inside the same preset bounds as MC.

### Temperature-dependent viscosity

Three transport models share one conductivity closure `k = mu(T) cp / Pr`:

- constant: `mu = mu_ref` (the V1 model; identical arithmetic);
- Sutherland: `mu = mu_ref (T / T_ref)^{3/2} (T_ref + S) / (T + S)`;
- power law: `mu = mu_ref (T / T_ref)^n`.

`mu` is evaluated at the face-average temperature of the two reconstructed transport states in
the viscous flux, at the cell temperature in the cylindrical hoop-stress source, and at the cell
temperature in the diffusive timestep bound. Reference and Sutherland temperatures are
nondimensionalised by `T0` with the rest of the thermodynamic state. Verification uses steady
compressible Couette flow, whose exact solution exists for any `mu(T)`: the energy equation gives
`T(u) = T_w + Pr u (U - u) / (2 cp)` independently of the viscosity law, and the momentum equation
`mu(T(u)) du/dy = tau` gives `y(u)` by quadrature.

### Wall quantities

At every wall face the production viscous flux function is evaluated on the final physical
snapshot. With `n_out` the outward normal of the fluid domain, the viscous traction on the fluid
is `tau . n_out`, and the conductive heat flux leaving the fluid into the wall is
`-k grad(T) . n_out`; for a moving Cartesian wall the work `tau . u_wall` carried by the viscous
energy flux is subtracted so that the reported quantity is conduction alone. The wall gradient
is the solver's own one-sided mirrored-ghost estimate, i.e. exactly the flux the discretisation
applies. Because the scheme conserves energy discretely, the steady Couette wall flux equals half
the discrete dissipation `mu U^2 / (2H)` to 2e-5 at 48 cells, independent of the truncation
error of the gradient itself. Integrated loads apply the full-revolution factor `2 pi` on
axisymmetric meshes. The integrated wall heat transfer rate is
algebraically the wall viscous/thermal energy flux of the four-equation balance, and the tests
check that identity to 1e-10.
