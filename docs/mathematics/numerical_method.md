# Numerical method

Cell-centred conservative finite volumes store four separate contiguous state planes (SoA). Primitive MUSCL reconstruction uses selectable minmod, van Leer or MC slopes; MC is the tested baseline. End-cell slopes vanish at transmissive boundaries.

Face states are rotated into normal/tangential coordinates for HLLC. Davis minimum/maximum acoustic wave bounds bracket the contact estimate. Nonfinite or nonphysical intermediate states trigger conservative HLLE using the same wave bounds. Invalid reconstructed states revert that face to first order and increment a diagnostic; the conservative solution itself is never clipped.

SSP-RK2 computes U1=U+dt L(U), then Unew=(U+U1+dt L(U1))/2. The 1D timestep is CFL*dx/max(|u|+a), capped by requested end time; CFL must be <=0.5. Every stage is checked for finite, positive density, pressure and temperature. Invalid stages throw with cell and iteration context. Fluxes are stored once per face, ensuring cancellation between neighbouring cells.

The exact Sod reference independently solves the pressure wave-curve equation by 100 bisections and samples analytical shock, contact and rarefaction regions. The reference is test-only. At shocks, second-order global convergence is not expected.
