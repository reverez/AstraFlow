# Governing equations

The implemented 1D reference solves dU/dt + dF/dx = 0, with U = (rho, rho u, rho v, rho E). F = (rho u, rho u²+p, rho u v, u(rho E+p)). Transverse momentum is retained for reuse in multidimensional face rotations.

The calorically perfect ideal gas uses p=(gamma-1)[rho E-rho(u²+v²)/2], T=p/(rho R), a=sqrt(gamma p/rho). Gas properties are configurable. Sod uses nondimensional reference density, pressure and length of one; its R=1. Floats and doubles share templated mathematical functions.

Axisymmetric and viscous extensions are subsequent milestones; this document will be extended when implemented.
