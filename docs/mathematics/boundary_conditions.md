# Boundary conditions

Boundary functions are shared by CPU and CUDA. Every RK stage updates boundary states from that stage's interior primitive state.

- Axis: rho,u,p even; v odd. The axis face has zero cylindrical area.
- Inviscid wall: reflect the face-normal velocity and retain tangential velocity, rho and p; the symmetric Riemann problem enforces impermeability.
- Reservoir inlet: prescribed stagnation P0,T0 and axial inflow. Retain the outgoing interior acoustic invariant J-=u-2a/(gamma-1). Solve J-=a(M)[M-2/(gamma-1)] on 0<=M<1, with a(M)²=gamma R T0/[1+(gamma-1)M²/2]. Reconstruct static temperature, pressure and density from isentropic relations. Bisection uses 52 iterations in FP64 and 28 in FP32. This model is intended for subsonic reservoir-fed inflow; it does not model arbitrary inlet backflow or supersonic injection.
- Outlet: for axial Mach >=1, extrapolate the outgoing primitive state. For subsonic flow prescribe static back pressure, retain entropy p/rho^gamma and J+=u+2a/(gamma-1), and derive density and axial velocity. Radial velocity is extrapolated.
- Periodic Cartesian verification meshes connect opposite cells; paired boundary-face fluxes use the same states, preserving conservation.

The 1D Euler reference uses transmissive boundaries with zero end-cell slopes. Viscous nozzle walls reflect both velocity components about the stationary wall velocity (no slip); density/pressure reflect evenly. The heat gradient is projected tangentially at adiabatic walls so its normal flux is exactly zero, including sloping faces. Verification additionally supports moving Cartesian upper walls and prescribed wall temperature: ghost velocity and temperature reflect about their specified face values. Production nozzle examples use stationary, adiabatic walls.
