# Axisymmetric finite volumes

The domain is the meridional (x,r) plane without swirl. The stored conservative state is (rho, rho u, rho v, rho E). Integrals use volume V=integral r dx dr and face area A=integral r ds, per radian. Physical volume/area for engineering integrals multiplies by 2 pi.

The inviscid update is dU/dt = -sum_faces(F_n A)/V + (0,0,p B/V,0), where B is the unweighted meridional cell area. This is the cylindrical radial-momentum pressure source. It is required even though radial weighting is already present in the face areas.

Vertices follow r(x,eta)=eta R(x). R is constant in the chamber and piecewise cubic Hermite (smoothstep) in contraction/expansion, with zero derivative at joins. Faces connect sampled contour vertices by straight segments. For a cell spanning x0 to x1 and eta0 to eta1, with R0,R1 at its axial ends:

- V = (eta1²-eta0²) dx (R0²+R0 R1+R1²)/6.
- B = dx (eta1-eta0)(R0+R1)/2.
- Axial-face A=(r_upper²-r_lower²)/2.
- Sloping radial-face A=(r0+r1) length/2, with the exact segment normal.

These formulas integrate the piecewise-linear mesh exactly. They obey sum A n_x=0 and sum A n_r=B, so a constant pressure cancels its geometric source discretely. The axis face has zero area; the first cell centre is at positive r. No source evaluates 1/r at the axis itself. The centreline reflects radial velocity and retains density, pressure and axial velocity.

Logical-coordinate MUSCL slopes reconstruct to the corresponding face; cell-centre radial positions are arithmetic interpolation points. A 50-step stationary-pressure test checks the full flux/source cancellation, not just the mesh identities.

With viscosity, use div(u)=u_x+v_r+v/r; tau_xx=2 mu u_x-(2/3)mu div(u), tau_rr=2 mu v_r-(2/3)mu div(u), tau_xr=mu(u_r+v_x), and tau_theta_theta=2 mu v/r-(2/3)mu div(u). The radial source becomes (p-tau_theta_theta) B/V. The conservative face flux is convective flux minus viscous flux. Its momentum components are tau dot n; its energy component is velocity dot (tau dot n) + k grad(T) dot n.

At r=0 the face evaluation uses the regular odd-velocity limit v/r -> v_r, while the zero face area removes its integrated flux. Cell sources evaluate v/r only at positive cell centres. Cartesian tests explicitly omit cylindrical hoop terms.
