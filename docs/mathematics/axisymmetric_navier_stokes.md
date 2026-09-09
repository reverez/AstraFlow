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

The viscous extension will add the azimuthal normal-stress source at M5.
