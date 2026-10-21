# V1.2 rocket closure protocol (pre-registered)

Written and committed before any V1.2 scientific run was launched. Nothing in this section is
edited after the results are known; measured outcomes are reported separately in
[rocket_closure_v12.md](rocket_closure_v12.md).

## Why a new protocol

V1.1 left the rocket case scientifically unresolved. The evidence collected there is accepted as
the starting point:

- [C2b](residual_closure_c2b.md) showed that the 256×64 MC residual stall is not a conservation or
  assembly defect. The true inlet/outlet numerical mass-flux mismatch is ~1e-7, every four-equation
  balance assembles to roundoff, and the dominant residual cells cycle in sign hundreds of times in
  step with MC limiter branch switches at near-axis cells.
- The V1 interior-cell mass estimate (first/last interior cell values integrated over inlet/exit
  areas) is a first-order quadrature of a different quantity from the reconstructed boundary flux.
  Its 1e-3 gate failed at 1.18e-3 on 256×64 while the numerical boundary flux balanced to 1.3e-7.
- [C2c](steady_acceleration_c2c.md) rejected local pseudo-time because two solutions stopped at a
  1e-6 residual differed by 2.8e-6 in Linf against a preset 1e-8 field bound. Both endpoints carry
  finite residuals, so that bound tested something the stopping rule cannot deliver; the report
  says so itself.

The V1.1 gates are not modified retroactively and no V1.1 result is re-labelled. V1.2 is a new
study with its own criteria, justified by that evidence, and the criteria below are fixed now.

## Hypothesis

H1: the medium-grid stall is a limiter-induced limit cycle of the non-smooth MC slope limiter.
Replacing MC by the smooth, symmetric van Albada limiter, with every other ingredient unchanged
(HLLC/HLLE, corrected viscous gradients, cylindrical source, boundary states, SSP-RK2, CFL 0.4,
FP64), lets 256×64 and 512×128 reach the unchanged 1e-6 residual target.

H2: local pseudo-time and physical-time integration converge to the same discrete root. Evidence
for H2 is that the conservative-field difference between the two endpoints shrinks as the shared
residual tolerance is tightened; a fixed field bound at a fixed tolerance is not used again.

## Production scheme under study

`examples/rocket_nozzle/v12_closure.json`. Physical case unchanged from V1.1: chamber/throat/exit
radii 0.006/0.0045/0.00675 m, lengths 0.03/0.06/0.12 m, P0 = 2 MPa, T0 = 2800 K, ambient
101325 Pa, γ = 1.4, R = 287 J/(kg K), constant μ = 4e-5 Pa s, Pr = 0.72, stationary adiabatic
no-slip wall. Numerics: MUSCL with the **van Albada** limiter
σ = ab(a+b)/(a²+b²) for ab > 0, else 0; HLLC with HLLE fallback; SSP-RK2; CFL 0.4; CUDA FP64.

Van Albada is accepted as a production limiter only because it passes the same software gates as
MC before this study starts: Sod L1 bounds (0.006/0.008/0.005), entropy-wave observed order > 1.5,
limiter unit properties (sign change gives zero, symmetric, equal slopes preserved), and CPU/CUDA
parity bounds. These results are in the verification report, not here.

## Convergence definition

Unchanged from V1.1 except for the mass gate:

- Sampling: first completed step and every 20 iterations; at least 5,000 iterations; final
  2,000-iteration window of 101 samples.
- Every sampled mass, axial-momentum, radial-momentum and energy residual ≤ 1e-6 (nondimensional
  RMS time derivatives at the second SSP-RK2 stage, as before).
- Relative range (max−min)/|mean| ≤ 5e-4 over the window for mass flow, thrust, Isp, throat Mach,
  exit Mach and exit pressure.
- **Mass gate: numerical boundary mass-flux mismatch** |ṁ_out − ṁ_in| / max(|ṁ_in|, |ṁ_out|)
  ≤ **1e-6** at every sampled point, where ṁ_in and ṁ_out are sums of the actual reconstructed
  inlet/outlet HLLC mass fluxes evaluated at the same stage as the residual statistic. The interior
  estimate is still recorded and reported, but it is no longer a gate.
- Station relative mass-flow spread ≤ 2e-3 at every sampled point (unchanged).

The gate tolerance 1e-6 is chosen from the V1.1 accepted coarse case, whose final numerical
mismatch was 9.5e-10, and from the rejected medium case, whose value was 1.3e-7 while its
residual was 1e-4. A converged state must therefore sit three decades below the failed one.

## Grid sequence and caps

| Grid | Iteration cap | Reason |
|---|---:|---|
| 128×32 | 150,000 | V1.1 cap; MC converged at 82,480 |
| 256×64 | 400,000 | V1.1 used 250,000 and stalled; the cap allows for the halved timestep plus the slower late decay observed on the coarse grid |
| 512×128 | 1,000,000 | Four times the coarse timestep count plus margin |

Grids run sequentially with `scripts/run_scientific_closure.py`; a grid that is not
`steady_converged` stops the sequence and no finer grid is launched. Caps are not raised after a
failure. Each grid starts from the quasi-1D initial state, not from a prolonged coarser solution,
so that every grid is an independent integration.

## Acceptance

V1.2 is **accepted** only if all of the following hold:

1. All three grids terminate `steady_converged` under the definition above with zero HLLC and
   reconstruction fallbacks.
2. An independent CPU FP64 integration of 128×32 reaches the same termination iteration and the
   maximum nondimensional conservative-field difference from the CUDA run is ≤ 1e-10.
3. Richardson analysis on the three accepted grids is valid (monotone, decreasing differences) for
   at least mass flow, thrust and exit Mach, and the fine-grid GCI for each is reported. No GCI
   threshold is imposed; the uncertainty is the result.
4. Consecutive-grid relative differences meet the V1.1 goals: ≤ 0.5% for mass flow, thrust, Isp,
   throat Mach, exit Mach and ≤ 1% for exit pressure, on the medium/fine pair.

Partial outcomes are reported as partial. If 256×64 stalls under van Albada, H1 is rejected and
the study records that; it does not switch limiters again inside the same study.

## Pseudo-time root study (secondary, 128×32 only)

For residual tolerances 1e-6, 3e-7 and 1e-7, run physical-time SSP-RK2 and local pseudo-time RK2
(steady CFL 0.4) with identical configuration otherwise, each capped at 400,000 iterations. For
each tolerance record the nondimensional conservative-field Linf difference between the two
endpoints, their wall times and iteration counts.

Pre-registered reading of H2: supported if the Linf difference decreases monotonically with the
tolerance and the ratio Linf(1e-6)/Linf(1e-7) is ≥ 3; refuted if the difference does not decrease.
If any run hits its cap the study is inconclusive at that tolerance and says so. No field bound at
a single tolerance is used. The pseudo-time mode remains an acceleration option, not the accepted
study integrator, regardless of outcome.

## Preservation

Every run keeps its config, engineering and residual histories, final VTK field, final true
balance, provenance (binary SHA256, source commit) and summary. Nothing is overwritten;
`run_scientific_closure.py` reuses a completed attempt only after independent validation and
checksum match. The analysis step reads preserved outputs only.
