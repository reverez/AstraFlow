# Developer context

> **Purpose:** give a future contributor or coding agent enough context to continue AstraFlow safely without replaying the entire project history.

Read this file first, then follow links to evidence before making numerical changes.

## Repository state at this documentation cut

- Repository: reverez/AstraFlow
- Stable release baseline: annotated tag v1.0.0 at ecd0337
- Remote scientific-closure branch at documentation cut: v1.1-scientific-closure at 6cd80f2
- This documentation work is isolated on docs/future-development-roadmap
- main must remain untouched until scientific review explicitly approves a merge.

### Known local-only work

The project owner reported a later local commit, 72a6ec2, that was **not pushed** when this documentation branch was created. That experiment implemented optional local pseudo-time RK2 while keeping physical-time SSP-RK2 as the default.

Reported 128×32 pseudo-time result:

| Gate | Result |
|---|---:|
| Maximum residual | 6.792e-7 — pass |
| Engineering stability | pass |
| Conservative-field Linf <= 1e-8 | 2.794e-6 — fail |
| Maximum primary engineering relative difference | 4.764e-9 — pass |
| Numerical boundary mass mismatch | 2.661e-9 — pass |
| Fallbacks | zero |
| Wall time | 44.502 s vs 64.689 s physical-time |
| Scientific status | **not accepted** |

Before continuing development, reconcile that local commit with GitHub rather than reimplementing it or allowing branches to diverge silently.

## What AstraFlow is

AstraFlow is a C++20/CUDA scientific-computing codebase for **2D axisymmetric compressible Euler/Navier–Stokes nozzle flow**. It has:

- a deterministic CPU reference backend;
- a CUDA production backend;
- shared finite-volume mathematics;
- MUSCL reconstruction with HLLC/HLLE fluxes;
- viscous stress and heat conduction;
- SSP-RK2 physical-time integration;
- parameterised nozzle geometry;
- engineering integrals;
- CLI and Dear ImGui/ImPlot GUI;
- JSON/CSV/VTK output;
- canonical verification and measured CPU/GPU parity.

The code is **AI-assisted, architect-directed**. Generated implementation is subordinate to mathematical design, tests, measured evidence, and review.

## Scientific status

The accepted 128×32 rocket case reaches steady closure at 82,480 iterations in CPU and CUDA FP64. The CUDA run has final maximum residual 8.56e-7, zero fallbacks, and an extremely stable engineering window.

The refined 256×64 case remains the scientific blocker. At 250,000 iterations its maximum residual remained about 1.14e-4; conservative coarse-to-fine initialization reduced the energy residual to about 4.57e-5, but not to the unchanged 1e-6 gate.

C2b diagnostics established:

- no demonstrated CPU/CUDA defect;
- finite-volume boundary/source assembly agrees to roundoff (maximum relative identity error about 6.88e-15);
- actual numerical inlet/outlet mass-flux mismatch is about 1e-7, much tighter than the legacy interior-cell estimate;
- strong reconstruction-sensitive signed residual cycling exists near the chamber/contraction join and throat;
- MC limiter branches switch frequently in neighboring cells;
- Van Leer changes oscillation amplitude, but altered reconstructions remain diagnostic only;
- a permanent limit cycle has **not** been proven.

Therefore:

> **V1.1 has baseline steady convergence but does not yet have refined-grid scientific acceptance, rocket grid independence, Richardson extrapolation, or GCI.**

## Non-negotiable invariants

Future development must preserve these unless a reviewed defect proves one wrong:

1. **Do not weaken convergence tolerances to make a case pass.**
2. **Do not use a visually stable flow field as evidence of steady convergence.**
3. **Do not treat the interior-cell inlet/outlet estimate as the same quantity as the true finite-volume boundary flux balance.**
4. **Physical-time SSP-RK2 remains the reference transient integrator.**
5. **Any pseudo-time method is a steady root finder only; pseudo-time is never physical simulated time.**
6. **CPU and CUDA must continue to instantiate the same mathematical operator where designed to do so.**
7. **Performance optimisation must be followed by scientific/parity regression.**
8. **No Phase 2 result may be presented as implemented or validated before its own gates pass.**
9. **No experimental rocket-engine validation claim is permitted from numerical verification alone.**

## Immediate continuation order

~~~mermaid
flowchart TD
    A["Reconcile local 72a6ec2 with remote branch"] --> B["Root-equivalence closure at 128×32"]
    B --> C{"Physical and pseudo solvers approach same discrete root?"}
    C -->|No| D["Diagnose root / reconstruction dependence<br/>STOP refined-grid study"]
    C -->|Yes| E["Retry 256×64 steady closure"]
    E --> F{"Residual + engineering gates pass?"}
    F -->|No| G["Investigate acceleration / nonlinear steady strategy<br/>without changing R(U)=0"]
    F -->|Yes| H["Run 512×128"]
    H --> I["Rocket grid comparison"]
    I --> J["Richardson / GCI where mathematically valid"]
    J --> K["V1.1 scientific acceptance review"]
    K --> L["Begin Phase 2"]
~~~

## Next scientific experiment

The next experiment should determine whether the physical-time and pseudo-time 128×32 terminal states are finite-tolerance approximations to the same root.

Required evidence:

- evaluate the same production residual R(U) directly on both saved states;
- cross-relax each state with the opposite integrator;
- run a prospective tighter residual experiment (for example 1e-8) with criteria fixed before results are inspected;
- localise the existing 2.794e-6 conservative Linf discrepancy;
- track field distance as residual decreases;
- stop before 256×64 if root equivalence is not established.

Do not reinterpret the historical failed 1e-8 field-equivalence gate after the fact. Preserve it as a failed experiment and create a new prospective protocol.

## Where to look

| Need | Source |
|---|---|
| Current code architecture | [../architecture/overview.md](../architecture/overview.md) |
| Future architecture | [../architecture/evolution.md](../architecture/evolution.md) |
| V1 evidence | [../verification/verification_report.md](../verification/verification_report.md) |
| V1.1 closure | [../verification/rocket_scientific_closure.md](../verification/rocket_scientific_closure.md) |
| Residual investigation | [../verification/residual_closure_c2b.md](../verification/residual_closure_c2b.md) |
| Performance methodology | [../performance/benchmark_report.md](../performance/benchmark_report.md) |
| Future releases | [../roadmap/roadmap.md](../roadmap/roadmap.md) |
| Phase 2 implementation | [../roadmap/phase2_plan.md](../roadmap/phase2_plan.md) |
| Engineering process | [workflow.md](workflow.md) |
| Durable decisions | [decision_log.md](decision_log.md) |

## Anti-patterns for future agents

Do not:

- perform broad rewrites before understanding verification coverage;
- introduce turbulence, chemistry, plume physics, AMR, or new Riemann solvers during V1.1 closure;
- rerun expensive studies after documentation-only edits;
- delete or overwrite failed scientific runs;
- report pseudo-time iteration count as physical time;
- hard-code expected speedup;
- add GPU optimisations without parity checks;
- turn the GUI into the owner of physics;
- create an alternative solver stack “just in case”.

Prefer targeted, evidence-driven changes with one scientific question per milestone.
