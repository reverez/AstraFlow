# AstraFlow documentation

This directory is the engineering and scientific source map for AstraFlow. It separates **implemented evidence**, **current architecture**, and **future intent** so contributors and coding agents do not confuse roadmap items with completed functionality.

## Current state

| Area | State |
|---|---|
| V1 implementation | Complete and preserved by v1.0.0 |
| CPU/CUDA verification | Complete for V1 canonical cases |
| V1.1 128×32 rocket steady closure | Accepted |
| V1.1 refined-grid closure | Unresolved at 256×64 |
| Rocket grid independence / GCI | Not established |
| Phase 2 physics | Not started |
| Main development rule | Do not weaken scientific gates to advance the roadmap |

The remote v1.1-scientific-closure branch was at 6cd80f2 when this documentation branch was created. The developer context records known local-only work that must be reconciled before future implementation continues.

## Documentation map

### Start here

- [Developer context](development/developer_context.md) — compact state handoff, invariants, blockers, and next actions.
- [Roadmap](roadmap/roadmap.md) — release sequence, dependencies, scope boundaries, and acceptance gates.
- [Phase 2 plan](roadmap/phase2_plan.md) — detailed implementation plan for the next major physics/product phase.
- [Development workflow](development/workflow.md) — branch, testing, review, evidence, and release discipline.
- [Decision log](development/decision_log.md) — durable architectural and scientific decisions.

### Current implementation

- [Architecture overview](architecture/overview.md)
- [Architecture evolution](architecture/evolution.md)
- [Environment](environment.md)
- [Phase 1 status](PHASE1_STATUS.md)

### Mathematics

- [Governing equations](mathematics/governing_equations.md)
- [Numerical method](mathematics/numerical_method.md)
- [Axisymmetric Navier–Stokes](mathematics/axisymmetric_navier_stokes.md)
- [Boundary conditions](mathematics/boundary_conditions.md)

### Evidence

- [V1 verification report](verification/verification_report.md)
- [V1.1 rocket scientific closure](verification/rocket_scientific_closure.md)
- [V1.1 C2b residual investigation](verification/residual_closure_c2b.md)
- [Performance report](performance/benchmark_report.md)

## Source-of-truth hierarchy

When documents disagree, use this order:

1. **Measured artifacts and verification reports** for what has actually been demonstrated.
2. **Checked-in mathematics and current architecture docs** for how the present solver is defined.
3. **Developer context and decision log** for current project intent and invariants.
4. **Roadmap documents** for planned work only.

A roadmap statement is never evidence that a feature exists.

## Project lifecycle

~~~mermaid
flowchart LR
    V10["v1.0<br/>Verified solver platform"] --> V11["v1.1<br/>Scientific closure"]
    V11 -->|closure gates pass| V12["v1.x<br/>Reliability + reproducibility"]
    V12 --> V20["v2.0<br/>External-flow CFD"]
    V20 --> V21["v2.x<br/>Thermophysical fidelity"]
    V21 --> V22["v2.x<br/>Design studies + optimisation"]
    V22 --> V30["v3.x<br/>Turbulence / reacting flow / CHT"]

    V11 -. "blocked today: refined-grid residual closure" .-> V11
~~~

The roadmap is intentionally staged. Higher-fidelity physics is added only after the lower-level numerical evidence remains intact.
