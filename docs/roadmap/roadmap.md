# AstraFlow development roadmap

This roadmap defines **sequence and acceptance**, not promises or completion claims. The project advances only when each scientific dependency is sufficiently resolved.

## Roadmap at a glance

~~~mermaid
flowchart TD
    P1["v1.0 — Verified internal-flow platform<br/>COMPLETE"] --> C["v1.1 — Rocket scientific closure<br/>ACTIVE / BLOCKED"]
    C --> H["v1.x — Reliability + reproducibility hardening"]
    H --> E["v2.0 — External nozzle/plume CFD"]
    E --> T["v2.x — Variable thermodynamics + transport"]
    T --> W["v2.x — Wall loads + design studies"]
    W --> O["v2.x — Nozzle optimisation"]
    O --> F["v3.x — Turbulence / species / chemistry / CHT"]

    C -. "must establish accepted refined-grid path" .-> C
~~~

## Release map

| Release | Objective | Major evidence required | Status |
|---|---|---|---|
| v1.0 | Verified CPU/CUDA internal nozzle CFD platform | canonical verification, parity, GUI, benchmark | complete |
| v1.1 | Scientifically close the rocket steady solution | root convergence, refined grids, uncertainty where valid | active |
| v1.x | Harden reproducibility and scientific UX | restart/provenance, robust report generation, regression assets | planned |
| v2.0 | Extend domain through nozzle exit into ambient | external boundaries, plume shock structures, conservation, grid study | planned |
| v2.x | Improve thermophysical fidelity | variable cp/gamma, viscosity, conductivity with verified limiting mode | planned |
| v2.x | Add wall/load analysis and design studies | wall pressure/shear/heat-flux evidence, controlled parameter sweeps | planned |
| v2.x | Add automated nozzle optimisation | deterministic objective/constraints, reproducible search | planned |
| v3.x | Turbulence, reacting species, thermal coupling | model-specific verification/reference suites | deferred |

## Gate 0 — Finish V1.1 scientific closure

This is the highest-priority dependency.

### Required before V1.1 acceptance

- reconcile and review the local pseudo-time experiment;
- establish whether pseudo-time and physical-time approach the same 128×32 discrete root;
- obtain an accepted 256×64 solution without weakening the spatial operator or convergence criteria;
- run 512×128 only after medium-grid acceptance;
- compare global engineering quantities and spatial profiles;
- perform Richardson extrapolation/GCI only where the convergence sequence satisfies the assumptions;
- document any quantity for which GCI is invalid rather than manufacturing an order.

### Exit statement

The desired defensible claim is:

> For the documented ideal-gas axisymmetric rocket-nozzle configuration, AstraFlow reaches a reproducible steady discrete solution under explicit residual and engineering-stability criteria, with mesh-sensitivity quantified across accepted systematically refined grids.

This remains numerical verification, not experimental validation.

## Gate 1 — V1.x reliability and reproducibility

Keep this release small. It exists to make future science safer.

Candidate work:

- robust restart/checkpoint format with explicit schema version;
- provenance manifest linking executable/config/source hashes;
- automated comparison/report generation;
- CI separation of fast unit tests from expensive scientific workflows;
- optional GPU diagnostic workflow when suitable runners/tools are available;
- documented migration policy for result schemas.

Do not turn this into a feature backlog.

## Phase 2 — Non-reacting high-fidelity nozzle + plume

Phase 2 should be an incremental physics extension, not a rewrite.

Primary themes:

1. external axisymmetric plume domain;
2. robust characteristic/far-field boundaries;
3. variable thermodynamics/transport behind interfaces;
4. wall pressure, shear and heat-flux diagnostics;
5. richer scientific visualisation and probes;
6. evidence-backed CUDA optimisation on the larger domain.

See [phase2_plan.md](phase2_plan.md).

### Phase 2 success criteria

A v2.0 release should demonstrate:

- nozzle-to-ambient flow in one conservative domain;
- stable under-expanded and over-expanded cases within the implemented model;
- no unexplained mass/energy imbalance;
- canonical/independent verification for every new model component;
- mesh sensitivity of selected plume and nozzle observables;
- CPU/CUDA parity at representative scale;
- measured GPU performance on the expanded workload;
- explicit limitations around axisymmetry, idealisation, turbulence, and chemistry.

## Phase 2 design studies

After v2.0 physics is trustworthy:

- controlled sweeps over area ratio, back pressure, chamber conditions, and contour parameters;
- wall-load and performance maps;
- sensitivity analysis;
- geometry optimisation with explicit objectives and constraints.

Optimisation must never become a substitute for solver verification.

## Phase 3 — Advanced fidelity

Treat each item as a separate scientific program:

### Turbulence

Start with a RANS model only after laminar/external-flow infrastructure is stable. Require canonical wall/jet benchmarks and near-wall resolution policy. LES is later still.

### Multi-species and reacting flow

Requires:

- species conservation equations;
- thermochemical property database;
- reaction mechanism representation;
- stiff chemistry integration strategy;
- independent reference comparisons;
- new positivity and conservation requirements.

### Thermal coupling

Wall conduction, conjugate heat transfer, regenerative cooling, and ablation require solid-domain models and interface energy conservation.

### 3D / AMR / multi-GPU

These are architectural scale changes, not prerequisites for a strong axisymmetric research platform.

## Dependency graph

~~~mermaid
graph TD
    A["V1.1 root + mesh closure"] --> B["Reproducible steady workflow"]
    B --> C["Extended axisymmetric domain"]
    C --> D["Far-field / plume boundaries"]
    D --> E["External plume verification"]
    E --> F["Thermophysical abstraction"]
    F --> G["Variable properties"]
    E --> H["Wall loads"]
    G --> H
    H --> I["Design sweeps"]
    I --> J["Optimisation"]

    E --> K["Turbulence tier"]
    G --> L["Species / chemistry tier"]
    H --> M["Conjugate thermal tier"]
~~~

## Scope control

A milestone may be promoted only if:

- its prerequisite evidence exists;
- the implementation has a single clear scientific purpose;
- acceptance criteria are set before the final experiment;
- failure can stop the roadmap without forcing a redesign of unrelated modules.

The roadmap should evolve through reviewed decisions, not accumulation of features.
