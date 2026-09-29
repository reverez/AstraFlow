# Development workflow

AstraFlow uses a **verification-gated scientific development process**. The goal is not simply to keep the build green; each change must preserve the meaning of the numerical results.

## Change lifecycle

~~~mermaid
flowchart LR
    Q["Question / defect / milestone"] --> S["Specification"]
    S --> I["Targeted implementation"]
    I --> T["Unit + numerical tests"]
    T --> E["Measured experiment"]
    E --> R["Scientific / architectural review"]
    R -->|pass| C["Commit + push"]
    R -->|fail| D["Preserve evidence<br/>refine diagnosis"]
    D --> S
~~~

A failed experiment is a valid result. Do not erase it by changing the gate after seeing the outcome.

## Branch strategy

| Branch / tag | Role |
|---|---|
| main | accepted public baseline only |
| v1.0.0 | immutable accepted Phase 1 tag |
| v1.1-scientific-closure | scientific-closure development |
| docs/* | isolated documentation work |
| future feature/* | one bounded implementation milestone |
| future experiment/* | non-production scientific diagnostics |

Rules:

- branch from a known accepted parent;
- do not force-push accepted/shared scientific history;
- commit only coherent milestones;
- never mix unrelated numerical changes and documentation cleanup in the same scientific commit;
- preserve raw/compact evidence according to repository policy;
- merge to main only after the milestone's acceptance document is complete.

## Work classification

### A. Documentation-only

Examples: roadmap, explanatory diagram, typo, link correction.

Minimum gate:

- links resolve;
- no measured claim is changed without checking its source;
- no roadmap item is written as completed functionality.

### B. Software-only, numerically neutral

Examples: CLI ergonomics, serialization plumbing, GUI layout, build scripts.

Minimum gate:

- targeted unit tests;
- CPU test suite;
- CUDA tests if GPU-facing code is touched;
- output/schema compatibility review where applicable.

### C. Numerical implementation

Examples: residual assembly, reconstruction, boundary condition, time integrator, thermodynamics.

Minimum gate:

- unit tests for the new mathematics;
- canonical verification affected by the change;
- CPU/CUDA parity where both backends apply;
- conservation/finite-state checks;
- before/after evidence;
- documented acceptance bounds fixed before final measurement.

### D. Performance change

Minimum gate:

- same scientific state or stated tolerance;
- CPU/CUDA parity unchanged;
- measured benchmark with warmup and defined timing region;
- profiler evidence for the claimed bottleneck;
- no speedup claim from one unstable sample.

### E. New physics

Minimum gate:

- mathematical model documented first;
- limiting case reduces to an already verified model when possible;
- independent analytical/reference case;
- grid/refinement evidence appropriate to the model;
- explicit model limitations;
- no broad fidelity claim from qualitative images.

## Acceptance matrix

| Change | Unit | CPU | CUDA | Analytical/reference | Grid study | Performance |
|---|---:|---:|---:|---:|---:|---:|
| Docs | optional | — | — | — | — | — |
| I/O / GUI | yes | yes | as touched | — | — | — |
| Numerical core | yes | yes | yes | yes | as relevant | regression |
| Thermodynamics | yes | yes | yes | yes | as relevant | regression |
| External plume | yes | yes | yes | yes | required | measured |
| Turbulence | yes | yes | yes | benchmark/reference | required | measured |
| Reacting flow | yes | yes | yes | independent reference | required | measured |

“Required” does not mean every commit runs the entire expensive suite; it means the milestone cannot be accepted without that evidence.

## Scientific evidence policy

Each accepted scientific milestone should leave:

- checked-in configuration;
- exact backend and precision;
- termination reason;
- convergence criteria;
- compact numerical result tables;
- reproduction command;
- raw-data location or retained checksum;
- figures generated from actual outputs;
- limitations and failed prerequisites.

Do not report:

- a transient quantity as a converged steady result;
- grid independence from a single mesh;
- GCI from an invalid/non-monotonic sequence;
- GPU speedup against an unspecified CPU baseline;
- “validation” when only verification has been performed.

## Commit guidance

Recommended prefixes:

- feat: production capability
- fix: demonstrated defect correction
- verify: scientific verification or closure evidence
- perf: evidence-backed performance change
- docs: documentation only
- test: test infrastructure

Scientific commit messages should describe the result, not the aspiration.

Good:

~~~text
verify: establish 128x32 steady rocket closure
perf: fuse residual diagnostics after profile evidence
docs: define phase 2 external-flow roadmap
~~~

Avoid vague messages such as “fix everything”, “improve CFD”, or “make GPU faster”.

## AI-assisted implementation protocol

AI coding agents may implement substantial code, but they are not the scientific authority.

Before implementation, give the agent:

1. exact milestone;
2. current branch/commit;
3. mathematical invariants;
4. files/evidence that are authoritative;
5. prohibited scope;
6. test gate;
7. stop conditions.

After implementation, review:

- whether the code implements the stated equations rather than an easier substitute;
- whether tolerances were changed;
- whether failed runs were hidden or overwritten;
- whether results are measured;
- whether README claims match verification documents.

## Release gate

~~~mermaid
flowchart TD
    A["Feature complete"] --> B["Targeted tests"]
    B --> C["Full required regression"]
    C --> D["Scientific evidence complete"]
    D --> E["Docs + limitations consistent"]
    E --> F{"Architectural review"}
    F -->|pass| G["Tag / merge"]
    F -->|fail| H["Keep branch; preserve evidence"]
~~~

A release is a statement about evidence, not merely source-code completeness.
