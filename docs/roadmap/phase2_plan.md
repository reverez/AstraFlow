# Phase 2 implementation plan

> **Target:** evolve AstraFlow from a verified internal axisymmetric nozzle solver into a verified **non-reacting nozzle-to-plume CFD research platform** while preserving V1 reference behavior.

Phase 2 begins only after the V1.1 scientific-closure decision is resolved sufficiently to prevent new physics from masking unresolved steady-state behavior.

## Design principles

- extend, do not rewrite, the verified finite-volume core;
- keep the constant-gamma ideal-gas model as a reference mode;
- preserve CPU/CUDA parity architecture;
- isolate thermodynamics, transport, boundaries, geometry, and analysis behind explicit interfaces;
- add one physical capability per milestone with its own verification;
- retain headless reproducibility as the primary scientific interface;
- treat GUI changes as consumers of solver state, never owners of physics.

## Target architecture

~~~mermaid
flowchart LR
    CFG["Configuration"] --> DOM["Domain / multi-block mesh"]
    DOM --> OP["Finite-volume residual R(U)"]
    TH["Thermodynamics model"] --> OP
    TR["Transport model"] --> OP
    BC["Boundary models"] --> OP
    OP --> CPU["CPU backend"]
    OP --> GPU["CUDA backend"]
    CPU --> SIM["Simulation / steady solver"]
    GPU --> SIM
    SIM --> ANA["Engineering + plume + wall analysis"]
    SIM --> IO["JSON / CSV / VTK"]
    SIM --> GUI["GUI snapshots"]
    ANA --> REP["Verification / reports"]

    REF["V1 constant-gamma reference"] -. regression .-> TH
    REF -. regression .-> TR
~~~

## P2-M0 — Entry gate and baseline lock

### Goal

Begin Phase 2 from an explicitly accepted scientific baseline.

### Required

- V1.1 decision recorded;
- main/release tag points to the accepted baseline;
- all canonical CPU/CUDA tests green;
- a representative steady rocket result is reproducible from a checked-in config;
- result schema/provenance policy documented.

### Stop condition

Do not begin plume work if the accepted baseline cannot be reproduced.

---

## P2-M1 — Physics/model abstraction

### Goal

Prepare for variable thermodynamics and transport without changing existing results.

Introduce narrow conceptual interfaces such as:

~~~text
ThermoModel
  pressure(U)
  temperature(U)
  sound_speed(U)
  cp(T)
  gamma(T)

TransportModel
  viscosity(T)
  conductivity(T)
~~~

Exact API shape should follow the codebase rather than this pseudocode.

### Requirements

- existing calorically perfect gas becomes one concrete model;
- existing constant viscosity/Prandtl behavior becomes one concrete transport model;
- CPU/CUDA-friendly data representation; avoid virtual dispatch inside hot device loops if it harms compilation/performance;
- no GUI-specific physics paths.

### Acceptance

Existing V1 verification and parity must remain numerically equivalent within preset regression tolerance.

---

## P2-M2 — Extended axisymmetric domain

### Goal

Represent the internal nozzle and an external ambient region in one conservative domain.

Preferred progression:

1. retain structured topology where possible;
2. introduce multi-block structured representation if one block cannot maintain quality;
3. defer a general unstructured rewrite unless evidence shows it is required.

Conceptual domain:

~~~text
 chamber        throat        nozzle exit             far field
+---------+       |       /----------------\
|         +-------+------/                  \________
|   internal flow |     /        plume                 -> x
|         +-------+------\                  /--------
+---------+       |       \----------------/
                         ambient domain
~~~

### Required geometry evidence

- positive cell volumes;
- face closure/orientation;
- block-interface conservation if multi-block;
- axis treatment;
- mesh-quality metrics;
- regression of the original internal-only domain.

### Stop condition

No plume simulation until uniform-flow preservation across the extended domain and any block interfaces passes.

---

## P2-M3 — External-flow boundary conditions

### Goal

Support a finite ambient domain without artificial reflections dominating the plume.

Boundary types likely required:

- axis symmetry;
- solid nozzle wall;
- reservoir inlet;
- internal supersonic/subsonic flow transition already supported;
- far-field/outflow characteristic treatment appropriate to the local Mach regime.

### Verification strategy

Use small cases before a rocket plume:

- uniform free-stream preservation;
- supersonic outflow with no upstream contamination;
- prescribed subsonic far-field relaxation;
- shock/expansion cases with known reference behavior where geometry permits;
- global finite-volume balance.

Do not accept a far-field condition based on a visually plausible plume alone.

---

## P2-M4 — Non-reacting plume solver

### Goal

Simulate the nozzle exit jet and ambient interaction using the existing conservative compressible operator.

Study cases:

- near ideally expanded;
- under-expanded;
- over-expanded.

Expected structures may include expansion fans, compression waves and shock-cell behavior, but the documentation must distinguish model-supported observations from experimental validation.

### Outputs

Add plume-oriented analysis:

- centerline Mach/pressure/temperature;
- plume width proxy;
- shock-cell locations where objectively detectable;
- axial momentum flux;
- domain mass/energy balance;
- nozzle thrust accounting consistent with the chosen control surface.

### Acceptance

- stable finite fields;
- no unexplained fallback growth;
- boundary/global conservation;
- mesh sensitivity on selected plume observables;
- CPU/CUDA parity on a representative case;
- far-field sensitivity study demonstrating boundaries are sufficiently remote or nonreflecting for the reported quantities.

---

## P2-M5 — Variable thermodynamics and transport

### Goal

Move beyond a constant-gamma, constant-property gas without immediately introducing chemistry.

Candidate progression:

1. Sutherland-law viscosity;
2. temperature-dependent conductivity / Prandtl treatment;
3. temperature-dependent cp(T);
4. derive h(T), e(T), a(T), and gamma(T) consistently.

### Numerical implications

Primitive recovery becomes more complex when e(T) is nonlinear. Implement a bounded, tested inversion strategy and ensure CPU/CUDA results match.

### Verification

- constant-property mode reproduces the existing solution;
- unit tests across temperature range;
- thermodynamic consistency identities;
- independent reference table or tool comparison;
- nozzle comparison showing the effect of the model without presenting it as experimental truth.

Reacting chemistry remains out of scope.

---

## P2-M6 — Wall loads and thermal diagnostics

### Goal

Make the internal nozzle solution useful for engineering analysis.

Add, where mathematically supported:

- wall static pressure;
- wall shear stress;
- adiabatic heat flux (zero by model at an adiabatic wall);
- non-adiabatic prescribed-temperature wall heat flux if that boundary mode is introduced;
- integrated axial wall force;
- distributions in normalized and SI coordinates.

### Acceptance

Use canonical viscous/thermal cases to verify gradients and integrated loads before publishing rocket values.

---

## P2-M7 — Scientific visualisation and UX

The GUI should expose Phase 2 results without becoming a new solver.

Candidate additions:

- external plume field rendering;
- centerline and wall probes;
- line plots with export;
- selectable control surfaces;
- shock-gradient overlay;
- run metadata and convergence state;
- comparison of two saved runs.

Avoid decorative 3D rendering unless it communicates new scientific information.

---

## P2-M8 — Performance engineering

The larger external domain changes the performance balance.

Required process:

~~~mermaid
flowchart LR
    M["Measure"] --> P["Profile"]
    P --> B["Identify bottleneck"]
    B --> O["One optimisation"]
    O --> V["Parity + verification"]
    V --> M
~~~

Candidates only after profiling:

- face-kernel fusion/splitting;
- memory-layout changes;
- reduction frequency;
- snapshot cadence;
- block-interface scheduling;
- CUDA graph capture if launch overhead becomes relevant.

Do not optimise by changing the numerical method.

---

## P2-M9 — Phase 2 scientific acceptance

A v2.0 release candidate must include:

### Software

- CPU build;
- CUDA build;
- CLI;
- GUI smoke;
- serialization and restart/provenance;
- no sanitizer regressions in CPU-supported code.

### Numerical

- all V1 canonical verification;
- extended-domain uniform preservation;
- far-field boundary tests;
- external plume reference cases;
- global conservation;
- CPU/CUDA parity;
- mesh/far-field sensitivity.

### Documentation

- governing equations and new thermodynamics;
- external boundary conditions;
- mesh/domain design;
- plume verification report;
- performance report;
- known limitations;
- reproduction commands.

## Explicit Phase 2 exclusions

Unless separately promoted through review:

- combustion chemistry;
- species transport;
- RANS/LES;
- full 3D;
- AMR;
- multi-GPU;
- regenerative cooling;
- ablation;
- structural mechanics.

These belong to later fidelity tiers.

## Post-v2.0 design and optimisation

After the non-reacting plume platform is accepted:

1. deterministic parameter sweep;
2. sensitivity analysis;
3. objective/constraint definition;
4. geometry search;
5. independent re-evaluation of candidate optima on finer grids.

Candidate objectives:

- thrust;
- specific impulse within model;
- exit-pressure matching;
- plume compactness proxies;
- wall-load constraints.

Optimisation results must carry the uncertainty of the underlying CFD model.
