# Decision log

This is a lightweight architectural decision record. It captures decisions that future contributors should treat as defaults until new evidence justifies a reviewed change.

## D-001 — Preserve the verified V1 baseline

**Status:** accepted.

v1.0.0 is the immutable Phase 1 reference. Later research branches may evolve architecture, but comparisons must retain a reproducible path back to the accepted solver.

**Consequence:** do not rewrite or move the tag; do not use later experimental behavior to retroactively redefine V1 evidence.

## D-002 — CPU reference + CUDA production backends

**Status:** accepted.

The CPU implementation remains a deterministic reference and CI-capable backend. CUDA is the high-performance backend. Shared mathematical functions should be reused where this improves parity without coupling execution/memory architecture.

**Consequence:** GPU correctness is demonstrated by numerical evidence, not by visual plausibility.

## D-003 — Scientific closure precedes higher-fidelity claims

**Status:** active.

The 128×32 rocket case is steady-converged; the 256×64 refined-grid case is not yet accepted. Rocket grid independence and GCI remain unresolved.

**Consequence:** Phase 2 may be designed in documentation, but major new physics should not obscure or replace the V1.1 closure question.

## D-004 — Distinguish conservation metrics

**Status:** accepted.

The legacy engineering inlet/outlet estimate samples interior cell states. It is not identical to the finite-volume numerical boundary-flux balance.

**Consequence:** preserve the legacy metric for compatibility, but report the true numerical flux balance explicitly. Do not silently substitute one metric for another in historical acceptance criteria.

## D-005 — Pseudo-time is a steady solver, not physical time

**Status:** provisional pending root-equivalence acceptance.

Local pseudo-time stepping may accelerate convergence to R(U)=0, but it must not be reported as physical transient evolution. The physical SSP-RK2 path remains the transient reference.

**Consequence:** pseudo-time must carry explicit metadata and must establish root equivalence before it is used for the refined grid study.

## D-006 — External plume before reacting chemistry

**Status:** planned.

The next major physics expansion should first extend the non-reacting axisymmetric solver outside the nozzle, where the existing compressible finite-volume machinery can be exercised against new far-field/shock structures.

**Consequence:** chemistry is not the first Phase 2 deliverable. It belongs after geometry/domain, far-field boundaries, thermophysical abstraction, and external-flow verification are mature.

## D-007 — Introduce thermophysical interfaces before high-fidelity models

**Status:** planned.

The current calorically perfect ideal gas remains the reference model. Future variable-property thermodynamics and transport should enter through narrow model interfaces rather than conditionals spread across flux kernels and GUI code.

**Consequence:** constant-gamma behavior must remain a verified limiting mode.

## D-008 — Prefer staged mesh evolution

**Status:** planned.

Phase 2 should first attempt an axisymmetric structured or multi-block structured external domain compatible with the current solver architecture. A full unstructured-mesh rewrite is not justified merely to render a plume.

**Consequence:** change mesh representation only when the required boundary topology or quality cannot be represented robustly by the staged approach.

## D-009 — Optimisation follows trustworthy simulation

**Status:** planned.

Automated nozzle design is downstream of solver closure and Phase 2 physics verification.

**Consequence:** first build deterministic parameter sweeps and sensitivity evidence; only then introduce search/optimisation algorithms.

## D-010 — Advanced turbulence/reacting flow is a separate fidelity tier

**Status:** planned for later major releases.

RANS/LES, species transport, chemistry, conjugate heat transfer, and regenerative cooling create new equations, closures, stiffness, reference-data requirements, and numerical risks.

**Consequence:** they are not small Phase 2 extras. They require their own verification programs and release gates.
