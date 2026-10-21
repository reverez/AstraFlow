# V1.3 continuation protocol (pre-registered)

Written on 2026-10-21 while the V1.2 256×64 run was at iteration 216,180, before its outcome
under the 400,000-iteration cap was known. Nothing here is edited after the V1.2 or V1.3 results
are in; results go to [rocket_closure_v12.md](rocket_closure_v12.md).

## Why

The V1.2 protocol fixed iteration caps before any run (150k / 400k / 1M) and states that caps are
not raised after a failure. The first 216k iterations of the V1.2 256×64 van Albada run show a
monotone exponential decay of the energy residual with no plateau: least-squares e-folding
lengths of 119,848 iterations (100k–150k) and 135,679 iterations (150k–200k), energy residual
1.1e-4 at 61k, 5.0e-5 at 150k, 3.4e-5 at 200k. Extrapolating that rate puts the 1e-6 target near
iteration 620k–680k, beyond the V1.2 cap. The V1.2 cap was therefore set too low for this
scheme, which is a budgeting error, not evidence of a limit cycle. V1.2 is completed exactly as
written; V1.3 is a separate study whose only change is the budget.

## Hypothesis

H1': the V1.2 medium and fine van Albada integrations converge to the unchanged 1e-6 residual
gate when given a budget consistent with their measured decay rate.

## Design

- **Medium (256×64):** continue the V1.2 trajectory from its preserved final state. The
  preserved `final_state.vts` (17 significant digits, physical units) is converted to a schema-1
  state file and supplied through `--initial-state` with an identical configuration. This is the
  same physical-time integration with a round-off-level restart perturbation (re-nondimensionalised
  double precision), not a different experiment; the convergence monitor restarts its 5,000-
  iteration minimum and 2,000-iteration window, which only makes acceptance harder. Budget for the
  continuation: **700,000 further iterations** (total ≤ 1.1M), i.e. about three e-folding lengths
  beyond the extrapolated crossing.
- **Fine (512×128):** fresh integration from the quasi-1D initial state as in V1.2, cap
  **2,500,000 iterations**. The timestep halves again and the decay length in iterations is
  expected to roughly double; the cap is 2.5× the medium total.
- **Coarse (128×32):** the V1.2 result (49,480 iterations, `steady_converged`) is reused unchanged.
- Convergence definition, gates, scheme and configuration are exactly those of V1.2.

## Acceptance

Identical to V1.2 items 1–4: three accepted grids with zero fallbacks, the CPU FP64 cross-check
on 128×32 (already part of the V1.2 campaign), valid Richardson analysis with reported GCI, and
the 0.5% / 1% consecutive-grid goals on the medium/fine pair. If the medium continuation reaches
its cap without closure, V1.3 stops there and reports it; no further budget is granted inside V1.3.

## Preservation

Continuation outputs live in `runs/v13-closure/<grid>/continuation-NNN` or `attempt-NNN` with
provenance (source run checksums, binary SHA256, commit). The V1.2 directories are never written
to. The analysis manifest lists the coarse V1.2 run, the medium continuation and the fine run.
