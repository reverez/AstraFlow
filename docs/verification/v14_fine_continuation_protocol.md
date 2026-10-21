# V1.4 fine-grid continuation protocol (pre-registered)

Written on 2026-10-21 after the V1.3 512×128 run reached its 2,500,000-iteration cap
(`iteration_limit`, energy residual 2.17e-5, zero fallbacks, 4,861 s wall) and before any
further fine-grid iteration was run. Nothing here is edited after the V1.4 result is known.

## Measured basis

Energy residual of the V1.3 fine run: 1.28e-4 at 0.5M, 7.74e-5 at 1.0M, 4.87e-5 at 1.5M,
3.20e-5 at 2.0M, 2.17e-5 at 2.5M iterations. Least-squares e-folding lengths: 1.10M (1.0–1.5M),
1.21M (1.5–2.0M), 1.27M (2.0–2.5M) iterations. The decay is monotone with no plateau but the
e-folding length is still growing slowly. Reaching 1e-6 needs 3.1 more e-folds; extrapolating the
last measured rate puts the crossing near 6.2–6.4M total iterations.

## Design

Continue the preserved V1.3 trajectory through `--initial-state` (schema-1 state rebuilt from the
17-digit `final_state.vts`, identical configuration, round-off-level restart perturbation), as in
V1.3 for the medium grid. Budget: **4,500,000 further iterations** (total ≤ 7.0M), which is the
extrapolated requirement plus about one e-folding of margin for the observed slow growth of the
decay length. Physical-time SSP-RK2, van Albada, CUDA FP64, every gate unchanged. If the cap is
reached without closure, V1.4 stops and reports; no further budget is granted inside V1.4.

## Acceptance

With three accepted grids: valid Richardson analysis and reported fine-grid GCI for at least mass
flow, thrust and exit Mach, and the 0.5% / 1% consecutive-grid goals on the medium/fine pair.
Otherwise the study is reported as partial with the measured decay record.
