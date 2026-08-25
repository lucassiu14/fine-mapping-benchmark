# Iteration 006 — what to revert when the iteration is retired

**Iterations 005 and 006 are exploratory. The simulator, the benchmark and the
analysis all go back to their standing form afterwards.** Iteration 004 remains
the project's standard and none of its files are touched by either iteration —
verified: every `iter004_*` script is byte-identical to the version that
produced the standing results.

Every temporary edit is wrapped in a banner so it can be found mechanically:

```
# >>> ITERATION 006 ONLY - TEMPORARY (see iteration-006-REVERT.md) <<<
...
# >>> END ITERATION 006 TEMPORARY BLOCK <<<
```

Find them all with:

```bash
grep -rn "ITERATION 006 ONLY\|ITER-006 (temp" R/ scripts/
```

Iteration 005's own list is in [iteration-005-REVERT.md](iteration-005-REVERT.md);
revert 006 first, then 005, since 006 builds on 005's plumbing.

## What to remove

| # | file | what to remove | notes |
|---|---|---|---|
| 1 | `R/simulate_phenotypes.R` | the ITERATION 006 block: `.w_for_r2()`, `.orthogonal_nonlinear()`, `.representability_log_weights()` | self-contained; delete the whole banner block |
| 2 | `R/simulate_phenotypes.R` | `target_r2` / `nl_depth` args on `simulate_phenotypes()` and `select_causal_variants()`, and the `is.finite(target_r2)` branch | the branch returns early, so removing it restores the Iteration 005 path untouched |
| 3 | `R/simulate_phenotypes.R` | `target_r2`, `realised_r2`, `mixture_w` from the stored `truth` list | 3 fields, one site |
| 4 | `R/run_simulation.R` | `target_r2` / `nl_depth` pass-through | 2 sites |
| 5 | `scripts/hpc/generate_params_grid_iter006.R` | delete the file | separate generator; the standing one is untouched |
| 6 | `R/wrapper_polyfun_est_heldout.R` | delete the file, and its entry in the method registry | held-out-prior variant, added to measure the circularity |
| 7 | `scripts/hpc/run_benchmark_job.R` | the ITER-006 method set and the `target_r2` / `nl_depth` grid plumbing | |
| 8 | `scripts/analysis/iter006_*.R` | delete all of them | a fork of the analysis path, as Iteration 005 did |

## What NOT to revert

These were written during Iteration 005/006 but are **not** iteration-specific
and should stay:

- `scripts/hpc/check_grid_columns.R` and its call in the submitter — asserts a
  generated grid carries every column the worker reads. Iteration 005 lost a
  240-task array to a missing `enrichment_values` column.
- `scripts/hpc/check_pkg_current.R`, its call in the submitter, and the
  `unknown_args` guard in `run_benchmark_job.R` — the worker prefers the
  INSTALLED package over the source tree, so a `git pull` touching `R/` does not
  reach the compute nodes until the package is reinstalled. Iteration 005 lost a
  second 240-task array to exactly this.
- The `FMB_GRID_GENERATOR` override in `submit_benchmark_pbs.sh` — lets an
  alternative generator be swapped in without weakening the unconditional
  regeneration that prevents stale-grid reuse.

## Why Iteration 006 exists (so the revert is an informed decision)

Iteration 005 asked whether the LassoNet prior in Functional BEATRICE earns its
complexity when the annotation → causality relationship is not the log-linear
form every competitor assumes. It could not answer, for two reasons found only
afterwards:

1. **The design factor was never controlled.** Linear-representability — the
   share of the true log-weight vector a log-linear-in-A prior can explain, and
   therefore the ceiling on every competitor — was sampled at 1.00 (x4),
   0.33–0.63 (x5) and 0.001 (x1). The single lowest cell was the only one where
   the LassoNet won.
2. **The binary `cooccur` arm was not one relationship.** It used a raw product
   of marks, and annotation proportions are drawn `runif(0.01, 0.30)` per
   region, so realised representability ranged roughly 0.05–0.65 *within* the
   arm.

Iteration 006 makes representability an explicit designed factor via a mixture
of a linear part and an orthogonal nonlinear part, with a closed-form weight
`w = sqrt(1-r) / (sqrt(r) + sqrt(1-r))`. The binary nonlinear part is a
**centred** product, which is orthogonal to all main effects at any annotation
proportion (measured R^2 = 0.0002 across proportions 0.05–0.50, against
0.168–0.801 for the raw product). None of this belongs in the standing
simulator: it exists to answer one question.
