# Iteration 005 — REVERT CHECKLIST

Iteration 005 is a **one-off, small-scale** study of how annotation-aware methods
behave when the annotation→causality relationship departs from the log-linear
form they all assume. **Nothing here belongs in the standing benchmark.** This
file is the complete list of what to undo once the iteration is reported.

Every change is inside an explicit banner:

```
# >>> ITERATION 005 ONLY - TEMPORARY. REMOVE WHEN THAT ITERATION IS DONE. <<<
...
# >>> END ITERATION 005 TEMPORARY BLOCK <<<
```

or carries an inline `ITER-005 (temporary)` / `ITERATION 005 (TEMPORARY)` comment.

Find them all with:

```bash
grep -rn "ITER-005\|ITERATION 005" R/ scripts/ docs/
```

## Checklist

| # | file | what to remove | notes |
|---|---|---|---|
| 1 | `R/simulate_phenotypes.R` | the banner block: `.causal_log_weights()`, `.calibrate_log_weights()`, `.concentration_target()` | self-contained; delete the whole block |
| 2 | `R/simulate_phenotypes.R` | `relationship` and `n_informative` args of `select_causal_variants()`, and the `else` branch that uses them | the `if (identical(relationship,"additive") && is.null(n_informative))` branch IS the original code — keep that, drop the `else` |
| 3 | `R/simulate_phenotypes.R` | `causal_probs` from the returned list and from `genotypes[[i]]$truth` | 3 sites |
| 4 | `R/wrapper_polyfun_oracle.R` | the `causal_probs` argument, the stored-probs branch, and the pass-through in `run_polyfun_oracle_region()` | restores `exp(A' log gamma)` reconstruction as the only path |
| 5 | `R/run_simulation.R` | `relationship` / `n_informative` pass-through | |
| 6 | `scripts/hpc/generate_params_grid_iter005.R` | delete the file | separate generator; the iteration-004 one is untouched |
| 7 | `scripts/hpc/run_benchmark_job.R` | the ITER-005 method-set and relationship blocks | |
| 8 | `scripts/hpc/submit_benchmark_pbs.sh` | the `export FMB_ITER005_METHODS` block inside the PBS heredoc | leave the `FMB_GRID_GENERATOR` override in place — see below |
| 9 | `scripts/analysis/iter004_lib.R` | the generalised `variance_components()` — **KEEP THIS ONE** | see below |

## What NOT to revert

**The `FMB_GRID_GENERATOR` override in `submit_benchmark_pbs.sh` stays.** It is
not Iteration-005-specific. The submitter regenerates `params_grid.csv`
unconditionally, which is the safety that stops a stale grid being silently
reused; the override lets an alternative *generator* be swapped in without
weakening that. Only the `FMB_ITER005_METHODS` export (row 8) is temporary.

**The generalised `variance_components()` is an improvement, not a hack.** The
original estimator required exactly two regions per size class and recovered
σ²_u by differencing them with a bias correction. The generalisation handles any
number of regions per class and reduces to the original when there are two. It
is strictly more general and better conditioned; keep it.

## Why the calibration exists (do not silently drop it if the arms are reused)

The five relationships produce very different enrichment strengths at the same
nominal fold — measured top-decile probability share at fold 10.8 ran 0.839
(additive/binary) down to 0.269 (nonmono/binary), and on continuous annotations
`cooccur` saturated at 1.000. Without calibration a difference between arms is
partly a difference in strength rather than shape. `enrichment_fold` therefore
no longer denotes a fold in this iteration; it indexes a concentration ladder.

## Verification after reverting

```bash
Rscript scripts/analysis/test_iter004.R      # 37 tests must still pass
grep -rn "ITER-005\|ITERATION 005" R/ scripts/   # must return nothing
```
