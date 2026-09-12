# Iteration 006 — REVERT CHECKLIST

Iteration 006 is exploratory, like Iteration 005. Once both are reported the
simulator returns to the standard design (Iteration 004's), which neither
iteration changes. This file lists everything Iteration 006 added.

Every change is inside a banner

```
# >>> ITERATION 006 ONLY - TEMPORARY. REMOVE WHEN THAT ITERATION IS DONE. <<<
...
# >>> END ITERATION 006 TEMPORARY BLOCK <<<
```

or tagged `# ITER-006 (temp)`. Find them all with:

```bash
grep -rn "ITERATION 006\|ITER-006" R/ scripts/ docs/
```

## Checklist

| # | file | what to remove |
|---|---|---|
| 1 | `R/simulate_phenotypes.R` | the ITERATION 006 block: `.ITER006_RELATIONSHIPS`, the two weight vectors, `.iter006_neighbour_avg()`, `.iter006_score()` |
| 2 | `R/simulate_phenotypes.R` | the two `# ITER-006 (temp)` lines at the top of `.causal_log_weights()` |
| 3 | `scripts/hpc/generate_params_grid_iter006.R` | delete the file |
| 4 | `scripts/analysis/test_iter006_relationships.R` | delete the file |
| 5 | `scripts/hpc/run_benchmark_job.R` | the ITERATION 006 block (`ITER006_METHODS`) and the `FMB_ITER006_METHODS` selector tagged `# ITER-006 (temp)` |
| 6 | `scripts/hpc/submit_benchmark_pbs.sh` | the `export FMB_ITER006_METHODS` block inside the PBS heredoc |
| 7 | `docs/autoresearch/iteration-006.md` and this file | archive with the report |

## Revert it together with Iteration 005

Row 2's hook sits inside Iteration 005's `.causal_log_weights()`, which
Iteration 005's checklist deletes whole. Reverting 005 alone would take the hook
with it and leave row 1 as dead code; reverting 006 alone leaves Iteration 005's
relationships working exactly as they were. Row 5 also builds on Iteration 005's
`ITER005_METHODS`, so reverting 005 alone would leave the worker referring to a
list that no longer exists.

Apart from its method-set selector (rows 5-6), Iteration 006 runs through
Iteration 005's temporary paths - the grid's `relationship` and `n_informative`
columns and the stored `causal_probs` that `polyfun_oracle` reads - so Iteration
005's checklist removes those. `FMB_GRID_GENERATOR` is permanent (see
iteration-005-REVERT.md).

## After reverting

R/ changes, so reinstall the package on the cluster before the next submission;
`check_pkg_current.R` blocks the submitter until you do.

```bash
grep -rn "ITERATION 006\|ITER-006" R/ scripts/     # must return nothing
Rscript scripts/analysis/test_iter004.R            # must still pass
```
