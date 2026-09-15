# Iteration 007 — REVERT CHECKLIST

Iteration 007 is exploratory, like Iterations 005 and 006. It changes nothing in
`R/` and nothing in Iteration 004's grid. It adds a separate grid generator and a
method-set selector, and both go once the iteration is reported. This file lists
everything Iteration 007 added.

Every change is inside a banner

```
# >>> ITERATION 007 ONLY - TEMPORARY (see iteration-007-REVERT.md) <<<
...
# >>> END ITERATION 007 TEMPORARY BLOCK <<<
```

or is a file named for the iteration. Find them all with:

```bash
grep -rn "ITERATION 007\|iter007" R/ scripts/ docs/
```

## Checklist

| # | file | what to remove |
|---|---|---|
| 1 | `scripts/hpc/generate_params_grid_iter007.R` | delete the file |
| 2 | `scripts/hpc/run_benchmark_job.R` | the ITERATION 007 block: `ITER007_METHODS` and the `FMB_ITER007_METHODS` selector |
| 3 | `scripts/hpc/submit_benchmark_pbs.sh` | the `export FMB_ITER007_METHODS` block inside the PBS heredoc |
| 4 | `scripts/hpc/check_iter007.sh` | delete the file |
| 5 | `scripts/analysis/submit_iter007_collect.sh` | delete the file |
| 6 | `docs/autoresearch/iteration-007.md` and this file | archive with the report |

Rows for the analysis scripts are added when those are written.

## Independent of Iterations 005 and 006

The Iteration 007 block defines its own method list rather than building on
`ITER005_METHODS`. Its grid has no `relationship` column, so it runs through the
standard simulator. It can be reverted before, after or together with 005 and 006.

## After reverting

No `R/` file changes, so no reinstall is needed.

```bash
grep -rn "ITERATION 007\|iter007" R/ scripts/    # must return nothing
```
