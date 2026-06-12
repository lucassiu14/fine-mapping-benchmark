# HPC benchmark workflow

A SLURM job array for running the full fine-mapping benchmark. Each task is
one *experimental condition* (model × region size × annotation regime) and
sweeps `S × phi × n_iter` internally. Per-task outputs are stitched into a
single long-form data frame at the end.

## Prerequisites (once per cluster)

1. **Get the code** and restore dependencies:
   ```bash
   git clone https://github.com/lucassiu14/fine-mapping-benchmark.git
   cd fine-mapping-benchmark
   R -e 'renv::restore()'
   ```
2. **Install the package** into the renv library so jobs can `library(fmbenchmark)`:
   ```bash
   R -e 'install.packages(".", repos = NULL, type = "source")'
   ```
3. **Reference data** (one-time, ~150 MB):
   ```bash
   Rscript inst/scripts/prepare_vcfs.R
   ```
   For genome-wide simulation, use `prepare_gwfm_vcfs.R` instead (~400 MB).
4. **(Optional) External methods** — only needed if you're including them in
   `METHODS`:
   - FINEMAP: `Rscript -e 'fmbenchmark::setup_finemap()'` (auto-downloads)
   - PAINTOR: `conda install -c bioconda paintor` (or compile from source)
   - Funmap / BEATRICE / Functional BEATRICE: create the Python env from
     `environment.yml` and pass `python = "/path/to/env/bin/python"` via
     `METHOD_ARGS` in `run_benchmark_job.R`.

## Smoke test (once per cluster, before submitting)

Catches dependency / path issues without burning array time. Finishes in a
minute or two:

```bash
Rscript scripts/hpc/smoke_test.R
```

You want to see `PASSED` and non-zero AUPRC for every Tier-1 method.

## Submit the array

```bash
bash scripts/hpc/submit_benchmark.sh
```

This:
1. Generates `scripts/hpc/params_grid.csv` (one row per job) if absent.
2. Submits a SLURM array with one task per row.
3. Writes logs to `logs/benchmark/slurm-<jobid>_<arrayid>.{out,err}` and
   per-job outputs to `results/benchmark/job_<id>_<label>/`.

**Before submitting**, edit the configuration block at the top of
`submit_benchmark.sh` to set:
- `SLURM_PARTITION` (required — your cluster's partition name)
- `SLURM_ACCOUNT` (if your cluster requires account codes)
- `SLURM_TIME`, `SLURM_MEM` (per-task resource limits)
- `SLURM_CONDA_ENV` (optional — for Tier 3 Python methods)

## Collect results

After the array finishes:

```bash
Rscript scripts/hpc/collect_results.R
```

Writes to `results/benchmark/`:
- `combined_evaluation.rds` / `.csv` — one row per (method × stratum × stratum
  value), joined to its job's parameter row. This is what you load to build
  paper figures.
- `run_summary.csv` — per-job completion status + runtime.

## Customising the experiment

| Want to | Edit | What changes |
|---|---|---|
| Add/remove jobs (different `p`, `model`, annotation) | `generate_params_grid.R` (`AXES`) | New rows in `params_grid.csv` → new SLURM tasks |
| Sweep different `S` / `phi` / `n_iter` per job | `generate_params_grid.R` (`WITHIN_JOB`) | Each task's inner grid |
| Include Tier 2/3 methods | `run_benchmark_job.R` (`METHODS`, `METHOD_ARGS`) | Each task runs the larger method list |
| Change SLURM resources | `submit_benchmark.sh` (config block) | `--time`, `--mem`, `--partition`, etc. |
| LD-mismatch experiment | `generate_params_grid.R` (`ADD_LD_MISMATCH <- TRUE`) | Doubles the grid: matched + mismatched LD per condition |

## File map

| File | Role |
|---|---|
| `generate_params_grid.R` | Produces `params_grid.csv` |
| `params_grid.csv` | One row per SLURM array task (generated, not committed) |
| `run_benchmark_job.R` | Array worker — runs the pipeline for one row |
| `smoke_test.R` | Pre-submit sanity check (tiny params, ~1 min) |
| `submit_benchmark.sh` | SBATCH wrapper that submits the array |
| `collect_results.R` | Stitches per-job outputs into a combined data frame |
