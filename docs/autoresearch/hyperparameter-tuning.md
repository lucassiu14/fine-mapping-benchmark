# Hyperparameter tuning for Functional BEATRICE (Optuna)

Replaces the abandoned Iteration 003 "Track A" knob grid (`fb_l1hi`, `fb_prreg5`,
…), which hand-picked ten configurations and evaluated each over the full
benchmark. That was both expensive and unprincipled. This is a proper search:
Optuna proposes configurations, each is scored on a small fixed tuning set, and
the study accumulates across jobs.

## What tuning optimises

Iteration 002 established that Functional BEATRICE's weakness is **calibration,
not ranking** (mass ratio 9–11× and high-PIP reliability 0.10–0.14 under
reference-panel LD). Tuning therefore scores both, using the *same* formulas as
`scripts/hpc/collect_results.R` — the objective is computed by the package's own
`evaluate_methods()`, not a re-implementation:

| metric | meaning |
|---|---|
| `ap` | step-integral average precision (ranking) |
| `max_fdr_violation_n20` | max over thresholds of FDR(t) − (1−t), ≥20 selected |
| `total_mass_ratio` | ΣPIP / #causal (1.0 = calibrated) |
| `hi_pip_reliab` | fraction of PIP ≥ 0.9 calls truly causal |

Two objective modes:

- **`multi` (default)** — multi-objective: maximise `ap`, minimise
  `max_fdr_violation_n20`. Returns a **Pareto front**, so no arbitrary weighting
  between accuracy and calibration is imposed. Optuna does not support pruning
  for multi-objective studies, so every trial runs in full.
- **`scalar`** — `ap − penalty × max_fdr_violation_n20` (default penalty 0.5).
  A single number, which **enables the MedianPruner**: bad configurations are
  abandoned after the first scenario or two. Faster, at the cost of fixing the
  accuracy/calibration trade-off up front.

All four metrics are stored on every trial as user attributes, so either
objective can be re-derived from a finished study without re-running trials.

**One study per stratum.** The never-pool rule applies: a prior tuned on binary
annotations under in-sample LD has no reason to be optimal for continuous
annotations under a reference panel. `STUDY` names the stratum and the tuning sim
is built for that stratum only.

## Search space

Every knob `beatrice_annot.py` actually consumes:

| parameter | range | scale |
|---|---|---|
| `lambda_l1` | 1e-4 – 1.0 | log |
| `prior_regularisation` | 0.1 – 50 | log |
| `sigma_sq` | 5e-3 – 0.5 | log |
| `hierarchy_M` | 1 – 100 | log |
| `n_caus` | 1 – 10 | int |
| `sparse_concrete` | 10 – 300 (step 10) | int |
| `gamma_coverage` | 0.8 – 0.99 (opt-in `--tune-gamma-coverage`) | float |

`max_iter` is deliberately **not** tuned: it trades compute for convergence, so
including it would let the search buy objective value with runtime and make
trials incomparable. Fixed at **2000** (the `run_functional_beatrice()` default).

> Note: the Iteration 002 benchmark ran Functional BEATRICE at `max_iter = 1500`
> (set in `scripts/hpc/run_benchmark_job.R`). Hyperparameters tuned at 2000 are
> optimal *for 2000*, so the validation run should use `max_iter = 2000` too —
> otherwise the tuned configuration is being judged under a different model.

## How it runs continually and survives walltime

- **Resumable.** The study lives in an Optuna `JournalStorage` file. Every worker
  opens it with `load_if_exists=True`, so **re-submitting continues the same
  study** — completed trials are remembered and the sampler keeps its history.
  A worker killed mid-trial leaves a stale `RUNNING` trial, which the sampler
  ignores (only `COMPLETE` trials inform it); no cleanup needed.
- **Walltime-aware.** The submit script computes
  `timeout = walltime − TRIAL_MARGIN` (default 30 min) and passes it to
  `study.optimize()`. A worker stops starting new trials at the timeout and exits
  cleanly, rather than being killed mid-write.
- **Parallel.** `WORKERS` array elements share one study through the storage
  file. `JournalStorage` is the backend Optuna recommends for networked
  filesystems — plain SQLite is prone to lock errors on RDS/GPFS.
- **Efficient.** Each trial is scored on a small fixed tuning simulation, one
  scenario at a time; in `scalar` mode the running mean is reported after each
  scenario and clearly-bad trials are pruned early.

Verified locally before any cluster use: resume (6 → 11 → 13 trials across
re-invocations), the timeout guard (worker self-stopped on schedule), pruning
(4 of 12 trials pruned), and 4 concurrent workers writing 20 trials to one study
with no lost or corrupted writes.

## Running it

One-off: install Optuna into the cluster Python venv. The venv python is
dynamically linked against the module's libpython, so **load the Python module
first** or it fails with `error while loading shared libraries:
libpython3.12.so.1.0`.

```bash
module load Python/3.12.3-GCCcore-13.3.0
~/tools/fmpy-venv/bin/python -m pip install optuna
```

(The submit script loads both the R and Python modules inside the job, so this
is only needed for the interactive install and the `--report` command.)

Build the fixed tuning set for the stratum (once per stratum):

```bash
module load R/4.5.2-gfbf-2025b
Rscript scripts/tuning/make_tuning_sim.R --out tuning/sim_binary_ref500.rds \
    --annotations binary --enrichment 10.8 --n_ref 500 --model sparse --regions 4
```

Submit the workers:

```bash
STUDY=fb_binary_ref500 SIM=tuning/sim_binary_ref500.rds \
WORKERS=10 PBS_WALLTIME=24:00:00 \
bash scripts/tuning/submit_optuna_pbs.sh
```

**To continue after walltime, run that exact command again.** Nothing is
recomputed; the study picks up where it left off.

Check progress at any time (no cluster job needed — reads the storage file):

```bash
~/tools/fmpy-venv/bin/python scripts/tuning/optuna_fb.py \
    --sim tuning/sim_binary_ref500.rds --study fb_binary_ref500 \
    --storage tuning/fb_binary_ref500.journal --report
```

## Files

| file | role |
|---|---|
| `scripts/tuning/make_tuning_sim.R` | builds/caches the fixed per-stratum tuning simulation |
| `scripts/tuning/fb_objective.R` | scores one config on one scenario; audited metrics |
| `scripts/tuning/optuna_fb.py` | Optuna driver: storage, resume, timeout, pruning, search space |
| `scripts/tuning/submit_optuna_pbs.sh` | PBS worker array; re-run to continue |

## After tuning

The Pareto front (or best scalar trial) gives a configuration, not a conclusion.
Validate the chosen config on the **full** Iteration 002 grid as a normal
supplemental method run, so it is judged on held-out scenarios rather than the
small tuning set it was selected on.
