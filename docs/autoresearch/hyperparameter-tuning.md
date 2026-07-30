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

**`multi` degenerated on `fb_binary_ref500` (observed 2026-07-30, 1586 trials).**
The Pareto front collapsed to a *single distinct parameter set* — the two front
entries, trials #1356 and #1566, carry identical parameters (NSGA-II re-evaluated
its own optimum). Every point on the front had `viol = 0.0000`, which implies no
trial anywhere in the study beat AP = 0.8089 at *any* violation level: a trial
that did would be non-dominated and would appear on the front. So the second
objective exerted no trade-off pressure at all and the study was effectively
single-objective on AP already.

That is not free. Because `multi` disables the MedianPruner, all 1586 trials ran
all four scenarios, where `scalar` would have abandoned most of them after the
first — up to ~4× the cost per trial for an objective that returned no
information on this stratum. **Prefer `scalar` on any stratum where violation
saturates at 0.** `--report` now prints the fraction of trials at exactly zero
violation and flags the degenerate case, so this is visible before a long run
rather than after it.

**Mass ratio is deliberately absent from both objective modes**, and the v1
optimum has `mass = 2.00` — total PIP mass twice the true causal count. Since
Iteration 003's central Functional BEATRICE result was mass inflation (3.06 →
2.70 via the joint prior), tuning purely on AP can return a configuration that
wins on ranking while remaining globally over-confident. `reliab = 1.00` and
`viol = 0` say the *top* of the ranking is clean, so the excess mass sits in the
tail, where neither objective looks. `--report` now prints the mass ratio of the
top 5 trials by AP so this stays in view; folding it into the `scalar` penalty is
an open option, not yet taken.

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
| `sigma_sq` | **1e-4** – 0.5 | log |
| `hierarchy_M` | 1 – 100 | log |
| `n_caus` | 1 – 10 | int |
| `sparse_concrete` | **1 – 300** | int, **log** |
| `gamma_coverage` | 0.8 – 0.99 (opt-in `--tune-gamma-coverage`) | float |

**Bounds revision, 2026-07-30** (bolded rows above). Study `fb_binary_ref500`
v1 reached 1586 completed trials, and its optimum sat *on* two of the original
bounds:

| parameter | v1 optimum | v1 floor | position in range |
|---|---|---|---|
| `sigma_sq` | 0.005323 | 5e-3 | 1.4% into a 2-decade log range |
| `sparse_concrete` | 10 | 10 (step 10) | exactly on the floor |

An optimum on a bound means the *bound* chose the value, not the objective, so
both floors were opened: `sigma_sq` to 1e-4, and `sparse_concrete` to the
package's own declared `lower_bound=1` (`beatrice_annot.py`). `sparse_concrete`
also switched from `step=10` to log-uniform, which puts sampling density at the
low end where the optimum lives rather than spreading it evenly over 10–300.
BEATRICE's README examples all use `--sparse_concrete 10`, so v1 did rediscover
the authors' value — but standing on a floor is not evidence that it is optimal.

`lambda_l1`, `prior_regularisation` and `n_caus` all found interior optima in v1
(0.1223, 8.61, 4), so their ranges are doing their job and are unchanged.

**`hierarchy_M` is under suspicion.** Its v1 optimum, 10.15, is almost exactly
the geometric centre of its [1, 100] log range — which is what a parameter the
objective does not respond to looks like. `--report` now prints an fANOVA param
importance for AP; if `hierarchy_M` comes out near zero it should be *fixed*
rather than tuned, since tuning an inert knob spends search budget for no
return. This is a hypothesis to test against the existing trials, not yet a
finding.

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

Then submit — **one command does everything**. If the tuning sim does not exist
yet it is built as its own job first, and the worker array is submitted with
`-W depend=afterok:<build job>` so it starts automatically once the build
succeeds (and never runs against a missing or partial sim). If the sim already
exists the build is skipped and the workers start immediately — which is also
what happens on every resume.

```bash
STUDY=fb_binary_ref500 SIM=tuning/sim_binary_ref500.rds \
WORKERS=20 PBS_WALLTIME=72:00:00 \
bash scripts/tuning/submit_optuna_pbs.sh
```

The stratum used to build the sim is set by `ANNOTATIONS` / `ENRICHMENT` /
`N_REF` / `MODEL` / `P_CAUSAL` (defaults: binary, 10.8, 500, sparse) and must
describe the same stratum `STUDY` names. To build the sim by hand instead:

```bash
module load R/4.5.2-gfbf-2025b
Rscript scripts/tuning/make_tuning_sim.R --out tuning/sim_binary_ref500.rds \
    --annotations binary --enrichment 10.8 --n_ref 500 --model sparse --regions 4
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
