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

**The `multi` front collapsed on `fb_binary_ref500` (final: 3102 trials,
2026-08-02).** The Pareto front is a *single distinct parameter set*, re-found by
the sampler over and over — the top five trials by AP (#2024, #2143, #2174,
#2236, #2241) all carry AP = 0.8091 to four decimals. Every front point has
`viol = 0.0000`.

The reason is **not** that violation is degenerate: 51% of trials (1577/3102)
have nonzero violation, so it is a genuinely varying quantity. It is that
violation never *conflicts* with AP — configurations that rank well also control
FDR well, so the frontier has no trade-off to trace and the multi-objective
formulation returns one point.

That costs something. `multi` disables the MedianPruner, so all 3102 trials ran
all four scenarios where `scalar` would have abandoned most after the first —
several times the cost per trial for a second objective that produced a
single-point front. **Prefer `scalar` unless the front actually spreads.**
`--report` prints the zero-violation fraction and the front, so this is
checkable early rather than after 72h.

### AP is saturated; calibration is the remaining headroom

v1 went from 1586 → 3102 trials and moved AP from **0.8089 → 0.8091**. 1516
additional trials bought +0.0002, with the sampler re-finding one point. AP on
this stratum is exhausted.

Meanwhile **`total_mass_ratio` sat at ~2.0 in every top trial** and is in neither
objective mode. Total PIP mass is twice the true causal count. Since Iteration
003's central Functional BEATRICE result was mass inflation (3.06 → 2.70 via the
joint prior), an AP-only search returns a configuration that wins on ranking
while staying globally over-confident — and, now that AP has saturated, returns
*nothing else*. `reliab = 1.00` and `viol = 0` locate the problem: the top of the
ranking is clean, so the excess mass is in the tail where neither objective
looks.

`--mass-penalty W` (scalar mode) subtracts `W × |total_mass_ratio − 1|` from the
score. It defaults to **0**, reproducing the previous objective exactly. `|·|` is
symmetric because both directions are miscalibration: above 1 is over-confident
(Iteration 002's 9–11× failure mode), below 1 leaves detectable mass unused. The
in-loop pruning score uses the identical formula, so the MedianPruner ranks
trials on the quantity actually being optimised.

### Importance-driven pruning of the search space

`--report`'s fANOVA importances for AP on v1's 3102 trials:

| parameter | importance |
|---|---|
| `sparse_concrete` | 0.617 |
| `sigma_sq` | 0.340 |
| `prior_regularisation` | 0.031 |
| `n_caus` | 0.008 |
| `hierarchy_M` | 0.002 |
| `lambda_l1` | 0.001 |

The two parameters that were **pinned at their bounds** are the two that carry
**96% of the variance in AP** — direct confirmation that v1's bounds, not its
objective, were the binding constraint, and that its 3102 trials were largely
spent on knobs explaining 4%.

`hierarchy_M` (0.002) is confirmed inert, as its geometric-centre optimum
predicted. `lambda_l1` (0.001) is inert too — its interior optimum was flat-
direction noise, not a located optimum.

`--fix NAME=VALUE ...` removes parameters from the search space and holds them
fixed, so budget concentrates on what moves the objective:

```bash
FIX="hierarchy_M=10.15 lambda_l1=0.1223" ...
```

> **Do not use AP importance to fix `n_caus`.** Its AP importance is 0.008, but
> it is the *prior number of causal variants* — the knob most likely to control
> `total_mass_ratio`, which is the metric v2 exists to improve and which the AP
> importance analysis cannot see. Low importance for one objective is not
> evidence of irrelevance to another. `prior_regularisation` (0.031) is small but
> non-zero and also stays in.

**One study per stratum.** The never-pool rule applies: a prior tuned on binary
annotations under in-sample LD has no reason to be optimal for continuous
annotations under a reference panel. `STUDY` names the stratum and the tuning sim
is built for that stratum only.


## Grid design is a user decision (2026-08-02)

The v1 grid (S ∈ {1,3}, phi ∈ {0.05,0.2}, p ≤ 500, n = 1000, 4 regions) was chosen by
Claude for trial cost, without checking that the reduction preserved the failure being
tuned against. Measured on Iteration 003 in its own stratum, it contained **9% of the
high-PIP calls and 4% of the high-PIP errors** — the study reported `hi_pip_reliab = 1.00`
where the benchmark measures 0.52, and 3102 trials moved AP by 0.0002.

**Do not set, thin or extend a grid axis without agreement, including "just for tuning".**
`submit_optuna_pbs.sh` now aborts when `SCENARIOS != |S| x |phi|`, because a trial scores
scenarios `1..SCENARIOS` and a mismatch silently drops cells.

### Agreed v2 grid — study `fb_insample_v2`

| axis | value | rationale |
|---|---|---|
| S | 1, 2, 3, 5, 10 | the dominant axis in-sample: `ece_hi` 0.189 → 0.033 and mass ratio 3.66 → 0.51 across it |
| phi | 0.05, 0.1, 0.2, 0.4 | 0.0075 dropped (4 high-PIP calls, pure noise); `ece_hi` is near-flat 0.082–0.097 over the rest |
| n | 5000 | matches BEATRICE's own paper; the old n = 1000 put p = 1000 regions at p = n, where in-sample LD is rank deficient |
| region sizes | 500, 1000, 2000 | one per scenario, **rotating**; cost is quadratic in p and p = 2000 alone is 76% of a 3-region trial |
| LD | in-sample | the regime FB is competitive in; ref500 deferred |
| model | sparse | sparse_inf out of scope by decision |

Rotation does not alias: `run_simulation` builds scenarios with `expand.grid(S, phi, iter)`
so S cycles with period 5 while the rotation has period 3, and gcd(5,3) = 1.

### Field practice, for reference

| study | causals/locus | per-locus h² | n | variants/region | LD |
|---|---|---|---|---|---|
| BEATRICE | 1, 4, 8, 12 | ω² = 0.1–0.8 | 5,000 | 1,000 | in-sample |
| CARMA | — | — | 10,000 | 1,500–4,000 | reference |
| SuSiE 2.0 | — | — | 1,000 | 5,000 | in-sample |
| ours (v2) | 1, 2, 3, 5, 10 | 0.05–0.4 | 5,000 | 500–2,000 | in-sample |

Open gap, deferred to Iteration 004: our phi tops out at 0.4 where BEATRICE's own paper
goes to 0.8, so we test our base method below the signal range its authors used.

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
