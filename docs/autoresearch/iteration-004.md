# Iteration 004 — corrected grid, complete method set, tuned Functional BEATRICE

**Date launched:** 2026-08-06
**PBS array:** `3599396`, 450 tasks × 72h, queue `v1_small72a`
**Scratch root:** `$EPHEMERAL/fmbench_iter004`
**Status:** running.

> **Iteration 004 is not comparable to Iterations 001–003.** The sample size, region
> sizes, heritability grid, `p_causal` grid and LD regime all changed. It replaces them
> as the benchmark rather than extending them. Anything quoting an iteration-003 number
> alongside an iteration-004 number is wrong.

---

## 1. Why this iteration exists

Iteration 003 left four things unresolved, and each is addressed here.

**Three methods had never produced a single usable result.** CARMA was in the registry
from Iteration 001 but was never installed, so it returned 100% NA for three iterations
without ever raising an error. FiniMOM and FINEMAP-inf were written and unit-tested
during the Iteration 003 method audit but never run on the grid.

**PAINTOR's results were invalid.** It was run per-locus, but PAINTOR's EM estimates
annotation enrichment *across all loci supplied in one `-input` file* — that is the
method's entire mechanism. Running it one locus at a time asks it to estimate a
genome-wide enrichment from a single region. It never errored, so this was
mis-configuration rather than failure, and every PAINTOR verdict in Iterations 001–003
must be discarded.

**Functional BEATRICE was running on hand-picked hyperparameters.** Iteration 003's
"Track A" attacked this with a grid of ten hand-chosen knob settings, which was both
expensive and unprincipled; it was abandoned in favour of a proper search.

**The grid had a design flaw.** Regions ran to p = 1000 against n = 1000, so the largest
cells had p = n, where the in-sample LD matrix is rank-deficient. This is a regime the
literature deliberately avoids (BEATRICE: n = 5000, m = 1000; CARMA: n = 10,000,
1,500–4,000 variants). Our hardest cells may have been hard for reasons unrelated to
fine-mapping.

---

## 2. The simulation grid

All grid parameters were chosen by the user; see
[hyperparameter-tuning.md](hyperparameter-tuning.md#grid-design-is-a-user-decision-2026-08-02)
for the standing rule that grid design is not Claude's to set.

| variable | Iteration 004 | was (001–003) | why |
|---|---|---|---|
| **n** | **5,000** | 1,000 | p/n peaked at 1.0, where in-sample LD is rank-deficient; now peaks at 0.4 |
| **region sizes** | **500, 750, 1000, 1500, 2000** (2 regions each) | 100, 200, 400, 500, 1000 | as above; also brings the range toward BEATRICE (1,000) and CARMA (1,500–4,000) |
| **phi** | **0.05, 0.1, 0.2, 0.4, 0.6** | 0.0075, 0.05, 0.1, 0.2, 0.4 | 0.0075 produced 6 high-PIP calls across an entire Iteration 003 stratum — no calibration information. 0.6 reaches BEATRICE's own upper bound (ω² = 0.1–0.8) |
| **p_causal** | **0.2, 0.4, 0.6, 0.8** | 0.5, 0.7, 0.9, 1.0 | extends toward the mostly-infinitesimal end; 1.0 removed (see below) |
| **LD** | **in-sample only** | in-sample, ref500, ref200 | deliberate narrowing of scope (see §6) |
| S | 1, 2, 3, 5, 10 | unchanged | |
| annotations | none; binary; continuous | unchanged | |
| enrichment | 2.7, 5.4, 8.1, 10.8 | unchanged | |
| fixed | 10 annotations, 10 regions/scenario, 10 replicates | unchanged | |

**45 cells × (5 S × 5 phi × 10 iterations) = 11,250 scenarios**, against Iteration 003's
135 cells and 33,750 scenarios.

Per-scenario cost rose 5.5× — Σp² over the ten regions went from 2.92×10⁶ to 1.61×10⁷ —
against a third as many scenarios, so total compute is roughly **1.8× Iteration 003**
before accounting for the three methods that had never run.

### Why `p_causal = 1.0` was removed

`p_causal` is the proportion of genetic variance carried by the sparse component, so the
infinitesimal component receives `(1 - p_causal) × phi`
([R/simulate_phenotypes.R:744-752](../../R/simulate_phenotypes.R)):

```r
g_inf_scaled <- g_inf * sqrt((1 - p_causal) * phi / var_g_inf)
```

At `p_causal = 1` this is exactly zero, so the cell **is** the sparse model and duplicates
the separate `model = "sparse"` arm. Verified further: `model = "sparse"` routes to
`generate_phenotype_sparse()`, which takes no `p_causal` argument at all, so the
`p_causal = 0.5` passed at [run_simulation.R:408](../../R/run_simulation.R) is inert on
that path — it is *not* silently introducing a 50% infinitesimal background.

The two are statistically identical but not bit-identical: the `sparse_inf` path still
draws the infinitesimal effects before scaling them to zero, so it consumes different
random numbers.

---

## 3. Methods — 19, for the first time

| role | methods |
|---|---|
| **Focus** | `functional_beatrice`, `fb_pooled`, `fb_xregion`, with `beatrice` as the controlled pair |
| **Ceiling** | `polyfun_oracle` (receives the simulator's true prior weights) |
| **Annotation-informed** | `funmap`, `paintor`, `polyfun_est`, `polyfun_ldsc`, `sbayesrc` |
| **Summary-statistic baselines** | `susie`, `susie_inf`, `finemap`, `finemap_inf`, `abf`, `marginal_z` |
| **Sparse / robust** | `carma`, `finimom`, `sparsepro` |
| **Dropped** | `caviar` — enumeration is combinatorial in region size; already too slow at p = 150, and this grid runs to p = 2000. Wrapper retained. |

### First appearances

- **`carma`** — installed at last. Its `LinkingTo: RcppGSL` needs the GSL module at build
  time, which is why three iterations of 100% NA went unnoticed. Fidelity fixes applied
  during the audit: `num.causal` 5 → 10 (CARMA's own default, and the grid's max S is 10),
  and a per-fit temp directory for `output.labels`, which previously defaulted to `'.'`
  and wrote four files per locus into the working directory — under a 500-task array that
  is the project root, with a race.
- **`finimom`** — fidelity fix: `standardize = FALSE`. FiniMOM's default `TRUE` rescales
  beta/se by `sqrt(2f(1-f))` to convert *per-allele* effects onto the standardised scale,
  but `simulate_genotypes(standardise = TRUE)` means our `beta_hat` is already per-SD.
  Ranking is unaffected (beta and se scale together, so z is preserved); what was wrong is
  the effect-size scale the prior is calibrated against.
- **`finemap_inf`** — now runs the authors' default `susieinf,finemapinf` pipeline rather
  than `finemapinf` alone, which cold-starts at σ² = 1, τ² = 0. Costs one extra SuSiE-inf
  fit per region. Note the authors scope both SuSiE-inf and FINEMAP-inf to **in-sample
  LD**, so this grid is finally their intended regime.
- **`paintor`** — `run_paintor_scenario_setup()` pools all of a scenario's regions into a
  single `-input` file and one PAINTOR call, z-fingerprint keyed.

### A method that legitimately produces nothing

**`funmap` returns 100% NA on the `none` arm, by construction** — it requires
annotations. Verified against Iteration 003: 0 usable / 45,554 NA on `none`, against
182,500 / 0 on binary and 182,442 / 0 on continuous. SuSiE, needing no annotations, is
fine on all three.

The `none` arm is 5 of 45 rows (11.1%), which is why funmap shows ~12% overall failures
in QC while every other method sits at 0.0%. This is correct behaviour and costs the
analysis nothing, because results are never pooled across annotation type.
`qc_run.R` reports a `%failed(annot)` column excluding `anNone` so this does not raise a
false alarm on every run.

**For the paper:** funmap's row on the no-annotation stratum is *empty*, not *poor*. A
reader must not infer a performance claim from an absence.

### A method that produces nothing *without saying so*

**`finemap_inf` returns an all-NA posterior at phi = 0.6 while reporting no error.**
Found after the fitting array completed: all 11,250 `results.rds` were present, but only
10,861 `evaluation.rds`. All 389 missing evaluations carried one identical message —
`'vec' must be sorted non-decreasingly and not contain NAs`, from `findInterval()`.

Decoding the scenario indices (`expand.grid(S, phi, iter)`, S fastest) put every failure
on phi level 5, the level added this iteration. The mechanism, confirmed by direct
inspection of `results.rds`:

| arm | S=1, phi=0.6 | S=10, phi=0.6 | `err_fits` |
|---|---|---|---|
| sparse | 100% NA | works | **0 — silent** |
| sparse_inf | works | 100% NA | **10 — declared** |

`finemap_inf` breaks at phi = 0.6 in *both* models, at opposite ends of S, and the two
failures report themselves differently. `evaluate_methods()` skips fits carrying an
`error` but has no defence against a fit that silently returns all-NA, so **the 389
evaluation failures are exactly the silent ones**. In the sparse arm that is S ∈ {1,2,3,5}
across all 10 iterations plus S=10 in 3 of them — 43 per row × 9 rows ≈ 387, plus
sparse_inf stragglers. The declared failures crashed nothing and were never a problem.

Two consequences.

**No compute was lost.** `results.rds` holds the raw PIPs and is written *before*
evaluation, which runs inside a `tryCatch`. The Iteration 004 analysis re-derives every
metric from `results.rds` + `sim.rds` and never reads `evaluation.rds`, so the 389 gaps
cost nothing.

**The collect had the same latent fault.** Its guard was `all(is.na(pip))`, which returns
`FALSE` for a *partially* non-finite vector. Verified against the real functions:
`.compute_ap_exact()` returns NA via a length-recycled comparison with only a warning, and
`.rank_set_metrics()` throws. These NAs turned out to be whole-vector, so the old code
would have survived this dataset — but the guard is now `!all(is.finite(pip))`, catching
`Inf` from `exp()` overflow and the `NaN` from `Inf/sum(Inf)`, with a `tryCatch` so one
pathological fit cannot kill a 250-scenario task.

Such a fit is **marked failed, not imputed**. Its total mass and the ranking of the
affected variants are undefined, so no metric derived from it is trustworthy; imputing the
non-finite entries to 0 would invent posterior mass and would flatter precisely the method
that overflowed. A new `bad_pip` column, scoped to fits that did *not* report an error,
makes the loss countable per method × phi and keeps a silent failure distinguishable from
an honest one like funmap's.

**For the paper:** this is a robustness result about `finemap_inf`, not a reason to drop
phi = 0.6 — the level is valid for the other 17 methods. But `finemap_inf`'s variance
decomposition must not simply skip the level where it collapses: dropping phi = 0.6 from
its own phi main effect removes the level it fails hardest on and makes it look *more*
phi-robust than it is. Report the failure rate as a result alongside the metric.

---

## 4. Functional BEATRICE hyperparameters

Adopted from Optuna study `fb_insample_v2`, **trial #401**, selected from 962 completed
trials. Applied to `functional_beatrice`; `fb_pooled` and `fb_xregion` inherit.

```
n_caus = 10          sigma_sq = 0.0899        sparse_concrete = 51
max_iter = 2000      prior_regularisation = 0.1255
lambda_l1 = 0.1223   hierarchy_M = 10.15      (both held fixed during the search)
```

On the tuning simulation this gives **AP 0.7702, mass ratio 1.73, high-PIP reliability
1.00, FDR@0.95 = 0.0000, ece_hi 0.244** — against the Pareto front's best-AP point at
AP 0.8109 with mass ratio **8.62** and reliability 0.90. We take a 5% AP loss for a 5×
reduction in posterior-mass inflation and zero incorrect high-PIP calls.

`max_iter` was raised 1500 → 2000 because the study ran every trial at 2000; running the
adopted values at 1500 would judge the configuration under a different model. **BEATRICE
uses the same `max_iter`, `n_caus`, `sigma_sq` and `sparse_concrete`**, so the only
difference between the two is the prior — BEATRICE fixes p₀ = 1/p uniform, FB learns
p₀ = 𝒢(v; ψ). Running FB tuned and BEATRICE at defaults would confound the annotation
prior with the tuning and make the ablation uninterpretable.

> **Caveat to state in the paper:** those values were selected against *FB's* objective,
> so they may be suboptimal for BEATRICE and could understate it. The alternative —
> BEATRICE at its published defaults — confounds the prior with the tuning, which is the
> worse of the two problems for an ablation.

The full tuning history, including why study v1 failed and what the `sparse_concrete`
sweep showed, is in [hyperparameter-tuning.md](hyperparameter-tuning.md).

---

## 5. Two bugs found at submission, both silent

Both were caught by the `DESIGN:` echo added to the submit script, and both would have
invalidated the run without erroring.

**1. The stale grid.** `submit_benchmark_pbs.sh` regenerated `params_grid.csv` only when
the file was *missing*. The CSV is untracked but persists on disk between iterations, so
the first Iteration 004 submission silently reused Iteration 003's file — 135 rows,
n = 1000, regions 100–1000, three LD levels. An entire array went to the queue against
the wrong simulation design. The grid is now regenerated on **every** submit
(deterministic, sub-second), and the design is echoed before `qsub`:

```
DESIGN: n=5000   regions={500,500,750,750,1000,1000,1500,1500,2000,2000}
        S={1,2,3,5,10}   phi={0.05,0.1,0.2,0.4,0.6}
        LD: in-sample
        p_causal: 0.2, 0.4, 0.6, 0.8   annotations: binary, continuous, none
```

**2. The incomplete method list.** The default `METHODS` vector held 15 of the 19
registered methods, omitting `fb_pooled`, `fb_xregion`, `finimom` and `finemap_inf`. The
first two are the cross-region joint-prior models — this project's headline contribution —
and would have had to be requested by hand through `FMB_METHODS` to run at all. A default
that silently drops the thing being tested.

---

## 6. Scope — what Iteration 004 cannot support

**No claim about reference-panel LD.** Iteration 003 found the in-sample/reference split
to be the factor separating methods most: under ref500, `fb_xregion` put 59% of its
excess PIP mass in the 0.5–0.9 band (SuSiE: 56% *below* 0.1), with 0.678 nominal against
0.128 observed. Iteration 004 removes that axis entirely. Every claim must be stated as
in-sample, as SuSiE-inf's authors do for the same reason.

**No comparison against Iterations 001–003.** Four grid axes changed.

**Calibration claims must be stratified by S.** Iteration 003's in-sample `fb_xregion`
mass ratio of 1.00 is an artefact of averaging 3.66 at S = 1 against 0.51 at S = 10. An
aggregate near 1.0 can conceal over-confidence at low S cancelling under-confidence at
high S. Never report a pooled mass ratio.

---

## 7. Running and monitoring

```bash
FMB_SCRATCH=$EPHEMERAL/fmbench_iter004 bash scripts/hpc/submit_benchmark_pbs.sh
```

450 tasks (45 rows × 10), 25 scenarios per task, 72h, 32gb.

**Tasks may exceed the walltime, and that is fine.** Checkpointing is per-scenario: a
walltime kill costs at most the single scenario in flight, and the resume check skips
completed scenarios *before* paying the simulation load. Re-run the identical submit
command to continue. Plan on two or three passes.

**Redundant simulation on the first pass.** All 10 tasks of a row start together, find no
cached `sim.rds`, and each simulate the same panel. Writes are atomic (tmp + rename) so
this is safe, just wasteful — roughly 900 core-hours, about 10% of the run. Avoidable by
submitting `ARRAY_RANGE=1-441:10` first (one task per row) to build the caches, then the
full range.

### QC

```bash
Rscript scripts/hpc/qc_run.R $EPHEMERAL/fmbench_iter004/results/benchmark --sample 300
```

Reports per-method fits and failures, and flags three cases by name: **MISSING ENTIRELY**,
**100% FAILED** (the CARMA mode — registered, running, returning nothing usable), and
**PARTIAL FAILURES >5%**. It also reports completion, per-row progress, scenarios that
wrote `results.rds` but no `evaluation.rds`, and unreadable files.

`--method NAME` breaks one method's failures down by region size with the distinct error
messages — a failure rate concentrated at p = 2000 *biases* what survives rather than
merely thinning it.

**First QC pass (2,608 of 11,250 scenarios):** all 45 rows present, all 19 methods
producing fits, 18 at exactly 0.0% failure, funmap at 12.3% raw / 0.0% annotated
(structural, §3). No log contained Error, Killed or walltime. `results.rds` and
`evaluation.rds` counts were identical — no scenario lost its evaluation.

---

## 8. Next — the analysis

The analysis is specified in
[variable-importance-analysis.md](variable-importance-analysis.md), which was written
anticipating exactly this run: it assumes in-sample LD only and quotes 45 jobs, 11,250
scenarios and 5,625 cells per method, all of which Iteration 004 delivers.

**Two things must happen before it can be executed.**

1. **Its factor levels are stale.** §2.1 and §2.3 still list the Iteration 003 grid
   (`p_causal` 0.5–1.0, `phi` from 0.0075, `region_size` 100–1000, `n = 1000`). The
   *structure* is unchanged, which is why the cell counts still hold, but every level
   value and every worked example needs updating.

2. **None of its §5 code prerequisites are implemented.** `.compute_ap_exact`,
   `.rank_set_metrics`, `sum_pip_sel`, `sum_pip_sq`, `pip_cal_breaks`, per-fit AP, and
   the L1/L2 artefacts (`combined_fit_metrics.rds`, `combined_replicate_metrics.rds`,
   `region_idx`) are all absent. Without `region_idx` in particular, σ²_u is
   unidentifiable and `p_causal` and `enrichment_fold` have **no valid error term at all**.

> **⚠ Deadline.** The §5 changes operate on stored PIPs and require no refitting — but
> `results.rds` lives in `$EPHEMERAL`, which auto-purges. Apply the changes and re-collect
> before that happens, or the raw PIPs are gone and the analysis becomes un-re-derivable.

Iteration 004 also resolves three of the §12 data-validity constraints in that spec:
`carma` is installed, `finimom` and `finemap_inf` have run, and `paintor` is corrected.
Those rows should be updated once the run completes.

---

## 9. Reproducing this iteration

| what | where |
|---|---|
| grid definition | [scripts/hpc/generate_params_grid.R](../../scripts/hpc/generate_params_grid.R) |
| method args, FB hyperparameters | [scripts/hpc/run_benchmark_job.R](../../scripts/hpc/run_benchmark_job.R) |
| submission | [scripts/hpc/submit_benchmark_pbs.sh](../../scripts/hpc/submit_benchmark_pbs.sh) |
| QC | [scripts/hpc/qc_run.R](../../scripts/hpc/qc_run.R) |
| tuning study + history | [hyperparameter-tuning.md](hyperparameter-tuning.md), `tuning/fb_insample_v2.journal` |
| tuning sweep tooling | [scripts/tuning/sweep_from_study.py](../../scripts/tuning/sweep_from_study.py) |
| method install recipes | [adding-methods.md](adding-methods.md) |
| wrapper audit | [iteration-003.md](iteration-003.md) |

---

## 11. Outcome — the run that produced the reported results

**Status: complete.** 45 rows, 11,250 scenarios, 112,500 fits per method, 18
methods, in-sample LD. All validity checks pass; all 18 analysed methods have
exactly 5,625 populated cells; declared-failure and silent-NA rates are 0.0%
throughout.

### What it took to get here

The arm covering rows 10–45 (`sparse_inf`) had to be fitted **three times**.

1. **First attempt** — the toolchain lived behind `~/tools`, a symlink onto
   `$EPHEMERAL`. Imperial's rolling 30-day purge ate it *mid-run*: sourcing the
   virtualenv's `activate` kept succeeding while imports inside it stopped
   resolving, so seven Python-backed methods recorded 38–50% failed fits for
   hours without the job stopping. Discarded.
2. **Second attempt** — toolchain rebuilt on project storage, but the package
   list was written from memory. BEATRICE also needs `imageio`, `seaborn` and
   `tqdm`; FINEMAP-inf needs `bgzip`. Every fit died on `ModuleNotFoundError`
   after a check that had passed. Discarded.
3. **Third attempt** — dependencies derived from the tool sources by AST walk
   and actually imported rather than merely located; `setuptools` pinned below
   81 (81+ removed `pkg_resources`); `fine-mapping-inf`'s C extension built with
   `--no-build-isolation`. All 1,800 tasks reported `preflight ok`.

Each failure is now a guard: `check_toolchain.sh`, the per-task preflight, the
`bad_pip` column, and `verify_and_analyse.sh`'s gates. Full detail in
`docs/reports/iter004_process_record.pdf` §11.

### Headline results

| | sparse / cont | sparse_inf / cont |
|---|---|---|
| `polyfun_est` | **0.8061** (2nd) | 0.5480, FDR@.9 0.005 → **0.35** |
| `fb_xregion` | 0.7933 (3rd) | **0.6521 (2nd)**, FDR@.9 0.049 → 0.068 |

**A PolyFun-style per-SNP prior collapses under misspecification; the
jointly-fitted variational prior does not.** `fb_xregion` rises from third to
second of eighteen, behind only the oracle.

**The ablation holds in all four annotated strata** (t = 17–82) and is *exactly*
zero in both unannotated ones over 6,250 paired replicates.

**Region sharing improves calibration as well as accuracy.** Adding the
annotation prior costs calibration (`beatrice` → `functional_beatrice` takes
FDR@0.9 from 0.020 to 0.100 on continuous annotations); sharing that prior
across regions partly repairs it (`fb_xregion` 0.049).

### Caveats that must travel with these numbers

- **Four methods are reimplementations**, not the published software:
  `polyfun_{oracle,est,ldsc}` and `sbayesrc`. Findings about them are findings
  about the *modelling strategy*.
- **`sbayesrc`'s calibration metrics measure a different quantity.** Its PIP
  counts assignment to any non-spike component and its smallest slab has
  variance 5e-5, so inconsequential effects count fully. Ranking metrics remain
  valid; calibration ones are excluded from comparison.
- **The infinitesimal arm scores the sparse component only.** Non-causal
  variants carry real effects that count as false discoveries, so absolute FDR
  is not comparable across arms — degradation between arms is.
- **Nominal region sizes are approximate above 1000** (the 2000 class has median
  1590 actual variants), so the swept range is 3.2× not 4×.
- **AP is not monotone in phi** — it dips at 0.6 for `sbayesrc`, `sparsepro` and
  `finimom` while `susie_inf` improves. Real and method-specific.

### Deliverables

| file | contents |
|---|---|
| `docs/reports/iter004_full_analysis.pdf` | 28pp, 17 figures, all six strata |
| `docs/reports/iter004_process_record.pdf` | 9pp, simulation → analysis, with the incident log |
| `docs/reports/iter004_sparse_arm_analysis.pdf` | 24pp, the sparse arm alone (superseded) |
