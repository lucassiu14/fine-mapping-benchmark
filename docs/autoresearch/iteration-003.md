# Iteration 003 — improving (Functional) BEATRICE

**Motivation (from the Iteration 002 analysis).** Functional BEATRICE's problem is
*calibration*, not ranking: it is badly over-confident (broken in 52–73 % of
annotated scenarios; high-PIP reliability down to 0.29 on continuous; total PIP
mass ratio 2.8–4.6), and its per-locus annotation learning *helps* only on binary
annotations with strong enrichment under in-sample LD, and *harms* elsewhere. The
root cause is that the `LassoNetPrior` is fit from a single locus — an underpowered
estimation problem that fits noise. Plain BEATRICE is also over-confident (26–39 %
broken, mass ratio ~2.1). So the ideas below target **over-confidence** and
**per-locus annotation over-fitting**.

All variants are evaluated on the **same Iteration 002 grid** (135 rows × 250
scenarios), reusing the cached row `sim.rds` files, and compared against the 14
baselines already collected.

---

## Track A — hyperparameter variants (DROPPED 2026-08 — do not resubmit)

Originally 10 methods (`fb_l1hi`, `fb_l1vhi`, `fb_prreg5`, `fb_prreg20`, `fb_ncaus2`,
`fb_concrete`, `fb_sigma_hi`, `fb_reg_combo`, `beatrice_ncaus2`, `beatrice_sigma_hi`)
that reused the FB/BEATRICE wrappers with a single knob changed each. **Abandoned:**
a hand-picked knob grid is the wrong tool for tuning — a principled hyperparameter
search (e.g. Optuna) in a dedicated study is the right approach, and these were
also ~4/5 of the per-scenario compute cost for little expected gain. The methods
have been removed from `.FM_REGISTRY` and the worker `METHOD_ARGS`; any partial
Track A results already on disk (rows ~1–81 of the interrupted supplemental run)
are simply ignored at collect/analysis time. **Only the Track B joint-prior models
are carried forward.**

## Track B — model changes (BUILT + LOCALLY VALIDATED 2026-07-22; ready to submit)

Status: `fb_pooled` (idea #1, linear head) and `fb_xregion` (idea #2, LassoNet
head) are implemented, registered, and validated end-to-end through `run_methods`.
Ideas #3–#4 (below) remain future work. R-side wiring:
`R/wrapper_fb_joint.R` (scenario_setup hooks + thin per-region lookup, z-fingerprint
keyed), `.FM_REGISTRY` entries `fb_pooled`/`fb_xregion`, and `METHOD_ARGS` in
`run_benchmark_job.R` (FB base args; joint driven by the scenario_setup hook).

Local validation (torch 2.13 venv, 3-region fixture through `run_methods`):
- both methods return per-region PIPs bit-identical to the direct Python trainer,
  `joint_fallback = FALSE` on every region (z-fingerprint cache hits, no silent
  degradation);
- `none`-arm (no annotations) falls back to plain BEATRICE, tagged
  `joint_fallback = TRUE`, no error;
- realistic full-scenario timing (10 regions, sizes 2×{100,200,400,500,1000},
  max_iter 500): **29 s wall, 585 MB peak RSS** → ≈ **87 s/scenario at max_iter 1500**,
  i.e. about one `functional_beatrice`-equivalent per scenario per joint method.
  Well inside 72 h walltime and node memory.

### (superseded design note) Track B — model changes

These require new torch code in `BEATRICE_annot_sparse/` and will be built and
locally validated before they go to the cluster (shipping untested model code has
cost us multi-day HPC runs before). The `finemapper` base loop already accepts a
per-SNP prior `p_0` and there is a `run_<method>_scenario_setup()` hook that pools
across regions (used by sbayesrc / polyfun_ldsc) — that is where these attach.

### Design decision (2026-07-22): cross-region prior training is JOINT, not sequential

**What the fork does today (rejected as the target).** The existing
`--prior_weights` / `--return_weights` scaffolding implements a *sequential
warm-start chain*: train the `LassoNetPrior` on region 1, save its weights, load
them to initialise region 2, keep training, save, region 3, … The prior network is
passed hand-to-hand. This is online / continual learning and is **order-dependent**
(a different region ordering gives a different prior) and suffers **catastrophic
forgetting** (by the last region, region 1's signal has decayed; late regions
dominate). It does not estimate the genome-wide enrichment — it approximates it with
the tail of a walk.

**What we are building instead (the user's proposal).** There is **one** shared
prior network `φ`, optimised against **all** regions **simultaneously**. Each
training step:

- for every region `r`, compute its ELBO using the *same* shared prior
  `p_0^(r) = f_φ(v_r)`;
- update `φ` **once** with the gradient **summed over all regions**;
- update each region's *own* finemapping posterior `ψ_r` (which SNPs are causal in
  region `r` — these stay region-specific and separate).

Objective (two-level / hierarchical):

```
L(φ, {ψ_r}) = Σ_r  ELBO_r( ψ_r ; p_0^(r) = f_φ(v_r) )
```

- `φ` — the shared prior head (annotations → per-SNP prior). ONE copy, sees every
  region's evidence at every step. Order-independent, no forgetting.
- `ψ_r` — region-specific variational finemapping posterior (one per region).
- Only the prior `p_0 = f_φ(v)` is shared across regions; the likelihood
  (`Z_r`, `LD_r`) is per-region.

This is the PolyFun / fgwas / TORUS pooling idea (their EM M-step pools all regions)
done by joint SGD instead of hard EM alternation — cleaner, no E/M split.

**Gradient aggregation over regions: EQUAL WEIGHT PER REGION** (decided 2026-07-22,
user-confirmed). `φ`'s gradient is the *unweighted* sum/mean over regions, i.e. each
locus is one exchangeable draw of the shared enrichment relationship. We do **not**
weight by region size (#SNPs) — size-weighting would let large loci dominate the
enrichment estimate, which contradicts the "annotations consistent across regions"
assumption we are exploiting.

**This unifies ideas #1 and #2.** They differ *only* in the form of the shared prior
head `f_φ`; the joint cross-region training harness is identical for both:

| idea | shared prior head `f_φ` |
|---|---|
| #1 (flagship) | simple pooled map: linear/logistic on annotations (funmap/PolyFun-style global coefficients) |
| #2 | the existing `LassoNetPrior` (neural, feature-sparse) |

Plan: build **one** multi-region training loop (instantiate the prior head once, hold
all regions' `(Z, LD, v)` + per-region posteriors together, equal-weight gradient on
the head), with two swappable heads. The current `prior_weights` chain does NOT give
this — it is the sequential scheme — so the multi-region loop is the genuinely new
code.

### Implementation + local validation (2026-07-22)

New code lives in the fork, leaving the tested single-region path untouched:
- `BEATRICE_annot_sparse/scripts/joint_trainer.py` — `run_joint()` (the joint loop),
  `_region_elbo`/`_abf` (per-region ELBO **copied verbatim** from
  `finemapper_lassonet.train`/`.abf` so joint vs single-region share identical
  likelihood/KL/reg terms — only the optimisation is cross-region), and
  `LinearPrior` (idea #1's shared logistic head).
- `BEATRICE_annot_sparse/beatrice_joint.py` — CLI. Takes a `--manifest` TSV
  (`z,LD,annot,target,N` per region) + `--prior_head lassonet|linear`, trains ONE
  shared head across all regions, writes per-region `pip.csv`/`credible_set.txt`
  via the unchanged `gen_cred`.

Scaling note: the objective is the **mean** over regions `(1/R)Σ_r ELBO_r` (equal
weight per region). Under Adam's per-parameter normalisation the global `1/R` washes
out, so each region posterior `ψ_r` gets the same ~lr per-step update as
single-region FB, and the shared prior `φ` gets the mean of the R per-region
gradients — results are directly comparable to single-region FB.

Validated locally on a 3-region synthetic fixture (venv: torch 2.13, p≈40, 1–2
causals, 500 iters, ~6 s):
- **LassoNet head (idea #2)** localised every causal as the top hit (region PIPs
  0.63 / 0.71+0.67 / 0.90) and — vs the single-locus FB on the same region0 —
  **cut a false positive from 0.58 → 0.12 and mass ratio from 3.06 → 2.70**. This is
  the intended effect: sharing the prior across regions regularises the annotation
  map and reduces the single-locus over-fitting that breaks FB's calibration.
- **Linear head (idea #1)** also localised all causals but was *more* over-confident
  here (Σpip 3.4–3.75) — expected, since with no L1 it uses all 5 noise annotations
  and can over-fit the map on only 3 regions. Real-grid behaviour (10 regions, 1500
  iters, swept enrichment) is the actual test.
- Both shared heads learned the enrichment in the correct direction (mean
  logit-contrast gap enriched−noise: +0.19 LassoNet, +0.14 linear).

1. **Pooled genome-wide annotation prior** (user idea #1; flagship). The joint
   scheme above with a **simple linear/logistic** shared head `f_φ`. Assumes shared
   annotation effects across regions. → the quantified "Idea C" target (learned
   S-LDSC captures only 3–24 % of the oracle gap). A cheaper fallback if the joint
   loop is delayed: precompute a pooled per-SNP prior in a `scenario_setup` hook
   (polyfun-style regression) and feed it as fixed `p_0` to the base loop — but the
   joint version is the real deliverable.
2. **Cross-region LassoNet — JOINT** (user idea #2). The same joint scheme with the
   **`LassoNetPrior`** as the shared head, trained on all regions at once (NOT the
   sequential `prior_weights` chain). Equal weight per region.
3. **Confidence-gated annotation prior** (mine). Scale the LassoNet's
   deviation-from-uniform by estimated annotation informativeness, so FB smoothly
   reduces to BEATRICE when annotations are uninformative (continuous / weak
   enrichment / noisy LD) — removing the "harmful" region of FB's behaviour.
4. **Cross-region consistency regulariser** (mine). Fit per-locus LassoNets but
   penalise their annotation weights for differing across loci — a soft version of
   pooling that operationalises "annotations are consistent for causal SNPs across
   regions".

---

## How to run (supplemental, reuses cached sims)

Run only the two joint-prior models against the cached Iteration 002 sims. This
resumes: scenarios already carrying an `evaluation_supp.rds` (from the interrupted
12-method run) are skipped, so only the not-yet-done scenarios are computed, and
now only for the 2 joint methods (~1/5 the per-scenario cost of the old 12-method
job).

```bash
export FMB_SCRATCH=$EPHEMERAL/fmbench_iter002          # SAME root -> reuse sim.rds + evaluation_supp.rds
export FMB_METHODS="fb_pooled,fb_xregion"
export FMB_SCENARIOS_PER_TASK=25                       # only 2 methods/scenario now -> larger chunk is fine
bash scripts/hpc/submit_benchmark_pbs.sh
```

Reinstalling the package first is optional here — `fb_pooled`/`fb_xregion` are
already registered; a reinstall is only needed to pick up the Track A *removal*
above, which does not affect these two methods.

- The worker reuses `job_*/sim.rds` if present (deterministic seed, so if the
  ephemeral cache was purged it re-simulates identical data — one-off cost).
- Supplemental mode writes `results_supp.rds` / `evaluation_supp.rds`; the 14
  baselines are **not** re-run. `collect_results.R` overlays the new methods.
- **Joint methods (`fb_pooled`/`fb_xregion`) run ONE Python process per scenario**
  over all 10 regions (the scenario_setup hook), not per region — measured ≈ 87 s
  and < 0.6 GB per scenario per method (max_iter 1500), ≈ one `functional_beatrice`
  per scenario. On the `none` arm they fall back to plain BEATRICE (tagged
  `additional$joint_fallback = TRUE`).
- chunk=5 → 6,750 tasks, each ≈ 5 scenarios × 12 methods. Fire a `ARRAY_RANGE=1-2`
  canary first and confirm `fb_pooled`/`fb_xregion` produced non-fallback results on
  the annotated arms before releasing the full array.
- Collect with the same root, then the notebook / split report picks up the new
  method names automatically.

### Analysis note — never-pool still applies

`fb_pooled` / `fb_xregion` are **new methods**, compared against the existing
baselines within each cell. The never-pool rule is unchanged: do NOT aggregate
their metrics across model (sparse vs sparse_inf) or annotation type
(none/binary/continuous). The `joint_fallback` flag lets analysis separate genuine
joint results (annotated arms) from the BEATRICE-equivalent `none`-arm fallback.

---

# RESULTS (collected 2026-07-28, analysed 2026-07-29)

Full report: **`results/iteration-003-joint-prior-analysis.pdf`** (10 pp).

## Provenance — what was run and where it lives

| item | location |
|---|---|
| Array | PBS `3449477[]`, 1350 tasks (chunk=25), `FMB_METHODS="fb_pooled,fb_xregion"` |
| Scratch root | `$EPHEMERAL/fmbench_iter002` (SAME root as Iteration 002 → reused `job_*/sim.rds`) |
| Per-scenario outputs | `job_*/scenario_*/{results_supp.rds, evaluation_supp.rds}` — 33,750/33,750 complete |
| Collect job | PBS `3470549`, `FMB_DUMP_FDR_CURVES=1` (first iteration to dump the true FDR curves) |
| Collected locally | `results/iter003/` — `combined_scenario_metrics.{rds,csv}` (367,025 cells), `combined_pip_calibration.rds`, `combined_fdr_curves.rds`, `combined_evaluation.{rds,csv}`, `run_summary.csv` |
| Analysis scripts | scratchpad (`i3_analyse.R`, `i3_figs*.R`, `i3_tabs.R`, `build_pdf_i3.py`) |

⚠️ Ephemeral auto-purges ~30 days from 2026-07-23. The `results/iter003/` copy is the durable one.

**Run history note.** The first submission (`3391109`, 12 methods incl. Track A) reached ~55 % (18,693/33,750)
over 5 days before Track A was dropped; it was cancelled and re-submitted as the 2-method job above, which
resumed on the existing checkpoints and finished in hours. Track A metrics survive in `evaluation_supp.rds`
for rows ~1–81 and are simply ignored at analysis time (9,890 cells each vs 16,875 for the real methods).

## Validity

- **Full coverage**: 16,875 scenario cells for each of `fb_pooled` / `fb_xregion` — identical to every baseline.
- **Fallback correct**: 100.0 % identical to `functional_beatrice` on the `none` arm (where a joint annotation
  prior is meaningless); differ in >99 % of annotated cells → real cross-region training, not silent degradation.
- NA rates 0.0 % for AUPRC / FDR-violation / mass ratio; 26 % for high-PIP reliability (cells where nothing
  reached the top PIP bin — the usual abstention signal).

## Findings

1. **The diagnosis was right and the fix works, partially.** Sharing the prior across regions removes roughly
   half of FB's excess posterior mass (continuous / n_ref=200: **11.45 → 5.04**; binary: 9.02 → 5.41) and
   improves FDR control at every threshold.
2. **Systematic, not an averaging artefact** — paired within-cell, the joint models improve mass calibration in
   **68–94 %** of cells. `fb_xregion` is above 50 % on essentially every metric × LD × annotation combination.
3. **The LassoNet head beats the linear head (idea #2 > idea #1).** `fb_xregion` wins AUPRC in **62–97 %** of
   cells vs `fb_pooled`. `fb_pooled` improves calibration but **loses ranking on binary annotations** (AUPRC
   better in only 6–23 % of cells) — the linear head is too weak. **Recommend dropping `fb_pooled`.**
4. **No ranking cost for `fb_xregion`** — matches or beats per-locus FB, and is the best deployable method on
   continuous annotations under in-sample LD (AUPRC 0.539).
5. **Better annotation use** — enrichment slope 0.0085/fold vs FB's 0.0070 (binary, in-sample); both → ~0 under
   panel LD while the oracle keeps 0.017–0.022.
6. **NOT sufficient.** Still worse calibrated than plain BEATRICE under panel LD (5–7× vs ~3.5×), and beats
   finemap in only 2–10 % of cells off in-sample LD (<2 % on FDR anywhere). **LD mis-specification, not the
   annotation prior, is the binding constraint** — per-locus over-fitting was a real cause, not the only one.
7. **Residual inflation is concentrated at low S** (S=1: fb_xregion 6.37, FB 9.96, finemap 2.23), implicating
   the fixed causal-count prior `n_caus=5` rather than the annotation prior → directly motivates the Optuna
   search over `n_caus` / `sigma_sq` / `prior_regularisation`.

## Next steps

1. Tune `n_caus` + effect-variance prior via the Optuna harness (running) — targets finding 7.
2. Adopt `fb_xregion` as the Functional BEATRICE default; drop `fb_pooled`.
3. Re-prioritise ideas #3 (confidence-gated prior) and #4 (cross-region consistency regulariser) *behind*
   attacking LD mis-specification, which binds for every annotation-using method.

---

# METHOD AUDIT + NEW METHODS (2026-07-29/30)

## PAINTOR was being run incorrectly — its Iterations 001–003 numbers are invalid

PAINTOR's EM estimates annotation enrichment across **all** loci in its `-input`
file; that is the method's central mechanism and how it is run on real data. Our
wrapper ran it **one locus at a time**, and on the annotated arms it *did* receive
annotations, so it was estimating enrichment from a single locus in every fit.

It never errored — it produced numbers — so this is not a crash but a
mis-configuration that handicapped the method. **Every PAINTOR result and verdict
in the Iteration 001–003 reports should be treated as invalid** (FDR 0.35
in-sample, "worst in the field", its variance-decomposition row, its entries in
the guidance table). Group-level claims averaging over the six
functionally-informed methods are mildly affected too.

Ironically this is the same pathology we diagnosed in functional BEATRICE —
learning an annotation→prior map from one locus. The difference: FB does it *by
design* (so that finding stands), whereas we *imposed* it on PAINTOR.

Fixed by `run_paintor_scenario_setup()`: all regions → one `input.txt` → one
PAINTOR call → per-region results keyed by z-fingerprint, with per-locus
fallback and an `additional$pooled_across_loci` tag.

## Wrapper audit — all 16 checked

| finding | status |
|---|---|
| **polyfun_oracle** π reconstruction | **verified identical** to the simulator's causal-selection probability (`exp(A %*% log(enrichment))`, normalised). Our ceiling is trustworthy — every "captured fraction" claim rests on this. |
| **ABF** | matches Wakefield's closed form to 8e-17 |
| **FINEMAP** | real `beta_hat`/`se`/`maf` passed; rsid matching for PIPs and credible sets |
| **SparsePro** | `.cs` holds rsids, correctly mapped via `match()` |
| **susie / susie_inf** | `n` passed to `susie_rss`; `unmappable_effects = "inf"` with a version guard |
| **polyfun_est/ldsc** | `.tau_to_prior_weights` normalises to sum 1 — correct for susieR's prior *probability* semantics |
| **funmap** | **BUG FIXED** — read annotations only from `region_pheno`, which drops them under `simulate_gwfm_data`. Latent only; iters 001–003 used the locus pipeline, so results unaffected. |
| **sbayesrc** | sampler core **verified correct** (log-BF algebra, residual update, `beta_hat = z/sqrt(n)` scale). But see below. |

### sbayesrc: a latent 15× failure, now guarded

sbayesrc has two annotation paths. Measured on a benchmark-like locus (p=300,
n=1000, φ=0.05, S=3), same data and seed:

| path | mass ratio |
|---|---|
| `pooled_gamma` from `scenario_setup` (intended) | **3.6–3.8** |
| in-loop refit fallback | **55** |

The fallback inflates posterior mass ~15× **even under the true LD**, so the
inflation is caused by the annotation refit, not LD. `scenario_setup` returns
`list()` in four situations and pooled_gamma then silently became NULL. Now warns.

**Iterations 002/003 took the good path** — confirmed by reproducing the observed
in-sample mass ratio (~3–6). Those results are unaffected.

## Four methods added; three kept

Selected against the deep-dive finding that **LD mis-specification is the binding
constraint**.

| method | status |
|---|---|
| **CARMA** | installed and **PASSES** (sum PIP 2.03 for S=2, both causals at 1.00). Registered since Iteration 001 but **never installed** → 100% NA for three iterations. Cause: `LinkingTo: RcppGSL` needs the GSL module at build time. |
| **FiniMOM** | installed, **PASSES** (2.00) |
| **FINEMAP-inf** | installed, **PASSES** (2.53) |
| ~~CAVIAR~~ | **dropped on cost** — enumeration is combinatorial in region size; already too slow at p=150 and the grid runs to p=1000. Wrapper retained. |

### Fidelity fixes — deviations the planted-causal test could NOT catch

All three still recovered the causals; they change the configuration, not the
ability to run. Found by reading each package's source against our call sites.

1. **CARMA wrote four files per locus into the current working directory**
   (`output.labels` defaults to `'.'`). Under the array that is the project root
   with ~500 concurrent tasks sharing one label → repo pollution **and a race**.
   Now a per-fit temp dir.
2. **CARMA `num.causal` was 5**; CARMA's default is 10 and the grid's max S is 10,
   so every S=10 scenario was under-modelled. Set to 10.
3. **FiniMOM was double-standardising.** Its `standardize=TRUE` default rescales
   beta/se by `sqrt(2f(1-f))` to move *per-allele* effects onto the standardised
   scale — but `simulate_genotypes(standardise=TRUE)` means our `beta_hat` is
   already per-SD. Now `standardize=FALSE`. Ranking was unaffected (both beta and
   se scale together, so z is preserved); the *prior's* effect-size scale was
   wrong, which matters most at weak signal. Measured impact on a strong locus:
   0.007 — which is exactly why testing could not find it.
4. **FINEMAP-inf ran a non-default configuration.** `--method finemapinf` alone
   takes the cold-start branch (σ²=1, τ²=0); the authors' default
   `susieinf,finemapinf` has FINEMAP-inf **inherit SuSiE-inf's** τ²/σ². Now runs
   the default pipeline.

Verified clean: CARMA's `outlier.switch = TRUE` (the LD-discrepancy detection we
added it for) is on by default. FiniMOM's `insampleLD` is auto-detected from the
presence of `X_ref`.

**Recorded deviation:** CARMA accepts annotations via `w.list`; we pass none, so
it runs annotation-free even on the annotated arms.

## Not yet run

The canary was deliberately **not** fired — the Optuna tuning has the queue.
Untested at scale: PAINTOR's pooling, CARMA's runtime on p=1000 regions,
FiniMOM's `insampleLD` against a real reference-panel sim, FINEMAP-inf's
credible-set branch, and full `run_methods` → `evaluate` → `collect` integration.

```bash
export FMB_SCRATCH=$EPHEMERAL/fmbench_iter004     # NEW root
export FMB_METHODS="carma,finimom,finemap_inf,paintor"
export FMB_SCENARIOS_PER_TASK=25
ARRAY_RANGE=1-2 bash scripts/hpc/submit_benchmark_pbs.sh
```
