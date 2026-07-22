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

## Track A — hyperparameter variants (READY to submit)

Each reuses the *tested* `run_functional_beatrice_region` / `run_beatrice_region`
wrapper with one knob changed (`.FM_REGISTRY` entries + worker `METHOD_ARGS`).
Verified locally: every variant dispatches and forwards its overridden flag to
`beatrice_annot.py` (capture test — e.g. `fb_l1hi` → `--lambda_l1 0.1`). No new
model code, so these are safe to run now.

| method | change vs FB/BEATRICE default | targets |
|---|---|---|
| `fb_l1hi` | `lambda_l1` 0.01 → 0.1 (feature sparsity ×10) | ignore noisy annotations |
| `fb_l1vhi` | `lambda_l1` → 0.5 | strong feature selection |
| `fb_prreg5` | `prior_regularisation` 1 → 5 (shrink prior → uniform) | harm-under-noise |
| `fb_prreg20` | `prior_regularisation` → 20 | strong shrink toward BEATRICE |
| `fb_ncaus2` | `n_caus` 5 → 2 | mass inflation at low S |
| `fb_concrete` | `sparse_concrete` 50 → 200 (more Concrete samples) | softer, less over-confident PIPs |
| `fb_sigma_hi` | `sigma_sq` 0.05 → 0.2 | effect-variance / calibration |
| `fb_reg_combo` | `lambda_l1` 0.1 + `prior_regularisation` 5 + `n_caus` 3 | combined "regularised FB" |
| `beatrice_ncaus2` | BEATRICE `n_caus` 5 → 2 | BEATRICE mass at low S |
| `beatrice_sigma_hi` | BEATRICE `sigma_sq` → 0.2 | BEATRICE calibration |

Note: the 8 `fb_*` variants only differ from each other on the **annotated** arms;
on the `none` arm they reduce to BEATRICE (no annotations → the LassoNet prior is
inert), so their `none`-arm results are expected to duplicate `beatrice`.

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

Register the variants (package change → reinstall on the cluster), then run only
the new methods against the cached Iteration 002 sims:

Per the user's choice (b) — build the model changes first — the Track B joint
methods (`fb_pooled`, `fb_xregion`) go in the **same array** as the Track A
hyperparameter variants:

```bash
git pull && Rscript -e 'install.packages(".", repos=NULL, type="source")'
export FMB_SCRATCH=$EPHEMERAL/fmbench_iter002          # SAME root -> reuse sim.rds + write *_supp
# Track B (model changes) + Track A (hyperparameter screen) in one array:
export FMB_METHODS="fb_pooled,fb_xregion,fb_l1hi,fb_l1vhi,fb_prreg5,fb_prreg20,fb_ncaus2,fb_concrete,fb_sigma_hi,fb_reg_combo,beatrice_ncaus2,beatrice_sigma_hi"
export FMB_SCENARIOS_PER_TASK=5                        # ~500 FB-family fits/scenario-chunk -> smaller chunk
bash scripts/hpc/submit_benchmark_pbs.sh
# To run ONLY the flagship model changes first: FMB_METHODS="fb_pooled,fb_xregion"
```

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
