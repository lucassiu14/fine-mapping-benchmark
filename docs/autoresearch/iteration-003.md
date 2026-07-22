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

## Track B — model changes (STAGED: need local dev + test before HPC)

These require new torch code in `BEATRICE_annot_sparse/` and will be built and
locally validated before they go to the cluster (shipping untested model code has
cost us multi-day HPC runs before). The `finemapper` base loop already accepts a
per-SNP prior `p_0` and there is a `run_<method>_scenario_setup()` hook that pools
across regions (used by sbayesrc / polyfun_ldsc) — that is where these attach.

1. **Pooled genome-wide annotation prior** (user idea #1; flagship). Learn one
   annotation→prior mapping across *all* regions (assuming shared annotation
   effects), then apply it per locus. Cleanest first version: compute a pooled
   per-SNP prior in a `scenario_setup` hook (reuse the polyfun-style pooled
   annotation regression) and feed it as `p_0` to the BEATRICE base loop — no
   change to the torch training. → the quantified "Idea C" target (learned S-LDSC
   captures only 3–24 % of the oracle gap).
2. **Cross-region LassoNet EM loop** (user idea #2). Keep the LassoNet, but share
   it across regions via an outer loop: alternate (a) per-region variational
   inference given the current shared prior, (b) one gradient step on the shared
   LassoNet from all regions' posteriors.
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

```bash
git pull && Rscript -e 'install.packages(".", repos=NULL, type="source")'
export FMB_SCRATCH=$EPHEMERAL/fmbench_iter002          # SAME root -> reuse sim.rds + write *_supp
export FMB_METHODS="fb_l1hi,fb_l1vhi,fb_prreg5,fb_prreg20,fb_ncaus2,fb_concrete,fb_sigma_hi,fb_reg_combo,beatrice_ncaus2,beatrice_sigma_hi"
export FMB_SCENARIOS_PER_TASK=5                        # 10 slow FB-family fits/scenario -> smaller chunk
bash scripts/hpc/submit_benchmark_pbs.sh
```

- The worker reuses `job_*/sim.rds` if present (deterministic seed, so if the
  ephemeral cache was purged it re-simulates identical data — one-off cost).
- Supplemental mode writes `results_supp.rds` / `evaluation_supp.rds`; the 14
  baselines are **not** re-run. `collect_results.R` overlays the new methods.
- chunk=5 → 6,750 tasks (each 5 scenarios × 10 variants × 10 regions ≈ 500
  FB-family fits, well inside 72 h). Fire a `ARRAY_RANGE=1-2` canary first.
- Collect with the same root, then the notebook / split report picks up the new
  method names automatically.
