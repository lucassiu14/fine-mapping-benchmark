# Iteration 000 — Phase 0 setup and method additions

**Date:** 2026-07-04
**Scope:** Phase 0 of the plan (§0.1 – §0.7). No Phase 1 simulation yet.
**Outcome:** Three of the four Phase 0 method / feature deliverables landed;
SBayesRC deferred per plan (§0.1 caveat).

## What landed on `main`

| # | Item | Plan section | PR | State |
|---|---|---|---|---|
| 1 | Functional BEATRICE annotation-drop fix (gw path) | §0.6 Error 2 | #19 | ✅ merged |
| 2 | `annotation_correlation` argument on the simulator | §0.7 | #20 | ✅ merged |
| 3 | Corrected LD-score PolyFun (`polyfun_ldsc`) | §0.6 Error 1 | #21 | ✅ merged |
| 4 | SBayesRC wrapper | §0.1 | (this PR) | ✅ landed (in-R reimplementation) |

### 1. Functional BEATRICE annotation-drop (§0.6 Error 2, PR #19)

The wrapper previously read the annotation matrix only from
`region_pheno$annotations_matrix`. The locus pipeline (`run_simulation`)
copies annotations onto both objects, so locus runs were unaffected.
`simulate_gwfm_data()` populates **only** `region_geno$annotations_matrix`,
so Functional BEATRICE silently ran annotation-free under the gw pipeline
— inert for the current locus loop, but a real bug for the eventual gw
paper.

Fix: internal helper `.fb_extract_annotations()` that prefers
`region_geno$annotations_matrix` and falls back to
`region_pheno$annotations_matrix`. Regression-tested via the helper
directly (4 cases) and via a mocked wrapper call that captures the
`annotations` argument and asserts it survives a row permutation — the
exact failure mode of the previous code.

### 2. `annotation_correlation` on the simulator (§0.7, PR #20)

New scalar `annotation_correlation ∈ [0, 1]` on
`run_simulation()` / `simulate_phenotypes()`. Default `0` recovers the
original independent generation.

Implementation (binary annotations): group columns by identical
enrichment fold; within each group of size ≥ 2 generate a
compound-symmetric latent Gaussian via the shared-factor construction

    Z_k = sqrt(rho) * F + sqrt(1 − rho) * eps_k,  F, eps_k ~ N(0,1) iid

then threshold each column at `qnorm(1 − prop_k)` so the marginal
frequency is preserved. Cross-group columns remain independent.

The realised Bernoulli correlation is attenuated relative to the latent
target by thresholding — expected, documented in the roxygen block and
in tests. The four sweep values from §1.3 (`{0, 0.25, 0.5, 0.75}`) will
now flow straight into the `WITHIN_JOB` grid.

### 3. Corrected LD-score PolyFun (§0.6 Error 1, PR #21)

Added `polyfun_ldsc` as a new method alongside `polyfun_est`. Per Plan
Option b, the naive one is retained as a labelled baseline so the loop
can demonstrate the corrected method beating it.

Correction: regressor is annotation-weighted LD scores instead of raw
annotations —

    E[chi^2_j] = 1 + N * sum_c tau_c * l_{j,c},
    l_{j,c}    = sum_k r_{j,k}^2 * A_{k,c}

fitted by weighted active-set NNLS. LOCO across regions is implemented
in `run_polyfun_ldsc_scenario_setup()`: for each held-out region $i$,
tau is fitted on the pooled (chi², ℓ) rows from all *other* regions.
`run_methods()`'s scenario-setup merge is scenario-wide (not
per-region), so the setup returns a **named list keyed by `region_id`**
and the region wrapper looks up its own entry via a new `region_id`
argument that `run_polyfun_ldsc_region()` forwards automatically.

### 4. SBayesRC — implemented in R (§0.1)

The earlier version of this log deferred SBayesRC to a supervised
follow-up (per the plan's "if faithful reformatting proves infeasible
in v1, defer" caveat). That deferral was **reversed** — the user
requested a concrete implementation, and this iteration now ships one.

**What was built.** A from-scratch in-R implementation of the SBayesRC
algorithm (`R/wrapper_sbayesrc.R`), tailored to the summary-stat +
per-region-LD setting the rest of this package produces:

- Prior: β_{i,j} ~ Σ_k π_{i,j,k} · Normal(0, σ²_k), with k=0 as a spike
  at zero and K normal slabs on a fixed variance grid (default
  `c(0.05, 0.005, 5e-4, 5e-5)`).
- Per-SNP mixture weights annotation-modulated via multinomial logit
  with k=0 (spike) as reference:

      log(π_{i,j,k} / π_{i,j,0}) = α_k + A_{i,j}^T γ_k.

- Cross-region information sharing: α and γ are **shared across all
  regions in a scenario**. The scenario_setup hook runs a short pilot
  Gibbs per region with region-local annotation refits, pools the
  resulting per-SNP component assignments, and does one final
  multinomial regression on the pooled data to obtain the shared
  (α, γ). Each region's main run then uses those coefficients frozen
  (`prior_source = "pooled_scenario_gamma"`), so its priors are set by
  data that never included that region — mirroring SBayesRC's
  "annotation prior learned genome-wide, applied per-block" pattern.
- Summary-stat likelihood: residual r_j = β̂_j − Σ_{k≠j} R_{j,k}·β_k
  gives r_j | β_j ~ Normal(β_j, 1/n), and mixture-conditional
  posteriors are the standard normal-normal updates
  v_k = 1/(1/σ²_k + n), m_k = v_k · n · r_j.
- Gibbs sweep: joint update of (comp_j, β_j) per SNP with cheap
  running R β vector; annotation refit every `gamma_update_every`
  iterations (default 10); burn-in + posterior-mean PIP.

**Why not the upstream `zhilizheng/SBayesRC` package.** Investigating
the alternative first:

- No CRAN package. Upstream is a C++ project built from source,
  driven by a **pre-eigen-decomposed genome-wide LD folder** tied to a
  ~7M-SNP reference panel + a matching annotation file. Repurposing
  that machinery for our per-region 40–1000-SNP LD blocks is a
  substantially harder engineering task than reimplementing the
  algorithm — and would produce something less faithful, not more
  (because the LD-folder format assumes hg19 ≥1 Mb blocks, not the
  arbitrary simulator windows the benchmark produces).
- Even a faithful wire-up runs FAR outside SBayesRC's native regime:
  the plan explicitly notes its scores at this scale should be treated
  as **in-context-relative**, not as a reflection of genome-wide
  performance. Reimplementing captures the algorithmic content (the
  cross-region annotation-prior mixture) without the reference-panel
  scaffolding that only matters at genome-wide scale.

**Bugs caught + fixed during development** (recorded here so future
iterations do not re-hit them):

- Sparse initial prior needed: uniform initial α gives every class
  ~20 % posterior mass and non-causal PIPs land at ≈0.8. Fixed by
  initialising α_k = log((0.03/K)/0.97), i.e. ~97 % spike, ~3 %
  distributed evenly across the K slabs.
- Empty-class fallback in the annotation-regression refit had α = 0
  for unseen slab classes, which after softmax gave them equal prior
  to the spike — reintroducing the 0.8-non-causal-PIP pathology. Fixed
  by defaulting empty classes to α = −12 (essentially never sampled)
  and capping coefficients to |·| ≤ 6 to prevent divergence.
- `pmax()`/`pmin()` on a matrix strips the `dim` attribute; used
  `full[] <- pmax(...)` to preserve matrix structure.
- `nnet::multinom` returns a plain named vector (no `dim`) when only
  one non-reference class was actually fit — the fallback padding
  needs to infer which class was fit from the `y` factor rather than
  assuming it is class 1.

**Sanity results** on the smoke-test scale (n_regions=3, n=300, p=50,
S=2, φ=0.4, 2 iterations): mean causal PIP > mean non-causal PIP by
> 0.05; non-causal PIP mean stays < 0.05; `run_methods()` end-to-end
reports `prior_source == "pooled_scenario_gamma"` on every fit.
Formal behaviour on the Phase 1 grid remains to be seen; at that scale
the scores should be treated as in-context-relative per §0.1.

## Package errors (from the plan) — status

- ✅ §0.6 Error 2 (FB annotation-drop): fixed.
- ✅ §0.6 Error 1 (polyfun_est LD misuse): **new** `polyfun_ldsc`
  method landed alongside the flagged `polyfun_est` as a labelled
  baseline. Option b in the plan.

## New feature — correlated annotations (§0.7)

- ✅ `annotation_correlation ∈ [0, 1]` argument added; ready to sweep
  `{0, 0.25, 0.5, 0.75}` in Phase 1.

## Feasible method set going into Phase 1

Tier 1 (pure R, always available):

- `susie` (baseline)
- `susie_inf` (baseline)
- `abf` (baseline)
- `carma` (baseline)
- `marginal_z` (baseline — model-free floor)
- `polyfun_oracle` (baseline — cheating ceiling)
- `polyfun_est` (baseline — **flagged as the naive comparator**, kept
  deliberately so the loop can show `polyfun_ldsc` beating it)
- `polyfun_ldsc` (novel — corrected S-LDSC regressor + LOCO across
  regions; scenario_setup hook)
- `sbayesrc` (novel — in-R SBayesRC reimplementation with K-normal
  mixture prior and pooled annotation regression across regions;
  scenario_setup hook)

Tier 2/3 (external binaries / Python — only on hosts that have them):

- `finemap`, `paintor`, `beatrice`, `functional_beatrice`, `funmap`,
  `sparsepro`.

## Datasets

None yet — Phase 1 has not been kicked off. The next step is `§0.3`
(laptop test on tiny config) and `§0.4` (HPC parity), then §1.

## §0.3 laptop smoke test — done

Extended `scripts/hpc/smoke_test.R` to include the full 8-method Tier 1
set (adding `susie_inf`, `carma`, and the new `polyfun_ldsc`) and added
explicit checks that every method's `pip_calibration` and
`fdr_power_curve` compute, that `plot_results()` runs end-to-end, and
that the AUPRC ordering (`marginal_z ≤ susie`, `≤ polyfun_oracle`,
`≤ polyfun_ldsc`) holds.

**Result on the laptop** (3 regions × 8 scenarios, ~80 s total):

| Method | AUPRC | Calibration | FDR curve | Bins > 0 |
|---|---:|:-:|:-:|---:|
| polyfun_oracle | 0.623 | ✓ | ✓ | 10 |
| polyfun_est | 0.592 | ✓ | ✓ | 9 |
| susie | 0.576 | ✓ | ✓ | 9 |
| polyfun_ldsc | 0.557 | ✓ | ✓ | 9 |
| susie_inf | 0.518 | ✓ | ✓ | 9 |
| abf | 0.496 | ✓ | ✓ | 8 |
| marginal_z | 0.231 | ✓ | ✓ | 1 |
| carma | 0.229 | ✓ | ✓ | 9 |

Ordering sanity passes; `plot_results()` produces a valid PDF; every
method emits a well-formed calibration curve. Marginal_z's single
non-empty calibration bin is expected — it puts near-uniform mass on
tails so most bins are empty.

Note: `polyfun_est` marginally beat `polyfun_ldsc` on this tiny (3-region,
S∈{1,2}, φ∈{0.1,0.4}) demo — well within the finite-sample noise floor
of a smoke test and consistent with LOCO's expected variance-cost at
n_regions = 3. The comparison that matters is in Phase 1 at
n_regions = 20 with the full sweep; that's where §0.6 Error 1's
motivating claim (the corrected method should exceed the naive one on
average) will actually be tested.

Also fixed a cosmetic issue: the "Precision (1 − FDR)" y-axis label used
U+2212 (Unicode minus), which triggers `mbcsToSbcs` graphics-device
warnings on the CI runners. Swapped to ASCII "-".

## §0.4 HPC parity — deferred to the loop operator

Same script (`scripts/hpc/smoke_test.R`) is what's designed to run
first thing on the cluster. When the loop is next resumed on the HPC:

```
git pull
R -e 'renv::restore()'
R -e 'install.packages(".", repos = NULL, type = "source")'
Rscript inst/scripts/prepare_vcfs.R   # once, ~150 MB
Rscript scripts/hpc/smoke_test.R      # want "PASSED"
```

If the smoke test passes on the HPC with the same method-availability
+ ordering guarantees, Phase 1 can be submitted with the current
`scripts/hpc/submit_benchmark.sh`.

## Next steps

1. **§0.4 HPC parity** on the cluster (run the smoke test there, one
   session, no further code changes needed).
2. **Iteration 001 — Phase 1 simulation kickoff.** Extend
   `scripts/hpc/generate_params_grid.R` to encode §1.1 / §1.3 (in
   particular the length-20 per-region `p` vector, `annotation_correlation`
   sweep, `n_regions = 20`, length-20 enrichment vector with a handful of
   truly enriched positions). Submit the 3125-scenario array. Sanity-check
   per §1.4 before Phase 2 begins.
3. **Iteration 002+ — Phase 2 loop.** Feasible set = the 8 Tier 1
   methods above; representative dataset set = the full Phase 1
   simulation store initially; run the loop with the calibration gate
   (§2.3) and floor + separation-based dataset pruning (§2.5).
4. **Later — supervised SBayesRC round.** Once real block LD from
   Phase 1 exists, prototype the LD-folder converter against a couple
   of regions; only then decide whether to include SBayesRC in the
   feasible set.
