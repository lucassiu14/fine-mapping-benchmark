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
| 4 | SBayesRC wrapper | §0.1 | — | ⏭ deferred (see below) |

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

### 4. SBayesRC — deferred to a supervised iteration (§0.1)

**Decision:** defer, per the plan's own caveat: *"If faithful
reformatting proves infeasible in v1, defer SBayesRC rather than
feeding it malformed input — note the deferral in the iteration log."*

**Why deferred (checked, not speculative):**

- **No CRAN package.** `available.packages()` returns nothing for
  `SBayesRC` / `sbayesrc`. The upstream is
  [`zhilizheng/SBayesRC`](https://github.com/zhilizheng/SBayesRC), a
  C++ project intended to be built from source and driven via a
  particular LD-folder input format.
- **No GCTB binary present locally** (`which gctb` → not found), so the
  compile-and-install stack would need to be introduced before any
  wrapper work.
- **The LD-folder format is tuned for millions of common SNPs across
  hg19 chromosomal blocks** (≥1 Mb, eigen-decomposed reference).
  Reformatting per-region 40–1000-SNP correlation matrices into that
  format faithfully is genuine engineering — the plan explicitly calls
  this out as *"not a drop-in wrapper"* and *"budget real engineering
  time"*.
- The plan also warns that even when faithfully wired up, SBayesRC on
  our scale runs *far outside its native regime* and its scores must
  be treated as **in-context-relative**, not a reflection of its
  real-world genome-wide performance.

Given the install / reformatting overhead vs. what the score would
mean at this scale, and the plan's explicit permission, this is left
to a **supervised follow-up iteration** rather than shipped now with
malformed input.

**When to revisit:** after Phase 1 simulation data exists and the loop
has a stable feasible set — at that point the value of SBayesRC as a
comparator is clearer and the LD-folder converter can be prototyped
against real block LD from the plan's `n_regions = 20`, `p = 100–1000`
grid.

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

Tier 2/3 (external binaries / Python — only on hosts that have them):

- `finemap`, `paintor`, `beatrice`, `functional_beatrice`, `funmap`,
  `sparsepro`.

## Datasets

None yet — Phase 1 has not been kicked off. The next step is `§0.3`
(laptop test on tiny config) and `§0.4` (HPC parity), then §1.

## Next steps

1. **Iteration 001 — Phase 0.3 / 0.4 smoke test.** Run the whole chain
   (`simulate → run_methods → evaluate_methods → plot_results`) with
   the current feasible set on a tiny config, confirm calibration curve
   + AUPRC + FDR compute, then repeat on the HPC.
2. **Iteration 002 — Phase 1 simulation kickoff.** Full 3125-scenario
   grid on the HPC (§1.1 / §1.3), sanity-checked per §1.4 before Phase
   2 begins.
3. **Later — supervised SBayesRC round.** Once real block LD from
   Phase 1 exists, prototype the LD-folder converter against a couple
   of regions; only then decide whether to include SBayesRC in the
   feasible set.
