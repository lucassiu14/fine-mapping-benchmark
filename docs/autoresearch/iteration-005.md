# Iteration 005 — when the annotation→causality relationship is not what the methods assume

**Status:** ready to submit.
**Scope:** deliberately small and one-off. See
[iteration-005-REVERT.md](iteration-005-REVERT.md) for everything to undo afterwards.

> **Not comparable to Iteration 004.** Different grid, different method set, and
> `enrichment_fold` means something different here (see *Calibration* below).

---

## 1. The question

Every annotation-aware method in the benchmark assumes a **linear or log-linear**
map from annotations to the per-SNP prior:

| method | prior model |
|---|---|
| `polyfun_est` | χ² regressed linearly on A |
| `polyfun_ldsc` | linear in annotation LD scores |
| `paintor` | logistic-linear |
| `sbayesrc` | multinomial logit |
| `fb_pooled` | **linear** head |
| `functional_beatrice`, `fb_xregion` | **LassoNet** — hidden layers, *can* be nonlinear |

So relationships a linear model cannot express are the sharpest available test of
whether the LassoNet architecture earns its complexity. If it does,
`fb_xregion` should separate from `fb_pooled` and from the PolyFun-style family
as the relationship departs further from log-linear.

## 2. Design

**12 rows** = 6 relationships × 2 annotation types.

| relationship | form | what it breaks |
|---|---|---|
| `additive` | `λ Σ A_k` | nothing — control, and the bridge to Iteration 004 |
| `cooccur` | `λ Σ_k A_k A_{k+1}` | needs two marks *together*; linear-in-A sees only the marginal leak |
| `nonmono` | `λ exp(−(count−2)²/2)` | peaks at *intermediate* load; a monotone model finds almost no slope |
| `mixed` | `λ(A₁+A₂+A₃ − A₄−A₅)` | two annotations **deplete**; `polyfun_est` clamps negative coefficients to zero and `polyfun_ldsc` constrains non-negativity, so neither can express it |
| `threshold` | `λ·1{count ≥ 3}` | a step; the mildest departure, partly capturable by a slope |
| `null` | — | annotations independent of causality, but still handed to every method |

**10 annotations per region, of which only the first 5 are informative.** The
other 5 are independent noise, present in the data every method receives but
never entering selection. Annotation *selection* is therefore a constant
challenge in every arm, and only the functional form varies.

**Grid** — depth over breadth, because Iteration 004 found only 4.7% of cells had
a statistically decidable winner at 10 iterations:

| | |
|---|---|
| S | 1, 3 |
| φ | 0.1, 0.4 |
| iterations | **25** (2.5× Iteration 004) |
| regions | 10, all at p = 1000 |
| n | 5000, in-sample LD |
| enrichment | 5.4 (fixed — see below) |

**100 scenarios/row × 12 rows = 1,200 scenarios.**

**Methods (11):** `susie`, `beatrice` (annotation-blind controls);
`polyfun_oracle`, `polyfun_est`, `polyfun_ldsc`, `paintor`, `sbayesrc`, `funmap`
(annotation-aware); `functional_beatrice`, `fb_pooled`, `fb_xregion` (the focus).

## 3. Calibration — and why `enrichment_fold` no longer means a fold

The five relationships produce **very different enrichment strengths at the same
nominal fold**. Measured over 30 replicate regions at p = 1000, top-decile
probability share at fold 10.8:

| | binary | continuous |
|---|---|---|
| additive | 0.839 | 0.999 |
| cooccur | 0.793 | **1.000** |
| nonmono | 0.269 | 0.333 |
| mixed | 0.723 | 0.999 |
| threshold | 0.285 | 0.517 |

Continuous annotations **saturate** — `cooccur` puts all selection probability on
a tenth of the variants, making causal location near-deterministic — and
`nonmono`/`threshold` are three times weaker than `additive`. Left uncorrected, a
difference between arms would be partly a difference in *strength* rather than in
*shape*, which is the only thing this iteration measures.

`.calibrate_log_weights()` bisects a scalar multiplier so every relationship hits
a common top-decile target, taken from what `additive`/binary produces at that
fold. **Verified: all twelve arms land within 0.001 of target.** One structural
exception — `nonmono`/binary saturates near 0.69, because binary counts admit
only six distinct weights — which is why **fold 5.4 is used**: it is the value at
which all twelve arms match (0.585–0.605).

Consequence: `enrichment_fold` indexes a concentration ladder, not a fold. Since
every row uses 5.4, this is invisible within the iteration, but it matters if the
arms are ever reused with varying strength.

## 4. Two things that had to be fixed first

**`polyfun_oracle` would have stopped being an oracle.** It reconstructed
`π_j ∝ exp(A′log γ)`, correct only for the additive form. Correlation between
that reconstruction and the true selection probabilities under the new
relationships: **0.297** (cooccur), **0.093** (nonmono), **0.064** (mixed),
**0.438** (threshold). The simulator now stores the exact probabilities on
`truth` and the oracle reads them — verified to match to machine precision
(max |diff| ≤ 2.2e-15).

**`variance_components()` would have returned NULL.** It assumed exactly two
region draws per size class; Iteration 005 uses ten regions of one size. It is
now generalised — at three or more draws the between-draw variance is estimated
directly, which is better conditioned, and at exactly two the original path runs
unchanged. A hardcoded `N_ITER = 10` was also found and fixed: at 25 iterations
it over-subtracted the sampling term and biased σ²_u low by 23%. Verified
unbiased over 200 replicates (0.0402 against a true 0.0400, z = +0.64).
**This fix should be kept, not reverted.**

## 5. Running it

```bash
module load R/4.5.2-gfbf-2025b
cd ~/fine-mapping-benchmark
bash scripts/hpc/check_toolchain.sh            # must print RESULT: PASS

Rscript scripts/hpc/generate_params_grid_iter005.R \
  scripts/hpc/params_grid.csv

FMB_ITER005_METHODS=1 FMB_SCENARIOS_PER_TASK=5 \
  PBS_QUEUE=v1_small72a PBS_WALLTIME=72:00:00 \
  FMB_SCRATCH=$EPHEMERAL/fmbench_iter005 \
  bash scripts/hpc/submit_benchmark_pbs.sh
```

240 tasks (12 rows × 20). At ~1.5 h per scenario for the reduced method set at
p = 1000, each task is roughly 7–8 h — comfortably inside walltime.

**Note the submitter regenerates `params_grid.csv` from Iteration 004's
generator on every run.** Either point `RSCRIPT`-generated output at the
iteration-005 generator first (as above, overwriting the same filename), or run
with `FMB_SKIP_GRID=1` if that switch exists. Check the echoed DESIGN block
before the array queues — it must show `S={1,3}` and `phi={0.1,0.4}`.

## 6. What to look for

The prediction, stated in advance so it can be wrong: **`fb_xregion` and
`functional_beatrice` should hold up across `cooccur` and `nonmono` where
`fb_pooled`, `polyfun_est`, `polyfun_ldsc` and `paintor` degrade**, because only
the LassoNet head can represent a non-linear map. On `mixed`, the PolyFun pair
should fail for a different and more specific reason — they clamp away the
depleting half by construction. On `additive` all annotation-aware methods should
be comparable, and on `null` all of them should collapse to their
annotation-blind counterparts.

If `fb_xregion` does *not* separate on `cooccur`/`nonmono`, that is evidence the
LassoNet's capacity is not being used, and the simpler linear head would be the
better default.
