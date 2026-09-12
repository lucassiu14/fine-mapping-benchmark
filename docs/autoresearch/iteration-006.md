# Iteration 006 — the relationships you specify, on continuous annotations

**Status:** ready to submit.
**Scope:** exploratory and TEMPORARY, like Iteration 005. The simulator returns to
the standard design afterwards — see [iteration-006-REVERT.md](iteration-006-REVERT.md).
Iteration 004 remains the standard and is untouched.

> **Not comparable to Iteration 004** (different grid, method set and strength scale),
> and **not directly comparable to Iteration 005** either: same grid and methods, but
> continuous annotations only and a weaker strength target (0.349 against 0.605).

---

## 1. The question

Iteration 005's: how do annotation-aware methods behave when the annotation →
causality relationship is not the log-linear form they assume? This time the
relationships are specified directly, and three of them make a variant's causal
probability depend on its neighbours' annotations.

## 2. Design

**8 rows** = 8 relationships × continuous annotations. Each region carries 10
annotation tracks, of which the first 5 are informative. A relationship turns
the informative tracks of a variant into a score, summed over the five; the
selection weight is exp(c · score), with c set by strength matching (§3).

| relationship | score contributed by each informative track *a* |
|---|---|
| `null` | none — annotations handed to every method, never used in selection |
| `additive` | *a* — the log-linear control, the form every method assumes |
| `wavg2` | weighted average of *a* over the variant and 2 neighbours each side, weights 1/3, 2/3, 1, 2/3, 1/3 |
| `valthresh` | *a* · 1{*a* > 1} (annotations are N(0,1), so the cut-off is 1 SD) |
| `valthresh_wavg2` | `valthresh` applied to each variant first, then `wavg2` |
| `square` | *a*² — symmetric, so strongly negative values enrich too |
| `cubic` | *a*³ |
| `cosine10` | weighted average over ±10 variants, weight e^(−\|d\|/5) · (1 + cos(2πd/10))/2: a variant 10 steps away counts again, as if brought close by the helix, but less than the variant itself; one 5 steps away counts not at all |

Distance is in **variant steps**. A region keeps its variants in genomic order
(`simulate_genotypes.R` sorts the sampled indices) but not their base-pair
positions. Near a region's ends, the weights of the neighbours that exist are
renormalised.

**Grid** — Iteration 005's:

| | |
|---|---|
| S | 1, 3 |
| φ | 0.1, 0.4 |
| iterations | 25 |
| regions | 10, all at p = 1000 |
| n | 5000, in-sample LD, sparse model |

**100 scenarios/row × 8 rows = 800 scenarios.**

**Methods (9)** — Iteration 005's set without `polyfun_est` and `fb_pooled`, which
the LSR leaves out of every figure: `susie`, `beatrice` (annotation-blind
controls); `polyfun_oracle`, `polyfun_ldsc`, `paintor`, `sbayesrc`, `funmap`;
`functional_beatrice`, `fb_xregion`. Selected with `FMB_ITER006_METHODS=1`.
`polyfun_oracle` reads the exact selection probabilities the simulator stores, so
it stays a true ceiling under every relationship.

## 3. Strength

Every non-null arm is scaled so the **top 10% of variants hold 34.9% of the prior
probability** — Iteration 005's lower concentration level, chosen over the 60.5%
it ran at. `select_causal_variants()` looks the target up from the informative
tracks' enrichment value, so the grid's value 2.7 selects it; `enrichment_fold`
indexes that ladder, not a fold.

Without matching the arms would differ mainly in strength. Scoring × log 5.4 with
no scaling, measured on this design:

| | top 10% hold | single most likely variant holds |
|---|---|---|
| additive | 99% | 41% |
| wavg2 | 70% | 7% |
| valthresh | 97% | 44% |
| valthresh_wavg2 | 49% | 4% |
| square | 100% | 84% |
| cubic | 100% | 98% |
| cosine10 | 49% | 3% |

## 4. How much of each relationship any method can learn

Every method's prior maps a variant's **own** annotations to its own prior
probability (checked for the LassoNet in `trainer_annot.py`: each row of the
annotation matrix passes through the network separately). Neighbours'
annotations reach no method. Share of the log selection weight explained, over
20 simulated regions:

| | by a log-linear prior | by any function of the variant's own annotations |
|---|---|---|
| additive | 1.00 | 1.00 |
| wavg2 | 0.48 | 0.49 |
| valthresh | 0.48 | 0.94 |
| valthresh_wavg2 | 0.22 | 0.46 |
| square | 0.01 | 1.00 |
| cubic | 0.61 | 0.99 |
| cosine10 | 0.38 | 0.39 |

The second column is a ceiling for every method, LassoNet included; the first is
the ceiling for the log-linear priors. `square` is the cleanest test of whether a
flexible prior earns its complexity. In the three neighbourhood arms about half
the signal or more is out of every method's reach. Both columns are
scale-invariant, so strength matching does not change them.

## 5. Checks

- `scripts/analysis/test_iter006_relationships.R` — 34 tests, all passing:
  - the weight vectors
  - the neighbour average against a brute-force loop, including both ends of a region
  - every score against an independent brute-force implementation
  - that only the informative tracks enter
  - that Iteration 005's relationships are unchanged
  - that every arm reaches the strength target at both 0.349 and 0.605
  - that `null` is uniform
- End to end through `run_simulation()`, with the arguments the worker builds and
  the bundled VCF at p = 400. Every arm reaches 0.344–0.354 in both regions (`null`
  0.100), and the stored `causal_probs` equal a recomputation from the region's own
  annotations.

## 6. Running it

```bash
cd ~/fine-mapping-benchmark && git pull
module load R/4.5.2-gfbf-2025b GSL/2.8-GCC-14.3.0
R -e 'install.packages(".", repos = NULL, type = "source")'   # R/ changed
bash scripts/hpc/check_toolchain.sh                           # must print RESULT: PASS
Rscript scripts/analysis/test_iter006_relationships.R          # must print 34 passed, 0 failed

FMB_GRID_GENERATOR=$PWD/scripts/hpc/generate_params_grid_iter006.R \
FMB_ITER006_METHODS=1 FMB_SCENARIOS_PER_TASK=5 \
PBS_QUEUE=v1_small72a PBS_WALLTIME=72:00:00 \
FMB_SCRATCH=$EPHEMERAL/fmbench_iter006 \
  bash scripts/hpc/submit_benchmark_pbs.sh
```

Check the echoed design before the array queues. It must show `S={1,3}`,
`phi={0.1,0.4}`, annotations `continuous`, the eight relationships, and
`20 tasks/row x 8 rows = 160 array tasks`.

Iteration 005 planned ~1.5 h per scenario for its 11 methods. Dropping `fb_pooled`
(about 11% of those methods' runtime in Iteration 004) and `polyfun_est` (seconds)
brings that down, so a task of 5 scenarios takes under 7.5 h - well inside 72 h -
and the run under 1,200 CPU-hours.

The job logs print `[iter006] reduced method set: 9 methods` and, because the
relationships run through Iteration 005's code, `[iter005] relationship=...`.

## 7. Afterwards — before ephemeral deletes anything

Results land on `$EPHEMERAL`, which deletes files after 30 days. Once every task
has finished, copy the PIPs and the runtimes/importances into home space:

```bash
BENCH_ROOT=$EPHEMERAL/fmbench_iter006/results/benchmark FLOOR=0 WALLTIME=12:00:00 \
  OUT_DIR=$PWD/results/iter006/piptail bash scripts/hpc/submit_extract_pip_tail.sh
BENCH_ROOT=$EPHEMERAL/fmbench_iter006/results/benchmark SKIP_COUNT=1 \
  OUT_DIR=$PWD/results/iter006/aux bash scripts/hpc/submit_extract_aux.sh
```

The analysis reuses Iteration 005's pipeline, which stratifies by relationship
from the grid. `submit_iter005_analysis.sh` still needs pointing at this
iteration before it is used: its generator path, its output folder, and a
completeness count that walks the whole results tree from the login node.
