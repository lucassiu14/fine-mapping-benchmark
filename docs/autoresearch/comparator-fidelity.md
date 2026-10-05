# Comparator fidelity: polyfun_ldsc and sbayesrc rewritten to the papers

Branch `fix/comparator-fidelity`, 2026-09-16.

## Why

`polyfun_ldsc` and `sbayesrc` were presented in the LSR as PolyFun and SBayesRC,
but neither followed its paper. A step-by-step comparison found:

**polyfun_ldsc**
- It used unregularised weighted NNLS on all other regions. It had no L2 penalty,
  no odd/even split, and no 20-bin re-estimation.
- It used the S-LDSC regression intercept as a baseline per-SNP heritability,
  which is an error: the intercept is confounding, not heritability.
- Its regression weights and LD scores were not LDSC's.

**sbayesrc**
- It used full LD where the paper uses a low-rank model, and a fixed residual variance.
- It used fixed absolute slab variances instead of γ_k σ²_g.
- It used a softmax link instead of the stick-breaking probit.
- It fitted the annotation model once, from heuristic pilot labels, instead of
  sampling it in every iteration.
- It ran 300 iterations with 150 burn-in; the paper uses 3,000 with 1,000 burn-in.

## What changed

- **`R/wrapper_polyfun_ldsc.R`**: PolyFun's five steps from the paper's Online Methods,
  followed by SuSiE with L = 10 and the HESS prior variance.
  - Where the paper is silent, it follows the released PolyFun code (MIT).
  - Where the paper and the code disagree, it follows the paper (the user's decision):
    - the step-2 parity;
    - λ is chosen within each parity set.
  - Every step is sourced in the file header.
  - New dependency: `Ckmeans.1d.dp` (Suggests), the package PolyFun itself calls.
- **`R/wrapper_sbayesrc.R`**: SBayesRC from the paper and its Supplementary Note only.
  No SBayesRC code is used (GPL-3), so the package stays MIT (the user's decision).
  - The implementation includes:
    - the low-rank model;
    - the mixture prior γ·σ²_g;
    - the stick-breaking probit annotation model (Albert–Chib);
    - the per-block residual variance;
    - pseudo-validation tuning of ρ;
    - 3,000 iterations with 1,000 burn-in.
  - The file header lists each point the paper leaves open and the choice made.
  - The biggest of those: the flat prior on the probit intercepts is taken as uniform
    on [−8, 8], because a flat prior on the whole line diverges under complete separation.
- **`R/wrapper_polyfun_oracle.R`**: still passes the simulation's true causal probabilities
  (a ceiling, not a PolyFun reimplementation), but now fine-maps through `run_polyfun_ldsc()`,
  so it uses exactly PolyFun's SuSiE settings (L = 10, HESS prior variance). The gap between
  the two methods is then the cost of estimating the prior and nothing else.
- **`scripts/hpc/run_benchmark_job.R`**: `sbayesrc = list()`, so the paper's settings apply.
- **Tests**: sections 14b and 14c of `tests/testthat/test-comprehensive.R` were rewritten.

## Existing results were NOT rerun

Every stored result was produced by the old implementations, at commit `11258e4` or earlier
for both files. That covers iterations 004 to 007, the λ sweep, and the LSR text and figures.
Nothing under `results/` was touched.

## Checks

- **SBayesRC**, on data simulated from its own model:
  - σ²_g was recovered (1.060 against 1.065).
  - With identifiable components, the causal-versus-null probit (μ₂, α₂) was recovered.
  - PIPs were calibrated: variants with PIP > 0.9 had a non-zero effect 97% of the time.
  - Script: `scratchpad/fidelity/check_sbrc_model_ident2.R`.
- **PolyFun**, end to end on a simulated scenario:
  - It ran with a stand-in k-median, because `Ckmeans.1d.dp` could not be compiled locally
    until the Xcode licence is accepted.
  - HESS returned zero heritability for some regions. PolyFun errors there, so those fits fail.

## SparsePro and CARMA now use annotations

Both methods can use annotations, but the package ran each region on its own
without them. When a scenario carries annotations, `run_methods()` now runs each
method's own annotation workflow on all regions of the scenario together, through
a `run_<method>_scenario_setup()` hook as PAINTOR does. Without annotations,
nothing changes.

**SparsePro** (`R/wrapper_sparsepro.R`). The workflow is two-step and follows
the upstream README and `sparsepro_zld.py`:

1. One call with every region listed in `--zld`, plus `--anno anno --pthres 1e-5`.
   The script:
   - fine-maps each region under a uniform prior;
   - pools the PIPs across the regions;
   - runs a G-test per annotation;
   - writes joint weights `W` for the annotations with p < pthres.
2. A second call with `--zld <step1>.h2 --anno anno --aW <step1>.W1e-05`. This
   fine-maps each region again under the prior softmax(A W).

Details:
- If no annotation passes the G-test, no weight file is written and the step-1
  result is kept. This matches the paper ("SparsePro+ was not implemented").
- `prior = "all"` uses `W1.0` instead.
- The enrichment estimator is a relative risk of 0/1 annotation status, so
  non-binary annotations are not used. The reason is recorded in
  `additional$annotation_note`.
- SparsePro is CC BY-NC-ND 4.0, so the script is called and never copied.

**Bug fix in the same file.** SparsePro joins the variants of a credible set with
`/`, but the `.cs` parser only split on commas and whitespace. As a result:
- every credible set with more than one variant was dropped;
- `cs_pip` was NA for those sets.

The stored SparsePro credible-set metrics (Iterations 004 to 007) were computed
with that bug. In `results/iter004/native_credible_sets.rds`, every SparsePro set
has exactly one variant (`mean_size` = 1 in all six strata), which is the bug's
signature. Those native coverage, power and size figures cover singleton sets
only.

**CARMA** (`R/wrapper_carma.R`). One `CARMA::CARMA()` call per scenario with all
regions, with `w.list[[r]] = cbind(1, A_r)` as in the CARMA vignette.
- CARMA fits logit Pr(causal) = w'θ inside its EM algorithm. Its M-step pools
  all loci in the call; the paper fits CARMA "at the chromosome level".
- CARMA's defaults are kept: `EM.dist = "Logistic"` and `input.alpha = 0`
  (ridge).
- The other settings are unchanged from the per-region wrapper, including
  `lambda = 1/sqrt(p)`.

**Reproducing earlier runs.** `method_args = list(sparsepro = list(use_annotations
= FALSE))` restores the earlier behaviour; the same applies to `carma`.
`scripts/hpc/run_benchmark_job.R` was not changed, so re-running an annotated
iteration with it would now use the annotations.

**Checks.**
- Tests: sections 10b and 18b of `tests/testthat/test-comprehensive.R`. Stand-ins
  for CARMA and the script check the arguments, file formats, pooling and
  fallbacks; separate tests run the real tools.
- End to end on the real tools: SparsePro in
  `scratchpad/annot/e2e.R` and `e2e2.R`, CARMA in `e2e.R` only.
  - Example: 10 regions, 3 binary annotations, and only the first annotation
    enriched (γ = 30).
  - The G-test selected only that annotation (p = 3.7e-7, W = 1.76).
  - Mean AP was 0.648 with annotations, against 0.516 for the step-1 fit without
    them.

---

## Results: Iteration 008 (2026-10-05)

The rewritten comparators were re-run over Iteration 004's exact grid, into
`$EPHEMERAL/fmbench/results/iter008_fidelity`, with `susie` carried along as a
reproduction control. Tables in `results/iter008_fidelity`.

**The simulations are Iteration 004's.** A canary over rows 1-2 compared 5,000 `susie`
fits: `n_variants` and per-fit AP identical on every one. The Iteration 006 edit to
`simulate_phenotypes.R` is inert for this grid, so the new comparator numbers can be
placed directly beside Iteration 004's for every method that was NOT re-run.

### Average precision by stratum

Mean over cells, 95% interval. Strata are never pooled.

| stratum | method | Iteration 008 | Iteration 004 |
|---|---|---|---|
| sparse, none | polyfun_ldsc / _oracle | 0.708 ± 0.023 | 0.682 ± 0.023 |
| | sbayesrc | 0.542 ± 0.023 | 0.484 ± 0.026 |
| sparse, binary | polyfun_ldsc | 0.737 ± 0.012 | 0.773 ± 0.012 |
| | polyfun_oracle | 0.857 ± 0.009 | 0.843 ± 0.010 |
| | sbayesrc | 0.679 ± 0.010 | 0.509 ± 0.013 |
| sparse, continuous | polyfun_ldsc | 0.724 ± 0.015 | 0.791 ± 0.013 |
| | polyfun_oracle | 0.937 ± 0.006 | 0.911 ± 0.007 |
| | sbayesrc | 0.793 ± 0.009 | 0.619 ± 0.016 |
| sparse_inf, none | polyfun_ldsc / _oracle | 0.536 ± 0.018 | 0.521 ± 0.017 |
| | sbayesrc | 0.451 ± 0.015 | 0.384 ± 0.016 |
| sparse_inf, binary | polyfun_ldsc | 0.535 ± 0.009 | 0.544 ± 0.009 |
| | polyfun_oracle | 0.640 ± 0.008 | 0.631 ± 0.008 |
| | sbayesrc | 0.502 ± 0.009 | 0.381 ± 0.007 |
| sparse_inf, continuous | polyfun_ldsc | 0.538 ± 0.010 | 0.560 ± 0.009 |
| | polyfun_oracle | 0.774 ± 0.007 | 0.753 ± 0.007 |
| | sbayesrc | 0.581 ± 0.009 | 0.464 ± 0.008 |

`sbayesrc` gains everywhere, by 0.06 in the unannotated arm and by up to 0.17 under
annotations. `polyfun_ldsc` is unchanged where there are no annotations to estimate from
and loses 0.01 to 0.07 where there are. `polyfun_oracle` gains slightly: it shares the new
SuSiE settings but is handed the true prior, so only the fine-mapping half changed for it.

### What this does to the headline

With comparators running their published algorithms, and `fb_xregion` as Iteration 004
measured it on the same simulations, paired on fits usable for all three:

| stratum | fb_xregion | polyfun_ldsc | sbayesrc | fb_xregion vs best |
|---|---|---|---|---|
| sparse, binary | 0.7549 | 0.7373 | 0.6867 | +0.018 (SE 0.0026, t 6.7) |
| sparse, continuous | 0.7948 | 0.7228 | 0.8003 | -0.006 (t -0.9) |
| sparse_inf, binary | 0.6113 | 0.5339 | 0.5093 | +0.077 (t 59.6) |
| sparse_inf, continuous | 0.6567 | 0.5379 | 0.5873 | +0.069 (t 17.7) |

Binary annotations under the sparse model are the one annotated stratum the LSR concedes
("under binary annotations polyfun_ldsc leads, 0.773 against 0.748"). It no longer does.
The sparse continuous stratum stays a tie but against `sbayesrc` rather than
`polyfun_ldsc`; by enrichment fold, `fb_xregion` wins at 2.7 (+0.056) and `sbayesrc` edges
it at 5.4, 8.1 and 10.8 by 0.015 to 0.024, none of it significant.

`fb_xregion`'s distance to the oracle prior is now shorter than `polyfun_ldsc`'s in every
annotated stratum; under binary annotations it was the longer of the two before (0.095
against 0.070, now 0.102 against 0.120).

### Read these with three qualifications

1. **Paper-faithful is not software-verified.** Neither comparator has been checked against
   the released program. `sbayesrc` cannot be: it is distributed for genome-wide analysis
   with pre-computed low-rank LD from a fixed panel, which does not apply to simulated
   regions. `polyfun_ldsc` could be, and that is the one test that would settle the drop
   below beyond argument.
2. **`polyfun_ldsc` got worse, which favours our own method**, so it deserves scepticism.
   Two things argue it is real: `polyfun_oracle`, which shares the identical new
   fine-mapping path, got slightly BETTER, so the loss is entirely in prior estimation; and
   the paper's prior estimation pays two costs the old ad-hoc NNLS did not - the odd/even
   parity split halves the regions available to each fit, and the 20-bin Ckmedian step
   discretises what was a continuous per-SNP prior.
3. **New failures.** The modified HESS returns zero heritability on some loci and PolyFun
   errors there, so `polyfun_ldsc` and `polyfun_oracle` lose about 4.5% of fits where
   Iteration 004 lost none. They concentrate at SMALL regions, not large: ~5.7% at p=500,
   7.4% at p=1000, 0.6% at p=2000, because the estimator uses the SNPs below the 0.005
   quantile of the locus p-values and so gets 2-3 of them in a small region against 10 in a
   large one. Any comparison with Iteration 004 must be paired on surviving fits.

### Also stale now

`validity_checks.txt` reports `[FAIL] polyfun_ldsc == susie on the none arm`. That equality
held for Iteration 004's plain SuSiE and is false for PolyFun's L = 10 with the HESS prior
variance. The invariant that holds now is `polyfun_ldsc == polyfun_oracle` on that arm:
586 of 625 cells exactly, maximum difference 0.0025, the residue because the modified HESS
averages over 100 random maximal independent sets.

The LSR's p26 implementation check - "Both PolyFun variants and susie coincide exactly at
0.682, as they must" - is superseded for the same reason. The two PolyFun variants now
coincide at 0.708 and `susie` sits at 0.682.
