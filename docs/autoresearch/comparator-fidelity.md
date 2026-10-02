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
