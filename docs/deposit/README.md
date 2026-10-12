# Main simulation: fine-mapping benchmark results

Results of the main simulation study from `fine-mapping-benchmark`, in which
nineteen statistical fine-mapping methods were run on identical simulated data
and scored with identical metrics.

Code: https://github.com/lucassiu14/fine-mapping-benchmark

Every file here is an R serialised object (`.rds`), readable with
`readRDS("path")`. The two entry points are also provided as CSV.

---

## Start here

| If you want | Open |
|---|---|
| One row per method per simulation cell, with every metric | `evaluation-metrics/per-cell.csv` (or `.rds`) |
| What each parameter row of the grid contains | `simulation-grid.csv` |
| What every column means | `DATA-DICTIONARY.md` |

---

## The simulation

Genotypes were simulated from 1000 Genomes Phase 3 haplotypes with `sim1000G`,
retaining variants with minor allele frequency of at least 0.01, and each region
was standardised so that the in-sample linkage disequilibrium matrix is a
correlation matrix. Every method received the same marginal summary statistics,
the same LD matrix and, where applicable, the same annotation matrix.

**The grid.** 45 parameter rows, each simulated 250 times, giving 11,250
scenarios. A scenario holds ten independent regions, so each method was fitted
to 112,500 regions; across nineteen methods the study contains 2,031,230 fits.

| Factor | Levels |
|---|---|
| Causal variants per region, `S` | 1, 2, 3, 5, 10 |
| Variance explained, `phi` (h²) | 0.05, 0.1, 0.2, 0.4, 0.6 |
| Variants per region, `region_size` | 500, 750, 1000, 1500, 2000 — two regions at each |
| Replicate iterations | 10 per (S, phi) combination |
| Regions per scenario | 10 |
| Sample size, `n` | 5000 |
| Generative model | `sparse`; `sparse_inf` at four values of `p_causal` (0.2, 0.4, 0.6, 0.8) |
| Annotation type | `none`; `binary`; `continuous` |
| Enrichment fold | 2.7, 5.4, 8.1, 10.8 (annotated rows only) |
| Annotation tracks | 10, of which the first 5 are enriched and 5 are inert |
| LD given to methods | in-sample throughout |

The 45 rows are the product of five model configurations (`sparse`, plus
`sparse_inf` at four values of `p_causal`) and nine annotation regimes (`none`,
plus `binary` and `continuous` at four enrichment folds each).

**Realised region sizes are smaller than requested** where a source genotype
file held fewer variants after filtering. Use `n_variants`, which records what
each fit actually saw; `region_size` is the level requested.

## Methods

`abf`, `beatrice`, `carma`, `fb_pooled`, `fb_xregion`, `finemap`, `finemap_inf`,
`finimom`, `functional_beatrice`, `funmap`, `marginal_z`, `paintor`,
`polyfun_est`, `polyfun_ldsc`, `polyfun_oracle`, `sbayesrc`, `sparsepro`,
`susie`, `susie_inf`.

`carma` covers only part of the grid, having been stopped on cost. `funmap`
requires annotations and is absent from the unannotated rows by construction.
`polyfun_oracle` is handed the simulation's true per-variant causal
probabilities and is a ceiling, not a deployable method.

---

## How this folder is organised

Every file is named for the parameter row it belongs to, and that name is also
the value of the `job_dir` column inside the tables, so any file joins to
`simulation-grid.csv` on `job_dir`. For example
`job_002_sparse_anBinary_e2.7_refInsample` is row 2: the sparse model, binary
annotations, enrichment fold 2.7, in-sample LD.

```
simulation-grid.csv                  the 45 parameter rows
DATA-DICTIONARY.md                   every column of every file

evaluation-metrics/
  per-cell.rds / .csv                one row per cell per method  (101,635 rows)
  per-replicate.rds                  one row per iteration per method
  per-fit.rds                        one row per region per method  (2,031,230)
  false-discovery-rate-by-threshold.rds
  calibration-by-pip-band.rds
  calibration-and-fdr-per-iteration.rds

credible-sets/
  summary-by-method.rds              coverage, power and size by method
  per-fit/<row>.rds                  the sets each method reported, per fit

posterior-inclusion-probabilities/
  <row>.rds                          the PIPs each method assigned (see below)

runtime-and-annotation-importance/
  <row>.rds                          seconds per fit; LassoNet annotation weights
```

## Three units of analysis

Metrics are reported at three levels, and the distinction matters.

- **fit** — one region, one method. `per-fit.rds`.
- **replicate** (or draw) — one iteration of one parameter row: the ten regions
  simulated together. `per-replicate.rds`.
- **cell** — `job_dir` × `S` × `phi` × `region_size`, which is 20 fits (ten
  iterations × the two regions of that size). `per-cell.rds`.

Rates such as false discovery rate are formed from counts pooled to at least the
region level, never from a single fit, because a rate on one or two selections is
not estimable. `calibration-and-fdr-per-iteration.rds` holds the alternative
estimator, which forms a rate within each replicate and averages across
replicates.

## Do not pool across strata

There are six strata: two generative models (`model`) crossed with three
annotation regimes (`annotation_type`). Methods reorder between them — a method
that leads under the sparse model can trail under the sparse-plus-infinitesimal
one — so an average taken across strata can invert the conclusion. Every table
carries `model` and `annotation_type` for this reason.

---

## What is and is not here

**The posterior inclusion probabilities are a tail, not the full vectors.** Each
fit stores every PIP at or above **0.01**, with its variant index and whether
that variant was causal, plus `n_below` and `sum_pip_below` summarising
everything under the floor. That is about 5% of all PIPs. Calibration and false
discovery rate can be recomputed at any threshold at or above 0.01. Average
precision cannot be recomputed exactly, because the identity of causal variants
below the floor is not retained — use the stored `ap`.

**Not included:** the simulated genotypes, the per-variant truth beyond what the
tables carry, the complete PIP vectors, and the figures and sensitivity analyses
of the accompanying report, all of which are reproducible from the code.

## Things a reader should know before using these numbers

- **`polyfun_ldsc`, `polyfun_oracle` and `sbayesrc` are reimplementations**
  written from the published methods, not the released programs. `sbayesrc` is
  implemented from its paper alone; `polyfun_ldsc` follows the paper where the
  paper and the released code disagree. Neither has been validated against its
  original software.
- **`polyfun_ldsc` and `polyfun_oracle` fail about 4.5% of fits**, all of them
  because the modified HESS estimator returns zero heritability, which PolyFun
  treats as an error. The failures concentrate at *small* regions, not large:
  roughly 5.7% at 500 variants and 0.6% at 2000. `failed` and `bad_pip` in
  `per-fit.rds`, and `n_failed` in `per-cell.rds`, record them. Comparisons
  against methods that do not fail should be paired on surviving fits.
- **`sbayesrc` reports no credible sets.** Its published method does not define
  them, so it is absent from `credible-sets/` and its rows in
  `summary-by-method.rds` show `nosets` equal to `fits`.
- **`sparsepro`'s credible sets are singletons only.** The parser that read them
  did not split sets joined on `/`, so sets of more than one variant were
  dropped. Its native coverage, power and size describe single-variant sets.
- **`marginal_z` is not a fine-mapping method.** It ranks variants by the
  magnitude of their marginal z-score and is included as a floor.

## Licence

Data: CC BY 4.0. Code: MIT, at the repository above.
