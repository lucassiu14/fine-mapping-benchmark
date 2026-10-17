# Main simulation: fine-mapping benchmark results

Results of the main simulation study from `fine-mapping-benchmark`, in which
eighteen statistical fine-mapping methods were run on identical simulated data
and scored with identical metrics.

Code: https://github.com/lucassiu14/fine-mapping-benchmark

---

## Everything here is one fit

A **fit** is one method applied to one region of one iteration. Every file is at
that grain, and only that grain — one row per fit, or one row per thing *within*
a fit. There are no per-cell averages, false discovery rates or calibration
curves here, because all of them are aggregates of these files and which
aggregation was used is a choice the reader should make, not inherit. The
`DATA-DICTIONARY.md` states the rules the accompanying report used.

The single exception is `scenario-runtimes.rds`, named so you can see it is one.

## A fit is identified by the simulation, not by how it was run

| | |
|---|---|
| the parameter row | `model`, `p_causal`, `annotation_type`, `enrichment_fold` |
| the fit within it | `S`, `phi`, `iter`, `region_id`, `method` |

Those nine columns are unique across all 2,031,230 fits. The four row-level
factors name the file in the per-row folders rather than repeating down every
row of it, so `sparse-inf_rho0.2_continuous_fold10.8.rds` *is* the parameter row
it holds.

Nothing here records the job scheduler, the directory layout or the scenario
numbering those produced.

---

## The files

```
simulation-grid.csv          the 45 parameter rows, and which file each is
fit-metrics.rds / .csv.gz    one row per fit: every scalar result
posterior-inclusion-probabilities/<row>.rds
                             one row per retained variant
credible-sets/<row>.rds      one row per variant per reported set
scenario-runtimes.rds        one row per scenario per method (see below)
DATA-DICTIONARY.md           every column of every file
```

Files are R serialised objects, read with `readRDS("path")`. `fit-metrics` is
also gzipped CSV, and `simulation-grid.csv` is plain CSV, so the two entry points
need no R.

**Start with `fit-metrics`.** It carries average precision, credible-set
measures, the threshold and calibration-band counts, runtime, and the per-fit
counts the two long tables must be read against. Most questions are answered by
grouping it.

**`scenario-runtimes.rds`** holds work done once per scenario rather than per
region. Methods that train one annotation prior across all ten regions before
fitting any of them do nearly all their work there, so their cost appears only
here: `runtime_seconds` in `fit-metrics` is missing entirely for `paintor` and
for about 89% of `fb_xregion` and `fb_pooled` fits. Add the two files to compare
methods on cost; using `fit-metrics` alone makes the cross-region methods look
nearly free.

---

## The simulation

Genotypes were simulated from 1000 Genomes Phase 3 haplotypes with `sim1000G`,
keeping variants with minor allele frequency of at least 0.01, and each region
was standardised so that the in-sample linkage disequilibrium matrix is a
correlation matrix. Every method received the same marginal summary statistics,
the same LD matrix and, where applicable, the same annotation matrix. LD was
in-sample throughout: no method was given a reference panel.

45 parameter rows, each simulated 250 times, giving 11,250 scenarios. A scenario
holds ten independent regions, so each method was fitted to 112,500 regions;
across eighteen methods, 1,918,730 fits.

| Factor | Levels |
|---|---|
| Causal variants per region, `S` | 1, 2, 3, 5, 10 |
| Variance explained, `phi` (h²) | 0.05, 0.1, 0.2, 0.4, 0.6 |
| Variants per region, `region_size` | 500, 750, 1000, 1500, 2000 — two regions at each |
| Iterations, `iter` | 10 per (S, phi) combination |
| Regions per scenario | 10 |
| Sample size | 5000 |
| `model` | `sparse`; `sparse_inf` |
| `p_causal` | 0.2, 0.4, 0.6, 0.8 — share of genetic variance in the sparse component; `NA` under `sparse` |
| `annotation_type` | `none`; `binary`; `continuous` |
| `enrichment_fold` | 2.7, 5.4, 8.1, 10.8; `NA` where there are no annotations |
| Annotation tracks | 10, of which the first 5 are enriched and 5 are inert |

The 45 rows are the product of five model configurations (`sparse`, and
`sparse_inf` at four values of `p_causal`) and nine annotation regimes (`none`,
and `binary` and `continuous` at four enrichment folds each).

**Realised region sizes are smaller than requested** where a source genotype file
held fewer variants after filtering. `n_variants` is what each fit actually saw;
`region_size` is the level requested. Use `n_variants`.

## Methods

`abf`, `beatrice`, `carma`, `fb_pooled`, `fb_xregion`, `finemap`, `finemap_inf`,
`finimom`, `functional_beatrice`, `funmap`, `marginal_z`, `paintor`,
`polyfun_ldsc`, `polyfun_oracle`, `sbayesrc`, `sparsepro`, `susie`,
`susie_inf`.

`carma` covers only part of the grid, having been stopped on cost. `funmap`
requires annotations and is absent from the unannotated rows by construction.
`polyfun_oracle` is handed the simulation's true per-variant causal probabilities
and is a ceiling, not a deployable method. `marginal_z` ranks variants by the
magnitude of their marginal z-score and is a floor, not a fine-mapping method.

## Do not pool across strata

Six strata: two generative models crossed with three annotation regimes. Methods
reorder between them — a method that leads under the sparse model can trail under
the sparse-plus-infinitesimal one — so averaging across strata can invert the
conclusion.

---

## What to know before using these numbers

- **The posterior inclusion probabilities are a tail, plus every causal
  variant.** Two things are kept: every PIP at or above **0.01**, and every
  causal variant whatever its PIP. The second matters — a causal variant a
  method scored near zero is precisely what a power calculation must count, and
  dropping it would silently inflate recall at every threshold. So about 2.3% of
  the rows sit *below* the floor, and `pip` can be as low as 0. Roughly 196
  million rows in all.
  `fit-metrics` carries `pip_n_below_floor` and `pip_mass_below_floor` for
  everything that did not meet the floor, so calibration and false discovery rate
  can be recomputed at any threshold at or above 0.01, and the lowest band can be
  completed exactly. **Average precision cannot be recomputed**: the non-causal
  variants below the floor are summarised rather than listed, so the rank of a
  low-scoring causal variant among them is unknown. Use the stored `ap`.
- **`polyfun_ldsc`, `polyfun_oracle` and `sbayesrc` are reimplementations**
  written from the published methods, not the released programs. `sbayesrc` is
  implemented from its paper alone; `polyfun_ldsc` follows the paper where paper
  and released code disagree. Neither has been validated against its original
  software.
- **`polyfun_ldsc` and `polyfun_oracle` fail about 4.5% of fits**, all because
  the modified HESS estimator returns zero heritability, which PolyFun treats as
  an error. The failures concentrate at *small* regions: roughly 5.7% at 500
  variants and 0.6% at 2000. `failed` and `bad_pip` mark them. Comparisons with
  methods that do not fail should be paired on surviving fits.
- **`sbayesrc` reports no credible sets.** Its published method does not define
  them, so it has no rows in `credible-sets/` and `native_n_sets` is 0 for every
  one of its fits. Note that `native_sets_reported` is nonetheless `TRUE` for it:
  that column records whether the method's output carried a credible-set field at
  all, not whether it ever populated one. `native_n_sets` is the column to use.
- **`sparsepro`'s credible sets are singletons only.** The parser that read them
  did not split sets joined on `/`, so sets of more than one variant were lost.
- **The uniform credible-set columns are not each method's own sets.**
  `set_size_95` and the rest come from accumulating ranked PIPs to 0.95 of the
  mass, applied identically to every method, which is the like-for-like
  comparison. `credible-sets/` holds what each method itself reported, under its
  own convention, which differs between methods and is not comparable without
  care.

## Licence

Data: CC BY 4.0. Code: MIT, at the repository above.
