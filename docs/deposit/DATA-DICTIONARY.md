# Data dictionary

Every file is an R serialised object, read with `readRDS("path")`. See
`README.md` for the simulation design and for the caveats that bear on these
numbers.

## Keys shared by every table

| column | meaning |
|---|---|
| `job_dir` | the parameter row, e.g. `job_002_sparse_anBinary_e2.7_refInsample`. Joins to `simulation-grid.csv`, and is the filename in the per-row folders. |
| `scenario_id` | 1–250 within a row. One scenario is ten regions simulated together. |
| `region_id` | 1–10 within a scenario. |
| `method` | one of the nineteen methods. |
| `S` | causal variants simulated in the region: 1, 2, 3, 5 or 10. |
| `phi` | proportion of phenotypic variance explained (h²): 0.05–0.6. |
| `region_size` | the region size *requested*: 500, 750, 1000, 1500 or 2000. |
| `n_variants` | the region size *realised*, after MAF filtering. Use this one. |
| `model` | `sparse` or `sparse_inf`. |
| `annotation_type` | `none`, `binary` or `continuous`. |
| `p_causal` / `pc` | share of genetic variance in the sparse component (`sparse_inf` only). |
| `enrichment_fold` / `enrich` | annotation enrichment, 2.7–10.8. |
| `n_ref` | LD reference panel size; `NA` throughout, since this study used in-sample LD. |

---

## `simulation-grid.csv`

One row per parameter row (45). `job_dir` plus the factor settings for that row:
`model`, `annotation_type`, `p_causal`, `enrichment_fold`, `annotation_correlation`,
`n_annotations`, `n_regions`, `n` (sample size), `n_iter`, and the swept vectors
`S_values`, `phi_values`, `p_values`, `enrichment_values` (pipe-separated).

---

## `evaluation-metrics/per-fit.rds` — one row per region per method (2,031,230)

| column | meaning |
|---|---|
| `iter`, `region_idx` | replicate index within the row; region index within the row's panel. |
| `failed`, `bad_pip` | the fit errored, or returned non-finite PIPs. Exclude both before averaging. |
| `ap` | average precision, computed on this fit's own ranking. |
| `n_causal` | causal variants in this region (equals `S`). |
| `sum_pip_total` | total posterior mass the method assigned. |
| `prec_at_1`, `prec_at_S`, `prec_at_2S` | precision among the top 1, top `S`, top 2`S` variants. |
| `mass_above_50`, `mean_pip_topS` | mass on variants with PIP > 0.5; mean PIP of the top `S`. |
| `set_size_5`, `set_size_95` | size of the credible set built by accumulating ranked PIPs to 0.5 / 0.95 mass. This is a *uniform* rule applied to every method, not the method's own convention. |
| `set_hit_5`, `set_hit_95` | whether that set contains a causal variant (0/1). |
| `set_prec_5`, `set_prec_95` | proportion of the set that is causal. |
| `set_reached_5`, `set_reached_95` | whether the cumulative mass reached the level at all. `FALSE` means the whole region was assigned. |
| `tp_at_T`, `nsel_at_T`, `sum_pip_sel_at_T` | at PIP threshold `T` ∈ {50, 80, 90, 95, 99} (percent): causal variants selected, variants selected, and their summed PIP. **Counts, not rates** — form rates after pooling. |
| `n_band_B`, `c_band_B`, `sum_pip_band_B` | per calibration band `B`: variants, causal variants, summed PIP. Bands are frozen at `lo` = [0, 0.1), `mid` = [0.1, 0.5), `hi` = [0.5, 0.9), `top` = [0.9, 1]. |

## `evaluation-metrics/per-replicate.rds`

The same counts aggregated to one iteration of one parameter row — the ten
regions simulated together. Count columns are summed; scalar columns are meaned.

## `evaluation-metrics/per-cell.rds` / `.csv` — one row per cell per method (101,635)

A cell is `job_dir` × `S` × `phi` × `region_size`, which is 20 fits. Count
columns are summed into the cell; scalar columns (`ap`, `prec_at_*`,
`set_size_*`, `set_prec_*`, `mass_above_50`, `mean_pip_topS`) are meaned.
`set_hit_*` are Bernoulli per fit and summed, so coverage is `set_hit_95 / 20`.

Beyond the per-fit columns above:

| column | meaning |
|---|---|
| `n_fits`, `n_failed` | fits in the cell, and how many failed. |
| `ap_se` | standard error of `ap` across the cell's fits. |
| `fdr_at_T` | false discovery rate at threshold `T`, formed from the cell's pooled counts: `(nsel - tp) / nsel`. |
| `exp_fp_at_T` | false positives the method's own PIPs imply: `nsel - sum_pip_sel`. |
| `honesty_ratio_at_T` | observed false positives over expected. 1 is honest; above 1 is overconfident. |
| `total_mass_ratio` | summed PIP over the number of causal variants. |
| `rel`, `res`, `unc`, `bs`, `bss` | Murphy decomposition of the Brier score from the band counts: reliability `Σ n(p̄−ȳ)²/N`, resolution `Σ n(ȳ−ȳ_tot)²/N`, uncertainty `ȳ_tot(1−ȳ_tot)`, Brier score `rel − res + unc`, and skill `1 − bs/unc`. Lower `rel` is better calibrated; higher `res` is more discriminating. |
| `ece_hi` | expected calibration error over the `mid`, `hi` and `top` bands only. Plain ECE is dominated by the ~94% of variants in the lowest band and mis-ranks methods. |
| `hi_pip_reliab`, `hi_pip_n` | observed causal fraction among variants in the `top` band, and how many there were. |
| `power` | proportion of causal variants captured by the 0.95-mass credible set. Formed at fit level as `set_prec_95 × set_size_95 / n_causal` and then averaged, because the mean of a product is not the product of the means. |

## `evaluation-metrics/false-discovery-rate-by-threshold.rds` (807,720 rows)

One row per cell per method per threshold `t` ∈ {0.05, 0.1, 0.2, 0.5, 0.8, 0.9,
0.95, 0.99}, with `nsel`, `tp`, `sum_pip_sel` and `fdr`. A finer threshold grid
than the five frozen into the per-fit table, recomputed from the stored PIPs.

## `evaluation-metrics/calibration-by-pip-band.rds` (1,048 rows)

One row per stratum per method per band, over ten bands: [0, 0.01), [0.01, 0.05),
[0.05, 0.1), [0.1, 0.2), [0.2, 0.5), [0.5, 0.8), [0.8, 0.9), [0.9, 0.95),
[0.95, 0.99), [0.99, 1].

| column | meaning |
|---|---|
| `x`, `y` | mean PIP assigned in the band, and the proportion of those variants that were causal. A perfectly calibrated method has `y = x`; `y < x` is overconfident. |
| `se`, `lo`, `hi` | standard error and interval on `y`. |
| `x_pooled`, `y_pooled` | the same quantities pooling all variants rather than averaging over cells. |
| `n_total`, `cells_used`, `cells_dropped` | variants behind the estimate; cells contributing; cells with nothing in the band. |
| `reportable` | whether the band met the minimum variant and cell counts. **Filter on this.** |

## `evaluation-metrics/calibration-and-fdr-per-iteration.rds`

The alternative estimator the accompanying report adopted: a rate formed *within
each replicate* and averaged across replicates, rather than pooled across a
whole cell. A list with `cal` and `fdr` (the estimates), `cal_draw_full` and
`fdr_draw_full` (per-draw), `units` (the underlying count matrices at cell and
draw level) and `meta` (the bands, thresholds and minimum counts used). Covers
the seven methods the report displays, not all nineteen.

Iterations holding no variants in a band, or making no selection at a threshold,
have no rate and are excluded rather than scored as zero; `units_used` and
`units_dropped` record this. Excluding them matters: scoring an abstention as
zero would reward a method for selecting nothing.

---

## `credible-sets/` — the methods' own sets

`summary-by-method.rds` — one row per stratum per method:

| column | meaning |
|---|---|
| `fits`, `nosets` | fits, and fits where the method reported no set at all. |
| `sets`, `sets_per_fit` | sets reported, and per fit. |
| `hits`, `coverage_per_set` | sets containing a causal variant, and the proportion. |
| `size`, `mean_size` | total and mean variants per set. |
| `power` | proportion of causal variants captured by any reported set. |
| `cap`, `S`, `pct_no_set` | the method's cap on set count where it has one; causal variants simulated; percentage of fits with no set. |

`per-fit/<row>.rds` — a list with `job_dir`, `fits`, and `rows`, one record per
fit: `job_dir`, `scenario_id`, `region_id`, `method`, `n_variants`, `reported`
(whether the method reports sets at all), `n_sets`, `set_sizes`, and `sets`, the
variant indices of each set.

These are each method's **own** convention, which differs between methods — one
set per single-effect component, or per selected causal count, or ranked mass to
a level — so coverage is not comparable across methods without care. The uniform
0.95-mass reconstruction in `per-fit.rds` (`set_*_95`) is the like-for-like
comparison.

---

## `posterior-inclusion-probabilities/<row>.rds`

A list with `job_dir`, `floor` (0.01), `source`, `fits`, and `tail`, one record
per fit:

| field | meaning |
|---|---|
| `idx` | variant indices of the retained PIPs, within the region. |
| `pip` | the PIPs themselves, all ≥ `floor`. |
| `is_causal` | whether each retained variant was causal. |
| `n_below`, `sum_pip_below` | how many variants fell below the floor and their summed mass. |
| `n_variants`, `n_causal` | region size realised, and causal variants simulated. |

About 5% of all PIPs are retained. Calibration and false discovery rate can be
recomputed at any threshold ≥ 0.01; average precision cannot, because which
below-floor variants were causal is not retained.

---

## `runtime-and-annotation-importance/<row>.rds`

A list with `job_dir`, `scenarios_read`, and two record lists:

- `runtime` — `job_dir`, `scenario_id`, `region_id`, `method`,
  `runtime_seconds`, `scope`, `n_total`, `n_failed`. `scope` distinguishes a
  per-region fit from work done once per scenario: methods that train one prior
  across all ten regions before fitting any of them do most of their work there,
  and charging it per region would misattribute it.
- `importance` — the LassoNet annotation weights the BEATRICE-family methods
  return, per scenario. Empty for every other method, which do not produce them.
