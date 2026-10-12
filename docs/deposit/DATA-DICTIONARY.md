# Data dictionary

Files are R serialised objects, read with `readRDS("path")`. See `README.md` for
the simulation design and the caveats that bear on these numbers.

## The key

Every file is at the grain of a **fit** — one method, one region, one iteration.

| column | meaning |
|---|---|
| `model` | `sparse`, or `sparse_inf` for sparse plus an infinitesimal background. |
| `p_causal` | share of genetic variance carried by the sparse component: 0.2, 0.4, 0.6, 0.8. `NA` under `sparse`, where all of it is. Written ρ in the accompanying report. |
| `annotation_type` | `none`, `binary` or `continuous`. |
| `enrichment_fold` | annotation enrichment: 2.7, 5.4, 8.1 or 10.8. `NA` where there are no annotations. |
| `S` | causal variants simulated in the region: 1, 2, 3, 5 or 10. |
| `phi` | proportion of phenotypic variance explained, h²: 0.05 to 0.6. |
| `iter` | replicate index, 1–10, within the (`S`, `phi`) combination. |
| `region_id` | which of the scenario's ten regions, 1–10. Regions 1 and 2 hold 500 variants, 3 and 4 hold 750, and so on. |
| `method` | one of the nineteen methods. |

The first four identify the parameter row and **name the file** in
`posterior-inclusion-probabilities/` and `credible-sets/`; the last five
identify the fit within it. Together they are unique over all 2,031,230 fits.

`region_id` is comparable across methods: the same `region_id` under the same
row, `S`, `phi` and `iter` is the same simulated region, which is what makes
paired comparisons between methods possible.

---

## `simulation-grid.csv` — 45 rows

`file` (the filename stem in the per-row folders), the four row factors, then
`annotation_correlation`, `n_annotations`, `n_regions`, `n` (sample size, 5000),
`n_iter`, and the swept vectors `S_values`, `phi_values`, `p_values`,
`enrichment_values`, pipe-separated.

---

## `fit-metrics.rds` / `.csv.gz` — one row per fit (2,031,230)

The key above, then `region_size` (the size requested; `n_variants` is what the
fit saw), then:

**Did it work**

| column | meaning |
|---|---|
| `failed`, `bad_pip` | the fit errored, or returned non-finite PIPs. Exclude both before averaging anything. |

**Ranking**

| column | meaning |
|---|---|
| `ap` | average precision on this fit's own ranking. Not recomputable from the deposited PIPs — see the README. |
| `prec_at_1`, `prec_at_S`, `prec_at_2S` | precision among the top 1, top `S`, top 2`S` variants. |
| `mass_above_50`, `mean_pip_topS` | posterior mass on variants with PIP > 0.5; mean PIP of the top `S`. |
| `sum_pip_total` | total posterior mass the method assigned. |
| `n_variants`, `n_causal` | variants in the region as fitted, and causal variants simulated (equals `S`). |

**Credible sets, uniform rule.** Accumulating ranked PIPs to 0.5 or 0.95 of the
mass, applied identically to every method. Not the method's own sets.

| column | meaning |
|---|---|
| `set_size_5`, `set_size_95` | variants in the set. |
| `set_hit_5`, `set_hit_95` | whether it contains a causal variant (0/1). |
| `set_prec_5`, `set_prec_95` | proportion of the set that is causal. |
| `set_reached_5`, `set_reached_95` | whether the cumulative mass reached the level. `FALSE` means the whole region was assigned. |

**Counts at thresholds and in bands.** These are **counts, not rates**. A rate
formed on one or two selections is not estimable, so form rates only after
pooling fits — the report pools to at least the region level.

| column | meaning |
|---|---|
| `tp_at_T`, `nsel_at_T`, `sum_pip_sel_at_T` | at PIP threshold `T` ∈ {50, 80, 90, 95, 99} percent: causal variants selected, variants selected, their summed PIP. False discovery rate is `(nsel − tp) / nsel` on pooled counts; `nsel − sum_pip_sel` is the false positives the method's own PIPs imply. |
| `n_band_B`, `c_band_B`, `sum_pip_band_B` | per band `B`: variants, causal variants, summed PIP. Bands are `lo` = [0, 0.1), `mid` = [0.1, 0.5), `hi` = [0.5, 0.9), `top` = [0.9, 1]. A calibration curve is `sum_pip_band_B / n_band_B` against `c_band_B / n_band_B`, pooled over whatever set of fits you choose. |

Written out, those 27 columns are `tp_at_50`, `nsel_at_50`, `sum_pip_sel_at_50`
and the same three at `_80`, `_90`, `_95` and `_99`; and `n_band_lo`,
`c_band_lo`, `sum_pip_band_lo` and the same three at `_mid`, `_hi` and `_top`.

**Cost**

| column | meaning |
|---|---|
| `runtime_seconds` | wall-clock for this fit. `NA` where the method does its work per scenario rather than per region: entirely so for `paintor`, and for about 89% of `fb_xregion` and `fb_pooled` fits. Their cost is in `scenario-runtimes.rds`. |

**Reading the two long tables against this one**

| column | meaning |
|---|---|
| `pip_n_retained` | rows this fit has in `posterior-inclusion-probabilities/`. |
| `pip_n_below_floor`, `pip_mass_below_floor` | variants that did not meet the 0.01 floor, and their summed PIP. **These overlap with the retained rows**: a causal variant scored below the floor is retained *and* counted here, so `pip_n_retained + pip_n_below_floor` exceeds `n_variants` by the number of such variants. To complete the lowest calibration band exactly, take `n` = `pip_n_below_floor`, `summed PIP` = `pip_mass_below_floor`, and count the causal ones directly as the retained rows with `is_causal == 1` and `pip < 0.01`. |
| `native_sets_reported` | whether the method's output carried a credible-set field at all — **not** whether it populated one. It is `TRUE` for `sbayesrc`, which never reports a set. Use `native_n_sets`. |
| `native_n_sets` | sets it reported for this fit. `0` means it reported none, which is why the fit has no rows in `credible-sets/` — distinguishing that from a missing fit. |

---

## `posterior-inclusion-probabilities/<row>.rds` — one row per retained variant

| column | meaning |
|---|---|
| `S`, `phi`, `iter`, `region_id`, `method` | the fit. The parameter row is the filename. |
| `variant` | the variant's index within the region, 1 to `n_variants`. Comparable across methods for the same region. |
| `pip` | the posterior inclusion probability. Usually ≥ 0.01 — see the retention rule below — and `NA` on the few fits that returned non-finite output, which `bad_pip` marks. |
| `is_causal` | whether that variant was simulated causal (0/1). |

**The retention rule.** Two sets of variants are kept: every variant with a PIP
at or above 0.01, and **every causal variant whatever its PIP**. The second is
deliberate — a causal variant a method scored near zero is exactly what a power
calculation has to count, and dropping it would inflate recall at every
threshold. So about 2.3% of rows carry a `pip` below the floor, and all of those
are causal. Around 196 million rows across the 45 files, about 5% of every PIP
the study produced.

Because causal variants are always retained, `sum(is_causal)` for a fit equals
its `n_causal` — except on fits with non-finite output, where the PIP is `NA`.

## `credible-sets/<row>.rds` — one row per variant per reported set

| column | meaning |
|---|---|
| `S`, `phi`, `iter`, `region_id`, `method` | the fit. |
| `set` | the set's label as the method gave it, e.g. `L1`, `L2`. Variants sharing a label are one set. |
| `variant` | the variant's index within the region, as above. |

These are each method's **own** convention — one set per single-effect
component, or per selected causal count, or ranked mass to a level — so coverage
and size are not comparable across methods without care. For a like-for-like
comparison use the `set_*` columns in `fit-metrics`.

A fit with no sets has no rows here; `native_n_sets` is 0 for it.

## `scenario-runtimes.rds` — one row per scenario per method

The four row factors, then `S`, `phi`, `iter`, `method`, and:

| column | meaning |
|---|---|
| `runtime_seconds` | wall-clock for work done once for the whole scenario rather than per region. |
| `n_total`, `n_failed` | fits attempted and failed in that scenario. |

The one file not at fit grain. Methods that train one annotation prior across all
ten regions before fitting any of them book nearly all their cost here, so
`runtime_seconds` in `fit-metrics` understates them badly on its own. Summing the
two gives a method's true cost for a scenario.
