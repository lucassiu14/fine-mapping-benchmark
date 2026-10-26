# Evaluation Methodology

This document describes every step of the evaluation pipeline implemented in
`R/evaluate.R` and `R/evaluate_per_region.R`: how raw method outputs are matched to
ground truth, how each metric is computed, how standard errors are derived, how
results are stratified, and what each plot shows.

---

## Contents

1. [Inputs](#inputs)
2. [Annotating fits with ground truth](#annotating-fits-with-ground-truth)
3. [Grouping: there are no fixed strata](#grouping-there-are-no-fixed-strata)
4. [Metric computation](#metric-computation)
   - [Power / FDR curve](#power--fdr-curve)
   - [Average precision](#average-precision-ap)
   - [PIP calibration](#pip-calibration)
   - [Credible-set metrics](#credible-set-metrics)
   - [Runtime](#runtime)
5. [Output structure](#output-structure)
6. [Saving to disk](#saving-to-disk)
7. [Quick reference — field names](#quick-reference--field-names)

---

## Inputs

`evaluate_methods()` takes two required arguments:

- **`simulation`** — the list returned by `run_simulation()`. The fields used
  are `simulation$scenarios` (a list of scenario objects, each containing a
  `regions` sub-list whose entries carry a `truth` object) and
  `simulation$params` (stored verbatim in the output for reproducibility).

- **`results`** — the list returned by `run_methods()`. The fields used are
  `results$methods_run` (the character vector of method names that were run) and
  `results[[method]]$results` (the flat list of per-fit output objects, one per
  scenario × region combination).

Optional arguments:

| Argument | Default | Meaning |
|---|---|---|
| `pip_thresholds` | `seq(0, 1, by = 0.005)` | Thresholds at which the power/FDR curve is evaluated (201 points) |
| `n_pip_cal_bins` | `10` | Number of equal-width bins for PIP calibration |
| `save` | `FALSE` | Write `evaluation.rds` and `evaluation_summary.csv` to `output_dir` |
| `output_dir` | `"results"` | Directory for saved outputs |
| `verbose` | `TRUE` | Print progress messages |

---

## Annotating fits with ground truth

Before any metric is computed, every fit object in `results[[method]]$results`
is augmented with two fields from the simulation ground truth.

For each fit `f`:

1. The corresponding scenario is looked up: `sc = simulation$scenarios[[f$scenario_id]]`.
2. The corresponding region truth object is looked up: `truth = sc$regions[[f$region_id]]$truth`.
3. Two fields are added to `f`:
   - `f$causal_indices` — integer vector of 1-based variant positions that are
     truly causal in this scenario × region. Length equals the number of causal
     variants `S` for this scenario.
   - `f$n_variants` — integer equal to `length(f$pip)`, i.e. the number of
     variants in the region.

This annotation is performed by `.annotate_fits_with_truth()` and is applied
to every fit before any metric is computed. Fits where `f$error` is non-NULL
(i.e. the method failed) are still annotated but are excluded from all metric
computations (see below).

---

## Grouping: there are no fixed strata

`evaluate_methods()` does not pre-compute a set of strata. It returns `global`
per method - the metrics over all of that method's fits - plus one
`per_region` data frame whose unit is **one method fitted to one region of one
iteration**.

The design factors travel with each row, so any grouping is a `aggregate()`
call away:

| Column | What it is |
|---|---|
| `method`, `scenario_id`, `region_id`, `iter` | The identity of the fit. `(method, scenario_id, region_id)` and `(method, S, phi, iter, region_id)` are both unique keys over the same rows, because each scenario is one `(S, phi, iter)` draw |
| `S` | Number of causal variants in the region |
| `phi` | Heritability of the region |
| `p_causal` | Sparse-component share (`NA` outside the `sparse_inf` model) |
| `region_size` | Nominal region length, the requested per-region `p`. Fine-mapping difficulty scales steeply with it, so it is a genuine axis and metrics should never be pooled across region lengths without saying so |
| `min_causal_maf`, `causal_maf_bin` | Lowest causal-variant MAF, and its bin: `rare` (<= 0.01), `low` (<= 0.05), `common`. `NA_character_` when the region has no causal variants |
| `true_annotation_type` | `none`, `binary` or `continuous`, from the simulation's parameters. Constant within one simulation; it earns its place when evaluation objects from differently-annotated runs are combined |
| `failed`, `bad_pip` | Whether the fit errored, and whether it returned non-finite PIPs |

A fixed set of strata computed inside the package would only be one analysis's
choice of axes. The cost of this freedom is that the aggregation rules are the
caller's to follow - see the next chapter.

---

## Metric computation

All metrics are implemented in `.metrics_for_fits()`. Fits where `f$error` is
non-NULL are counted in `n_failed` but excluded from all computations. If no
valid (non-failed) fits remain, all metrics are returned as `NA` or `NULL`.

### What is pooled across regions, and what is not

This distinction decides which numbers are meaningful, so it is worth stating
before the formulas:

> **Rates are pooled from counts. Per-region scalars are averaged.**

- **Average precision is a per-region quantity**, so `ap` is computed on each
  region separately and then averaged. The AP of two regions' merged ranked
  list is *not* the mean of their APs: merging lets one region's false
  positives outrank the other region's true positive, which describes an
  analysis nobody performs. The difference is not small and not uniform — on a
  three-method test fixture, the pooled estimator moved `susie` up by 0.065,
  `abf` down by 0.048 and `marginal_z` down by 0.195, so it reorders methods.
  `marginal_z` and `abf` rank variants identically *within* a region yet the
  pooled estimator separated them by a factor of 3.5, because it was partly
  measuring their PIP scale rather than their ranking.

- **FDR, power and calibration are also reported per region**, because
  fine-mapping is applied to a region: the question a practitioner asks is "I
  fine-mapped my locus — how trustworthy are these calls?", not "across every
  call made anywhere in the benchmark, what fraction were wrong". Those are
  different estimands, not a precise and an imprecise version of one. The
  call-weighted rate lets a region that makes many calls dominate; under the
  region-weighted one every region counts once.

  So `fdr`, `precision`, `power` and `frac_causal` are the **means over
  regions** of each region's own rate, and the call-weighted values sit beside
  them with a `_pooled` suffix. The underlying counts (`tp`, `fp`, `fn`, `n`,
  `n_causal`) are kept, so either rate can be re-formed at any granularity.

  **One convention has to change when averaging.** FDR is 0 by convention for
  a region that selects nothing. That is harmless when pooling — such a region
  contributes no calls — but under averaging it *rewards abstention*: a method
  selecting nothing anywhere would average to a perfect FDR of 0. So those
  regions are excluded from the mean, and `n_regions` reports how many
  contributed. Never read the rate without it: at a high PIP threshold, most
  regions select nothing. In a test fixture, `marginal_z` selected nothing
  anywhere at PIP 0.95 — the pooled convention scored it a flawless FDR of
  0.0000, while the region-weighted version correctly returned `NA` on 0
  regions.

- **Credible-set coverage, power and size are region-weighted too.**
  `cs_coverage` is the mean over regions of the share of *that region's* sets
  holding a causal variant; pooling over sets instead would weight a region
  reporting five sets five times as heavily as one reporting a single set, and
  `cs_power` pooled over causal variants weights a region with four causals
  four times as heavily as one with a single causal. Measured with S varying
  (1 and 4), `cs_power` came out 0.6500 pooled against 0.7812 region-weighted,
  and `cs_size_median` 2.50 against 3.72.

  A region that reported no sets has *undefined* coverage and set size — there
  is no set to be right or wrong about — so it is excluded from those means
  and counted in `cs_coverage_n`. Its power is 0, which is defined: none of its
  causal variants were captured.

**There is no `auprc`.** It and `ap_pooled_grid` estimated the same quantity as
`ap`, the area under the precision–recall curve, but on the pooled ranking —
and `auprc` additionally used the trapezoidal rule, whose linear interpolation
is not valid in PR space (Davis & Goadrich 2006; the segment between two PR
points is not achievable by any classifier). The three disagreed about which
method was best: on a four-method fixture `ap` ranked `susie_inf > susie`,
while both pooled estimators reversed it. Keeping them as alternatives meant
shipping two wrong answers alongside the right one, so they were removed and
`ap` is the single PR-curve summary.

### The per-region table

`per_region_metrics()` returns one row per (method, scenario, region): the
design factors, the per-region scalars, and the counts that rates are formed
from. The rules above are the caller's to follow - see the
[Quick reference](#quick-reference--field-names) for which column is in which
family.
Both are also reachable as `evaluate_methods(...)$per_region`.

A failed fit contributes an all-`NA` metric row rather than being dropped, so
that a method which cannot run on part of the grid stays visible — such absence
is structural, not missing at random. A fit whose PIP vector contains any
non-finite value is marked `bad_pip` and treated as failed, never repaired:
such a posterior has undefined total mass and an undefined ranking, and
imputing the non-finite entries to zero would invent posterior mass and flatter
precisely the method that overflowed.

### Power / FDR curve

**What is computed.** For a collection of valid fits, all PIP values and their
corresponding causal-variant labels are pooled into two vectors:

```
pip_vals  ← unlist of f$pip for every valid fit f
is_causal ← logical vector of the same length;
            is_causal[i] = TRUE if variant i (in the pooled ordering)
            has its 1-based position in f$causal_indices
```

For a PIP threshold `t`, every variant with `pip >= t` is declared a positive
call. The four confusion-matrix quantities are:

```
TP(t) = number of called variants that are truly causal
FP(t) = number of called variants that are not truly causal
FN(t) = number of truly causal variants not called
n_sel(t) = TP(t) + FP(t)   (total positives called)
total_causal = sum(is_causal)
```

Derived quantities at each threshold:

```
FDR(t)       = FP(t) / n_sel(t)      [0 when n_sel(t) = 0]
Power(t)     = TP(t) / total_causal
Precision(t) = TP(t) / n_sel(t)      [1 when n_sel(t) = 0]
Recall(t)    = TP(t) / total_causal  [identical to Power(t)]
```

**Algorithm.** Computing this naively at `T` thresholds requires O(N × T) time.
The implementation uses an O(N log N + T) approach:

1. Sort all variants by decreasing PIP. Let `sorted_pips` be the sorted vector
   and `cum_tp[k]` be the cumulative number of causal variants in the top-`k`
   entries (with `cum_tp[0] = 0`).
2. For each threshold `t`, the number of selected variants is
   `n_sel(t) = max { k : sorted_pips[k] >= t }`. Because `sorted_pips` is
   non-increasing, `-sorted_pips` is non-decreasing. `n_sel(t)` is found with
   `findInterval(-t, -sorted_pips)`, which is an O(log N) binary search.
3. `TP(t) = cum_tp[n_sel(t)]`, looked up in O(1).

**Output.** A `data.frame` with one row per threshold and columns:
`threshold`, `tp`, `fp`, `fn`, `fdr`, `power`, `precision`, `recall`.

After SE computation (see below), two additional columns are present:
`n_regions`, the `*_pooled` call-weighted variants, and the `tp`/`fp`/`fn`
totals.

**Boundary behaviour:**
- At `t = 0`: all variants are called, so `n_sel = N`, `TP = total_causal`,
  `FN = 0`, `Power = 1`, `Precision = total_causal / N`.
- At `t = 1`: only variants with `pip = 1.0` exactly are called. In practice
  this is usually zero or very few variants.
- When `n_sel(t) = 0`: `FDR = 0` and `Precision = 1` by convention (no false
  calls have been made).
- When `total_causal = 0` (no causal variants in the pooled set): all metrics
  are set to 0 and `Precision = 1`.

---

### Average precision (`ap`)

`ap` is the **mean over regions** of each region's own average precision,
computed by `compute_ap_exact()`:

```
AP = sum over distinct observed PIP values of:
     (recall_i - recall_{i-1}) x precision_i
```

Three choices in that definition matter, and each was a defect in the estimator
this replaced:

1. **The step integral, not the trapezoid.** Linear interpolation between two
   precision-recall points is not valid — the segment between them is not
   achievable by any classifier (Davis & Goadrich 2006). The trapezoid
   understated by roughly 29% in the weak-signal regime and roughly 0% at
   strong signal, so its error was signal-dependent rather than a constant
   offset.

2. **Evaluated at the distinct observed PIP values, not on a fixed grid.** A
   grid treats every variant inside one band as entering the selection
   together, so a causal variant is credited with the precision computed after
   every non-causal in its band has also been admitted. A causal variant at
   PIP 0.0032 with true rank 61 of 1000 has exact AP 1/61 = 0.0164 but grid AP
   1/1000 = 0.001. The error grows as PIPs compress toward zero, making it a
   heritability-dependent bias.

3. **Per region, then averaged.** The AP of several regions' merged ranked list
   is not the mean of their APs, because merging lets one region's false
   positives outrank another region's true positive.

Ties are scored as a block, which is the conservative convention: a method
assigning the same PIP to a causal and a non-causal variant gets the precision
of the whole block rather than the benefit of an arbitrary tie-break.

AP is preferred over AUROC for fine-mapping because the class imbalance is
extreme: in a typical region with 100-500 variants and S = 1-5 causal
variants, the positive class is under 5% of variants, and AUROC is insensitive
to performance in that regime.

A perfect method scores 1. A random PIP assignment scores about
`total_causal / N`, the prevalence of causal variants, which is typically
below 0.05.

---

### PIP calibration

PIP calibration measures whether a method's PIP values are well-calibrated:
a variant assigned PIP = 0.8 should be causal in approximately 80% of cases.

**Binning.** The unit interval [0, 1] is divided into `n_pip_cal_bins`
equal-width bins using breakpoints `seq(0, 1, length.out = n_pip_cal_bins + 1)`.
Variants are assigned to bins using `findInterval(..., rightmost.closed = TRUE)`,
so the final bin includes `pip = 1.0`.

**Per-bin quantities.** For each bin `b` (1-indexed):

```
n[b]           = number of variants in bin b
n_causal[b]    = number of causal variants in bin b
mean_pip[b]    = mean PIP of variants in bin b  (NA if n[b] = 0)
frac_causal[b] = n_causal[b] / n[b]            (NA if n[b] = 0)
```

**Calibration interpretation.** A well-calibrated method satisfies
`mean_pip[b] ≈ frac_causal[b]` for all bins — the expected causal probability
(mean PIP) equals the observed fraction of causal variants. Perfect calibration
corresponds to all points lying on the identity line `y = x`.

**Output.** A `data.frame` with one row per bin and columns:
`bin`, `bin_lower`, `bin_upper`, `bin_mid`, `n`, `n_causal`, `mean_pip`,
`frac_causal`.

`n_regions` gives the number of regions populating each bin, and
`mean_pip_pooled` / `frac_causal_pooled` the call-weighted values.

---

### Credible-set metrics

Credible-set (CS) metrics are computed by `.cs_metrics()` over the list of
valid fits. Each fit carries `f$credible_sets` (a list of integer vectors, each
giving the 1-based variant indices in one CS) and `f$causal_indices` (the true
causal positions).

**Per-CS quantities.** For each CS `cs` within each fit `f`:

```
cs_hit[cs]  = any(cs %in% f$causal_indices)   # TRUE if CS contains ≥1 causal
cs_size[cs] = length(cs)                       # number of variants in CS
```

**Per-causal-variant coverage.** For each fit `f`, the set of all variants
covered by any CS is `all_cs_variants = unique(unlist(f$credible_sets))`. Each
truly causal variant is recorded as either covered (`f$causal_indices %in%
all_cs_variants`) or not. Fits with no reported CSs contribute `FALSE` for
every causal variant.

**Aggregated CS metrics:**

| Metric | Formula | Interpretation |
|---|---|---|
| `cs_coverage` | `mean(cs_hit)` over all CSs from all fits | Proportion of reported CSs containing ≥1 true causal variant. Should be ≥ the nominal coverage level (e.g. 0.95 for 95% CSs). |
| `cs_power` | `mean(causal_captured)` over all causal variants from all fits | Proportion of true causal variants captured by at least one CS. |
| `cs_size_median` | `median(cs_size)` over all CSs | Median number of variants per CS. Smaller is better, provided coverage is maintained. |
| `cs_size_mean` | `mean(cs_size)` over all CSs | Mean number of variants per CS. |
| `n_cs_reported` | `length(cs_hit)` | Total number of CSs reported across all fits. |

If no CSs are reported across all fits, `cs_coverage`, `cs_power`,
`cs_size_median`, `cs_size_mean` are `NA` and `n_cs_reported = 0`.

**Note on ABF.** ABF returns exactly one CS per fit under its single-causal
assumption, so `n_cs_reported = n_valid_fits` always. When S > 1, ABF's CS
will contain a causal variant (cs_coverage remains high) but will miss the
remaining causal variants (cs_power is correspondingly lower than for methods
that support multiple CSs).

---

### Runtime

Runtime is measured per fit as `runtime_seconds` (wall-clock time, recorded by
each method wrapper using `proc.time()`). Aggregate statistics are computed
over valid (non-failed) fits only:

```
runtime_mean = mean(runtime_seconds)   over valid fits
runtime_sd   = sd(runtime_seconds)     over valid fits (NA if n_valid < 2)
```

`runtime_sd` is the SD over valid fits, and each region's own
`runtime_seconds` is in the per-region table.

---

## Output structure

`evaluate_methods()` returns a named list: one element per evaluated method,
one shared per-region table, and three metadata fields.

```
eval_out
├── methods_evaluated       character vector of method names
├── simulation_params       copy of simulation$params
├── pip_thresholds_used     the pip_thresholds vector used
├── per_region              data.frame, one row per (method, region, iteration)
└── <method>                one entry per method (e.g. "susie", "abf")
    └── global              the metrics over all of that method's fits
```

There are no `by_*` strata and no `_se` fields. Group `per_region` instead, and
take standard errors across iterations (see
[Grouping](#grouping-there-are-no-fixed-strata)).

### `global`

Every rate here is **region-weighted**: computed within each region, then
averaged over the regions where it is defined. Each ships with the count of
contributing regions, and with the call-weighted value under a `_pooled`
suffix.

| Field | Type | Description |
|---|---|---|
| `ap` | numeric | Mean over regions of the per-region average precision. `NA` when no valid fits. |
| `ap_n_fits` | integer | Regions contributing to `ap`. A failed fit and a region with no causal variants both contribute none. |
| `fdr_power_curve` | data.frame or NULL | One row per threshold. Region means in `fdr`, `precision`, `power`, `recall`; `n_regions` and `n_regions_with_causal`; call-weighted `fdr_pooled`, `precision_pooled`, `power_pooled`, `recall_pooled`; and the `tp`/`fp`/`fn` totals. **NULL when no valid fits** - a method that produced nothing gets no curve, not a zero FDR. |
| `pip_calibration` | data.frame or NULL | `n_pip_cal_bins` rows. `bin`, `bin_lower`, `bin_upper`, `bin_mid`; region means in `mean_pip` and `frac_causal`; `n_regions` per bin; call-weighted `mean_pip_pooled` and `frac_causal_pooled`; and the `n`/`n_causal` totals. NULL when no valid fits. |
| `cs_coverage` | numeric | Mean over regions of the share of that region's credible sets holding a causal variant. `NA` when no region reported a set. |
| `cs_coverage_n` | integer | Regions that reported at least one credible set. |
| `cs_power` | numeric | Mean over regions of the share of that region's causal variants captured by any set. |
| `cs_power_n` | integer | Regions with at least one causal variant. |
| `cs_size_median`, `cs_size_mean` | numeric | Means over regions of that region's own median and mean set size. |
| `cs_coverage_pooled`, `cs_power_pooled`, `cs_size_median_pooled`, `cs_size_mean_pooled` | numeric | The same quantities pooled over sets, and over causal variants for power. |
| `n_cs_reported` | integer | Total credible sets reported across all fits. |
| `runtime_mean` | numeric | Mean runtime in seconds over valid fits. |
| `runtime_sd` | numeric | SD of runtime over valid fits. `NA` if fewer than 2 valid fits. |
| `n_fits` | integer | Fits attempted, including failed. |
| `n_failed` | integer | Fits where `f$error` is non-NULL. |

---

## Saving to disk

When `save = TRUE`, two files are written to `output_dir`:

### `evaluation.rds`

The complete `eval_out` list as described above, serialised with `saveRDS()`.
Load with:

```r
ev <- readRDS("results/my_run/evaluation.rds")
ev$susie$global$ap
ev$susie$global$fdr_power_curve
```

### `evaluation_summary.csv`

A flat CSV with one row per method:

| Column | Description |
|---|---|
| `method` | Method name |
| `ap` | Mean over regions of the per-region average precision |
| `ap_n_fits` | Regions contributing to `ap` |
| `cs_coverage` | Region-weighted credible-set coverage |
| `cs_coverage_n` | Regions that reported at least one set |
| `cs_power` | Region-weighted credible-set power |
| `cs_power_n` | Regions with at least one causal variant |
| `cs_size_median`, `cs_size_mean` | Region-weighted set sizes |
| `n_cs_reported` | Total credible sets reported |
| `runtime_mean`, `runtime_sd` | Runtime in seconds over valid fits |
| `n_fits`, `n_failed` | Fits attempted, and fits that errored |

The curve and calibration tables are not in the CSV - they are data frames, and
live in `evaluation.rds` and in `per_region`.

---

## Quick reference — field names

### Inside `global`

```
$ap                      numeric scalar (mean over regions)
$ap_n_fits               integer scalar
$fdr_power_curve         data.frame(threshold, fdr, power, precision, recall,
                                    n_regions, n_regions_with_causal,
                                    fdr_pooled, precision_pooled,
                                    power_pooled, recall_pooled, tp, fp, fn)
$pip_calibration         data.frame(bin, bin_lower, bin_upper, bin_mid,
                                    mean_pip, frac_causal, n_regions,
                                    mean_pip_pooled, frac_causal_pooled,
                                    n, n_causal)
$cs_coverage             numeric scalar      $cs_coverage_n   integer scalar
$cs_power                numeric scalar      $cs_power_n      integer scalar
$cs_size_median          numeric scalar
$cs_size_mean            numeric scalar
$cs_coverage_pooled      numeric scalar
$cs_power_pooled         numeric scalar
$cs_size_median_pooled   numeric scalar
$cs_size_mean_pooled     numeric scalar
$n_cs_reported           integer scalar
$runtime_mean            numeric scalar
$runtime_sd              numeric scalar
$n_fits                  integer scalar
$n_failed                integer scalar
```

### Inside `per_region`

One row per (method, region, iteration). Three column families, each with its
own aggregation rule:

```
identity   method, scenario_id, region_id, iter
factors    S, phi, p_causal, region_size, min_causal_maf, causal_maf_bin,
           true_annotation_type
status     failed, bad_pip

COUNTS -> sum, then form the rate (call-weighted)
           nsel_at_*, tp_at_*, sum_pip_sel_at_*,
           n_band_*, c_band_*, sum_pip_band_*,
           n_cs, n_cs_with_causal, n_causal_captured, cs_size_sum,
           n_variants, n_causal, sum_pip_total

RATES  -> average over the regions where DEFINED (region-weighted), and
          report how many; NA where undefined, never 0
           fdr_at_*, precision_at_*, power_at_*,
           frac_causal_band_*, mean_pip_band_*,
           cs_coverage, cs_power

SCALARS -> average
           ap, max_pip, cs_size_mean, cs_size_median, runtime_seconds
```

### Working with it

```r
ev <- readRDS("results/my_run/evaluation.rds")
pr <- ev$per_region

# A method's headline AP, and the same thing from the table.
ev$susie$global$ap
mean(pr$ap[pr$method == "susie"], na.rm = TRUE)

# Any axis, or combination of axes.
aggregate(ap ~ method + S, data = pr, FUN = mean)
aggregate(ap ~ method + region_size + phi, data = pr, FUN = mean)
aggregate(ap ~ method + causal_maf_bin, data = pr, FUN = mean)

# A rate: average over the regions where it is defined, and count them.
v <- pr$fdr_at_095[pr$method == "susie"]
c(fdr = mean(v, na.rm = TRUE), n_regions = sum(!is.na(v)))

# A standard error: collapse to one value per iteration FIRST.
per_iter <- aggregate(ap ~ method + iter, data = pr, FUN = mean)
aggregate(ap ~ method, data = per_iter,
          FUN = function(v) sd(v) / sqrt(length(v)))

# The PR curve for one method.
df <- ev$susie$global$fdr_power_curve
plot(df$recall, df$precision, type = "l")
```
