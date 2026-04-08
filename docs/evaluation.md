# Evaluation Methodology

This document describes every step of the evaluation pipeline implemented in
`R/evaluate.R` and `R/plot_results.R`: how raw method outputs are matched to
ground truth, how each metric is computed, how standard errors are derived, how
results are stratified, and what each plot shows.

---

## Contents

1. [Inputs](#inputs)
2. [Annotating fits with ground truth](#annotating-fits-with-ground-truth)
3. [Stratification](#stratification)
4. [Metric computation](#metric-computation)
   - [Power / FDR curve](#power--fdr-curve)
   - [AUPRC](#auprc)
   - [PIP calibration](#pip-calibration)
   - [Credible-set metrics](#credible-set-metrics)
   - [Runtime](#runtime)
5. [Standard error computation](#standard-error-computation)
6. [Output structure](#output-structure)
7. [Saving to disk](#saving-to-disk)
8. [Plots](#plots)
   - [Global PR curve](#global-pr-curve)
   - [Global PIP calibration](#global-pip-calibration)
   - [Global summary table](#global-summary-table)
   - [Stratified PR curve grid](#stratified-pr-curve-grid)
   - [Stratified PIP calibration grid](#stratified-pip-calibration-grid)
   - [Metrics vs stratification variable](#metrics-vs-stratification-variable)
9. [Quick reference — field names](#quick-reference--field-names)

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

## Stratification

Metrics are computed four times per method, each over a different subset of fits:

| Stratum | Key in output | Fits included |
|---|---|---|
| **Global** | `$global` | All fits across all scenarios and regions |
| **By S** | `$by_S` | Fits grouped by the number of causal variants (`f$S`). One entry per unique `S` value. |
| **By phi** | `$by_phi` | Fits grouped by the proportion of variance explained (`f$phi`). One entry per unique `phi` value. |
| **By p_causal** | `$by_p_causal` | Fits grouped by the sparse-component proportion (`f$p_causal`). One entry per unique `p_causal` value. `NULL` for the `sparse` model. |

For each stratification variable, `.stratify_metrics()` identifies all unique
non-`NA` values of that variable across the annotated fits, then calls
`.metrics_with_se()` on the subset of fits matching each value.

When a stratification variable is `NULL` or `NA` on all fits (e.g. `p_causal`
for a `sparse` model), the stratum is returned as `NULL`.

**What "pooling" means:** within each stratum, all (PIP, is-causal) pairs from
all fits in that stratum are concatenated into a single vector before metric
computation. Metrics are therefore computed on this aggregate — they reflect the
average behaviour across regions and, within a stratification level, across all
values of the other parameters and all replicates.

---

## Metric computation

All metrics are implemented in `.metrics_for_fits()`. Fits where `f$error` is
non-NULL are counted in `n_failed` but excluded from all computations. If no
valid (non-failed) fits remain, all metrics are returned as `NA` or `NULL`.

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
`power_se` and `precision_se`.

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

### AUPRC

The area under the precision-recall curve is computed from the
`precision` and `recall` columns of the FDR/power curve data frame using
the **trapezoidal rule**:

```
AUPRC = sum over adjacent threshold pairs of:
        Δrecall × (precision_left + precision_right) / 2
```

Before summation, rows are sorted by increasing `recall`. At ties in `recall`,
the row with higher `precision` is placed first (so that precision is not
artificially depressed at the same recall level).

The result is a single scalar in [0, 1]. A perfect method would achieve
AUPRC = 1 (all causal variants receive the highest PIPs). A method equivalent
to a random PIP assignment would achieve AUPRC ≈ total_causal / N (the
prevalence of causal variants), which is typically very small (< 0.05).

AUPRC is preferred over AUROC for fine-mapping because the class imbalance
is extreme: in a typical region with 100–500 variants and S = 1–5 causal
variants, the positive class constitutes < 5% of variants, and AUROC is
insensitive to performance in this regime.

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

After SE computation, one additional column is present: `frac_causal_se`.

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

After SE computation, `runtime_mean_se` is also present (see below).

---

## Standard error computation

All SE fields are computed by `.metrics_with_se()`, which wraps
`.metrics_for_fits()`. The SE represents variability **across replicates**
(the `iter` dimension of the simulation design).

**Grouping.** The `iter` field on each fit records which independent replicate
it belongs to (1 through `n_iter`). Fits are grouped by their `iter` value.
Within a given stratum (e.g. all fits with `S = 1`), there is one group per
unique `iter` value found in that stratum.

**Per-replicate metrics.** For each replicate `r`, `.metrics_for_fits()` is
called on the subset of fits with `iter == r`. This produces one complete set
of metrics pooled across all regions and parameter values (S, phi, p_causal)
within that replicate and stratum.

**SE formula.** Let `x_1, ..., x_K` be the per-replicate values of a scalar
metric (some may be `NA`). Let `K*` = number of non-NA values. Then:

```
SE = sd(x_1, ..., x_K, na.rm = TRUE) / sqrt(K*)
```

If `K* < 2`, the SE is `NA` (cannot estimate from fewer than two replicates).

**Which fields have SE.** The following scalar metrics have a paired `_se` field:

| Metric field | SE field |
|---|---|
| `auprc` | `auprc_se` |
| `cs_coverage` | `cs_coverage_se` |
| `cs_power` | `cs_power_se` |
| `cs_size_median` | `cs_size_median_se` |
| `cs_size_mean` | `cs_size_mean_se` |
| `runtime_mean` | `runtime_mean_se` |

In addition, two data frames carry per-row SE fields:

- **`fdr_power_curve`** gains `power_se` and `precision_se`: the SE of power
  (recall) and precision at each threshold across replicates. These are computed
  by running `.metrics_for_fits()` per replicate and collecting the
  `power`/`precision` vectors into a matrix (rows = thresholds,
  columns = replicates), then applying the SE formula row-wise.

- **`pip_calibration`** gains `frac_causal_se`: the SE of `frac_causal` in
  each bin across replicates. Computed similarly as a matrix of per-replicate
  `frac_causal` vectors.

**Interpretation.** Because grouping is by `iter` within the stratum, the SE
reflects variability due to both random noise in the phenotype simulation and
stochastic variation in the fine-mapping algorithm (where applicable). It does
**not** reflect uncertainty from the choice of S or phi (those are design
parameters, not sources of noise). For the global stratum, per-replicate
pooling spans all S and phi values, so the global SE also captures
cross-parameter variability — it is larger than the within-parameter SE.

---

## Output structure

`evaluate_methods()` returns a named list with one element per evaluated method
plus three top-level metadata fields:

```
eval_out
├── methods_evaluated       character vector of method names
├── simulation_params       copy of simulation$params
├── pip_thresholds_used     the pip_thresholds vector used
└── <method>                one entry per method (e.g. "susie", "abf")
    ├── global              metrics pooled over all scenarios × regions
    ├── by_S                named list; one entry per unique S value
    │   ├── "1"             metrics for fits with S = 1
    │   ├── "2"             metrics for fits with S = 2
    │   └── ...
    ├── by_phi              named list; one entry per unique phi value
    │   ├── "0.1"
    │   └── ...
    └── by_p_causal         named list (sparse_inf model only); NULL otherwise
        ├── "0.2"
        └── ...
```

Each stratum object (global, or an entry in by_S / by_phi / by_p_causal)
contains:

| Field | Type | Description |
|---|---|---|
| `fdr_power_curve` | data.frame or NULL | 201 rows (one per threshold). Columns: `threshold`, `tp`, `fp`, `fn`, `fdr`, `power`, `precision`, `recall`, `power_se`, `precision_se`. NULL when no valid fits exist. |
| `auprc` | numeric | Area under precision-recall curve. NA when no valid fits. |
| `auprc_se` | numeric | SE of AUPRC across replicates. NA when n_iter < 2. |
| `pip_calibration` | data.frame or NULL | `n_pip_cal_bins` rows. Columns: `bin`, `bin_lower`, `bin_upper`, `bin_mid`, `n`, `n_causal`, `mean_pip`, `frac_causal`, `frac_causal_se`. NULL when no valid fits. |
| `cs_coverage` | numeric | Proportion of CSs containing ≥1 causal variant. NA if no CSs reported. |
| `cs_coverage_se` | numeric | SE of cs_coverage across replicates. |
| `cs_power` | numeric | Proportion of causal variants captured by any CS. |
| `cs_power_se` | numeric | SE of cs_power across replicates. |
| `cs_size_median` | numeric | Median variants per CS. NA if no CSs reported. |
| `cs_size_median_se` | numeric | SE of cs_size_median across replicates. |
| `cs_size_mean` | numeric | Mean variants per CS. NA if no CSs reported. |
| `cs_size_mean_se` | numeric | SE of cs_size_mean across replicates. |
| `n_cs_reported` | integer | Total CSs reported across all fits in stratum. |
| `runtime_mean` | numeric | Mean runtime in seconds over valid fits. |
| `runtime_mean_se` | numeric | SE of runtime_mean across replicates. |
| `runtime_sd` | numeric | SD of runtime over valid fits. NA if n_valid < 2. |
| `n_fits` | integer | Total fits attempted (including failed). |
| `n_failed` | integer | Fits where `f$error` is non-NULL. |

---

## Saving to disk

When `save = TRUE`, two files are written to `output_dir`:

### `evaluation.rds`

The complete `eval_out` list as described above, serialised with `saveRDS()`.
Load with:

```r
ev <- readRDS("results/my_run/evaluation.rds")
ev$susie$global$auprc
ev$susie$by_S[["1"]]$fdr_power_curve
```

### `evaluation_summary.csv`

A flat CSV with one row per method and the following columns. All `_se` columns
are NA when `n_iter < 2`.

| Column | Description |
|---|---|
| `method` | Method name |
| `auprc` | Global AUPRC |
| `auprc_se` | SE of AUPRC |
| `cs_coverage` | Global CS coverage |
| `cs_coverage_se` | SE of CS coverage |
| `cs_power` | Global CS power |
| `cs_power_se` | SE of CS power |
| `cs_size_median` | Global median CS size |
| `cs_size_median_se` | SE of median CS size |
| `cs_size_mean` | Global mean CS size |
| `cs_size_mean_se` | SE of mean CS size |
| `n_cs_reported` | Total CSs reported globally |
| `runtime_mean` | Mean runtime (s) |
| `runtime_mean_se` | SE of mean runtime |
| `runtime_sd` | SD of runtime |
| `n_fits` | Total fits attempted |
| `n_failed` | Fits that errored |

---

## Plots

`plot_results()` writes a multi-page PDF. The number of pages depends on the
number of methods with successful fits and the number of unique parameter values.
For the `sparse` model with M successful methods, K_S unique S values, and
K_phi unique phi values, the PDF has 3 + 3×K_S + 3×K_phi pages. For
`sparse_inf` with K_pc unique p_causal values, 3 additional pages are added.

Call:

```r
source("R/plot_results.R")
plot_results(eval_out, output_file = "results/my_run/results.pdf")
```

Arguments:

| Argument | Default | Meaning |
|---|---|---|
| `eval_out` | — | Output of `evaluate_methods()` |
| `output_file` | `"results/evaluation.pdf"` | Path of PDF to write |
| `methods` | `NULL` (all) | Character vector to restrict which methods are plotted |
| `verbose` | `TRUE` | Print progress messages |

Methods where all fits failed (i.e. `fdr_power_curve` is NULL) are excluded
from PR curve and calibration plots but are included in the summary table and
in the metrics line plots (where they appear as NA points).

---

### Global PR curve

**Page 1.** One panel. x-axis: recall (= power), y-axis: precision (= 1 − FDR).
Each method is a separate coloured line. The curve is derived from the global
`fdr_power_curve` data frame: each threshold gives one (recall, precision) point.

The curve starts at (recall=1, precision=total_causal/N) at threshold=0 (all
variants called) and moves towards (recall=0, precision=1) at threshold=1 (no
variants called). A method with strong discriminative power has a curve that
stays close to the top-right corner of the plot.

No error bands are shown on this plot. SE on precision and power is available
in the evaluation object (`precision_se`, `power_se` columns) but is omitted
from this figure for clarity.

---

### Global PIP calibration

**Page 2.** Faceted by method (up to 3 columns). x-axis: `mean_pip` (expected
causal probability), y-axis: `frac_causal` (observed fraction of causal
variants). Each point is one of the `n_pip_cal_bins` (default 10) PIP bins.
Error bars show ±1 SE (`frac_causal_se`) on the y-axis. The dashed diagonal
line is the identity `y = x`; a perfectly calibrated method has all points on
this line.

Bins with no variants (`n = 0`) are omitted (they carry `NA` for `mean_pip`
and `frac_causal`).

---

### Global summary table

**Page 3.** A formatted table with one row per method (including failed methods).
Columns: Method, AUPRC, CS Coverage, CS Power, Median CS Size, Runtime (s),
n_fits, n_failed. Where n_iter ≥ 2, values are formatted as `mean ± SE`.
Where n_iter < 2 or the method failed entirely, values are shown as `NA`.

---

### Stratified PR curve grid

**Pages 4, 7, (10).** One page per stratification variable (by S, by phi, by
p_causal). A `facet_grid(method ~ stratum_value)` layout: rows are methods
(only those with at least one successful fit), columns are the unique values of
the stratification variable.

Each panel shows the PR curve for that (method, parameter-value) combination,
pooled across all other parameters and all replicates. The curve is computed
from the corresponding stratum's `fdr_power_curve` object.

The layout for by_S with 4 methods and 3 S values produces a 4×3 grid of PR
curves (12 panels). Methods with no successful fits in a given stratum show an
empty panel.

---

### Stratified PIP calibration grid

**Pages 5, 8, (11).** Same `facet_grid(method ~ stratum_value)` layout as the
PR curve grids, but showing calibration plots. Each panel: x = `mean_pip`,
y = `frac_causal` with ±1 SE error bars, dashed identity line.

---

### Metrics vs stratification variable

**Pages 6, 9, (12).** One page per stratification variable. A `facet_wrap` of
5 panels (one per metric: AUPRC, CS Coverage, CS Power, Median CS Size,
Runtime). In each panel:

- x-axis: the stratification variable value (e.g. S = 1, 2, 3)
- y-axis: the metric value (scale is free across panels)
- one coloured line + points per method
- error bars: ±1 SE across replicates

All methods (including those that failed entirely) are included, but methods
with all-NA values produce invisible lines. The y-axis scale is free across
panels (each metric uses its own range).

This plot is the primary tool for assessing how each method's performance
changes as the simulation difficulty changes (more causal variants, lower PVE,
or smaller sparse proportion).

---

## Quick reference — field names

### Within a stratum object

```
$fdr_power_curve         data.frame(threshold, tp, fp, fn, fdr, power,
                                    precision, recall, power_se, precision_se)
$auprc                   numeric scalar
$auprc_se                numeric scalar (NA if n_iter < 2)
$pip_calibration         data.frame(bin, bin_lower, bin_upper, bin_mid,
                                    n, n_causal, mean_pip, frac_causal,
                                    frac_causal_se)
$cs_coverage             numeric scalar
$cs_coverage_se          numeric scalar
$cs_power                numeric scalar
$cs_power_se             numeric scalar
$cs_size_median          numeric scalar
$cs_size_median_se       numeric scalar
$cs_size_mean            numeric scalar
$cs_size_mean_se         numeric scalar
$n_cs_reported           integer scalar
$runtime_mean            numeric scalar
$runtime_mean_se         numeric scalar
$runtime_sd              numeric scalar
$n_fits                  integer scalar
$n_failed                integer scalar
```

### Accessing stratified results

```r
ev <- readRDS("results/my_run/evaluation.rds")

# Global metrics for SuSiE
ev$susie$global$auprc
ev$susie$global$auprc_se

# By S
ev$susie$by_S[["1"]]$auprc          # AUPRC at S = 1
ev$susie$by_S[["2"]]$cs_coverage    # CS coverage at S = 2

# By phi
ev$abf$by_phi[["0.2"]]$cs_power

# By p_causal (sparse_inf model only)
ev$susie_inf$by_p_causal[["0.4"]]$auprc

# PR curve at S = 1
df <- ev$susie$by_S[["1"]]$fdr_power_curve
plot(df$recall, df$precision, type = "l")

# Calibration globally
cal <- ev$abf$global$pip_calibration
plot(cal$mean_pip, cal$frac_causal)
abline(0, 1, lty = 2)
```
