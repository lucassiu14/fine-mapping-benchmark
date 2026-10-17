# =============================================================================
# R/evaluate_per_region.R
#
# PER-REGION EVALUATION: one row per (method, scenario, region), so a user can
# aggregate the benchmark along whatever axes they care about instead of being
# limited to the strata `evaluate_methods()` happens to return.
#
# THE AGGREGATION RULE, which this file exists to make possible:
#
#     Rates are pooled from COUNTS. Per-locus scalars are AVERAGED.
#
#   - Average precision is a PER-LOCUS quantity. The AP of two regions' merged
#     ranked list is not the mean of their APs, because merging lets one
#     region's false positives outrank the other region's true positive. A
#     pooled ranking describes an analysis nobody performs, so AP is computed
#     per region and averaged - at every level.
#
#   - FDR, power and calibration are properties of a POPULATION OF CALLS. A rate
#     computed on one region is unestimable, because a region typically carries
#     0-2 selections. So this table stores the COUNTS (`tp_at_*`, `nsel_at_*`,
#     `n_band_*`, `c_band_*`), and the rate is formed after summing them over
#     whatever set of regions the user is aggregating.
#
#   Never average rates. Never pool rankings. Both change the answer, not just
#   its precision.
#
# Two conveniences follow from that rule and are provided here:
# `aggregate_region_metrics()` applies it correctly for a chosen grouping, and
# `add_rate_metrics()` forms the rates once the counts are summed.
# =============================================================================


#' Exact rank-based average precision for one fit
#'
#' AP is a function of the ranking alone, so it is evaluated at the DISTINCT
#' OBSERVED PIP values rather than on a fixed threshold grid. A grid treats
#' every variant inside one band as entering the selection simultaneously, so a
#' causal variant is credited with the precision computed after every
#' non-causal in its band has also been admitted. A causal variant at PIP
#' 0.0032 with true rank 61 of 1000 has exact AP 1/61 = 0.0164 but grid AP
#' 1/1000 = 0.001.
#'
#' That error is zero at strong signal and grows as PIPs compress toward zero,
#' so on a grid it becomes a heritability-dependent bias in a metric whose
#' heritability dependence one is usually trying to measure. A constant bias
#' would be harmless; this one is not.
#'
#' Tied PIPs are scored as a block, which is the conservative convention: a
#' method that assigns the same PIP to a causal and a non-causal variant gets
#' the precision of the whole block, not the benefit of an arbitrary tie-break.
#'
#' @param pip_vals Numeric vector of PIPs for ONE fit.
#' @param is_causal Logical vector, the same length as `pip_vals`.
#' @return Scalar average precision, or `NA_real_` when the region has no
#'   causal variants.
#' @export
compute_ap_exact <- function(pip_vals, is_causal) {
  S <- sum(is_causal)
  if (S == 0L || length(pip_vals) == 0L) return(NA_real_)
  if (!all(is.finite(pip_vals))) return(NA_real_)
  ord  <- order(pip_vals, decreasing = TRUE)
  p    <- pip_vals[ord]
  y    <- as.integer(is_causal[ord])
  last <- c(p[-1L] != p[-length(p)], TRUE)   # TRUE at the end of each tied block
  k    <- which(last)                        # n selected at each distinct value
  tp   <- cumsum(y)[last]
  prec <- tp / k
  rec  <- tp / S
  sum(diff(c(0, rec)) * prec)                # step integral, not trapezoid
}


# The calibration bands. Bands rather than the full bin table: a single region
# has too few variants for a 10-bin table to mean anything, and the bands can be
# summed into whatever binning the user wants afterwards.
.PIP_BANDS <- list(lo  = c(0,   0.1),
                   mid = c(0.1, 0.5),
                   hi  = c(0.5, 0.9),
                   top = c(0.9, 1 + 1e-9))


# Everything computed from one (pip, truth) pair. This is the only place raw
# PIPs are touched; every column downstream is a function of these.
.region_metrics <- function(pip, causal_idx, n_var, credible_sets,
                            count_thresholds) {
  y <- logical(n_var)
  y[causal_idx] <- TRUE
  S <- sum(y)

  out <- list(
    ap            = compute_ap_exact(pip, y),
    n_variants    = n_var,
    n_causal      = S,
    sum_pip_total = sum(pip),
    max_pip       = if (n_var > 0L) max(pip) else NA_real_,
    top_pip_is_causal = if (n_var > 0L) y[which.max(pip)] else NA
  )

  # Counts at the decision thresholds, and THIS REGION's rates. The counts are
  # kept so a caller can also form the call-weighted rate over a set of
  # regions; the rates are what "how trustworthy are my calls in this region"
  # actually asks.
  #
  # fdr and precision are NA when the region selected nothing - NOT 0 and 1.
  # The convention that FDR is 0 for an empty selection is harmless when
  # pooling, since such a region contributes no calls, but it REWARDS
  # ABSTENTION under averaging: a method that selects nothing anywhere would
  # average to a perfect FDR. Those regions are excluded from the mean and
  # counted instead.
  for (t in count_thresholds) {
    tag  <- .thresh_tag(t)
    sel  <- pip >= t
    nsel <- sum(sel)
    tp   <- sum(y & sel)
    out[[paste0("nsel_at_",        tag)]] <- nsel
    out[[paste0("tp_at_",          tag)]] <- tp
    out[[paste0("sum_pip_sel_at_", tag)]] <- sum(pip[sel])
    out[[paste0("fdr_at_",         tag)]] <- if (nsel > 0L) (nsel - tp) / nsel
                                             else NA_real_
    out[[paste0("precision_at_",   tag)]] <- if (nsel > 0L) tp / nsel
                                             else NA_real_
    out[[paste0("power_at_",       tag)]] <- if (S > 0L) tp / S else NA_real_
  }

  # Calibration band counts and this region's observed frequency per band. A
  # band with no variant in this region is NA, not 0.
  for (nm in names(.PIP_BANDS)) {
    b   <- .PIP_BANDS[[nm]]
    sel <- pip >= b[1] & pip < b[2]
    nb  <- sum(sel)
    out[[paste0("n_band_",           nm)]] <- nb
    out[[paste0("c_band_",           nm)]] <- sum(y & sel)
    out[[paste0("sum_pip_band_",     nm)]] <- sum(pip[sel])
    out[[paste0("frac_causal_band_", nm)]] <- if (nb > 0L) sum(y & sel) / nb
                                              else NA_real_
    out[[paste0("mean_pip_band_",    nm)]] <- if (nb > 0L) sum(pip[sel]) / nb
                                              else NA_real_
  }

  # Credible sets. `cs_hits` is a Bernoulli per SET and `causal_captured` a
  # Bernoulli per CAUSAL VARIANT, so both are counts: coverage and power are
  # formed from them after summing, not by averaging per-region rates.
  if (length(credible_sets) == 0L) {
    out$n_cs              <- 0L
    out$n_cs_with_causal  <- 0L
    out$n_causal_captured <- 0L
    out$cs_size_sum       <- 0L
    out$cs_size_mean      <- NA_real_
    out$cs_size_median    <- NA_real_
    out$cs_coverage       <- NA_real_   # no sets: coverage undefined, not 0
    out$cs_power          <- if (S > 0L) 0 else NA_real_
  } else {
    sizes <- vapply(credible_sets, length, integer(1))
    hits  <- vapply(credible_sets, function(cs) any(cs %in% causal_idx), logical(1))
    out$n_cs              <- length(credible_sets)
    out$n_cs_with_causal  <- sum(hits)
    out$n_causal_captured <- sum(causal_idx %in% unique(unlist(credible_sets)))
    out$cs_size_sum       <- sum(sizes)
    out$cs_size_mean      <- mean(sizes)
    out$cs_size_median    <- stats::median(sizes)
    out$cs_coverage       <- mean(hits)
    out$cs_power          <- if (S > 0L) out$n_causal_captured / S else NA_real_
  }
  out
}


.thresh_tag <- function(t) gsub("[.]", "", format(t, trim = TRUE, scientific = FALSE))


#' Per-region evaluation metrics
#'
#' One row per (method, scenario, region) fit, with the design factors that
#' describe the region alongside the metrics computed from its PIPs. This is
#' the table to aggregate from: nothing in it is pooled across regions, so you
#' choose the axes.
#'
#' Average precision is already a finished per-region number, so average it.
#' FDR, power, coverage and calibration are not: they are rates, and a rate on
#' the 0-2 selections a single region carries is unestimable. Their ingredients
#' are stored here as counts (`nsel_at_*`, `tp_at_*`, `n_band_*`, `c_band_*`,
#' `n_cs`, `n_cs_with_causal`, `n_causal_captured`), to be summed over the
#' regions you are aggregating before the rate is formed.
#' [aggregate_region_metrics()] does this for you.
#'
#' A failed fit contributes an all-`NA` metric row rather than being dropped.
#' Its absence has to stay visible: a method that cannot run on part of the
#' grid is not missing at random, and silently dropping it changes the method
#' set per cell.
#'
#' @param simulation A simulation object from [run_simulation()].
#' @param results A results object from [run_methods()].
#' @param count_thresholds Numeric. PIP thresholds at which selection counts
#'   are stored. Default: 0.5, 0.8, 0.9, 0.95, 0.99.
#' @param verbose Logical. Report progress per method. Default: TRUE.
#'
#' @return A data frame with one row per fit. Design columns: `method`,
#'   `scenario_id`, `region_id`, `S`, `phi`, `p_causal`, `region_size`,
#'   `causal_maf_bin`, `min_causal_maf`, `true_annotation_type`, `iter`.
#'   Status columns: `failed`, `bad_pip`. Metric columns as described above,
#'   plus `runtime_seconds`.
#'
#' @seealso [aggregate_region_metrics()] to aggregate it correctly,
#'   [evaluate_methods()] for the pre-aggregated summary.
#' @export
per_region_metrics <- function(simulation,
                               results,
                               count_thresholds = c(0.5, 0.8, 0.9, 0.95, 0.99),
                               verbose          = TRUE) {
  stopifnot(
    "simulation must contain 'scenarios'" = !is.null(simulation$scenarios),
    "results must contain 'methods_run'"  = !is.null(results$methods_run),
    "count_thresholds must be in [0, 1]"  =
      is.numeric(count_thresholds) && all(count_thresholds >= 0 &
                                          count_thresholds <= 1)
  )
  count_thresholds <- sort(unique(count_thresholds))

  methods <- intersect(results$methods_run, names(results))
  if (length(methods) == 0L)
    stop("No recognised methods found in results.", call. = FALSE)

  rows <- list()
  k    <- 0L

  for (method in methods) {
    if (verbose) message(sprintf("Per-region metrics for %s...", method))
    fits <- .annotate_fits_with_truth(results[[method]]$results, simulation)

    for (f in fits) {
      n_var  <- length(f$pip)
      failed <- !is.null(f$error)
      # is.finite(), not is.na(): exp() overflow gives Inf and Inf/sum gives
      # NaN, so both must be caught. A partially non-finite PIP vector is the
      # dangerous case - it has undefined total mass and an undefined ranking,
      # so no metric from it is trustworthy. The fit is marked failed, not
      # repaired: imputing the non-finite entries to zero would invent
      # posterior mass and flatter precisely the method that overflowed.
      bad_pip <- !failed && (n_var == 0L || !all(is.finite(f$pip)))

      met <- if (failed || bad_pip) {
        proto <- .region_metrics(rep(0, max(n_var, 1L)), integer(0),
                                 max(n_var, 1L), list(), count_thresholds)
        lapply(proto, function(z) NA)
      } else {
        .region_metrics(f$pip, f$causal_indices, n_var,
                        f$credible_sets %||% list(), count_thresholds)
      }

      k <- k + 1L
      rows[[k]] <- c(
        list(method               = method,
             scenario_id          = f$scenario_id,
             region_id            = f$region_id,
             iter                 = f$iter %||% NA_integer_,
             S                    = f$S %||% NA_real_,
             phi                  = f$phi %||% NA_real_,
             p_causal             = f$p_causal %||% NA_real_,
             region_size          = f$region_size %||% NA_integer_,
             causal_maf_bin       = f$causal_maf_bin %||% NA_character_,
             min_causal_maf       = f$min_causal_maf %||% NA_real_,
             true_annotation_type = f$true_annotation_type %||% NA_character_,
             failed               = failed,
             bad_pip              = bad_pip,
             runtime_seconds      = f$runtime_seconds %||% NA_real_),
        met)
    }
  }

  if (k == 0L)
    stop("No fits found in results.", call. = FALSE)

  out <- do.call(rbind, lapply(rows, function(r)
    as.data.frame(r, stringsAsFactors = FALSE)))
  attr(out, "count_thresholds") <- count_thresholds
  if (verbose)
    message(sprintf("  %d region-fits across %d methods",
                    nrow(out), length(methods)))
  out
}


#' Form rate metrics from pooled counts
#'
#' Call this on a table whose count columns have already been summed over the
#' regions being aggregated - never on a single region, where a rate on 0-2
#' selections is unestimable.
#'
#' These are the CALL-WEIGHTED rates: a region that makes many calls dominates
#' them. They are written with a `_pooled` suffix, because the rate normally
#' reported is the region-weighted one - the mean over regions of each region's
#' own rate, which is what [aggregate_region_metrics()] puts in the unsuffixed
#' columns. Fine-mapping is applied to a region, so "how trustworthy are the
#' calls in a region" is usually the question; "what fraction of all calls made
#' anywhere were wrong" is a different one.
#'
#' @param d A data frame of summed counts, as [aggregate_region_metrics()]
#'   produces.
#' @param count_thresholds The thresholds whose counts `d` carries.
#' @return `d` with the call-weighted rate columns added: `fdr_at_*_pooled`,
#'   `precision_at_*_pooled`, `power_at_*_pooled`, `cs_coverage_pooled`,
#'   `cs_power_pooled`, `cs_size_mean_pooled`, and per-band
#'   `frac_causal_band_*_pooled` against `mean_pip_band_*_pooled`.
#' @export
add_rate_metrics <- function(d, count_thresholds = c(0.5, 0.8, 0.9, 0.95, 0.99)) {
  for (t in sort(unique(count_thresholds))) {
    tag  <- .thresh_tag(t)
    nsel <- d[[paste0("nsel_at_", tag)]]
    tp   <- d[[paste0("tp_at_",   tag)]]
    if (is.null(nsel) || is.null(tp)) next
    fp <- nsel - tp
    # FDR is 0 by convention when nothing is selected, which is harmless here -
    # a region with no calls contributes none to a call-weighted rate - but is
    # exactly why the region-weighted version must not use it. nsel is kept
    # beside the rate either way: never report one without the other.
    d[[paste0("fdr_at_",       tag, "_pooled")]] <-
      ifelse(nsel > 0, fp / nsel, 0)
    d[[paste0("precision_at_", tag, "_pooled")]] <-
      ifelse(nsel > 0, tp / nsel, NA_real_)
    d[[paste0("power_at_",     tag, "_pooled")]] <-
      ifelse(d$n_causal > 0, tp / d$n_causal, NA_real_)
  }

  if (!is.null(d$n_cs))
    d$cs_coverage_pooled <- ifelse(d$n_cs > 0, d$n_cs_with_causal / d$n_cs,
                                   NA_real_)
  if (!is.null(d$n_causal_captured))
    d$cs_power_pooled <- ifelse(d$n_causal > 0,
                                d$n_causal_captured / d$n_causal, NA_real_)
  if (!is.null(d$cs_size_sum))
    d$cs_size_mean_pooled <- ifelse(d$n_cs > 0, d$cs_size_sum / d$n_cs,
                                    NA_real_)

  for (nm in names(.PIP_BANDS)) {
    n <- d[[paste0("n_band_", nm)]]
    cc <- d[[paste0("c_band_", nm)]]
    sp <- d[[paste0("sum_pip_band_", nm)]]
    if (is.null(n)) next
    d[[paste0("frac_causal_band_", nm, "_pooled")]] <-
      ifelse(n > 0, cc / n, NA_real_)
    d[[paste0("mean_pip_band_",    nm, "_pooled")]] <-
      ifelse(n > 0, sp / n, NA_real_)
  }
  d
}


#' Aggregate per-region metrics along chosen axes
#'
#' Applies the rule the per-region table exists for: COUNT columns are summed,
#' per-region SCALAR columns are averaged, and the rates are formed from the
#' summed counts afterwards. Aggregating the table yourself with `aggregate()`
#' is fine as long as you follow that rule; averaging a rate or pooling a
#' ranking will give a different answer, not merely a noisier one.
#'
#' @param per_region A data frame from [per_region_metrics()].
#' @param by Character vector of column names to group by. `"method"` is always
#'   included. Default: `"method"` alone, i.e. one row per method.
#' @param na.rm Logical. Drop failed fits from the scalar averages. Counts from
#'   failed fits are `NA` and are dropped either way. Default: TRUE.
#'
#' @return A data frame with one row per group: the summed counts, the averaged
#'   scalars, the rates formed from the counts, and `n_fits` / `n_failed`.
#'
#' @examples
#' \dontrun{
#' pr <- per_region_metrics(sim, res)
#'
#' # One row per method.
#' aggregate_region_metrics(pr)
#'
#' # By anything the table carries - including combinations evaluate_methods()
#' # does not offer.
#' aggregate_region_metrics(pr, by = c("region_size", "phi"))
#' aggregate_region_metrics(pr, by = c("true_annotation_type", "S"))
#'
#' # Or aggregate it yourself; just sum the counts and average the scalars.
#' }
#' @export
aggregate_region_metrics <- function(per_region, by = "method", na.rm = TRUE) {
  stopifnot(is.data.frame(per_region))
  by <- unique(c("method", setdiff(by, "method")))
  missing_by <- setdiff(by, names(per_region))
  if (length(missing_by))
    stop("not a column of per_region: ", paste(missing_by, collapse = ", "),
         call. = FALSE)

  thresholds <- attr(per_region, "count_thresholds") %||%
    c(0.5, 0.8, 0.9, 0.95, 0.99)
  tags <- vapply(sort(unique(thresholds)), .thresh_tag, character(1))

  count_cols <- c(paste0("nsel_at_",        tags),
                  paste0("tp_at_",          tags),
                  paste0("sum_pip_sel_at_", tags),
                  paste0("n_band_",         names(.PIP_BANDS)),
                  paste0("c_band_",         names(.PIP_BANDS)),
                  paste0("sum_pip_band_",   names(.PIP_BANDS)),
                  "n_variants", "n_causal", "sum_pip_total",
                  "n_cs", "n_cs_with_causal", "n_causal_captured",
                  "cs_size_sum")
  # Per-region RATES are averaged over the regions where they are defined, and
  # the count of contributing regions is reported with each. A region that
  # selected nothing has an undefined FDR, not an FDR of zero.
  rate_cols <- c(paste0("fdr_at_",           tags),
                 paste0("precision_at_",     tags),
                 paste0("power_at_",         tags),
                 paste0("frac_causal_band_", names(.PIP_BANDS)),
                 paste0("mean_pip_band_",    names(.PIP_BANDS)),
                 "cs_coverage", "cs_power")
  scalar_cols <- c("ap", "max_pip", "cs_size_mean", "cs_size_median",
                   "runtime_seconds")

  count_cols  <- intersect(count_cols,  names(per_region))
  rate_cols   <- intersect(rate_cols,   names(per_region))
  scalar_cols <- intersect(scalar_cols, names(per_region))

  # NA-as-level before grouping: aggregate(), split() and table() all drop NA
  # groups silently, and p_causal is NA on every sparse row while
  # enrichment-related columns are NA on every unannotated row. Grouping on the
  # raw columns would delete those arms entirely.
  keys <- lapply(by, function(b) {
    v <- per_region[[b]]
    factor(ifelse(is.na(v), ".NA", as.character(v)))
  })
  names(keys) <- by
  grp <- interaction(keys, drop = TRUE, sep = "\r")

  pick <- function(cols, fn) {
    if (!length(cols)) return(NULL)
    m <- vapply(cols, function(cl)
      vapply(split(per_region[[cl]], grp), function(v) fn(v), numeric(1)),
      numeric(nlevels(grp)))
    if (nlevels(grp) == 1L) m <- matrix(m, nrow = 1L, dimnames = list(NULL, cols))
    as.data.frame(m)
  }

  defined_mean <- function(v) {
    v <- as.numeric(v)
    v <- v[is.finite(v)]
    if (!length(v)) NA_real_ else mean(v)
  }
  sums  <- pick(count_cols, function(v) sum(as.numeric(v), na.rm = TRUE))
  rates <- pick(rate_cols,  defined_mean)
  means <- pick(scalar_cols, function(v) {
    v <- as.numeric(v)
    if (na.rm) v <- v[!is.na(v)]
    if (!length(v)) NA_real_ else mean(v)
  })
  # How many regions each averaged rate rests on. Mandatory alongside it: at a
  # high PIP threshold most regions select nothing, so a mean FDR can sit on a
  # handful of regions out of thousands.
  rate_n <- pick(rate_cols, function(v) sum(is.finite(as.numeric(v))))
  if (!is.null(rate_n)) names(rate_n) <- paste0(names(rate_n), "_n")

  lev <- do.call(rbind, strsplit(levels(grp), "\r", fixed = TRUE))
  out <- as.data.frame(lev, stringsAsFactors = FALSE)
  names(out) <- by
  for (b in by) out[[b]][out[[b]] == ".NA"] <- NA

  out <- cbind(out, sums, rates, rate_n, means)
  out$n_fits   <- as.integer(table(grp))
  out$n_failed <- as.integer(vapply(split(per_region$failed | per_region$bad_pip,
                                          grp), sum, numeric(1)))
  # ap_n_fits, not n_fits, is AP's denominator: a failed fit and a region with
  # no causal variants both contribute NA.
  out$ap_n_fits <- as.integer(vapply(split(per_region$ap, grp),
                                     function(v) sum(is.finite(v)), numeric(1)))

  out <- add_rate_metrics(out, thresholds)
  rownames(out) <- NULL
  out
}
