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
# THE TABLE'S COLUMNS COME IN THREE FAMILIES, and the family decides how to
# aggregate:
#
#   COUNTS   nsel_at_*, tp_at_*, sum_pip_sel_at_*, n_band_*, c_band_*,
#            sum_pip_band_*, n_cs, n_cs_with_causal, n_causal_captured,
#            cs_size_sum, n_variants, n_causal, sum_pip_total
#            -> SUM them, then form the rate (the call-weighted rate).
#
#   RATES    fdr_at_*, precision_at_*, power_at_*, frac_causal_band_*,
#            mean_pip_band_*, cs_coverage, cs_power
#            -> AVERAGE them over the regions where they are DEFINED, and
#               report how many that was (the region-weighted rate). They are
#               NA where undefined - a region that selected nothing has no
#               FDR, and scoring it 0 would reward abstention.
#
#   SCALARS  ap, max_pip, cs_size_mean, cs_size_median, runtime_seconds
#            -> AVERAGE them.
#
# AND FOR A STANDARD ERROR, THE UNIT IS THE ITERATION, NOT THE REGION. Each
# region's genotypes, LD matrix and annotation matrix are drawn ONCE and reused
# by every iteration; only the causal-variant selection, the effect sizes and
# the phenotype noise are redrawn. So collapse to one value per `iter` first,
# then take sd/sqrt(n) across those. Measured on 112,500 fits across 5,625
# cells of a completed benchmark, an SE taken across region-rows instead had
# median 0.96 and IQR 0.83-1.08 of the iteration-based one - understating in
# 42% of cells and overstating in the rest, so it errs in no predictable
# direction.
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
#' The columns come in three families, and the family decides how to
#' aggregate. **Counts** (`nsel_at_*`, `tp_at_*`, `sum_pip_sel_at_*`,
#' `n_band_*`, `c_band_*`, `sum_pip_band_*`, `n_cs`, `n_cs_with_causal`,
#' `n_causal_captured`, `cs_size_sum`, `n_variants`, `n_causal`,
#' `sum_pip_total`) are SUMMED, and a rate formed from them afterwards is
#' call-weighted. **Rates** (`fdr_at_*`, `precision_at_*`, `power_at_*`,
#' `frac_causal_band_*`, `mean_pip_band_*`, `cs_coverage`, `cs_power`) are
#' AVERAGED over the regions where they are defined, giving the
#' region-weighted rate; they are `NA` where undefined, never 0, so count the
#' contributing regions alongside the mean. **Scalars** (`ap`, `max_pip`,
#' `cs_size_mean`, `cs_size_median`, `runtime_seconds`) are AVERAGED.
#'
#' For a standard error the unit is the ITERATION, not the region: a region's
#' genotypes, LD matrix and annotation matrix are drawn once and reused by
#' every iteration, so collapse to one value per `iter` before taking
#' `sd()/sqrt(n)`. Measured on 112,500 fits across 5,625 cells of a completed
#' benchmark, an SE taken across region-rows instead had median 0.96 and IQR
#' 0.83-1.08 of the iteration-based one - understating in 42% of cells and
#' overstating in the rest, so it errs in no predictable direction.
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
#' @seealso [evaluate_methods()], whose `global` entry is these metrics over
#'   every fit of a method.
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
