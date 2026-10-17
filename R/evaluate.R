# =============================================================================
# evaluate.R
#
# Evaluation metrics for the fine-mapping benchmark.
#
# Usage:
#   eval_out <- evaluate_methods(simulation, results)
#   eval_out <- evaluate_methods(simulation, results, save = TRUE, output_dir = "results/my_run")
#
# Output structure:
#   eval_out$<method>$global          — metrics pooled across all scenarios
#   eval_out$<method>$by_S            — named list, one entry per S value
#   eval_out$<method>$by_phi          — named list, one entry per phi value
#   eval_out$<method>$by_p_causal     — named list (sparse_inf only); NULL otherwise
#   eval_out$<method>$by_region_size  — named list, one entry per nominal region
#                                       length (e.g. "100","200","400","500",
#                                       "1000"). Region size is a real scenario
#                                       axis — fine-mapping difficulty scales with
#                                       it — so metrics must never be pooled
#                                       across lengths.
#   eval_out$<method>$by_causal_maf   — named list, keys "rare" / "low" / "common"
#                                       (bins on the minimum causal-variant MAF
#                                       per region; NULL if no fits have causal
#                                       MAFs available)
#   eval_out$<method>$by_true_annotation_type — named list, keys "none" /
#                                       "binary" / "continuous" / "user_supplied"
#                                       (true annotation regime used in
#                                       simulation; scenario-level, so typically
#                                       single-bin per evaluate_methods() call
#                                       and becomes meaningful only when eval
#                                       objects from differently-annotated sims
#                                       are merged)
#
# Each stratum contains:
# EVERY RATE IS REGION-WEIGHTED: computed within each region, then averaged
# over regions. Fine-mapping is applied to a region, so the question is "how
# trustworthy are the calls in a locus", not "what fraction of all calls made
# anywhere were wrong". Each averaged rate ships with the count of regions it
# rests on, and the call-weighted value under a _pooled suffix.
#
#   ap               — mean of the per-region average precisions
#   ap_n_fits        — regions contributing to it
#   fdr_power_curve  — data.frame: threshold, fdr, power, precision, recall
#                                  (region means), n_regions, the _pooled
#                                  variants, and the tp/fp/fn totals
#   pip_calibration  — data.frame: bin, bin_lower, bin_upper, bin_mid,
#                                  mean_pip, frac_causal (region means),
#                                  n_regions, the _pooled variants, n, n_causal
#   cs_coverage      — mean over regions of the share of that region's CSs
#                      holding ≥1 true causal variant
#   cs_power         — mean over regions of the share of that region's causal
#                      variants captured by any CS
#   cs_size_median   — mean over regions of the region's median CS size
#   cs_size_mean     — mean over regions of the region's mean CS size
#   cs_coverage_n    — regions reporting ≥1 CS (coverage and the sizes rest on
#   cs_power_n         these; power rests on regions with ≥1 causal variant)
#   n_cs_reported    — total CSs reported
#   runtime_mean     — mean runtime in seconds (successful fits only)
#   runtime_sd       — SD of runtime
#   n_fits           — total fits attempted
#   n_failed         — fits that errored
# =============================================================================


#' Evaluate fine-mapping methods against ground truth
#'
#' @param simulation  Output of \code{run_simulation()}.
#' @param results     Output of \code{run_methods()}.
#' @param pip_thresholds Numeric vector of PIP thresholds for the power/FDR
#'   curve. Default: seq(0, 1, by = 0.005).
#' @param n_pip_cal_bins Integer. Number of equal-width bins for PIP
#'   calibration. Default: 10.
#' @param per_region Logical. Also attach the per-region table from
#'   \code{\link{per_region_metrics}} as \code{$per_region}: one row per
#'   (method, scenario, region), nothing pooled, so it can be aggregated along
#'   axes this function does not offer. Default: TRUE.
#' @param save Logical. If TRUE, save the full evaluation object as
#'   \code{evaluation.rds} inside \code{output_dir}, and write a flat
#'   per-method summary table as \code{evaluation_summary.csv}.
#'   Default: FALSE.
#' @param output_dir Character. Directory in which to save results when
#'   \code{save = TRUE}. Created automatically if it does not exist.
#'   Default: \code{"results"}.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A named list with one element per method, each containing
#'   \code{global} and the \code{by_*} strata (see the file header for the
#'   fields in each), plus \code{$per_region}: a single data frame covering
#'   every method, one row per region-fit.
#'
#'   \strong{Every rate is region-weighted}: computed within each region, then
#'   averaged over regions, because fine-mapping is applied to a region. Each
#'   comes with the number of regions it rests on - \code{ap_n_fits},
#'   \code{cs_coverage_n}, \code{cs_power_n}, and \code{n_regions} on the
#'   curve and calibration tables - which must be read alongside it: at a high
#'   PIP threshold most regions select nothing. The call-weighted value, in
#'   which a region making many calls dominates, sits beside each under a
#'   \code{_pooled} suffix.
#'
#'   To aggregate on your own axes, work from \code{$per_region} with
#'   \code{\link{aggregate_region_metrics}}.
#'
#' @seealso \code{\link{per_region_metrics}},
#'   \code{\link{aggregate_region_metrics}}
#' @export
evaluate_methods <- function(simulation,
                             results,
                             pip_thresholds = seq(0, 1, by = 0.005),
                             n_pip_cal_bins = 10L,
                             per_region     = TRUE,
                             save           = FALSE,
                             output_dir     = "results",
                             verbose        = TRUE) {

  stopifnot(
    "simulation must contain 'scenarios'" = !is.null(simulation$scenarios),
    "results must contain 'methods_run'"  = !is.null(results$methods_run)
  )

  methods <- intersect(results$methods_run, names(results))
  if (length(methods) == 0) {
    stop("No recognised methods found in results.", call. = FALSE)
  }

  output <- list()

  for (method in methods) {
    if (verbose) message(sprintf("Evaluating %s...", method))

    fits_raw <- results[[method]]$results
    fits     <- .annotate_fits_with_truth(fits_raw, simulation)

    output[[method]] <- list(
      global      = .metrics_with_se(fits, pip_thresholds, n_pip_cal_bins),
      by_S                    = .stratify_metrics(fits, "S",           pip_thresholds, n_pip_cal_bins),
      by_phi                  = .stratify_metrics(fits, "phi",         pip_thresholds, n_pip_cal_bins),
      by_p_causal             = .stratify_metrics(fits, "p_causal",    pip_thresholds, n_pip_cal_bins),
      by_region_size          = .stratify_metrics(fits, "region_size", pip_thresholds, n_pip_cal_bins),
      by_causal_maf           = .stratify_metrics_by_maf(fits, pip_thresholds, n_pip_cal_bins),
      by_true_annotation_type = .stratify_metrics_by_annotation_type(fits, pip_thresholds, n_pip_cal_bins)
    )
  }

  # The per-region table: one row per (method, scenario, region), nothing pooled
  # across regions, so a caller can aggregate along axes this function does not
  # offer. See per_region_metrics() for the aggregation rule.
  output$per_region <- if (isTRUE(per_region))
    per_region_metrics(simulation, results, verbose = FALSE) else NULL

  output$methods_evaluated   <- methods
  output$simulation_params   <- simulation$params
  output$pip_thresholds_used <- pip_thresholds

  # --- Save to disk -----------------------------------------------------------

  if (save) {
    stopifnot(
      "output_dir must be a single character string" =
        is.character(output_dir) && length(output_dir) == 1
    )
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    # Full evaluation object
    rds_path <- file.path(output_dir, "evaluation.rds")
    saveRDS(output, file = rds_path)
    if (verbose) message(sprintf("  Saved evaluation object: %s", rds_path))

    # Flat summary table — one row per method, global metrics only
    summary_rows <- lapply(methods, function(m) {
      g <- output[[m]]$global
      data.frame(
        method              = m,
        ap                  = g$ap,
        ap_se               = g$ap_se,
        cs_coverage         = g$cs_coverage,
        cs_coverage_se      = g$cs_coverage_se,
        cs_power            = g$cs_power,
        cs_power_se         = g$cs_power_se,
        cs_size_median      = g$cs_size_median,
        cs_size_median_se   = g$cs_size_median_se,
        cs_size_mean        = g$cs_size_mean,
        cs_size_mean_se     = g$cs_size_mean_se,
        n_cs_reported       = g$n_cs_reported,
        runtime_mean        = g$runtime_mean,
        runtime_mean_se     = g$runtime_mean_se,
        runtime_sd          = g$runtime_sd,
        n_fits              = g$n_fits,
        n_failed            = g$n_failed,
        stringsAsFactors    = FALSE
      )
    })
    summary_df  <- do.call(rbind, summary_rows)
    csv_path    <- file.path(output_dir, "evaluation_summary.csv")
    write.csv(summary_df, file = csv_path, row.names = FALSE)
    if (verbose) message(sprintf("  Saved summary table:     %s", csv_path))
  }

  output
}


# =============================================================================
# Internal: annotate fits with ground truth
# =============================================================================

.annotate_fits_with_truth <- function(fits, simulation) {
  lapply(fits, function(f) {
    sc    <- simulation$scenarios[[f$scenario_id]]
    truth <- sc$regions[[f$region_id]]$truth
    f$causal_indices <- truth$causal_indices
    f$n_variants     <- length(f$pip)

    # Nominal region length (the requested per-region p, NOT the realised
    # post-MAF-filter count in params$p_actual). Fine-mapping difficulty scales
    # steeply with region size, so this is a genuine scenario axis: metrics must
    # be stratified by it, never pooled across region lengths. Indexed by
    # region_id; falls back to p_actual then the realised variant count when the
    # nominal vector is unavailable (e.g. legacy simulations).
    p_nominal <- simulation$params$p
    rid       <- f$region_id
    f$region_size <- if (!is.null(p_nominal) && length(p_nominal) >= rid &&
                         !is.na(p_nominal[rid])) {
      as.integer(p_nominal[rid])
    } else {
      p_actual <- simulation$params$p_actual
      if (!is.null(p_actual) && length(p_actual) >= rid && !is.na(p_actual[rid]))
        as.integer(p_actual[rid]) else length(f$pip)
    }

    # Attach causal-variant MAFs and a coarse stratification bin so the
    # evaluator can stratify performance by how rare the rarest causal
    # variant is in the region (the rarest causal is the bottleneck for
    # fine-mapping). MAFs live on genotypes (shared across scenarios
    # for the same region).
    geno_maf <- simulation$genotypes[[f$region_id]]$maf
    if (!is.null(geno_maf) && length(truth$causal_indices) > 0L) {
      f$causal_maf     <- geno_maf[truth$causal_indices]
      f$min_causal_maf <- min(f$causal_maf)
    } else {
      f$causal_maf     <- numeric(0)
      f$min_causal_maf <- NA_real_
    }
    f$causal_maf_bin <- .maf_bin(f$min_causal_maf)

    # True annotation type used by the simulation. Scenario-level (the
    # whole sim was generated with one annotation regime), but attached
    # to each fit so the evaluator can stratify across simulations that
    # have been merged into one results object.
    f$true_annotation_type <- if (!is.null(simulation$params$annotation_type)) {
      as.character(simulation$params$annotation_type)
    } else {
      NA_character_
    }
    f
  })
}


# =============================================================================
# Internal: bin a MAF value into "rare" / "low" / "common"
# =============================================================================

# Bins follow the standard human-genetics partition. `NA_real_` input
# (e.g. no causal variants in a region) maps to NA_character_ so the
# fit is excluded from MAF stratification rather than counted in any bin.
.maf_bin <- function(maf) {
  if (is.null(maf) || length(maf) == 0L || is.na(maf)) return(NA_character_)
  if (maf <= 0.01)  return("rare")
  if (maf <= 0.05)  return("low")
  "common"
}


# =============================================================================
# Internal: stratify fits by a scenario-level variable and compute metrics
# =============================================================================

.stratify_metrics <- function(fits, by_var, pip_thresholds, n_pip_cal_bins) {
  values <- sort(unique(sapply(fits, function(f) {
    v <- f[[by_var]]
    if (is.null(v)) NA_real_ else v
  })))
  values <- values[!is.na(values)]

  if (length(values) == 0) return(NULL)

  result <- lapply(values, function(v) {
    subset <- Filter(function(f) {
      fv <- f[[by_var]]
      !is.null(fv) && !is.na(fv) && fv == v
    }, fits)
    .metrics_with_se(subset, pip_thresholds, n_pip_cal_bins)
  })
  names(result) <- as.character(values)
  result
}


# =============================================================================
# Internal: stratify fits by causal-variant MAF bin
#
# Groups fits by `causal_maf_bin` (attached in .annotate_fits_with_truth) and
# computes metrics per bin. Bins are kept in the canonical order
# "rare" -> "low" -> "common" rather than alphabetical, so downstream plots
# read left-to-right by increasing MAF.
#
# Returns NULL when no fit has a valid bin assignment (e.g. simulation with
# no causal variants assigned to any region, or genotypes without an attached
# `maf` field).
# =============================================================================

.stratify_metrics_by_maf <- function(fits, pip_thresholds, n_pip_cal_bins) {
  bin_order <- c("rare", "low", "common")

  present <- unique(vapply(fits, function(f) {
    b <- f$causal_maf_bin
    if (is.null(b) || length(b) == 0L) NA_character_ else as.character(b[[1L]])
  }, character(1L)))
  present <- present[!is.na(present)]
  present <- bin_order[bin_order %in% present]

  if (length(present) == 0L) return(NULL)

  result <- lapply(present, function(b) {
    subset <- Filter(function(f) {
      bv <- f$causal_maf_bin
      !is.null(bv) && !is.na(bv) && bv == b
    }, fits)
    .metrics_with_se(subset, pip_thresholds, n_pip_cal_bins)
  })
  names(result) <- present
  result
}


# =============================================================================
# Internal: stratify fits by the true annotation type used in simulation
#
# Groups fits by `true_annotation_type` (attached in
# .annotate_fits_with_truth from simulation$params$annotation_type) and
# computes metrics per type. Kept in canonical order
# "none" -> "binary" -> "continuous" -> "user_supplied" so plots read from
# the misspecified annotation-null case through to richer annotation
# regimes.
#
# Within a single evaluate_methods() call this typically produces a
# single bin, because annotation_type is scenario-level. The axis
# becomes meaningful when evaluation objects from differently-annotated
# simulations are merged.
# =============================================================================

.stratify_metrics_by_annotation_type <- function(fits, pip_thresholds,
                                                  n_pip_cal_bins) {
  type_order <- c("none", "binary", "continuous", "user_supplied")

  present <- unique(vapply(fits, function(f) {
    v <- f$true_annotation_type
    if (is.null(v) || length(v) == 0L) NA_character_ else as.character(v[[1L]])
  }, character(1L)))
  present <- present[!is.na(present)]
  present <- c(type_order[type_order %in% present],
               setdiff(present, type_order))  # tolerate unknown labels

  if (length(present) == 0L) return(NULL)

  result <- lapply(present, function(t) {
    subset <- Filter(function(f) {
      tv <- f$true_annotation_type
      !is.null(tv) && !is.na(tv) && tv == t
    }, fits)
    .metrics_with_se(subset, pip_thresholds, n_pip_cal_bins)
  })
  names(result) <- present
  result
}


# =============================================================================
# Internal: compute all metrics for a collection of annotated fits
# =============================================================================

.metrics_for_fits <- function(fits, pip_thresholds, n_pip_cal_bins) {

  n_fits   <- length(fits)
  valid    <- Filter(function(f) is.null(f$error), fits)
  n_failed <- n_fits - length(valid)

  if (length(valid) == 0) {
    return(list(
      fdr_power_curve = NULL,
      ap              = NA_real_,
      ap_n_fits       = 0L,
      pip_calibration = NULL,
      cs_coverage     = NA_real_,
      cs_power        = NA_real_,
      cs_size_median  = NA_real_,
      cs_size_mean    = NA_real_,
      cs_coverage_n   = 0L,
      cs_power_n      = 0L,
      cs_coverage_pooled    = NA_real_,
      cs_power_pooled       = NA_real_,
      cs_size_median_pooled = NA_real_,
      cs_size_mean_pooled   = NA_real_,
      n_cs_reported   = NA_integer_,
      runtime_mean    = NA_real_,
      runtime_sd      = NA_real_,
      n_fits          = n_fits,
      n_failed        = n_failed
    ))
  }

  # Pool (pip, is_causal) pairs across all valid fits
  pip_vals  <- unlist(lapply(valid, function(f) f$pip))
  is_causal <- unlist(lapply(valid, function(f) {
    v <- logical(f$n_variants)
    v[f$causal_indices] <- TRUE
    v
  }))

  # Power / FDR curve, REGION-WEIGHTED: one curve per region, averaged
  # pointwise. The practitioner's question is about the region in front of
  # them, so the quantity to report is the FDR of a typical region, not the
  # fraction of all calls made anywhere that were wrong. The two differ
  # whenever regions make unequal numbers of calls, which they always do.
  #
  # `n_regions` travels with it: regions selecting nothing at a threshold are
  # excluded from the mean rather than scored FDR 0, because averaging that
  # convention would reward abstention.
  fdr_power <- .fdr_power_curve_per_region(valid, pip_thresholds)

  # The call-weighted curve, kept because it is what most of the literature
  # reports, and because its tp/fp/fn totals are counts a reader can
  # re-aggregate.
  fdr_power_pooled <- .fdr_power_curve(pip_vals, is_causal, pip_thresholds)
  if (!is.null(fdr_power)) {
    # tp/fp/fn are TOTALS over the stratum's regions, so they are counts and
    # belong with the pooled columns - they are what the pooled rates are
    # formed from, and are kept because a count is the one thing a reader can
    # re-aggregate any way they like.
    fdr_power$tp <- fdr_power_pooled$tp
    fdr_power$fp <- fdr_power_pooled$fp
    fdr_power$fn <- fdr_power_pooled$fn
    fdr_power$fdr_pooled       <- fdr_power_pooled$fdr
    fdr_power$precision_pooled <- fdr_power_pooled$precision
    fdr_power$power_pooled     <- fdr_power_pooled$power
    fdr_power$recall_pooled    <- fdr_power_pooled$recall
  } else {
    fdr_power <- fdr_power_pooled
  }

  # AP is a PER-REGION quantity and is averaged over regions, never read off
  # the pooled curve. The AP of several regions' merged ranked list is not the
  # mean of their APs: merging lets one region's false positives outrank
  # another region's true positive, which describes an analysis nobody
  # performs. Computing it on the pooled curve also reads it off the threshold
  # grid, which biases AP downward by more at low heritability than at high -
  # a bias in exactly the dependence one usually wants to measure.
  ap_per_fit <- vapply(valid, function(f) {
    y <- logical(f$n_variants)
    y[f$causal_indices] <- TRUE
    compute_ap_exact(f$pip, y)
  }, numeric(1))
  ap        <- if (any(is.finite(ap_per_fit)))
                 mean(ap_per_fit[is.finite(ap_per_fit)]) else NA_real_
  ap_n_fits <- sum(is.finite(ap_per_fit))


  # PIP calibration, REGION-WEIGHTED for the same reason: "are this region's
  # PIPs trustworthy" is the question being asked. The pooled columns are kept
  # beside it.
  pip_cal        <- .pip_calibration_per_region(valid, n_pip_cal_bins)
  pip_cal_pooled <- .pip_calibration(pip_vals, is_causal, n_pip_cal_bins)
  if (!is.null(pip_cal) && !is.null(pip_cal_pooled) &&
      nrow(pip_cal) == nrow(pip_cal_pooled)) {
    # n / n_causal are totals over the stratum's regions, hence counts; they
    # stay so the binning can be re-aggregated.
    pip_cal$bin_mid            <- pip_cal_pooled$bin_mid
    pip_cal$n                  <- pip_cal_pooled$n
    pip_cal$n_causal           <- pip_cal_pooled$n_causal
    pip_cal$mean_pip_pooled    <- pip_cal_pooled$mean_pip
    pip_cal$frac_causal_pooled <- pip_cal_pooled$frac_causal
  } else if (is.null(pip_cal)) {
    pip_cal <- pip_cal_pooled
  }

  # Credible set metrics
  cs <- .cs_metrics(valid)

  # Runtime
  rts <- Filter(function(x) !is.na(x),
                sapply(valid, function(f) f$runtime_seconds))

  list(
    fdr_power_curve = fdr_power,
    ap              = ap,
    ap_n_fits       = ap_n_fits,
    pip_calibration = pip_cal,
    cs_coverage     = cs$coverage,
    cs_power        = cs$power,
    cs_size_median  = cs$size_median,
    cs_size_mean    = cs$size_mean,
    cs_coverage_n   = cs$coverage_n,
    cs_power_n      = cs$power_n,
    cs_coverage_pooled    = cs$coverage_pooled,
    cs_power_pooled       = cs$power_pooled,
    cs_size_median_pooled = cs$size_median_pooled,
    cs_size_mean_pooled   = cs$size_mean_pooled,
    n_cs_reported   = cs$n_cs,
    runtime_mean    = if (length(rts) > 0) mean(rts) else NA_real_,
    runtime_sd      = if (length(rts) > 1) sd(rts)   else NA_real_,
    n_fits          = n_fits,
    n_failed        = n_failed
  )
}


# =============================================================================
# Internal: per-replicate SE wrapper around .metrics_for_fits
#
# Groups fits by iter value, computes .metrics_for_fits for each iter, then
# adds SE fields:
#   <metric>_se            — SE of each scalar metric across iters
#   fdr_power_curve$power_se, $precision_se — pointwise SE across iters
#   pip_calibration$frac_causal_se          — per-bin SE across iters
#
# When n_iter < 2 all SE fields are NA.
# =============================================================================

.metrics_with_se <- function(fits, pip_thresholds, n_pip_cal_bins) {

  main <- .metrics_for_fits(fits, pip_thresholds, n_pip_cal_bins)

  .se <- function(x) {
    n <- sum(!is.na(x))
    if (n > 1L) stats::sd(x, na.rm = TRUE) / sqrt(n) else NA_real_
  }

  scalar_fields <- c("ap", "cs_coverage", "cs_power",
                     "cs_size_median", "cs_size_mean", "runtime_mean")

  # Identify unique iter values present in these fits
  iter_vals <- sort(unique(vapply(fits, function(f) {
    v <- f$iter
    if (is.null(v) || length(v) == 0L) NA_integer_ else as.integer(v[[1L]])
  }, integer(1L))))
  iter_vals <- iter_vals[!is.na(iter_vals)]

  if (length(iter_vals) < 2L) {
    # Cannot compute SE — attach NA placeholders
    for (fld in scalar_fields)
      main[[paste0(fld, "_se")]] <- NA_real_
    if (!is.null(main$fdr_power_curve)) {
      main$fdr_power_curve$power_se     <- NA_real_
      main$fdr_power_curve$precision_se <- NA_real_
    }
    if (!is.null(main$pip_calibration))
      main$pip_calibration$frac_causal_se <- NA_real_
    return(main)
  }

  # Per-iter metrics (pooling across regions / S / phi within each iter)
  per_iter <- lapply(iter_vals, function(it) {
    sub <- Filter(function(f) {
      !is.null(f$iter) && !is.na(f$iter) && f$iter == it
    }, fits)
    .metrics_for_fits(sub, pip_thresholds, n_pip_cal_bins)
  })

  # Scalar SEs
  for (fld in scalar_fields) {
    vals <- vapply(per_iter, function(m) {
      v <- m[[fld]]; if (is.null(v) || length(v) == 0L) NA_real_ else as.numeric(v[[1L]])
    }, numeric(1L))
    main[[paste0(fld, "_se")]] <- .se(vals)
  }

  # Pointwise SE on fdr_power_curve
  if (!is.null(main$fdr_power_curve)) {
    n_t <- nrow(main$fdr_power_curve)
    pow_mat <- vapply(per_iter, function(m) {
      if (is.null(m$fdr_power_curve)) rep(NA_real_, n_t) else m$fdr_power_curve$power
    }, numeric(n_t))
    prec_mat <- vapply(per_iter, function(m) {
      if (is.null(m$fdr_power_curve)) rep(NA_real_, n_t) else m$fdr_power_curve$precision
    }, numeric(n_t))
    main$fdr_power_curve$power_se     <- apply(pow_mat,  1L, .se)
    main$fdr_power_curve$precision_se <- apply(prec_mat, 1L, .se)
  }

  # Per-bin SE on pip_calibration
  if (!is.null(main$pip_calibration)) {
    n_b <- nrow(main$pip_calibration)
    frac_mat <- vapply(per_iter, function(m) {
      if (is.null(m$pip_calibration)) rep(NA_real_, n_b) else m$pip_calibration$frac_causal
    }, numeric(n_b))
    main$pip_calibration$frac_causal_se <- apply(frac_mat, 1L, .se)
  }

  main
}


# =============================================================================
# Internal: power / FDR curve
#
# Uses a sort + cumsum approach: O(N log N + T) rather than O(N * T).
#
# This computes the CALL-WEIGHTED curve over whatever (pip, is_causal) pairs it
# is given. .fdr_power_curve_per_region() below calls it once per region and
# averages, which is the REGION-WEIGHTED curve - the one reported, because the
# practitioner's question is "I fine-mapped my region; how trustworthy are
# these calls?", not "across every call made anywhere in the benchmark, what
# fraction were wrong". A region that makes many calls dominates the
# call-weighted rate; under the region-weighted one every region counts once.
# =============================================================================

.fdr_power_curve <- function(pip_vals, is_causal, thresholds) {

  total_causal <- sum(is_causal)
  N            <- length(pip_vals)

  if (N == 0 || total_causal == 0) {
    return(data.frame(
      threshold = thresholds,
      tp = 0L, fp = 0L, fn = 0L,
      fdr = 0, power = 0, precision = 1, recall = 0
    ))
  }

  # Sort by decreasing PIP; build cumulative TP and selection counts
  ord          <- order(pip_vals, decreasing = TRUE)
  sorted_pips  <- pip_vals[ord]
  cum_tp       <- c(0L, cumsum(as.integer(is_causal[ord])))  # length N+1

  # For threshold t, n_selected = #{pip >= t}
  # sorted_pips is non-increasing => -sorted_pips is non-decreasing
  # findInterval(-t, -sorted_pips) = largest k such that -sorted_pips[k] <= -t
  #                                = largest k such that sorted_pips[k] >= t
  #                                = n_selected(t)
  neg_sorted <- -sorted_pips
  n_sel      <- findInterval(-thresholds, neg_sorted)
  n_sel      <- pmin(pmax(n_sel, 0L), N)

  tp <- cum_tp[n_sel + 1L]
  fp <- n_sel - tp
  fn <- total_causal - tp

  data.frame(
    threshold = thresholds,
    tp        = tp,
    fp        = fp,
    fn        = fn,
    fdr       = ifelse(n_sel > 0L, fp / n_sel,          0),
    power     = tp / total_causal,
    precision = ifelse(n_sel > 0L, tp / n_sel,          1),
    recall    = tp / total_causal
  )
}


# =============================================================================
# Internal: REGION-WEIGHTED power / FDR curve
#
# One curve per region, then averaged pointwise. Two conventions matter:
#
#   - FDR and precision are NA for a region that selected nothing at that
#     threshold, not 0 and 1. The "FDR is 0 when nothing is selected"
#     convention is harmless when pooling (it contributes no calls), but under
#     averaging it REWARDS ABSTENTION: a method selecting nothing anywhere
#     would average to a perfect FDR of 0. So such regions are excluded from
#     the mean and counted in `n_regions`, which must be read alongside it.
#
#   - power and recall are defined for any region with a causal variant, so
#     they average over all such regions.
# =============================================================================

.fdr_power_curve_per_region <- function(fits, thresholds) {
  n_t <- length(thresholds)
  if (!length(fits)) return(NULL)

  acc <- lapply(fits, function(f) {
    y <- logical(f$n_variants)
    y[f$causal_indices] <- TRUE
    if (sum(y) == 0L || !all(is.finite(f$pip))) return(NULL)
    cur <- .fdr_power_curve(f$pip, y, thresholds)
    n_sel <- cur$tp + cur$fp
    cur$fdr[n_sel == 0L]       <- NA_real_
    cur$precision[n_sel == 0L] <- NA_real_
    cur
  })
  acc <- Filter(Negate(is.null), acc)
  if (!length(acc)) return(NULL)

  grab <- function(col) vapply(acc, function(d) d[[col]], numeric(n_t))
  mrow <- function(m) rowMeans(m, na.rm = TRUE)
  nrow_def <- function(m) rowSums(!is.na(m))

  fdr_m  <- grab("fdr")
  prec_m <- grab("precision")
  out <- data.frame(
    threshold  = thresholds,
    # Region-weighted: the mean over regions, which is what gets reported.
    fdr        = mrow(fdr_m),
    precision  = mrow(prec_m),
    power      = mrow(grab("power")),
    recall     = mrow(grab("recall")),
    # How many regions the fdr/precision means rest on. Never report the rate
    # without this: at a high threshold most regions select nothing.
    n_regions          = nrow_def(fdr_m),
    n_regions_with_causal = length(acc)
  )
  out$fdr[out$n_regions == 0L]       <- NA_real_
  out$precision[out$n_regions == 0L] <- NA_real_
  out
}


# =============================================================================
# Internal: REGION-WEIGHTED PIP calibration
#
# Per-bin observed frequency computed within each region, then averaged over
# the regions that have any variant in that bin. `n_regions` per bin is part of
# the answer: the top bin is usually populated in only a minority of regions.
# =============================================================================

.pip_calibration_per_region <- function(fits, n_bins) {
  if (!length(fits)) return(NULL)
  acc <- lapply(fits, function(f) {
    y <- logical(f$n_variants)
    y[f$causal_indices] <- TRUE
    if (!all(is.finite(f$pip)) || !length(f$pip)) return(NULL)
    .pip_calibration(f$pip, y, n_bins)
  })
  acc <- Filter(function(d) !is.null(d) && nrow(d) > 0L, acc)
  if (!length(acc)) return(NULL)

  base <- acc[[1L]][, c("bin", "bin_lower", "bin_upper")]
  n_b  <- nrow(base)
  grab <- function(col) vapply(acc, function(d) d[[col]], numeric(n_b))

  frac <- grab("frac_causal")
  base$mean_pip    <- rowMeans(grab("mean_pip"),  na.rm = TRUE)
  base$frac_causal <- rowMeans(frac,              na.rm = TRUE)
  base$n_regions   <- rowSums(!is.na(frac))
  base$mean_pip[base$n_regions == 0L]    <- NA_real_
  base$frac_causal[base$n_regions == 0L] <- NA_real_
  base
}


# =============================================================================
# Internal: PIP calibration
#
# Bins all (pip, is_causal) pairs into n_bins equal-width intervals on [0,1].
# Returns mean_pip (expected) vs frac_causal (observed) per bin.
# =============================================================================

.pip_calibration <- function(pip_vals, is_causal, n_bins = 10L) {

  breaks  <- seq(0, 1, length.out = n_bins + 1L)
  bin_ids <- findInterval(pip_vals, breaks, rightmost.closed = TRUE)
  bin_ids <- pmin(pmax(bin_ids, 1L), n_bins)

  ns          <- tabulate(bin_ids, nbins = n_bins)
  n_causal    <- tapply(as.integer(is_causal), bin_ids, sum)
  n_causal    <- as.integer(n_causal[as.character(seq_len(n_bins))])
  n_causal[is.na(n_causal)] <- 0L

  # Mean PIP per bin (only over non-empty bins)
  mean_pip <- vapply(seq_len(n_bins), function(b) {
    idx <- bin_ids == b
    if (!any(idx)) NA_real_ else mean(pip_vals[idx])
  }, numeric(1))

  frac_causal <- ifelse(ns > 0L, n_causal / ns, NA_real_)

  data.frame(
    bin         = seq_len(n_bins),
    bin_lower   = breaks[-length(breaks)],
    bin_upper   = breaks[-1],
    bin_mid     = (breaks[-length(breaks)] + breaks[-1]) / 2,
    n           = ns,
    n_causal    = n_causal,
    mean_pip    = mean_pip,
    frac_causal = frac_causal
  )
}


# =============================================================================
# Internal: credible set metrics
#
# Coverage : proportion of reported CSs containing ≥1 true causal variant
# Power    : proportion of true causal variants captured by any CS
# Size     : median and mean number of variants per CS
# =============================================================================

.cs_metrics <- function(valid_fits) {

  # Pooled accumulators: one entry per credible SET, one per CAUSAL VARIANT,
  # across every region.
  cs_hits  <- logical(0)
  cs_sizes <- integer(0)
  causal_captured <- logical(0)

  # Region-weighted accumulators: one entry per REGION. These are the reported
  # values, because fine-mapping is applied to a region and the question is "of
  # the sets I got for my locus, how many hold a causal variant". Pooling over
  # sets weights a region reporting five sets five times as heavily as one
  # reporting a single set.
  r_cov <- numeric(0); r_pow <- numeric(0)
  r_sz_mean <- numeric(0); r_sz_med <- numeric(0)

  for (f in valid_fits) {
    n_causal <- length(f$causal_indices)

    if (length(f$credible_sets) == 0) {
      # No sets reported: coverage and set size are UNDEFINED for this region
      # rather than 0 - there is no set to be right or wrong about, and
      # scoring it 0 would punish abstention as hard as a wrong call. Power IS
      # defined: none of the causal variants were captured.
      causal_captured <- c(causal_captured, rep(FALSE, n_causal))
      if (n_causal > 0L) r_pow <- c(r_pow, 0)
      next
    }

    all_cs_variants <- unique(unlist(f$credible_sets))
    sizes <- vapply(f$credible_sets, length, integer(1))
    hits  <- vapply(f$credible_sets,
                    function(cs) any(cs %in% f$causal_indices), logical(1))

    cs_hits         <- c(cs_hits,  hits)
    cs_sizes        <- c(cs_sizes, sizes)
    causal_captured <- c(causal_captured,
                         f$causal_indices %in% all_cs_variants)

    r_cov     <- c(r_cov,     mean(hits))
    r_sz_mean <- c(r_sz_mean, mean(sizes))
    r_sz_med  <- c(r_sz_med,  stats::median(sizes))
    if (n_causal > 0L)
      r_pow <- c(r_pow, sum(f$causal_indices %in% all_cs_variants) / n_causal)
  }

  n_cs <- length(cs_hits)
  mn   <- function(v) if (length(v)) mean(v) else NA_real_

  list(
    # Region-weighted - reported.
    coverage    = mn(r_cov),
    power       = mn(r_pow),
    size_median = mn(r_sz_med),
    size_mean   = mn(r_sz_mean),
    # Regions each mean rests on. coverage and the sizes count only regions
    # that reported at least one set.
    coverage_n  = length(r_cov),
    power_n     = length(r_pow),
    # Call-weighted, over sets and over causal variants respectively.
    coverage_pooled    = if (n_cs > 0L) mean(cs_hits) else NA_real_,
    power_pooled       = if (length(causal_captured) > 0L)
                           mean(causal_captured) else NA_real_,
    size_median_pooled = if (n_cs > 0L) stats::median(cs_sizes) else NA_real_,
    size_mean_pooled   = if (n_cs > 0L) mean(cs_sizes) else NA_real_,
    n_cs        = n_cs
  )
}
