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
#   fdr_power_curve  — data.frame: threshold, tp, fp, fn, fdr, power, precision, recall
#   auprc            — scalar (area under precision-recall curve)
#   pip_calibration  — data.frame: bin, bin_lower, bin_upper, bin_mid, n, n_causal,
#                                  mean_pip, frac_causal
#   cs_coverage      — proportion of reported CSs containing ≥1 true causal variant
#   cs_power         — proportion of true causal variants captured by any CS
#   cs_size_median   — median variants per CS
#   cs_size_mean     — mean variants per CS
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
#' @param save Logical. If TRUE, save the full evaluation object as
#'   \code{evaluation.rds} inside \code{output_dir}, and write a flat
#'   per-method summary table as \code{evaluation_summary.csv}.
#'   Default: FALSE.
#' @param output_dir Character. Directory in which to save results when
#'   \code{save = TRUE}. Created automatically if it does not exist.
#'   Default: \code{"results"}.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A named list, one element per method, each containing
#'   \code{global}, \code{by_S}, \code{by_phi}, and \code{by_p_causal}
#'   sub-lists.  See file header for the fields in each stratum.
#'
#' @export
evaluate_methods <- function(simulation,
                             results,
                             pip_thresholds = seq(0, 1, by = 0.005),
                             n_pip_cal_bins = 10L,
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
      by_S                    = .stratify_metrics(fits, "S",        pip_thresholds, n_pip_cal_bins),
      by_phi                  = .stratify_metrics(fits, "phi",      pip_thresholds, n_pip_cal_bins),
      by_p_causal             = .stratify_metrics(fits, "p_causal", pip_thresholds, n_pip_cal_bins),
      by_causal_maf           = .stratify_metrics_by_maf(fits, pip_thresholds, n_pip_cal_bins),
      by_true_annotation_type = .stratify_metrics_by_annotation_type(fits, pip_thresholds, n_pip_cal_bins)
    )
  }

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
        auprc               = g$auprc,
        auprc_se            = g$auprc_se,
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
      auprc           = NA_real_,
      ap              = NA_real_,
      pip_calibration = NULL,
      cs_coverage     = NA_real_,
      cs_power        = NA_real_,
      cs_size_median  = NA_real_,
      cs_size_mean    = NA_real_,
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

  # Power / FDR curve
  fdr_power <- .fdr_power_curve(pip_vals, is_causal, pip_thresholds)

  # AUPRC (trapezoid, legacy) + AP (correct estimator)
  auprc <- .compute_auprc(fdr_power$precision, fdr_power$recall)
  ap    <- .compute_ap(fdr_power$precision, fdr_power$recall)

  # PIP calibration
  pip_cal <- .pip_calibration(pip_vals, is_causal, n_pip_cal_bins)

  # Credible set metrics
  cs <- .cs_metrics(valid)

  # Runtime
  rts <- Filter(function(x) !is.na(x),
                sapply(valid, function(f) f$runtime_seconds))

  list(
    fdr_power_curve = fdr_power,
    auprc           = auprc,
    ap              = ap,
    pip_calibration = pip_cal,
    cs_coverage     = cs$coverage,
    cs_power        = cs$power,
    cs_size_median  = cs$size_median,
    cs_size_mean    = cs$size_mean,
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

  scalar_fields <- c("auprc", "cs_coverage", "cs_power",
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
# Pools all variants across all fits in the supplied collection.
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
# Internal: AUPRC via the trapezoidal rule
# =============================================================================

.compute_auprc <- function(precision, recall) {
  if (length(recall) < 2) return(NA_real_)

  # Sort by increasing recall; at recall ties keep the higher precision
  ord   <- order(recall, -precision)
  prec  <- precision[ord]
  rec   <- recall[ord]

  # Trapezoidal rule: sum of trapezoids.
  # NOTE: linear interpolation in PR space is not valid (Davis & Goadrich
  # 2006) - the segment between two PR points is not achievable by any
  # classifier. This estimator is RETAINED for backwards comparability with
  # earlier runs, but .compute_ap() below is the correct one and is what
  # downstream analysis should use. Empirically the trapezoid understates by
  # ~29% in the weak-signal regime and ~0% at strong signal, i.e. the error
  # is signal-dependent, not a constant offset.
  sum(diff(rec) * (head(prec, -1) + tail(prec, -1)) / 2)
}


#' Average precision (the correct PR-curve summary)
#'
#' Step-function integral of precision over recall - no interpolation, so it
#' avoids the invalid linear interpolation the trapezoidal rule performs in
#' PR space. Computed from the same (precision, recall) grid as
#' \code{.compute_auprc}, which recovers the exact ranking-based AP to
#' within ~1-5%.
#'
#' @param precision,recall Numeric vectors from the FDR/power curve.
#' @return Scalar average precision, or NA if the curve is degenerate.
#' @keywords internal
.compute_ap <- function(precision, recall) {
  if (length(recall) < 2) return(NA_real_)
  ord  <- order(recall, -precision)
  prec <- precision[ord]
  rec  <- recall[ord]
  # sum over steps: (R_i - R_{i-1}) * P_i, with R_0 = 0
  sum(diff(c(0, rec)) * prec)
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

  cs_hits  <- logical(0)   # per-CS: did it contain a causal?
  cs_sizes <- integer(0)   # per-CS: how many variants?
  causal_captured <- logical(0)  # per-causal-variant: was it in any CS?

  for (f in valid_fits) {
    n_causal <- length(f$causal_indices)

    if (length(f$credible_sets) == 0) {
      # No CSs: no causal variants captured
      causal_captured <- c(causal_captured, rep(FALSE, n_causal))
      next
    }

    # Variants covered by any CS in this fit
    all_cs_variants <- unique(unlist(f$credible_sets))

    for (cs in f$credible_sets) {
      cs_hits  <- c(cs_hits,  any(cs %in% f$causal_indices))
      cs_sizes <- c(cs_sizes, length(cs))
    }

    causal_captured <- c(
      causal_captured,
      f$causal_indices %in% all_cs_variants
    )
  }

  n_cs <- length(cs_hits)

  list(
    coverage    = if (n_cs > 0L)                   mean(cs_hits)       else NA_real_,
    power       = if (length(causal_captured) > 0L) mean(causal_captured) else NA_real_,
    size_median = if (n_cs > 0L)                   median(cs_sizes)    else NA_real_,
    size_mean   = if (n_cs > 0L)                   mean(cs_sizes)      else NA_real_,
    n_cs        = n_cs
  )
}
