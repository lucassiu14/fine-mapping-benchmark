# =============================================================================
# R/evaluate_extras.R
#
# The per-fit metric functions specified in variable-importance-analysis.md §5.1,
# §5.2b and §4.4. Kept in a SEPARATE FILE so they can be reviewed and unit-tested
# before being wired into R/evaluate.R, which is load-bearing for every existing
# result.
#
# WIRING (see the spec's §5.1b). In `.metrics_for_fits`, AP must be computed PER
# FIT and averaged, not on the concatenated PIP vector of the two same-size
# regions:
#
#   ap_per_fit <- vapply(valid, function(f) {
#     y <- logical(f$n_variants); y[f$causal_indices] <- TRUE
#     .compute_ap_exact(f$pip, y)
#   }, numeric(1))
#   ap        <- mean(ap_per_fit, na.rm = TRUE)
#   ap_n_fits <- sum(!is.na(ap_per_fit))
#
# Leave the FDR and calibration COUNTS pooled across the two fits - that is
# correct (spec §3.2). Only the ranking metrics change.
#
# STATUS: written, never executed. Nothing here has been run against real data.
# Add the §5.8 tests and run them before trusting any output.
# =============================================================================


#' Exact rank-based average precision
#'
#' AP is a function of the ranking alone, so it must be evaluated at the DISTINCT
#' OBSERVED PIP values rather than on the fixed 0.005 grid the FDR curve needs.
#' The grid version treats every variant inside a 0.005-wide band as entering the
#' selection simultaneously, so a causal is credited with the precision computed
#' after every non-causal in its band has also been admitted. Worked example from
#' the spec: a causal at PIP 0.0032 with true rank 61 of 1000 has exact AP 1/61 =
#' 0.0164 but grid AP 1/1000 = 0.001.
#'
#' That error is zero at strong signal and grows as PIPs compress toward zero, so
#' it is a PHI-DEPENDENT BIAS in the estimator of a metric whose phi dependence
#' the analysis is trying to measure. A constant bias would be harmless.
#'
#' Tie handling is by construction: all members of a tied block share one
#' operating point. Do NOT implement this as order() + precision-at-causal-
#' positions, which would let order()'s index tie-break - genomic position -
#' decide AP for methods emitting many exactly-tied PIPs.
#'
#' Matches sklearn.metrics.average_precision_score.
#'
#' @param pip_vals Numeric vector of PIPs for ONE fit.
#' @param is_causal Logical vector, same length.
#' @return Scalar AP, or NA if there are no causal variants.
.compute_ap_exact <- function(pip_vals, is_causal) {
  S <- sum(is_causal)
  if (S == 0L || length(pip_vals) == 0L) return(NA_real_)
  ord  <- order(pip_vals, decreasing = TRUE)
  p    <- pip_vals[ord]
  y    <- as.integer(is_causal[ord])
  last <- c(p[-1L] != p[-length(p)], TRUE)   # TRUE at the end of each tied block
  k    <- which(last)                        # n_selected at each distinct value
  tp   <- cumsum(y)[last]
  prec <- tp / k
  rec  <- tp / S
  sum(diff(c(0, rec)) * prec)                # step integral, NOT trapezoid
}


#' Threshold-free per-fit ranking and set metrics (spec §5.2b: P4, P5, Q5)
#'
#' Per-locus quantities, so these are AVERAGED across fits, never pooled - with
#' one exception: `set_hit_*` is a Bernoulli per fit and is SUMMED, its cell-level
#' rate being the coverage.
#'
#' `prec_at_k` matters because methods' PIP scales differ by orders of magnitude,
#' so any fixed-threshold comparison partly measures calibration rather than
#' ranking. `set_size_*` / `set_hit_*` is the practitioner's actual question -
#' how many variants must I follow up, and will the causal be among them - and is
#' the volume/error pair the spec's rule R1 demands.
#'
#' @param pip Numeric vector of PIPs for one fit.
#' @param is_causal Logical vector, same length.
#' @param alpha Cumulative-mass targets for the uniform credible sets.
#' @return Named list of scalars.
.rank_set_metrics <- function(pip, is_causal, alpha = c(0.5, 0.95),
                              S = sum(is_causal)) {
  stopifnot(length(pip) == length(is_causal))
  ord <- order(pip, decreasing = TRUE)
  p   <- pip[ord]
  y   <- as.integer(is_causal[ord])
  cum <- cumsum(p)
  tot <- sum(p)

  out <- list(
    prec_at_1     = if (length(y)) y[1] else NA_real_,
    prec_at_S     = if (S > 0) mean(y[seq_len(min(S,      length(y)))]) else NA_real_,
    prec_at_2S    = if (S > 0) mean(y[seq_len(min(2L * S, length(y)))]) else NA_real_,
    # Q5 - sharpness, the mandatory companion to any calibration number.
    mass_above_50 = if (tot > 0) sum(p[p >= 0.5]) / tot else NA_real_,
    mean_pip_topS = if (S > 0) mean(p[seq_len(min(S, length(p)))]) else NA_real_
  )

  for (a in alpha) {
    reached <- tot >= a
    k <- if (reached) which(cum >= a)[1L] else length(p)
    tag <- sub("0\\.", "", format(a))
    out[[paste0("set_size_",    tag)]] <- k
    out[[paste0("set_hit_",     tag)]] <- as.integer(any(y[seq_len(k)] == 1L))
    out[[paste0("set_prec_",    tag)]] <- mean(y[seq_len(k)])
    # A method whose total mass never reaches alpha has k = the whole region.
    # That is a RESULT, not a missing value - record it so downstream code does
    # not silently treat a maximally-diffuse posterior as a large credible set.
    out[[paste0("set_reached_", tag)]] <- reached
  }
  out
}


#' Murphy decomposition of the Brier score (spec §4.4)
#'
#' Exists because both headline calibration metrics score a USELESS predictor
#' perfectly. A method emitting p_j = S/p for every variant has mass ratio 1.000,
#' ECE 0.000 and hi_pip_reliab NA. RES is what closes that hole: it is exactly
#' zero for the uniform predictor.
#'
#'   BS = REL - RES + UNC
#'   REL  calibration error   (lower better)
#'   RES  sharpness           (higher better; 0 for the uniform predictor)
#'   UNC  base rate           (a property of the cell, not the method)
#'
#' Raw Brier is NOT comparable across cells because UNC varies over two orders of
#' magnitude (ybar = S/p), hence the skill score BSS.
#'
#' Computable entirely from bin counts already stored - no collection change.
#'
#' @param n_b,c_b,sum_pip_b Per-bin counts: variants, causals, summed PIP.
#' @return Named list: rel, res, unc, bs, bss.
.murphy_decomposition <- function(n_b, c_b, sum_pip_b) {
  keep <- n_b > 0
  if (!any(keep)) {
    return(list(rel = NA_real_, res = NA_real_, unc = NA_real_,
                bs = NA_real_, bss = NA_real_))
  }
  n_b <- n_b[keep]; c_b <- c_b[keep]; sum_pip_b <- sum_pip_b[keep]
  N    <- sum(n_b)
  pbar <- sum_pip_b / n_b          # mean predicted PIP in the bin
  ybar <- c_b / n_b                # observed causal fraction in the bin
  ytot <- sum(c_b) / N             # base rate

  rel <- sum(n_b * (pbar - ybar)^2) / N
  res <- sum(n_b * (ybar - ytot)^2) / N
  unc <- ytot * (1 - ytot)
  bs  <- rel - res + unc
  list(rel = rel, res = res, unc = unc, bs = bs,
       bss = if (unc > 0) 1 - bs / unc else NA_real_)
}


#' Excess PIP mass by band (spec §4.3, Q2)
#'
#' `total_mass_ratio` says a method is over-confident; this says WHERE. Excess
#' spread thinly over nulls is harmless; excess concentrated in high-PIP calls is
#' confidently wrong. No existing scalar distinguishes these, and the distinction
#' decides whether a miscalibration matters.
#'
#' Measured on Iteration 003 (ref500): susie put 56% of its excess below PIP 0.1,
#' fb_xregion put 59% in the 0.5-0.9 band - at near-identical mass ratios.
#'
#' @param bin_lower Numeric vector of each bin's lower edge.
#' @param c_b,sum_pip_b Per-bin causal counts and summed PIP.
#' @return Named numeric vector of excess mass per band.
.excess_mass_by_band <- function(bin_lower, c_b, sum_pip_b) {
  bands <- list(lo = c(0, 0.1), mid = c(0.1, 0.5),
                hi = c(0.5, 0.9), top = c(0.9, 1.0 + 1e-9))
  vapply(bands, function(b) {
    sel <- bin_lower >= b[1] & bin_lower < b[2]
    if (!any(sel)) return(0)
    sum(sum_pip_b[sel]) - sum(c_b[sel])
  }, numeric(1))
}


#' Self-consistency of a method's own error claim (spec §4.2, P3)
#'
#' Requires no external calibration theory: it asks the method to state its own
#' expected false-positive count among its selections, then scores it.
#'
#'   exp_fp(t) = sum over selected of (1 - p_j) = nsel(t) - sum_pip_sel(t)
#'   ExFP(t)   = FP(t) - exp_fp(t)      "3.2 more false calls than it admits to"
#'   rho(t)    = FP(t) / exp_fp(t)      > 1 means more errors than admitted
#'
#' Zero for an abstainer by construction, which is honest - neither rewarded nor
#' punished, unlike raw FDR where selecting nothing wins every comparison.
.honesty <- function(fp, nsel, sum_pip_sel) {
  exp_fp <- nsel - sum_pip_sel
  list(exp_fp = exp_fp,
       ex_fp  = fp - exp_fp,
       rho    = ifelse(exp_fp > 0, fp / exp_fp, NA_real_))
}
