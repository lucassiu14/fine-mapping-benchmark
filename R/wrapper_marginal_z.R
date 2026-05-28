# =============================================================================
# marginal_z.R
#
# Wrapper for the marginal-z baseline.
#
# This is a model-free baseline that ranks variants by the magnitude of their
# marginal z-score and converts the ranking into normalised "PIPs":
#
#   pip_j = |z_j| / sum_k |z_k|
#
# The PIPs sum to 1, consistent with a single-causal-variant assumption. A
# credible set is constructed greedily: variants are sorted by descending PIP
# and accumulated until the cumulative PIP reaches `coverage`. This is the
# same CS construction used by ABF and PAINTOR for fairness across baselines.
#
# Why include it: it is the "are you better than the raw z-scores?" baseline.
# No real fine-mapping method should be appreciably worse than marginal_z on
# any setting; if it is, the method is mis-tuned or actively harmful. The gap
# between marginal_z and ABF is informative — it shows what shrinkage adds
# under a single-causal model.
#
# No external packages or binaries are required.
#
# This file provides:
#   - run_marginal_z()        : runs the baseline on a single region (explicit inputs)
#   - run_marginal_z_region() : adapter called by run_methods()
# =============================================================================


# =============================================================================
# Run marginal_z on a single region
# =============================================================================

#' Run the marginal-z baseline on a single region
#'
#' Computes PIPs as the normalised absolute z-scores
#' (\code{pip_j = |z_j| / sum_k |z_k|}) and derives a single credible set
#' greedily from the ranking until cumulative PIP reaches \code{coverage}.
#'
#' This is a model-free baseline — no LD matrix, sample size, or prior
#' assumption is used. Its purpose is to anchor reported AUPRC numbers: any
#' real fine-mapping method should outperform it on a non-degenerate setting.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param coverage Numeric. Coverage level for the credible set. Default: 0.95.
#'
#' @return A list with the standardised fine-mapping output:
#' \describe{
#'   \item{pip}{Numeric vector (length p). Normalised |z|-based PIPs that sum
#'     to 1.}
#'   \item{credible_sets}{List containing one integer vector: the indices of
#'     variants in the credible set (sorted ascending), constructed greedily
#'     from PIPs until cumulative PIP \eqn{\ge} \code{coverage}.}
#'   \item{method}{Character. Always \code{"marginal_z"}.}
#'   \item{input_type}{Character. Always \code{"summary"}.}
#'   \item{params}{List. Hyperparameters used (\code{coverage}).}
#'   \item{runtime_seconds}{Numeric. Wall-clock time in seconds (effectively
#'     zero — included for interface consistency with the other methods).}
#'   \item{additional}{Empty list (no method-specific outputs).}
#' }
#'
#' @export
run_marginal_z <- function(z, coverage = 0.95) {

  # --- Validate ---------------------------------------------------------------

  p <- length(z)

  stopifnot(
    "z must be a non-empty numeric vector" =
      is.numeric(z) && p >= 1,
    "coverage must be a single number in (0, 1)" =
      is.numeric(coverage) && length(coverage) == 1 &&
      coverage > 0 && coverage < 1
  )

  # --- Compute PIPs -----------------------------------------------------------

  start_time <- proc.time()

  abs_z <- abs(z)

  if (sum(abs_z) == 0) {
    # Degenerate case: all z-scores are zero. Spread PIP uniformly.
    pip <- rep(1 / p, p)
  } else {
    pip <- abs_z / sum(abs_z)
  }

  # --- Derive credible set ----------------------------------------------------
  # Greedy: sort variants by PIP descending, accumulate until >= coverage.

  ord        <- order(pip, decreasing = TRUE)
  cumulative <- cumsum(pip[ord])
  n_in_cs    <- which(cumulative >= coverage)[1]
  if (is.na(n_in_cs)) n_in_cs <- p   # safety: include all if PIPs don't sum to 1
  cs_indices <- sort(ord[seq_len(n_in_cs)])

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  # --- Return -----------------------------------------------------------------

  list(
    pip             = pip,
    credible_sets   = list(cs_indices),
    method          = "marginal_z",
    input_type      = "summary",
    params          = list(coverage = coverage),
    runtime_seconds = elapsed,
    additional      = list()
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run the marginal-z baseline on a single region from simulation data
#'
#' Thin adapter that extracts z-scores from the simulation's
#' \code{region_pheno} object and calls \code{\link{run_marginal_z}}.
#'
#' @param region_geno List. One element of \code{simulation$genotypes}. Not
#'   used by marginal_z but required for a consistent adapter signature.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z}.
#' @param ... Additional arguments passed to \code{\link{run_marginal_z}}
#'   (e.g. \code{coverage}).
#'
#' @return The output of \code{\link{run_marginal_z}}.
#' @export
run_marginal_z_region <- function(region_geno, region_pheno, ...) {
  run_marginal_z(z = region_pheno$z, ...)
}
