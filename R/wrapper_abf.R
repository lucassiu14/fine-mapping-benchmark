# =============================================================================
# abf.R
#
# Wrapper for ABF (Approximate Bayes Factor) fine-mapping.
#
# ABF uses the Wakefield (2009) approximation to compute a Bayes factor for
# each variant being causal, under the assumption that exactly one variant in
# the region is causal. PIPs are obtained by normalising the ABFs. A single
# credible set is derived by greedily accumulating the top variants until the
# cumulative PIP reaches the coverage threshold.
#
# No external packages or binaries are required.
#
# Reference:
#   Wakefield J (2009). Bayes factors for genome-wide association studies:
#   comparison with p-values. Genetic Epidemiology, 33(1), 79-86.
#
# This file provides:
#   - run_abf()         : runs ABF on a single region (explicit inputs)
#   - run_abf_region()  : adapter called by run_methods()
# =============================================================================


# =============================================================================
# Run ABF on a single region
# =============================================================================

#' Run ABF fine-mapping on a single region
#'
#' Computes Wakefield approximate Bayes factors and derives posterior inclusion
#' probabilities (PIPs) under a single-causal-variant model. A single credible
#' set is returned as the minimal set of top-ranked variants whose cumulative
#' PIP reaches \code{coverage}.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param se Numeric vector or NULL. Standard errors (length p). If NULL,
#'   set to 1 for all variants (equivalent to treating z-scores as
#'   standardised effect estimates).
#' @param prior_variance Numeric. Prior variance W on effect sizes under H1
#'   (variant is causal). Default: 0.04 (prior SD = 0.2), following the
#'   Wakefield convention for standardised effects.
#' @param coverage Numeric. Coverage level for the credible set. Default: 0.95.
#'
#' @return A list with the standardised fine-mapping output:
#' \describe{
#'   \item{pip}{Numeric vector (length p). Posterior inclusion probabilities
#'     under the single-causal-variant model.}
#'   \item{credible_sets}{List containing one integer vector: the indices of
#'     variants in the 95\% (or \code{coverage}) credible set, ordered by
#'     decreasing PIP.}
#'   \item{method}{Character. Always \code{"abf"}.}
#'   \item{input_type}{Character. Always \code{"summary"}.}
#'   \item{params}{List. Hyperparameters used.}
#'   \item{runtime_seconds}{Numeric. Wall-clock time in seconds.}
#'   \item{additional}{List:
#'     \describe{
#'       \item{log10_abf}{Numeric vector (length p). Log10 approximate Bayes
#'         factor per variant.}
#'     }
#'   }
#' }
#'
#' @export
run_abf <- function(z,
                    se             = NULL,
                    prior_variance = 0.04,
                    coverage       = 0.95) {

  # --- Validate ---------------------------------------------------------------

  p <- length(z)

  stopifnot(
    "prior_variance must be a single positive number" =
      is.numeric(prior_variance) && length(prior_variance) == 1 &&
      prior_variance > 0,
    "coverage must be a single number in (0, 1)" =
      is.numeric(coverage) && length(coverage) == 1 &&
      coverage > 0 && coverage < 1
  )

  if (is.null(se)) se <- rep(1.0, p)

  stopifnot(
    "se must be a numeric vector of the same length as z" =
      is.numeric(se) && length(se) == p,
    "se values must be positive" = all(se > 0)
  )

  # --- Compute ABFs -----------------------------------------------------------
  # Wakefield (2009) approximation:
  #   V   = se^2  (sampling variance)
  #   W   = prior_variance
  #   ABF = sqrt(V / (V + W)) * exp(z^2 / 2 * W / (V + W))

  start_time <- proc.time()

  V   <- se^2
  W   <- prior_variance
  r   <- W / (V + W)           # shrinkage factor

  log_abf    <- 0.5 * log(1 - r) + 0.5 * r * z^2
  log10_abf  <- log_abf / log(10)

  # Normalise to PIPs (in log space for numerical stability)
  log_abf_max <- max(log_abf)
  pip <- exp(log_abf - log_abf_max)
  pip <- pip / sum(pip)

  # --- Derive credible set ----------------------------------------------------
  # Greedy: sort variants by PIP descending, accumulate until >= coverage.

  ord <- order(pip, decreasing = TRUE)
  cumulative <- cumsum(pip[ord])
  n_in_cs    <- which(cumulative >= coverage)[1]
  cs_indices <- sort(ord[seq_len(n_in_cs)])

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  # --- Return -----------------------------------------------------------------

  list(
    pip             = pip,
    credible_sets   = list(cs_indices),
    method          = "abf",
    input_type      = "summary",
    params          = list(
      prior_variance = prior_variance,
      coverage       = coverage
    ),
    runtime_seconds = elapsed,
    additional      = list(
      log10_abf = log10_abf
    )
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run ABF on a single region from simulation data structures
#'
#' Thin adapter that extracts inputs from the simulation's \code{region_pheno}
#' object and calls \code{\link{run_abf}}.
#'
#' @param region_geno List. One element of \code{simulation$genotypes}.
#'   Not used by ABF but required for a consistent adapter signature.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z} and \code{se}.
#' @param ... Additional arguments passed to \code{\link{run_abf}}
#'   (e.g. \code{prior_variance}, \code{coverage}).
#'
#' @return The output of \code{\link{run_abf}}.
#' @export
run_abf_region <- function(region_geno, region_pheno, ...) {
  run_abf(
    z  = region_pheno$z,
    se = region_pheno$se,
    ...
  )
}
