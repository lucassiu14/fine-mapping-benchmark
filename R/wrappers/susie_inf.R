# =============================================================================
# susie_inf.R
#
# Wrapper for SuSiE-inf (SuSiE with an infinitesimal polygenic component).
#
# SuSiE-inf extends SuSiE by adding a genome-wide infinitesimal term to the
# model, making it more robust when there is polygenic background signal in
# the region. It is implemented in susieR via the unmappable_effects = "inf"
# argument, which requires susieR >= 0.15.0.
#
# This file provides:
#   - setup_susie_inf()        : checks susieR is available and new enough
#   - run_susie_inf()          : runs SuSiE-inf on a single region
#   - run_susie_inf_region()   : adapter called by run_methods()
#
# Standard output format:
#   pip              Numeric vector (length p). Posterior inclusion probabilities.
#   credible_sets    List of integer vectors (variant indices), or empty list.
#   method           Character. "susie_inf".
#   input_type       Character. "individual" or "summary".
#   params           List. Hyperparameters used.
#   runtime_seconds  Numeric. Wall-clock time.
#   additional       List. SuSiE-inf-specific outputs (see run_susie_inf() docs).
# =============================================================================


# =============================================================================
# Setup
# =============================================================================

#' Check susieR is available and supports SuSiE-inf
#'
#' Verifies that susieR is installed and is a version that supports the
#' \code{unmappable_effects} argument (>= 0.15.0). If susieR is present but
#' too old, installation instructions are printed.
#'
#' @return Invisible TRUE if the requirement is met.
#' @export
setup_susie_inf <- function() {

  if (!requireNamespace("susieR", quietly = TRUE)) {
    stop(
      "susieR is not installed. Install it with:\n",
      "  install.packages('susieR')\n",
      "or from GitHub:\n",
      "  remotes::install_github('stephenslab/susieR')",
      call. = FALSE
    )
  }

  ver <- utils::packageVersion("susieR")
  if (ver < "0.15.0") {
    stop(
      "SuSiE-inf requires susieR >= 0.15.0 (installed: ", ver, ").\n\n",
      "The unmappable_effects argument was added after v0.14.x.\n",
      "Update with:\n",
      "  renv::update('susieR')   # if using renv\n",
      "  install.packages('susieR')  # otherwise\n\n",
      "After updating, run renv::snapshot() to record the new version.",
      call. = FALSE
    )
  }

  # Confirm the argument actually exists (guard against unexpected API changes)
  rss_args <- names(formals(susieR::susie_rss))
  if (!"unmappable_effects" %in% rss_args) {
    stop(
      "susieR ", ver, " is installed but susie_rss() does not have an\n",
      "'unmappable_effects' argument. Please update to the latest susieR:\n",
      "  remotes::install_github('stephenslab/susieR')",
      call. = FALSE
    )
  }

  message("susieR ", ver, " is available and supports SuSiE-inf.")
  invisible(TRUE)
}


# =============================================================================
# Run SuSiE-inf on a single region
# =============================================================================

#' Run SuSiE-inf fine-mapping on a single region
#'
#' Calls \code{susieR::susie_rss} (or \code{susieR::susie} for individual-level
#' data) with \code{unmappable_effects = "inf"}, which adds a polygenic
#' infinitesimal component to the standard SuSiE model.
#'
#' Requires susieR >= 0.15.0. Run \code{setup_susie_inf()} to check.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param n Integer. Sample size.
#' @param X Optional matrix. Individual-level genotype data (n x p).
#'   If provided with \code{y}, individual-level mode is used.
#' @param y Optional numeric vector. Phenotype vector (length n).
#' @param L Integer. Maximum number of sparse single-effect components.
#'   Default: 10.
#' @param estimate_residual_variance Logical. Default: TRUE.
#' @param estimate_prior_variance Logical. Default: TRUE.
#' @param prior_variance Numeric or NULL. Prior effect variance. If NULL,
#'   set to \code{0.1 * var(y)} (individual-level) or 0.1 (summary stats).
#' @param coverage Numeric. Coverage level for credible sets. Default: 0.95.
#' @param min_abs_corr Numeric. Purity threshold for filtering credible sets.
#'   Default: 0.5.
#' @param max_iter Integer. Maximum IBSS iterations. Default: 100.
#' @param verbose_susie Logical. Print susieR's internal progress. Default: FALSE.
#'
#' @return A list with the standardised fine-mapping output:
#' \describe{
#'   \item{pip}{Numeric vector (length p). Posterior inclusion probabilities.}
#'   \item{credible_sets}{List of integer vectors, or empty list.}
#'   \item{method}{Character. Always \code{"susie_inf"}.}
#'   \item{input_type}{Character. \code{"individual"} or \code{"summary"}.}
#'   \item{params}{List. Hyperparameters used.}
#'   \item{runtime_seconds}{Numeric.}
#'   \item{additional}{List of SuSiE-inf-specific outputs:
#'     \describe{
#'       \item{alpha}{Matrix (L x p). Per-component posterior assignment
#'         probabilities.}
#'       \item{posterior_mean}{Numeric vector (length p). Posterior mean
#'         sparse effect size summed across components.}
#'       \item{lbf}{Numeric vector (length L). Log Bayes factor per component.}
#'       \item{cs_purity}{Data frame or list. Purity statistics per CS.}
#'       \item{converged}{Logical.}
#'       \item{elbo}{Numeric. Final ELBO.}
#'       \item{n_iter_run}{Integer. Number of iterations run.}
#'       \item{tau2}{Numeric. Estimated infinitesimal (polygenic) variance
#'         component.}
#'       \item{theta}{Numeric vector (length p) or NULL. Posterior means of
#'         the infinitesimal effects (BLUP estimates), if returned by susieR.}
#'     }
#'   }
#'   \item{error}{Character or NULL. Error message if the fit failed.}
#' }
#'
#' @export
run_susie_inf <- function(z = NULL,
                          LD = NULL,
                          n = NULL,
                          X = NULL,
                          y = NULL,
                          L = 10,
                          estimate_residual_variance = TRUE,
                          estimate_prior_variance = TRUE,
                          prior_variance = NULL,
                          coverage = 0.95,
                          min_abs_corr = 0.5,
                          max_iter = 100,
                          verbose_susie = FALSE) {

  # --- Validate inputs --------------------------------------------------------

  if (!requireNamespace("susieR", quietly = TRUE)) {
    stop("susieR is not installed. Run setup_susie_inf() first.", call. = FALSE)
  }

  ver <- utils::packageVersion("susieR")
  if (ver < "0.15.0") {
    stop(
      "SuSiE-inf requires susieR >= 0.15.0 (installed: ", ver, "). ",
      "Run setup_susie_inf() for upgrade instructions.",
      call. = FALSE
    )
  }

  use_individual <- !is.null(X) && !is.null(y)
  use_summary    <- !is.null(z) && !is.null(LD) && !is.null(n)

  if (!use_individual && !use_summary) {
    stop(
      "Provide either (X, y) for individual-level data, or (z, LD, n) for ",
      "summary statistics.",
      call. = FALSE
    )
  }

  if (use_individual && use_summary) {
    message("Both individual-level and summary data provided. ",
            "Using individual-level data.")
    use_summary <- FALSE
  }

  if (is.null(prior_variance)) {
    prior_variance <- if (use_individual) 0.1 * var(y) else 0.1
  }

  params <- list(
    L = L,
    estimate_residual_variance = estimate_residual_variance,
    estimate_prior_variance = estimate_prior_variance,
    prior_variance = prior_variance,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    max_iter = max_iter
  )

  p_snps     <- if (use_individual) ncol(X) else length(z)
  input_type <- if (use_individual) "individual" else "summary"

  # --- Run SuSiE-inf ----------------------------------------------------------

  start_time <- proc.time()

  fit <- tryCatch({
    if (use_individual) {
      susieR::susie(
        X = X,
        y = y,
        L = L,
        scaled_prior_variance = prior_variance / var(y),
        estimate_residual_variance = estimate_residual_variance,
        estimate_prior_variance = estimate_prior_variance,
        coverage = coverage,
        min_abs_corr = min_abs_corr,
        max_iter = max_iter,
        unmappable_effects = "inf",
        verbose = verbose_susie
      )
    } else {
      susieR::susie_rss(
        z = z,
        R = LD,
        n = n,
        L = L,
        prior_variance = prior_variance,
        estimate_residual_variance = estimate_residual_variance,
        estimate_prior_variance = estimate_prior_variance,
        coverage = coverage,
        min_abs_corr = min_abs_corr,
        max_iter = max_iter,
        unmappable_effects = "inf",
        verbose = verbose_susie
      )
    }
  }, error = function(e) {
    list(error = conditionMessage(e))
  })

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  # --- Handle errors ----------------------------------------------------------

  if (!is.null(fit$error)) {
    return(list(
      pip             = rep(NA_real_, p_snps),
      credible_sets   = list(),
      method          = "susie_inf",
      input_type      = input_type,
      params          = params,
      runtime_seconds = elapsed,
      additional      = list(
        alpha          = matrix(NA_real_, nrow = L, ncol = p_snps),
        posterior_mean = rep(NA_real_, p_snps),
        lbf            = rep(NA_real_, L),
        cs_purity      = list(),
        converged      = FALSE,
        elbo           = NA_real_,
        n_iter_run     = NA_integer_,
        tau2           = NA_real_,
        theta          = NULL
      ),
      error = fit$error
    ))
  }

  # --- Extract standard outputs -----------------------------------------------

  pip <- susieR::susie_get_pip(fit)

  cs_raw <- fit$sets
  if (is.null(cs_raw) || is.null(cs_raw$cs) || length(cs_raw$cs) == 0) {
    credible_sets <- list()
    cs_purity     <- list()
  } else {
    credible_sets <- cs_raw$cs
    cs_purity     <- cs_raw$purity
  }

  # --- Extract SuSiE-inf-specific outputs -------------------------------------

  additional <- list(
    alpha          = fit$alpha,
    posterior_mean = colSums(fit$alpha * fit$mu),
    lbf            = if (!is.null(fit$lbf)) fit$lbf else rep(NA_real_, L),
    cs_purity      = cs_purity,
    converged      = fit$converged,
    elbo           = fit$elbo[length(fit$elbo)],
    n_iter_run     = length(fit$elbo),
    tau2           = if (!is.null(fit$tau2)) fit$tau2 else NA_real_,
    theta          = fit$theta   # NULL if not returned
  )

  # --- Return -----------------------------------------------------------------

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = "susie_inf",
    input_type      = input_type,
    params          = params,
    runtime_seconds = elapsed,
    additional      = additional
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run SuSiE-inf on a single region from simulation data structures
#'
#' Thin adapter that extracts inputs from the simulation's \code{region_geno}
#' and \code{region_pheno} objects and calls \code{\link{run_susie_inf}}.
#'
#' @param region_geno List. One element of \code{simulation$genotypes}.
#' @param region_pheno List. One element of a scenario's \code{regions}.
#' @param use_individual Logical. If TRUE and X/y are available, use
#'   individual-level data. Default: FALSE.
#' @param ... Additional arguments passed to \code{\link{run_susie_inf}}.
#'
#' @return The output of \code{\link{run_susie_inf}}.
#' @export
run_susie_inf_region <- function(region_geno, region_pheno,
                                 use_individual = FALSE, ...) {
  if (use_individual &&
      !is.null(region_geno$X) &&
      !is.null(region_pheno$y)) {
    run_susie_inf(X = region_geno$X, y = region_pheno$y, ...)
  } else {
    run_susie_inf(z = region_pheno$z, LD = region_geno$LD, n = region_geno$n, ...)
  }
}
