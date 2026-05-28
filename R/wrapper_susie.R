# =============================================================================
# susie.R
#
# Wrapper for SuSiE (Sum of Single Effects) fine-mapping.
#
# This file provides:
#   - setup_susie()        : checks/installs the susieR package
#   - run_susie()          : runs SuSiE on a single region (explicit inputs)
#   - run_susie_region()   : adapter called by run_methods(); extracts inputs
#                            from simulation data structures and calls run_susie()
#
# Standard output format (shared across all methods):
#   pip              Numeric vector (length p). Posterior inclusion probabilities.
#   credible_sets    List of integer vectors (variant indices), or NULL.
#   method           Character. "susie".
#   input_type       Character. "individual" or "summary".
#   params           List. All hyperparameters used.
#   runtime_seconds  Numeric. Wall-clock time.
#   additional       List. SuSiE-specific outputs (see run_susie() docs).
# =============================================================================


# =============================================================================
# Setup
# =============================================================================

#' Check and install susieR
#'
#' Verifies that the susieR package is available. If not, attempts to install
#' it from CRAN.
#'
#' @return Invisible TRUE if susieR is available.
#' @export
setup_susie <- function() {
  if (!requireNamespace("susieR", quietly = TRUE)) {
    message("susieR not found. Attempting to install from CRAN...")
    utils::install.packages("susieR")
    if (!requireNamespace("susieR", quietly = TRUE)) {
      stop(
        "Failed to install susieR. Please install manually:\n",
        "  install.packages('susieR')\n",
        "or from GitHub:\n",
        "  remotes::install_github('stephenslab/susieR')",
        call. = FALSE
      )
    }
    message("susieR installed successfully.")
  } else {
    message("susieR is available (version ",
            utils::packageVersion("susieR"), ").")
  }
  invisible(TRUE)
}


# =============================================================================
# Run SuSiE on a single region
# =============================================================================

#' Run SuSiE fine-mapping on a single region
#'
#' Takes the data for one region (summary statistics + LD, or individual-level
#' data) and runs SuSiE, returning results in the standardised output format.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param n Integer. Sample size.
#' @param X Optional matrix. Individual-level genotype data (n x p).
#'   If provided along with \code{y}, SuSiE is run with individual-level
#'   data (\code{susie}) rather than summary statistics (\code{susie_rss}).
#' @param y Optional numeric vector. Phenotype vector (length n).
#' @param L Integer. Maximum number of single-effect components. Default: 10.
#' @param estimate_residual_variance Logical. Default: TRUE.
#' @param estimate_prior_variance Logical. Default: TRUE.
#' @param prior_variance Numeric or NULL. Prior effect variance. If NULL,
#'   set to \code{0.1 * var(y)} (individual-level) or 0.1 (summary stats).
#' @param coverage Numeric. Coverage level for credible sets. Default: 0.95.
#' @param min_abs_corr Numeric. Purity threshold for filtering credible sets.
#'   Default: 0.5.
#' @param max_iter Integer. Maximum IBSS iterations. Default: 100.
#' @param verbose_susie Logical. Print SuSiE's own progress messages.
#'   Default: FALSE.
#'
#' @return A list with the standardised fine-mapping output:
#' \describe{
#'   \item{pip}{Numeric vector (length p). Posterior inclusion probabilities.}
#'   \item{credible_sets}{List of integer vectors (variant indices per CS),
#'     or an empty list if no credible sets pass the purity filter.}
#'   \item{method}{Character. Always \code{"susie"}.}
#'   \item{input_type}{Character. \code{"individual"} or \code{"summary"}.}
#'   \item{params}{List. Hyperparameters used.}
#'   \item{runtime_seconds}{Numeric. Wall-clock time in seconds.}
#'   \item{additional}{List of SuSiE-specific outputs:
#'     \describe{
#'       \item{alpha}{Matrix (L x p). Per-component posterior assignment
#'         probabilities.}
#'       \item{posterior_mean}{Numeric vector (length p). Posterior mean
#'         effect size summed across components.}
#'       \item{lbf}{Numeric vector (length L). Log Bayes factor per
#'         single-effect component.}
#'       \item{cs_purity}{Data frame. Min/mean/median absolute correlation
#'         for each reported credible set.}
#'       \item{converged}{Logical. Whether IBSS converged.}
#'       \item{elbo}{Numeric. Final ELBO value.}
#'       \item{n_iter_run}{Integer. Number of IBSS iterations run.}
#'     }
#'   }
#'   \item{error}{Character or NULL. Error message if SuSiE failed.}
#' }
#'
#' @export
run_susie <- function(z = NULL,
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
    stop("susieR is not installed. Run setup_susie() first.", call. = FALSE)
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

  # --- Set prior variance default ---------------------------------------------

  if (is.null(prior_variance)) {
    prior_variance <- if (use_individual) 0.1 * var(y) else 0.1
  }

  # --- Collect params ---------------------------------------------------------

  params <- list(
    L = L,
    estimate_residual_variance = estimate_residual_variance,
    estimate_prior_variance = estimate_prior_variance,
    prior_variance = prior_variance,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    max_iter = max_iter
  )

  p_snps <- if (use_individual) ncol(X) else length(z)

  # --- Run SuSiE --------------------------------------------------------------

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
        verbose = verbose_susie
      )
    }
  }, error = function(e) {
    list(error = conditionMessage(e))
  })

  elapsed <- (proc.time() - start_time)["elapsed"]

  input_type <- if (use_individual) "individual" else "summary"

  # --- Handle errors ----------------------------------------------------------

  if (!is.null(fit$error)) {
    return(list(
      pip             = rep(NA_real_, p_snps),
      credible_sets   = list(),
      method          = "susie",
      input_type      = input_type,
      params          = params,
      runtime_seconds = as.numeric(elapsed),
      additional      = list(
        alpha          = matrix(NA_real_, nrow = L, ncol = p_snps),
        posterior_mean = rep(NA_real_, p_snps),
        lbf            = rep(NA_real_, L),
        cs_purity      = list(),
        converged      = FALSE,
        elbo           = NA_real_,
        n_iter_run     = NA_integer_
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

  # --- Extract SuSiE-specific outputs (go into $additional) -------------------

  additional <- list(
    alpha          = fit$alpha,
    posterior_mean = colSums(fit$alpha * fit$mu),
    lbf            = if (!is.null(fit$lbf)) fit$lbf else rep(NA_real_, L),
    cs_purity      = cs_purity,
    converged      = fit$converged,
    elbo           = fit$elbo[length(fit$elbo)],
    n_iter_run     = length(fit$elbo)
  )

  # --- Return -----------------------------------------------------------------

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = "susie",
    input_type      = input_type,
    params          = params,
    runtime_seconds = as.numeric(elapsed),
    additional      = additional
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run SuSiE on a single region from simulation data structures
#'
#' Thin adapter that extracts the appropriate inputs from the simulation's
#' \code{region_geno} and \code{region_pheno} objects and calls
#' \code{\link{run_susie}}. This is the function registered in the method
#' registry and called by \code{\link{run_methods}}.
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing at minimum \code{LD}, \code{n}, and optionally \code{X}.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing at minimum \code{z} and optionally \code{y}.
#' @param use_individual Logical. If TRUE and both \code{X} and \code{y} are
#'   available, use individual-level data instead of summary statistics.
#'   Default: FALSE.
#' @param ... Additional arguments passed to \code{\link{run_susie}}
#'   (e.g. \code{L}, \code{coverage}, \code{min_abs_corr}).
#'
#' @return The output of \code{\link{run_susie}}.
#' @export
run_susie_region <- function(region_geno, region_pheno,
                             use_individual = FALSE, ...) {
  if (use_individual &&
      !is.null(region_geno$X) &&
      !is.null(region_pheno$y)) {
    run_susie(X = region_geno$X, y = region_pheno$y, ...)
  } else {
    run_susie(z = region_pheno$z, LD = region_geno$LD, n = region_geno$n, ...)
  }
}
