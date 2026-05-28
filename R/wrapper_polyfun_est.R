# =============================================================================
# polyfun_est.R
#
# Wrapper for PolyFun-style fine-mapping with **estimated** per-SNP priors.
#
# Provides a fair, "no-cheating" annotation-aware comparator to Funmap,
# PAINTOR, Functional BEATRICE, etc. The methodology is:
#
#   1. Estimate per-annotation contributions tau_k by regressing the
#      observed chi-squared statistics on the annotation matrix:
#
#          chi2_j = tau_0 + sum_k tau_k * A_{j,k} + noise
#
#      Non-negative coefficients are kept (negative heritability
#      contributions are not physical). This is "LDSC-lite" — much
#      simpler than canonical stratified LDSC, which would also account
#      for LD scores and reference-panel calibration. Those steps are
#      unnecessary on simulated benchmark data and would tie us to
#      pretrained UKB infrastructure that does not apply here.
#
#   2. Convert estimated coefficients into per-SNP prior variances:
#
#          sigma2(j) = max(tau_0 + sum_k tau_k * A_{j,k}, floor)
#
#      Normalised so the resulting vector sums to 1, then passed to
#      susieR::susie_rss(..., prior_weights = ...).
#
# Compared to polyfun_oracle (which uses the simulator's true per-SNP
# causal probabilities), polyfun_est pays the realistic "cost of
# estimation". The gap between the two is informative.
#
# Pooling across regions
# ----------------------
# Per-region tau estimates are statistically noisy with small p and small
# m. When run_methods() drives this wrapper, it calls
# run_polyfun_est_scenario_setup() once per scenario, which fits a single
# tau on the concatenated chi^2 and annotation rows across all regions and
# threads it back into each region's call via the `pooled_tau` argument.
#
# When the wrapper is invoked manually with `pooled_tau = NULL` (the
# default), tau is estimated from the supplied single-region data only.
# When `pooled_tau` is supplied, it is used directly and no per-region fit
# happens.
#
# When annotations are absent (annotations = NULL), the wrapper falls back
# to a uniform prior and reports prior_source = "uniform_fallback".
#
# This file provides:
#   - setup_polyfun_est()                    : checks susieR is installed
#   - run_polyfun_est()                      : runs on a single region
#                                              (explicit inputs)
#   - run_polyfun_est_region()               : adapter called by run_methods()
#   - run_polyfun_est_scenario_setup()       : scenario-level hook that pools
#                                              tau across all regions
#
# Reference (inspiration, not reimplementation):
#   Weissbrod O et al. (2020). Nature Genetics 52, 1355-1363.
# =============================================================================


# =============================================================================
# Setup
# =============================================================================

#' Check that susieR is available for polyfun_est
#'
#' @return Invisible TRUE if susieR is available.
#' @export
setup_polyfun_est <- function() {
  if (!requireNamespace("susieR", quietly = TRUE)) {
    stop(
      "susieR is required for polyfun_est. Install with:\n",
      "  install.packages('susieR')\n",
      "or from GitHub:\n",
      "  remotes::install_github('stephenslab/susieR')",
      call. = FALSE
    )
  }
  message("susieR ", utils::packageVersion("susieR"),
          " is available for polyfun_est.")
  invisible(TRUE)
}


# =============================================================================
# Run polyfun_est on a single region
# =============================================================================

#' Run PolyFun-style fine-mapping on a single region
#'
#' Estimates per-annotation contributions by regressing the chi-squared
#' statistics on the annotation matrix, converts them into per-SNP prior
#' variances, and runs \code{susieR::susie_rss} with those priors as
#' \code{prior_weights}.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param n Integer. Sample size.
#' @param annotations Matrix or NULL. Functional annotation matrix (p x m).
#'   When NULL, the method runs with a uniform prior.
#' @param pooled_tau Numeric vector or NULL. Pre-computed per-annotation
#'   coefficients (length m + 1, where the first entry is the intercept).
#'   When supplied, the per-region LDSC-lite fit is skipped and these
#'   coefficients are used directly to build the priors. Typically supplied
#'   by \code{run_polyfun_est_scenario_setup()} via \code{run_methods()}.
#'   Default: NULL.
#' @param L Integer. Maximum number of single-effect components. Default: 10.
#' @param coverage Numeric. Coverage level for credible sets. Default: 0.95.
#' @param min_abs_corr Numeric. Purity threshold for filtering credible sets.
#'   Default: 0.5.
#' @param max_iter Integer. Maximum IBSS iterations. Default: 100.
#' @param estimate_residual_variance Logical. Default: TRUE.
#' @param estimate_prior_variance Logical. Default: TRUE.
#' @param verbose_susie Logical. Print SuSiE's own progress messages.
#'   Default: FALSE.
#'
#' @return A list with the standardised fine-mapping output. The
#'   \code{additional} sub-list contains:
#' \describe{
#'   \item{prior_source}{Character. \code{"estimated_per_region"} when tau was
#'     fit on this region alone, \code{"estimated_pooled"} when supplied via
#'     \code{pooled_tau}, or \code{"uniform_fallback"} when no annotations
#'     were available.}
#'   \item{prior_weights}{Numeric vector (length p) of per-SNP prior
#'     probabilities (sums to 1).}
#'   \item{tau_hat}{Numeric vector (length m + 1). Estimated intercept and
#'     per-annotation contributions used to build the priors. NULL when
#'     uniform fallback.}
#'   \item{alpha, posterior_mean, lbf, cs_purity, converged, elbo,
#'     n_iter_run}{Standard SuSiE outputs.}
#' }
#'
#' @export
run_polyfun_est <- function(z,
                             LD,
                             n,
                             annotations = NULL,
                             pooled_tau  = NULL,
                             L           = 10,
                             coverage    = 0.95,
                             min_abs_corr = 0.5,
                             max_iter    = 100,
                             estimate_residual_variance = TRUE,
                             estimate_prior_variance    = TRUE,
                             verbose_susie = FALSE) {

  # --- Validate ---------------------------------------------------------------

  p <- length(z)

  stopifnot(
    "LD must be a p x p matrix" =
      is.matrix(LD) && nrow(LD) == p && ncol(LD) == p,
    "n must be a positive integer" =
      is.numeric(n) && length(n) == 1 && n > 0,
    "coverage must be a single number in (0, 1)" =
      is.numeric(coverage) && length(coverage) == 1 &&
      coverage > 0 && coverage < 1
  )

  if (!requireNamespace("susieR", quietly = TRUE)) {
    return(.polyfun_est_error_result(
      p, L,
      .polyfun_est_params(L, coverage, min_abs_corr, max_iter,
                          estimate_residual_variance,
                          estimate_prior_variance,
                          prior_source = "missing_susier"),
      0, "susieR is not installed. Run setup_polyfun_est() first."
    ))
  }

  # --- Resolve per-SNP prior weights -----------------------------------------

  have_annotations <- !is.null(annotations) &&
                      is.matrix(annotations) &&
                      nrow(annotations) == p

  if (!have_annotations) {
    prior_weights <- rep(1 / p, p)
    prior_source  <- "uniform_fallback"
    tau_hat       <- NULL
  } else if (!is.null(pooled_tau)) {
    if (length(pooled_tau) != ncol(annotations) + 1L) {
      return(.polyfun_est_error_result(
        p, L,
        .polyfun_est_params(L, coverage, min_abs_corr, max_iter,
                            estimate_residual_variance,
                            estimate_prior_variance,
                            prior_source = "bad_pooled_tau"),
        0, sprintf(
          "pooled_tau has length %d but annotation matrix has %d columns; expected %d.",
          length(pooled_tau), ncol(annotations), ncol(annotations) + 1L
        )
      ))
    }
    tau_hat       <- as.numeric(pooled_tau)
    prior_weights <- .tau_to_prior_weights(annotations, tau_hat)
    prior_source  <- "estimated_pooled"
  } else {
    tau_hat       <- .fit_ldsc_lite(z^2, annotations)
    prior_weights <- .tau_to_prior_weights(annotations, tau_hat)
    prior_source  <- "estimated_per_region"
  }

  params <- .polyfun_est_params(L, coverage, min_abs_corr, max_iter,
                                 estimate_residual_variance,
                                 estimate_prior_variance,
                                 prior_source = prior_source)

  # --- Run SuSiE-RSS ----------------------------------------------------------

  start_time <- proc.time()

  fit <- tryCatch({
    susieR::susie_rss(
      z = z,
      R = LD,
      n = n,
      L = L,
      prior_weights = prior_weights,
      estimate_residual_variance = estimate_residual_variance,
      estimate_prior_variance    = estimate_prior_variance,
      coverage     = coverage,
      min_abs_corr = min_abs_corr,
      max_iter     = max_iter,
      verbose      = verbose_susie
    )
  }, error = function(e) {
    list(error = conditionMessage(e))
  })

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  if (!is.null(fit$error)) {
    return(.polyfun_est_error_result(p, L, params, elapsed, fit$error))
  }

  # --- Extract outputs -------------------------------------------------------

  pip <- susieR::susie_get_pip(fit)

  cs_raw <- fit$sets
  if (is.null(cs_raw) || is.null(cs_raw$cs) || length(cs_raw$cs) == 0) {
    credible_sets <- list()
    cs_purity     <- list()
  } else {
    credible_sets <- cs_raw$cs
    cs_purity     <- cs_raw$purity
  }

  additional <- list(
    prior_source   = prior_source,
    prior_weights  = prior_weights,
    tau_hat        = tau_hat,
    alpha          = fit$alpha,
    posterior_mean = colSums(fit$alpha * fit$mu),
    lbf            = if (!is.null(fit$lbf)) fit$lbf else rep(NA_real_, L),
    cs_purity      = cs_purity,
    converged      = fit$converged,
    elbo           = fit$elbo[length(fit$elbo)],
    n_iter_run     = length(fit$elbo)
  )

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = "polyfun_est",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = additional
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run polyfun_est on a single region from simulation data structures
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing \code{LD}, \code{n}, and optionally \code{annotations_matrix}.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z}. \code{annotations_matrix} may also be here.
#' @param pooled_tau Numeric vector or NULL. Pre-computed pooled tau.
#'   Typically supplied automatically by
#'   \code{run_polyfun_est_scenario_setup()}. Default: NULL.
#' @param ... Additional arguments passed to \code{\link{run_polyfun_est}}.
#'
#' @return The output of \code{\link{run_polyfun_est}}.
#' @export
run_polyfun_est_region <- function(region_geno, region_pheno,
                                    pooled_tau = NULL, ...) {

  A <- region_geno$annotations_matrix
  if (is.null(A)) A <- region_pheno$annotations_matrix

  run_polyfun_est(
    z           = region_pheno$z,
    LD          = region_geno$LD,
    n           = region_geno$n,
    annotations = A,
    pooled_tau  = pooled_tau,
    ...
  )
}


# =============================================================================
# Scenario-level setup hook
# =============================================================================

#' Fit a pooled tau across all regions of a scenario for polyfun_est
#'
#' Concatenates \code{z^2} values and annotation matrix rows across every
#' region in the scenario, fits a single LDSC-lite regression, and returns
#' the resulting coefficients in the format expected by
#' \code{run_polyfun_est}.
#'
#' Called automatically by \code{run_methods()} once per scenario when
#' polyfun_est is in the methods list. The returned list is merged into the
#' user's \code{method_args} for that scenario, so every region call gets
#' the same pooled coefficients.
#'
#' If annotation matrices are missing or inconsistent across regions
#' (different m), the function returns an empty list, in which case
#' \code{run_polyfun_est} falls back to per-region estimation.
#'
#' @param genotypes List. The simulation's \code{genotypes} list (all regions).
#' @param regions List. The scenario's per-region phenotype objects.
#' @param user_args List. The user's method_args for this method.
#'
#' @return A named list with entry \code{pooled_tau} (numeric vector of length
#'   m + 1) if pooling succeeded, otherwise an empty list.
#' @export
run_polyfun_est_scenario_setup <- function(genotypes, regions, user_args) {

  n_regions <- length(regions)
  if (n_regions == 0L) return(list())

  chi2_list <- vector("list", n_regions)
  A_list    <- vector("list", n_regions)
  ncols     <- integer(n_regions)

  for (i in seq_len(n_regions)) {
    A_i <- genotypes[[i]]$annotations_matrix
    if (is.null(A_i)) A_i <- regions[[i]]$annotations_matrix
    if (is.null(A_i) || !is.matrix(A_i)) return(list())

    z_i <- regions[[i]]$z
    if (is.null(z_i) || length(z_i) != nrow(A_i)) return(list())

    chi2_list[[i]] <- z_i^2
    A_list[[i]]    <- A_i
    ncols[i]       <- ncol(A_i)
  }

  if (length(unique(ncols)) != 1L) return(list())

  chi2_pooled <- unlist(chi2_list, use.names = FALSE)
  A_pooled    <- do.call(rbind, A_list)

  tau <- tryCatch(.fit_ldsc_lite(chi2_pooled, A_pooled),
                  error = function(e) NULL)
  if (is.null(tau)) return(list())

  list(pooled_tau = tau)
}


# =============================================================================
# Internal helpers
# =============================================================================

# LDSC-lite: regress chi-squared statistics on annotations.
# Returns a numeric vector of length m + 1: c(intercept, tau_1, ..., tau_m),
# clamped to be non-negative (negative per-annotation heritability is not
# physical and would push the per-SNP prior negative).
.fit_ldsc_lite <- function(chi2, A) {
  df <- data.frame(chi2 = as.numeric(chi2), as.data.frame(A))
  fit <- stats::lm(chi2 ~ ., data = df,
                    singular.ok = TRUE)
  coefs <- stats::coef(fit)
  coefs[is.na(coefs)] <- 0
  # Clamp negatives to zero (non-physical) but retain the intercept
  # baseline as a positive offset to keep priors strictly positive.
  coefs <- pmax(coefs, 0)
  as.numeric(coefs)
}

# Convert tau coefficients to a normalised per-SNP prior weight vector.
# sigma2(j) = tau_0 + sum_k tau_k * A_{j,k}; floored at a small positive
# value to avoid zero probabilities; normalised to sum to 1.
.tau_to_prior_weights <- function(A, tau) {
  p <- nrow(A)
  m <- ncol(A)
  stopifnot(length(tau) == m + 1L)
  intercept <- tau[1L]
  slope     <- tau[-1L]
  sigma2    <- intercept + as.numeric(A %*% slope)
  floor_val <- max(.Machine$double.eps, 1e-8 * mean(abs(sigma2)))
  sigma2    <- pmax(sigma2, floor_val)
  sigma2 / sum(sigma2)
}


# Standard parameter list for polyfun_est outputs.
.polyfun_est_params <- function(L, coverage, min_abs_corr, max_iter,
                                 estimate_residual_variance,
                                 estimate_prior_variance,
                                 prior_source) {
  list(
    L                          = L,
    coverage                   = coverage,
    min_abs_corr               = min_abs_corr,
    max_iter                   = max_iter,
    estimate_residual_variance = estimate_residual_variance,
    estimate_prior_variance    = estimate_prior_variance,
    prior_source               = prior_source
  )
}

# Standard error result.
.polyfun_est_error_result <- function(p, L, params, elapsed, error_msg) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "polyfun_est",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      prior_source   = params$prior_source %||% NA_character_,
      prior_weights  = rep(NA_real_, p),
      tau_hat        = NULL,
      alpha          = matrix(NA_real_, nrow = L, ncol = p),
      posterior_mean = rep(NA_real_, p),
      lbf            = rep(NA_real_, L),
      cs_purity      = list(),
      converged      = FALSE,
      elbo           = NA_real_,
      n_iter_run     = NA_integer_
    ),
    error           = error_msg
  )
}
