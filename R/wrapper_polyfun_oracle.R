# =============================================================================
# polyfun_oracle.R
#
# Wrapper for PolyFun-oracle: PolyFun-style fine-mapping with **true** per-SNP
# priors recovered from the simulator. Serves as the methodological upper
# bound for any annotation-aware method.
#
# Canonical PolyFun (Weissbrod 2020) is a two-stage method: (1) estimate
# per-SNP heritability h2(j) by stratified LDSC on a pre-computed UKB
# annotation matrix; (2) feed those priors into SuSiE. Stage 1 depends on
# real human SNP coordinates, reference LD scores tuned to a specific GWAS
# population, and a ~25 GB pre-computed annotation file — none of which apply
# to simulated benchmark data.
#
# polyfun_oracle skips Stage 1 entirely by using the **true** per-SNP causal
# probability that the simulator used to assign causal variants:
#
#   pi_j ∝ exp(A_j' log gamma)
#
# where A is the region's annotation matrix and gamma is the enrichment
# vector. Both are stored by the simulator (truth$enrichment and
# region_geno$annotations_matrix), so we can reconstruct pi_j exactly.
# These per-SNP priors are then passed to susieR::susie_rss() via
# prior_weights.
#
# Because the prior is exact, polyfun_oracle defines the ceiling: no
# annotation-aware method should outperform it on a non-degenerate setting.
# If one does, either the comparator is cheating or the simulator is leaking
# information through some other channel.
#
# The fine-mapping step is exactly that of polyfun_ldsc (PolyFun + SuSiE as
# in Weissbrod et al. 2020: L = 10, causal effect variance from the modified
# HESS estimator, residual variance estimated; see R/wrapper_polyfun_ldsc.R),
# so the two methods differ only in the prior and the gap between them is the
# cost of estimating it. polyfun_oracle remains a ceiling: it is given the
# simulator's true prior and is not a deployable method.
#
# When no annotations were simulated (annotations = "none"), there is no
# per-SNP signal to reconstruct, and the non-functionally-informed PolyFun +
# SuSiE (uniform prior) is run; prior_source = "uniform_fallback".
#
# This file provides:
#   - setup_polyfun_oracle()        : checks susieR is installed
#   - run_polyfun_oracle()          : runs the method on a single region
#                                      (explicit inputs)
#   - run_polyfun_oracle_region()   : adapter called by run_methods()
#
# Reference:
#   Weissbrod O et al. (2020). Functionally informed fine-mapping and
#   polygenic localization of complex trait heritability enrichment.
#   Nature Genetics, 52(12), 1355-1363.
#   https://doi.org/10.1038/s41588-020-00735-5
# =============================================================================


# =============================================================================
# Setup
# =============================================================================

#' Check that susieR is available for polyfun_oracle
#'
#' polyfun_oracle is built on top of \code{susieR::susie_rss}. This function
#' verifies the dependency.
#'
#' @return Invisible TRUE if susieR is available.
#' @export
setup_polyfun_oracle <- function() {
  if (!requireNamespace("susieR", quietly = TRUE)) {
    stop(
      "susieR is required for polyfun_oracle. Install with:\n",
      "  install.packages('susieR')\n",
      "or from GitHub:\n",
      "  remotes::install_github('stephenslab/susieR')",
      call. = FALSE
    )
  }
  message("susieR ", utils::packageVersion("susieR"),
          " is available for polyfun_oracle.")
  invisible(TRUE)
}


# =============================================================================
# Run polyfun_oracle on a single region
# =============================================================================

#' Run PolyFun-oracle fine-mapping on a single region
#'
#' Reconstructs the true per-SNP causal probability from the supplied
#' annotation matrix and enrichment vector, then runs
#' \code{susieR::susie_rss} with those probabilities as
#' \code{prior_weights}. Defines an upper bound for annotation-aware
#' fine-mapping on the benchmark.
#'
#' If \code{annotations} is \code{NULL} or \code{enrichment} is \code{NULL},
#' the method falls back to a uniform prior (PolyFun's non-functionally
#' informed mode) and records this in the \code{params$prior_source} field. This makes the
#' wrapper safe to call on \code{annotations = "none"} simulations.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param n Integer. Sample size.
#' @param annotations Matrix or NULL. Functional annotation matrix (p x m).
#'   When NULL, the method runs with a uniform prior.
#' @param enrichment Numeric vector or NULL. Per-annotation fold-enrichment
#'   used by the simulator (length m). When NULL, the method runs with a
#'   uniform prior.
#' @param causal_probs Numeric vector or NULL. The simulator's stored
#'   selection probabilities; used in preference to \code{enrichment}.
#' @param L Integer. Maximum number of single-effect components. Default: 10.
#' @param hess Logical. Set the causal effect variance by PolyFun's modified
#'   HESS estimator. Default: TRUE.
#' @param hess_iter Integer. HESS draws to average. Default: 100.
#' @param max_iter Integer. Maximum IBSS iterations. Default: 100.
#' @param ... Ignored (for wrapper compatibility).
#'
#' The fine-mapping step is \code{\link{run_polyfun_ldsc}}'s, so that
#' polyfun_oracle and polyfun_ldsc differ only in the prior.
#'
#' @return A list with the standardised fine-mapping output. The
#'   \code{additional} sub-list contains:
#' \describe{
#'   \item{prior_source}{Character. \code{"oracle"} when annotation priors
#'     were used, \code{"uniform_fallback"} otherwise.}
#'   \item{prior_weights}{Numeric vector (length p). The per-SNP prior
#'     probabilities actually passed to SuSiE (sums to 1).}
#'   \item{enrichment_used}{The \code{enrichment} vector that was used to
#'     compute the prior, or NULL if uniform fallback.}
#'   \item{alpha}{Matrix (L x p). Per-component posterior assignment
#'     probabilities from SuSiE.}
#'   \item{posterior_mean}{Numeric vector (length p). Posterior mean effect
#'     summed across components.}
#'   \item{lbf}{Numeric vector (length L). Log Bayes factor per component.}
#'   \item{cs_purity}{Data frame or list. Min / mean / median absolute
#'     correlation per credible set.}
#'   \item{converged}{Logical. Whether IBSS converged.}
#'   \item{elbo}{Numeric. Final ELBO value.}
#'   \item{n_iter_run}{Integer. Number of IBSS iterations run.}
#' }
#'
#' @export
run_polyfun_oracle <- function(z,
                                LD,
                                n,
                                annotations = NULL,
                                enrichment  = NULL,
                                causal_probs = NULL,   # ITER-005 (temporary)
                                L           = 10,
                                hess        = TRUE,
                                hess_iter   = 100L,
                                max_iter    = 100L,
                                ...) {

  # --- Validate ---------------------------------------------------------------

  p <- length(z)

  stopifnot(
    "LD must be a p x p matrix" =
      is.matrix(LD) && nrow(LD) == p && ncol(LD) == p,
    "n must be a positive integer" =
      is.numeric(n) && length(n) == 1 && n > 0
  )

  if (!requireNamespace("susieR", quietly = TRUE)) {
    return(.polyfun_oracle_error_result(
      p, .polyfun_oracle_params(L, hess, hess_iter, max_iter,
                                prior_source = "missing_susier"),
      0, "susieR is not installed. Run setup_polyfun_oracle() first."
    ))
  }

  # --- Reconstruct per-SNP priors --------------------------------------------
  # If we have annotations + enrichment, recover the simulator's per-SNP
  # causal probabilities. Otherwise fall back to uniform.

  # ITERATION 005 (TEMPORARY - see R/simulate_phenotypes.R banner). Prefer the
  # simulator's STORED selection probabilities when present. Reconstructing
  # pi_j = exp(A' log gamma) below is only correct for the additive
  # relationship; under the co-occurrence, non-monotone, mixed-direction or
  # threshold forms it produces the WRONG prior, and this method silently stops
  # being a ceiling while still being reported as one. Remove this branch when
  # Iteration 005 is done and the relationships are withdrawn.
  if (!is.null(annotations) && !is.null(causal_probs) && length(causal_probs) == p &&
      all(is.finite(causal_probs)) && sum(causal_probs) > 0) {
    prior_weights   <- causal_probs / sum(causal_probs)
    prior_source    <- "oracle_stored_probs"
    enrichment_used <- if (is.numeric(enrichment)) enrichment else NULL
  } else {

  use_oracle <- !is.null(annotations) && !is.null(enrichment) &&
                is.matrix(annotations) && nrow(annotations) == p &&
                is.numeric(enrichment) && length(enrichment) == ncol(annotations)

  if (use_oracle) {
    # Clamp enrichments to a tiny positive value to keep log well-defined
    enrich_safe <- pmax(as.numeric(enrichment), .Machine$double.eps)
    log_enrich  <- log(enrich_safe)
    log_w       <- as.numeric(annotations %*% log_enrich)
    # Numerical stability: subtract max before exp
    w           <- exp(log_w - max(log_w))
    prior_weights <- w / sum(w)
    prior_source  <- "oracle"
    enrichment_used <- enrich_safe
  } else {
    prior_weights <- rep(1 / p, p)
    prior_source  <- "uniform_fallback"
    enrichment_used <- NULL
  }
  }   # closes the ITERATION 005 stored-probs branch

  params <- .polyfun_oracle_params(L, hess, hess_iter, max_iter,
                                   prior_source = prior_source)

  # --- PolyFun's fine-mapping step with the true prior -------------------------

  start_time <- proc.time()
  fit <- tryCatch(
    run_polyfun_ldsc(
      z = z, LD = LD, n = n,
      prior = if (prior_source == "uniform_fallback") NULL else prior_weights,
      L = L, hess = hess, hess_iter = hess_iter, max_iter = max_iter),
    error = function(e) list(error = conditionMessage(e)))
  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  if (!is.null(fit$error)) {
    return(.polyfun_oracle_error_result(p, params, elapsed, fit$error))
  }

  fit$method          <- "polyfun_oracle"
  fit$params          <- params
  fit$runtime_seconds <- elapsed
  fit$additional$prior_source    <- prior_source
  fit$additional$prior_weights   <- prior_weights
  fit$additional$enrichment_used <- enrichment_used
  fit
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run polyfun_oracle on a single region from simulation data structures
#'
#' Thin adapter that pulls the annotation matrix from
#' \code{region_geno$annotations_matrix} and the enrichment vector from
#' \code{region_pheno$truth$enrichment}, then calls
#' \code{\link{run_polyfun_oracle}}. If either is missing, the method falls
#' back to a uniform prior automatically.
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing at minimum \code{LD}, \code{n}, and optionally
#'   \code{annotations_matrix}.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing at minimum \code{z} and (when annotations are simulated)
#'   \code{truth$enrichment}. \code{annotations_matrix} on this list is also
#'   accepted as a fallback location for the matrix.
#' @param ... Additional arguments passed to \code{\link{run_polyfun_oracle}}
#'   (e.g. \code{L}, \code{hess}).
#'
#' @return The output of \code{\link{run_polyfun_oracle}}.
#' @export
run_polyfun_oracle_region <- function(region_geno, region_pheno, ...) {

  # Annotation matrix is stored on either side depending on pipeline path:
  # - run_simulation() places it on region_geno (shared across scenarios)
  # - simulate_phenotypes() places a copy on region_pheno
  # Prefer the geno-side matrix; fall back to pheno-side.
  A <- region_geno$annotations_matrix
  if (is.null(A)) A <- region_pheno$annotations_matrix

  # Enrichment lives on truth and may be NULL when annotations = "none".
  enrich <- if (!is.null(region_pheno$truth))
    region_pheno$truth$enrichment else NULL

  # ITERATION 005 (TEMPORARY): the simulator stores the exact selection
  # probabilities it used. Under any non-additive relationship these are the
  # ONLY correct oracle prior; reconstructing exp(A' log gamma) would be wrong.
  cprobs <- if (!is.null(region_pheno$truth))
    region_pheno$truth$causal_probs else NULL

  run_polyfun_oracle(
    z            = region_pheno$z,
    LD           = region_geno$LD,
    n            = region_geno$n,
    annotations  = A,
    enrichment   = enrich,
    causal_probs = cprobs,
    ...
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

.polyfun_oracle_params <- function(L, hess, hess_iter, max_iter,
                                    prior_source) {
  list(
    L            = L,
    hess         = hess,
    hess_iter    = hess_iter,
    max_iter     = max_iter,
    prior_source = prior_source
  )
}

.polyfun_oracle_error_result <- function(p, params, elapsed, error_msg) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "polyfun_oracle",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      prior_source    = params$prior_source %||% NA_character_,
      prior_weights   = rep(NA_real_, p),
      enrichment_used = NULL
    ),
    error           = error_msg
  )
}
