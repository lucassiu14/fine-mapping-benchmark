# =============================================================================
# oracle_susie.R
#
# SuSiE given the **true** per-variant prior the simulator used to choose causal
# variants. The methodological ceiling for any annotation-aware method.
#
# This is not PolyFun and uses none of PolyFun's machinery. PolyFun's
# contribution is Stage 1: estimating a per-SNP prior from annotations by
# stratified LD score regression. That stage is exactly what an oracle replaces,
# so what remains is SuSiE with a known prior, run through run_susie() at the
# same settings the `susie` method uses.
#
# RENAMED AND RE-POINTED 2026-10-17. It was oracle_susie, and fine-mapped
# through run_polyfun_susie(), which fixes the causal effect variance to the
# modified HESS estimate rather than fitting it. That made it comparable to
# polyfun_susie but not to susie. The comparison that matters is the other one:
# susie and oracle_susie now differ in the prior and in nothing else, so the gap
# between them is what knowing the prior is worth. RESULTS STORED BEFORE THAT
# DATE ARE NOT COMPARABLE TO RESULTS FROM THIS VERSION.
#
# The prior
# ---------
# The simulator draws causal variants with probability
#
#   pi_j = exp(A_j' log gamma) / sum_i exp(A_i' log gamma)
#
# and stores that vector. This method uses the stored vector where present and
# reconstructs it from the annotation matrix and the enrichment vector
# otherwise. The reconstruction is correct only for the additive relationship;
# under the alternative annotation-causality forms the stored vector is the only
# correct source, which is why the simulator records it.
#
# Because the prior is exact, oracle_susie is a ceiling: no annotation-aware
# method on this backbone should beat it on a non-degenerate setting. If one
# does, either that comparator is cheating or the simulator is leaking
# information through another channel.
#
# It reads ground truth, so it is NOT deployable and NOT a fair comparator.
#
# Where no annotations were simulated there is no per-SNP signal to recover and
# a uniform prior is used (prior_source = "uniform_fallback"); the method is
# then identical to susie, which is a useful check on the plumbing.
#
# This file provides:
#   - setup_oracle_susie()        : dependency check
#   - run_oracle_susie()          : runs the method on a single region
#   - run_oracle_susie_region()   : adapter called by run_methods()
# =============================================================================

# =============================================================================
# Setup
# =============================================================================

#' Check that susieR is available for oracle_susie
#'
#' oracle_susie is built on top of \code{susieR::susie_rss}. This function
#' verifies the dependency.
#'
#' @return Invisible TRUE if susieR is available.
#' @export
setup_oracle_susie <- function() {
  if (!requireNamespace("susieR", quietly = TRUE)) {
    stop(
      "susieR is required for oracle_susie. Install with:\n",
      "  install.packages('susieR')\n",
      "or from GitHub:\n",
      "  remotes::install_github('stephenslab/susieR')",
      call. = FALSE
    )
  }
  message("susieR ", utils::packageVersion("susieR"),
          " is available for oracle_susie.")
  invisible(TRUE)
}


# =============================================================================
# Run oracle_susie on a single region
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
#' @param max_iter Integer. Maximum IBSS iterations. Default: 100.
#' @param ... Ignored (for wrapper compatibility).
#'
#' The fine-mapping step is \code{\link{run_susie}}'s, so that \code{susie}
#' and \code{oracle_susie} differ in the prior and in nothing else.
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
run_oracle_susie <- function(z,
                                LD,
                                n,
                                annotations = NULL,
                                enrichment  = NULL,
                                causal_probs = NULL,   # ITER-005 (temporary)
                                L           = 10,
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
    return(.oracle_susie_error_result(
      p, .oracle_susie_params(L, max_iter,
                                prior_source = "missing_susier"),
      0, "susieR is not installed. Run setup_oracle_susie() first."
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

  params <- .oracle_susie_params(L, max_iter,
                                   prior_source = prior_source)

  # --- SuSiE, with the true prior -------------------------------------------
  # Plain SuSiE as run_susie() configures it, differing from the `susie` method
  # in the prior weights and nothing else. That is the point of the comparison:
  # susie -> oracle_susie isolates what knowing the prior is worth, with the
  # fine-mapping model held fixed.

  start_time <- proc.time()
  fit <- tryCatch(
    run_susie(
      z = z, LD = LD, n = n, L = L, max_iter = max_iter,
      prior_weights = if (prior_source == "uniform_fallback") NULL else prior_weights),
    error = function(e) list(error = conditionMessage(e)))
  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  if (!is.null(fit$error)) {
    return(.oracle_susie_error_result(p, params, elapsed, fit$error))
  }

  fit$method          <- "oracle_susie"
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

#' Run oracle_susie on a single region from simulation data structures
#'
#' Thin adapter that pulls the annotation matrix from
#' \code{region_geno$annotations_matrix} and the enrichment vector from
#' \code{region_pheno$truth$enrichment}, then calls
#' \code{\link{run_oracle_susie}}. If either is missing, the method falls
#' back to a uniform prior automatically.
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing at minimum \code{LD}, \code{n}, and optionally
#'   \code{annotations_matrix}.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing at minimum \code{z} and (when annotations are simulated)
#'   \code{truth$enrichment}. \code{annotations_matrix} on this list is also
#'   accepted as a fallback location for the matrix.
#' @param ... Additional arguments passed to \code{\link{run_oracle_susie}}
#'   (e.g. \code{L}).
#'
#' @return The output of \code{\link{run_oracle_susie}}.
#' @export
run_oracle_susie_region <- function(region_geno, region_pheno, ...) {

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

  run_oracle_susie(
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

.oracle_susie_params <- function(L, max_iter,
                                    prior_source) {
  list(
    L            = L,
    max_iter     = max_iter,
    prior_source = prior_source
  )
}

.oracle_susie_error_result <- function(p, params, elapsed, error_msg) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "oracle_susie",
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
