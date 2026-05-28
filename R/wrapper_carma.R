# =============================================================================
# carma.R
#
# Wrapper for CARMA (Contextual Adaptive Robust Marginal Analysis,
# Yang et al. 2023, Nature Genetics) fine-mapping.
#
# CARMA is an R package that performs Bayesian fine-mapping from summary
# statistics while accounting for potential LD discrepancies between the
# GWAS sample and the reference panel, via an outlier-detection mechanism.
# It uses a spike-slab (or Cauchy/Hyper-g) prior on effect sizes and fits
# an EM algorithm over causal configurations.
#
# Unlike SuSiE, CARMA does not decompose signals into independent components;
# instead it returns a single global credible set and per-variant PIPs.
#
# This file provides:
#   - setup_carma()        : installs the package if needed and verifies it
#   - run_carma()          : runs CARMA on a single region (explicit inputs)
#   - run_carma_region()   : adapter called by run_methods()
#
# Standard output format:
#   pip              Numeric vector (length p). Marginal PIPs.
#   credible_sets    List containing one integer vector: the global credible
#                    set (1-based variant indices, sorted). Empty list if
#                    CARMA returns no variants.
#   method           Character. "carma".
#   input_type       Character. Always "summary".
#   params           List. Hyperparameters used.
#   runtime_seconds  Numeric. Wall-clock time.
#   additional       List. CARMA-specific outputs (see run_carma() docs).
#
# Reference:
#   Yang Z et al. (2023). CARMA is a new Bayesian model for fine-mapping in
#   genome-wide association meta-analyses. Nature Genetics, 55, 1057-1065.
#   https://doi.org/10.1038/s41588-023-01392-0
# =============================================================================


# =============================================================================
# Setup
# =============================================================================

#' Set up the CARMA R package
#'
#' Checks that CARMA is installed and loads it. If not installed, installs it
#' from GitHub using \code{remotes}.
#'
#' \preformatted{
#'   remotes::install_github("ZikunY/CARMA")
#' }
#'
#' @return Invisible TRUE if CARMA is available.
#' @export
setup_carma <- function() {
  if (!requireNamespace("CARMA", quietly = TRUE)) {
    message("CARMA not found. Installing from GitHub (ZikunY/CARMA)...")
    if (!requireNamespace("remotes", quietly = TRUE)) {
      stop(
        "The 'remotes' package is needed to install CARMA.\n",
        "Install it with: install.packages('remotes')",
        call. = FALSE
      )
    }
    remotes::install_github("ZikunY/CARMA", quiet = TRUE)
  }
  ver <- tryCatch(
    packageDescription("CARMA")[["Version"]],
    error = function(e) "?"
  )
  message("CARMA v", ver, " ready.")
  invisible(TRUE)
}


# =============================================================================
# Run CARMA on a single region
# =============================================================================

#' Run CARMA fine-mapping on a single region
#'
#' Calls \code{CARMA::CARMA()} on a single locus and returns results in the
#' standardised format.
#'
#' CARMA returns one global credible set per region (not one per causal
#' signal as SuSiE does). The credible set contains the minimal set of
#' variants whose joint posterior probability reaches \code{rho.index}.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param n Integer. Sample size. Currently unused by CARMA's core algorithm
#'   but retained for a consistent interface.
#' @param rho.index Numeric. Coverage threshold for the credible set.
#'   Default: 0.95.
#' @param num.causal Integer. Maximum number of causal variants to consider.
#'   Default: 10.
#' @param tau Numeric. Prior variance on effect sizes under the spike-slab
#'   prior. Default: 0.04 (prior SD = 0.2, matching ABF convention).
#' @param effect.size.prior Character. Prior distribution on effect sizes:
#'   \code{"Spike-slab"} (default), \code{"Cauchy"}, or \code{"Hyper-g"}.
#' @param outlier.switch Logical. Enable LD-discrepancy outlier detection.
#'   Default: TRUE.
#' @param all.iter Integer. Number of outer EM iterations. Default: 3.
#'
#' @return A list with the standardised fine-mapping output:
#' \describe{
#'   \item{pip}{Numeric vector (length p). Marginal posterior inclusion
#'     probabilities.}
#'   \item{credible_sets}{List containing one integer vector: the global
#'     credible set (1-based indices, sorted). Empty list if CARMA produced
#'     no credible set.}
#'   \item{method}{Character. Always \code{"carma"}.}
#'   \item{input_type}{Character. Always \code{"summary"}.}
#'   \item{params}{List. Hyperparameters used.}
#'   \item{runtime_seconds}{Numeric. Wall-clock time in seconds.}
#'   \item{additional}{List of CARMA-specific outputs:
#'     \describe{
#'       \item{outliers}{Data frame of detected outlier variants (variants
#'         with discrepancies between z-scores and the LD matrix). Zero rows
#'         if none detected or \code{outlier.switch = FALSE}.}
#'     }
#'   }
#'   \item{error}{Character or NULL. Error message if CARMA failed.}
#' }
#'
#' @export
run_carma <- function(z,
                      LD,
                      n                = NULL,
                      rho.index        = 0.95,
                      num.causal       = 10,
                      tau              = 0.04,
                      effect.size.prior = "Spike-slab",
                      outlier.switch   = TRUE,
                      all.iter         = 3) {

  # --- Validate ---------------------------------------------------------------

  p <- length(z)

  stopifnot(
    "LD must be a p x p matrix" =
      is.matrix(LD) && nrow(LD) == p && ncol(LD) == p,
    "rho.index must be in (0, 1)" =
      is.numeric(rho.index) && rho.index > 0 && rho.index < 1
  )

  if (!requireNamespace("CARMA", quietly = TRUE)) {
    return(.carma_error_result(
      p,
      .carma_params(rho.index, num.causal, tau, effect.size.prior,
                    outlier.switch, all.iter),
      0,
      "CARMA is not installed. Run setup_carma() to install it."
    ))
  }

  params <- .carma_params(rho.index, num.causal, tau, effect.size.prior,
                          outlier.switch, all.iter)

  # lambda = 1/sqrt(p) is the standard CARMA default for the logistic prior
  lambda <- 1 / sqrt(p)

  # --- Run CARMA --------------------------------------------------------------

  start_time <- proc.time()

  fit <- tryCatch({
    # capture.output suppresses CARMA's per-locus timing print statements
    invisible(utils::capture.output(
      res_list <- CARMA::CARMA(
        z.list            = list(z),
        ld.list           = list(LD),
        lambda.list       = list(lambda),
        rho.index         = rho.index,
        num.causal        = as.integer(num.causal),
        tau               = tau,
        effect.size.prior = effect.size.prior,
        outlier.switch    = outlier.switch,
        all.iter          = as.integer(all.iter),
        printing.log      = FALSE
      )
    ))
    res_list[[1]]
  }, error = function(e) {
    list(error = conditionMessage(e))
  })

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  # --- Handle errors ----------------------------------------------------------

  if (!is.null(fit$error)) {
    return(.carma_error_result(p, params, elapsed, fit$error))
  }

  # --- Extract PIPs -----------------------------------------------------------

  pip <- as.numeric(fit$PIPs)
  pip <- pmax(0, pmin(1, pip))

  # --- Extract credible set ---------------------------------------------------
  # CARMA returns one global credible set per locus.
  # Structure: fit[["Credible set"]][[2]] is a list of 1-based variant indices.

  credible_sets <- list()

  cs_raw <- fit[["Credible set"]]
  if (!is.null(cs_raw) && length(cs_raw) >= 2) {
    cs_indices <- sort(as.integer(unlist(cs_raw[[2]])))
    cs_indices <- cs_indices[!is.na(cs_indices) & cs_indices >= 1L & cs_indices <= p]
    if (length(cs_indices) > 0) {
      credible_sets <- list(cs_indices)
    }
  }

  # --- Extract outliers -------------------------------------------------------

  outliers <- fit$Outliers
  if (is.null(outliers)) outliers <- data.frame()

  # --- Return -----------------------------------------------------------------

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = "carma",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      outliers = outliers
    )
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run CARMA on a single region from simulation data structures
#'
#' Thin adapter that extracts the appropriate inputs from the simulation's
#' \code{region_geno} and \code{region_pheno} objects and calls
#' \code{\link{run_carma}}.
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing \code{LD} and \code{n}.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z}.
#' @param ... Additional arguments passed to \code{\link{run_carma}}
#'   (e.g. \code{rho.index}, \code{num.causal}, \code{tau}).
#'
#' @return The output of \code{\link{run_carma}}.
#' @export
run_carma_region <- function(region_geno, region_pheno, ...) {
  run_carma(
    z  = region_pheno$z,
    LD = region_geno$LD,
    n  = region_geno$n,
    ...
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

.carma_params <- function(rho.index, num.causal, tau, effect.size.prior,
                          outlier.switch, all.iter) {
  list(
    rho.index         = rho.index,
    num.causal        = num.causal,
    tau               = tau,
    effect.size.prior = effect.size.prior,
    outlier.switch    = outlier.switch,
    all.iter          = all.iter
  )
}

.carma_error_result <- function(p, params, elapsed, error_msg) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "carma",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(outliers = data.frame()),
    error           = error_msg
  )
}
