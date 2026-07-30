# =============================================================================
# wrapper_finimom.R
#
# Wrapper for FiniMOM (Karhunen et al., Bioinformatics 2023), fine-mapping from
# summary statistics using a non-local product inverse-moment prior on the
# non-zero effects and a beta-binomial prior on the model dimension. The
# non-local prior penalises effects near zero, which is designed to improve
# separation of multiple causal variants in high LD.
#
# FiniMOM is an R package (vkarhune/finimom). It is summary-statistic based:
# it needs beta, se, effect-allele frequency, an LD matrix and n.
#
# INSTALL NOTE (verified 2026-07-29). The package's src/Makevars pins
#   CXX_STD = CXX11
# but current RcppArmadillo requires C++14 or later, so a plain
# remotes::install_github("vkarhune/finimom") FAILS with
#   "*** C++14 compiler required; enable C++14 mode in your compiler"
# Setting CXX_STD in ~/.R/Makevars does NOT help - the package's own Makevars
# wins. Patch the source instead; see setup_finimom() for the exact recipe.
#
# This file provides:
#   - setup_finimom()      : checks/install-instructs the finimom package
#   - run_finimom()        : runs FiniMOM on a single region (explicit inputs)
#   - run_finimom_region() : adapter called by run_methods()
#
# Standard output format:
#   pip              Numeric vector (length p). PIPs, in the input variant order.
#   credible_sets    List of integer vectors (1-based). FiniMOM returns R-native
#                    1-based indices, so no index-base conversion is applied.
#   method           Character. "finimom".
#   input_type       Character. Always "summary".
#   params           List. Hyperparameters used.
#   runtime_seconds  Numeric.
#   additional       List. signals table and the resolved insampleLD flag.
#
# Reference:
#   Karhunen V, Launonen I, Jarvelin MR, Sebert S, Sillanpaa MJ (2023).
#   Genetic fine-mapping from summary data using a nonlocal prior improves the
#   detection of multiple causal variants. Bioinformatics, 39(7), btad396.
# =============================================================================


#' Check that the finimom package is available
#'
#' @return Invisible TRUE if finimom can be loaded.
#' @export
setup_finimom <- function() {
  if (requireNamespace("finimom", quietly = TRUE)) {
    message("finimom v", utils::packageDescription("finimom")[["Version"]], " ready.")
    return(invisible(TRUE))
  }
  stop(
    "The 'finimom' package is not installed.\n\n",
    "A plain remotes::install_github('vkarhune/finimom') FAILS because the\n",
    "package pins CXX_STD = CXX11 in src/Makevars while current\n",
    "RcppArmadillo requires C++14+. Patch the source and install locally:\n\n",
    "  git clone --depth 1 https://github.com/vkarhune/finimom.git\n",
    "  sed -i 's/^CXX_STD = CXX11/CXX_STD = CXX17/' finimom/src/Makevars\n",
    "  R -e \"install.packages('finimom', repos = NULL, type = 'source')\"\n\n",
    "(Editing ~/.R/Makevars does NOT work - the package's own Makevars wins.)",
    call. = FALSE
  )
}


#' Run FiniMOM fine-mapping on a single region
#'
#' @param z Numeric vector. Marginal z-scores (length p). Used to derive
#'   \code{beta}/\code{se} when those are not supplied.
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param n Integer. GWAS sample size.
#' @param beta,se Numeric vectors or NULL. Marginal effect estimates and their
#'   standard errors. When NULL they are derived as \code{se = 1/sqrt(n)} and
#'   \code{beta = z * se}, i.e. the per-standard-deviation scale the benchmark's
#'   simulator uses. Verified empirically to recover planted causals with
#'   FiniMOM's default \code{standardize = TRUE}.
#' @param eaf Numeric vector or NULL. Effect-allele frequencies. When NULL,
#'   0.3 is used for every variant. The benchmark passes the region's MAF, which
#'   equals the effect-allele frequency whenever the effect allele is the minor
#'   allele; FiniMOM uses it only for effect-size standardisation.
#' @param variant_ids Character vector or NULL. Unused by FiniMOM (indices are
#'   positional) but accepted for interface consistency.
#' @param insample_ld Logical or NULL. Whether \code{LD} is in-sample. NULL lets
#'   FiniMOM decide. \code{run_finimom_region()} sets this automatically.
#' @param maxsize Integer. Maximum model size. Default: 10.
#' @param niter,burnin Integer. MCMC iterations and burn-in. Defaults: 12500, 2500.
#' @param cs_level Numeric. Credible set level. Default: 0.95.
#' @param purity Numeric or NULL. Credible set purity threshold.
#' @param clump Logical. Clump highly-correlated variants. Default: TRUE.
#' @param check_ld Logical. Run FiniMOM's LD-discrepancy check. Default: FALSE.
#' @param seed Integer. RNG seed. Default: 456.
#'
#' @return A list in the benchmark's standard fine-mapping output format.
#' @export
run_finimom <- function(z,
                        LD,
                        n,
                        beta        = NULL,
                        se          = NULL,
                        eaf         = NULL,
                        variant_ids = NULL,
                        insample_ld = NULL,
                        maxsize     = 10L,
                        niter       = 12500L,
                        burnin      = 2500L,
                        cs_level    = 0.95,
                        purity      = NULL,
                        clump       = TRUE,
                        check_ld    = FALSE,
                        seed        = 456L) {

  p <- length(z)
  params <- list(maxsize = maxsize, niter = niter, burnin = burnin,
                 cs_level = cs_level, purity = purity, clump = clump,
                 check_ld = check_ld, seed = seed, insample_ld = insample_ld, n = n)

  if (!requireNamespace("finimom", quietly = TRUE)) {
    return(.finimom_error_result(p, params, 0,
      "finimom is not installed. Run setup_finimom() for the install recipe."))
  }

  stopifnot(
    "LD must be a p x p matrix" = is.matrix(LD) && nrow(LD) == p && ncol(LD) == p,
    "n must be a positive integer" = is.numeric(n) && length(n) == 1 && n > 0
  )

  # Derive beta/se on the per-SD scale when not supplied (see @param se).
  if (is.null(se))   se   <- rep(1 / sqrt(n), p)
  if (is.null(beta)) beta <- z * se
  if (is.null(eaf))  eaf  <- rep(0.3, p)
  # FiniMOM standardises using eaf; degenerate frequencies would divide by zero.
  eaf <- pmin(pmax(as.numeric(eaf), 1e-4), 1 - 1e-4)

  args <- list(beta = as.numeric(beta), se = as.numeric(se), eaf = eaf,
               R = LD, n = as.integer(n),
               pip = TRUE, cs = TRUE, cs_level = cs_level,
               maxsize = as.integer(maxsize),
               niter = as.integer(niter), burnin = as.integer(burnin),
               seed = as.integer(seed), clump = clump, check_ld = check_ld,
               verbose = FALSE)
  if (!is.null(insample_ld)) args$insampleLD <- isTRUE(insample_ld)
  if (!is.null(purity))      args$purity     <- purity

  start_time <- proc.time()
  fit <- tryCatch(do.call(finimom::finimom, args),
                  error = function(e) structure(conditionMessage(e), class = "finimom_error"))
  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  if (inherits(fit, "finimom_error")) {
    return(.finimom_error_result(p, params, elapsed, as.character(fit)))
  }
  if (is.null(fit$pip) || length(fit$pip) != p) {
    return(.finimom_error_result(
      p, params, elapsed,
      sprintf("finimom returned no usable pip vector (length %s, expected %d).",
              if (is.null(fit$pip)) "NULL" else length(fit$pip), p)))
  }

  pip <- pmax(0, pmin(1, as.numeric(fit$pip)))

  # FiniMOM's $sets are already 1-based R indices - do NOT add 1.
  credible_sets <- list()
  if (!is.null(fit$sets) && length(fit$sets) > 0L) {
    credible_sets <- lapply(fit$sets, function(s) sort(as.integer(s)))
    credible_sets <- credible_sets[lengths(credible_sets) > 0L]
  }

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = "finimom",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(signals = fit$signals,
                           insample_ld_used = insample_ld,
                           n_sets = length(credible_sets))
  )
}


#' Run FiniMOM on a single region from simulation data structures
#'
#' Passes the simulator's \code{beta_hat}/\code{se}/\code{maf} straight through,
#' and sets FiniMOM's \code{insampleLD} flag automatically: \code{run_simulation()}
#' attaches \code{X_ref} to a region exactly when a reference panel was drawn
#' (in which case \code{region_geno$LD} is \code{cor(X_ref)}), so the LD is
#' in-sample precisely when \code{X_ref} is absent.
#'
#' @param region_geno One element of \code{simulation$genotypes}.
#' @param region_pheno One element of a scenario's \code{regions}.
#' @param ... Passed to \code{\link{run_finimom}}.
#' @return The output of \code{\link{run_finimom}}.
#' @export
run_finimom_region <- function(region_geno, region_pheno, ...) {
  dots <- list(...)
  insample <- if (!is.null(dots$insample_ld)) dots$insample_ld
              else is.null(region_geno$X_ref)
  dots$insample_ld <- NULL
  do.call(run_finimom, c(list(
    z           = region_pheno$z,
    LD          = region_geno$LD,
    n           = region_geno$n,
    beta        = region_pheno$beta_hat,
    se          = region_pheno$se,
    eaf         = region_geno$maf,
    variant_ids = region_geno$variant_ids,
    insample_ld = insample
  ), dots))
}


# =============================================================================
# Internal helpers
# =============================================================================

.finimom_error_result <- function(p, params, elapsed, error_msg) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "finimom",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(signals = NULL, insample_ld_used = NULL, n_sets = 0L),
    error           = error_msg
  )
}
