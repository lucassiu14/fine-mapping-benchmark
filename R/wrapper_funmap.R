# =============================================================================
# funmap.R
#
# Wrapper for Funmap (Li et al., Bioinformatics 2025), a fine-mapping method
# that integrates functional annotations via a random effects model.
#
# Funmap is a Python package called via reticulate. Users must have:
#   - Python >= 3.6
#   - The funmap Python package installed (GitHub only, see setup_funmap())
#   - numpy, scipy, pandas, matplotlib
#
# Funmap REQUIRES functional annotations. If no annotations are available
# for a region, run_funmap_region() returns an informative error result
# rather than attempting to run the method.
#
# Reference:
#   Li Y et al. (2025). Funmap: integrating functional annotations for
#   fine-mapping. Bioinformatics, 41(1), btaf017.
#   https://doi.org/10.1093/bioinformatics/btaf017
#
# This file provides:
#   - setup_funmap()        : checks Python, reticulate, and the funmap module
#   - run_funmap()          : runs Funmap on a single region (explicit inputs)
#   - run_funmap_region()   : adapter called by run_methods()
# =============================================================================


# =============================================================================
# Setup
# =============================================================================

#' Check Python, reticulate, and the Funmap module are available
#'
#' Verifies that the \code{reticulate} R package is installed, Python is
#' accessible, and the \code{funmap} Python package can be imported. Prints
#' installation instructions for any missing component.
#'
#' To install the Funmap Python package:
#' \preformatted{
#'   # In a terminal:
#'   git clone https://github.com/LeeHITsz/Funmap.git
#'   cd Funmap
#'   pip install -r requirements.txt
#'   pip install .
#' }
#'
#' @return Invisible TRUE if all requirements are met.
#' @export
setup_funmap <- function(python = NULL) {

  # --- reticulate -------------------------------------------------------------

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      "The 'reticulate' R package is required to run Funmap.\n",
      "Install it with:\n",
      "  install.packages('reticulate')",
      call. = FALSE
    )
  }

  # --- Python -----------------------------------------------------------------

  if (!is.null(python)) {
    reticulate::use_python(python, required = TRUE)
  }

  if (!reticulate::py_available(initialize = TRUE)) {
    stop(
      "Python is not available. reticulate could not find a Python installation.\n\n",
      "Pass the path to your Python explicitly:\n",
      "  setup_funmap(python = '/opt/anaconda3/bin/python3')\n\n",
      "Or configure reticulate manually first:\n",
      "  reticulate::use_python('/path/to/python')\n",
      "  reticulate::use_condaenv('your-env')",
      call. = FALSE
    )
  }

  # --- funmap Python package --------------------------------------------------

  if (!reticulate::py_module_available("funmap")) {
    stop(
      "The 'funmap' Python package is not installed.\n\n",
      "Install it from GitHub:\n",
      "  pip install git+https://github.com/LeeHITsz/Funmap.git\n\n",
      "Make sure you are installing into the same Python that reticulate will use.\n",
      "If using a specific Python, pass it to setup_funmap():\n",
      "  setup_funmap(python = '/opt/anaconda3/bin/python3')",
      call. = FALSE
    )
  }

  py_ver <- reticulate::py_config()$version
  message("All Funmap requirements met (Python ", py_ver, ", funmap module found).")
  invisible(TRUE)
}


# =============================================================================
# Run Funmap on a single region
# =============================================================================

#' Run Funmap fine-mapping on a single region
#'
#' Calls the \code{FUNMAP()} Python function via reticulate, passing z-scores,
#' LD matrix, sample size, and a functional annotations matrix. Requires the
#' \code{funmap} Python package (see \code{\link{setup_funmap}}).
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param n Integer. Sample size.
#' @param annotations Matrix. Functional annotations matrix (p x K), where K
#'   is the number of annotations. Numeric values; binary or continuous
#'   annotations are both supported.
#' @param L Integer. Maximum number of causal components. Default: 10.
#' @param max_iter Integer. Maximum iterations for convergence. Default: 100.
#' @param tol Numeric. Convergence tolerance. Default: 5e-5.
#'
#' @return A list with the standardised fine-mapping output:
#' \describe{
#'   \item{pip}{Numeric vector (length p). Posterior inclusion probabilities.}
#'   \item{credible_sets}{List of integer vectors. One per identified causal
#'     component (e.g. L0, L1, ...). Indices are 1-based (R convention).
#'     Empty list if no credible sets pass the purity filter.}
#'   \item{method}{Character. Always \code{"funmap"}.}
#'   \item{input_type}{Character. Always \code{"summary"}.}
#'   \item{params}{List. Hyperparameters used.}
#'   \item{runtime_seconds}{Numeric.}
#'   \item{additional}{List of Funmap-specific outputs:
#'     \describe{
#'       \item{alpha}{Matrix (L x p). Posterior expectations of inclusion
#'         indicators per component.}
#'       \item{posterior_mean}{Numeric vector (length p). Posterior mean
#'         effect size summed across components.}
#'       \item{sigma2}{Numeric. Estimated residual variance.}
#'       \item{converged}{Logical.}
#'       \item{cs_purity}{Data frame or NULL. Min/mean/median absolute
#'         correlations per credible set.}
#'     }
#'   }
#'   \item{error}{Character or NULL. Error message if Funmap failed.}
#' }
#'
#' @export
run_funmap <- function(z,
                       LD,
                       n,
                       annotations,
                       L        = 10,
                       max_iter = 100,
                       tol      = 5e-5,
                       python   = NULL) {

  # --- Validate ---------------------------------------------------------------

  p <- length(z)

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(.funmap_error_result(p, L, list(L = L, max_iter = max_iter, tol = tol), 0,
      "reticulate is not installed. Run setup_funmap() for instructions."))
  }

  if (!is.null(python)) {
    reticulate::use_python(python, required = TRUE)
  }

  if (!reticulate::py_module_available("funmap")) {
    return(.funmap_error_result(p, L, list(L = L, max_iter = max_iter, tol = tol), 0,
      "funmap Python package not found. Run setup_funmap(python=...) for instructions."))
  }

  stopifnot(
    "LD must be a p x p matrix" =
      is.matrix(LD) && nrow(LD) == p && ncol(LD) == p,
    "annotations must be a numeric matrix with p rows" =
      is.matrix(annotations) && is.numeric(annotations) && nrow(annotations) == p,
    "n must be a positive integer" =
      is.numeric(n) && length(n) == 1 && n > 0
  )

  params <- list(
    L        = L,
    max_iter = max_iter,
    tol      = tol
  )

  # --- Import Python modules --------------------------------------------------

  funmap_mod <- reticulate::import("funmap")
  np         <- reticulate::import("numpy")

  # --- Convert inputs to numpy arrays -----------------------------------------

  z_np    <- np$array(z,           dtype = np$float64)
  LD_np   <- np$array(LD,          dtype = np$float64)
  A_np    <- np$array(annotations, dtype = np$float64)

  # --- Run Funmap -------------------------------------------------------------

  start_time <- proc.time()

  fit <- tryCatch({
    funmap_mod$FUNMAP(
      z        = z_np,
      R        = LD_np,
      A        = A_np,
      n        = as.integer(n),
      L        = as.integer(L),
      max_iter = as.integer(max_iter),
      tol      = tol,
      verbose  = FALSE
    )
  }, error = function(e) {
    list(error = conditionMessage(e))
  })

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  # --- Handle errors ----------------------------------------------------------

  # fit is an R list (with $error) only if the tryCatch above caught an error;
  # a successful FUNMAP() call returns a Python ResultFunmap object.
  if (is.list(fit) && !is.null(fit$error)) {
    return(.funmap_error_result(p, L, params, elapsed, fit$error))
  }

  # --- Extract standard outputs -----------------------------------------------

  pip <- as.numeric(reticulate::py_to_r(fit$pip))

  # Credible sets: result.sets['cs'] is a dict like {"L0": array([5,12]), ...}
  # Indices are 0-based in Python; convert to 1-based for R.
  credible_sets <- list()
  cs_purity     <- NULL

  sets_raw <- fit$sets
  if (!is.null(sets_raw)) {
    cs_dict <- sets_raw[["cs"]]
    if (!is.null(cs_dict) && length(cs_dict) > 0) {
      credible_sets <- lapply(cs_dict, function(idx_arr) {
        sort(as.integer(reticulate::py_to_r(idx_arr)) + 1L)
      })
      # purity is a pandas DataFrame; convert to R data frame
      purity_raw <- sets_raw[["purity"]]
      if (!is.null(purity_raw)) {
        cs_purity <- tryCatch(
          reticulate::py_to_r(purity_raw),
          error = function(e) NULL
        )
      }
    }
  }

  # --- Extract Funmap-specific outputs ----------------------------------------

  alpha_raw      <- reticulate::py_to_r(fit$alpha)   # L x p
  mu_raw         <- reticulate::py_to_r(fit$mu)       # L x p
  posterior_mean <- colSums(alpha_raw * mu_raw)

  converged <- tryCatch(
    as.logical(reticulate::py_to_r(fit$converged)),
    error = function(e) NA
  )

  sigma2 <- tryCatch(
    as.numeric(reticulate::py_to_r(fit$sigma2)),
    error = function(e) NA_real_
  )

  # --- Return -----------------------------------------------------------------

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = "funmap",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      alpha          = alpha_raw,
      posterior_mean = posterior_mean,
      sigma2         = sigma2,
      converged      = converged,
      cs_purity      = cs_purity
    )
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run Funmap on a single region from simulation data structures
#'
#' Extracts inputs from \code{region_geno} and \code{region_pheno} and calls
#' \code{\link{run_funmap}}. If no annotations matrix is present in
#' \code{region_pheno}, returns an informative error result without running
#' the method.
#'
#' Annotations are produced by \code{run_simulation()} when
#' \code{annotations = "binary"} or \code{annotations = "continuous"} is
#' specified. If \code{annotations = "none"} was used (the default), Funmap
#' cannot be run.
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing \code{LD} and \code{n}.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z} and \code{annotations_matrix}.
#' @param ... Additional arguments passed to \code{\link{run_funmap}}
#'   (e.g. \code{L}, \code{max_iter}, \code{tol}).
#'
#' @return The output of \code{\link{run_funmap}}, or an error result if no
#'   annotations are available.
#' @export
run_funmap_region <- function(region_geno, region_pheno, ...) {

  if (is.null(region_pheno$annotations_matrix)) {
    p <- length(region_pheno$z)
    return(.funmap_error_result(
      p     = p,
      L     = list(...)$L %||% 10L,
      params = list(...),
      elapsed = 0,
      error_msg = paste(
        "No annotations available for this region.",
        "Funmap requires functional annotations.",
        "Re-run simulation with annotations = 'binary' or annotations = 'continuous'."
      )
    ))
  }

  run_funmap(
    z           = region_pheno$z,
    LD          = region_geno$LD,
    n           = region_geno$n,
    annotations = region_pheno$annotations_matrix,
    ...
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

.funmap_error_result <- function(p, L, params, elapsed, error_msg) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "funmap",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      alpha          = matrix(NA_real_, nrow = L, ncol = p),
      posterior_mean = rep(NA_real_, p),
      sigma2         = NA_real_,
      converged      = FALSE,
      cs_purity      = NULL
    ),
    error = error_msg
  )
}

# %||% is defined in R/utils.R.
