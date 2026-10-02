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
#   - setup_carma()              : installs the package if needed and verifies it
#   - run_carma()                : runs CARMA on a single region without
#                                   annotations (explicit inputs)
#   - run_carma_annotated()      : runs CARMA with annotations over several
#                                   regions (explicit inputs)
#   - run_carma_scenario_setup() : per-scenario hook for run_methods()
#   - run_carma_region()         : adapter called by run_methods()
#
# Annotation mode
# ---------------
# When the simulated regions carry annotations, run_methods() fits CARMA once
# per scenario on all regions together, passing each region's annotations as
# cbind(1, A) in w.list, as the CARMA vignette does for several loci. CARMA
# links Pr(causal) to the annotations by a logistic model fitted inside its EM
# algorithm; the M-step regression pools all loci in the call (Yang et al.
# 2023 fit it "at the chromosome level"). CARMA's own defaults are kept for the
# annotation model (EM.dist = "Logistic", input.alpha = 0, i.e. ridge).
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

  params <- .carma_params(rho.index, num.causal, tau, effect.size.prior,
                          outlier.switch, all.iter)

  if (!.carma_available()) {
    return(.carma_error_result(
      p, params, 0,
      "CARMA is not installed. Run setup_carma() to install it."
    ))
  }

  # lambda = 1/sqrt(p) is the standard CARMA default for the logistic prior
  lambda <- 1 / sqrt(p)

  # Per-fit scratch directory for CARMA's side-effect output files (see the
  # output.labels note in .carma_call()).
  out_dir <- tempfile(pattern = "carma_out_")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # --- Run CARMA --------------------------------------------------------------

  start_time <- proc.time()

  fit <- tryCatch({
    .carma_call(
      z.list            = list(z),
      ld.list           = list(LD),
      lambda.list       = list(lambda),
      rho.index         = rho.index,
      num.causal        = as.integer(num.causal),
      tau               = tau,
      effect.size.prior = effect.size.prior,
      outlier.switch    = outlier.switch,
      all.iter          = as.integer(all.iter),
      output.labels     = out_dir
    )[[1]]
  }, error = function(e) {
    list(error = conditionMessage(e))
  })

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  # --- Handle errors ----------------------------------------------------------

  if (!is.null(fit$error)) {
    return(.carma_error_result(p, params, elapsed, fit$error))
  }

  # --- Return -----------------------------------------------------------------

  out <- .carma_extract(fit, p)
  list(
    pip             = out$pip,
    credible_sets   = out$credible_sets,
    method          = "carma",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      outliers        = out$outliers,
      annotation_mode = "none"
    )
  )
}


# =============================================================================
# Run CARMA with functional annotations over several regions
# =============================================================================

#' Run CARMA with functional annotations over several regions
#'
#' Calls \code{CARMA::CARMA()} once on all supplied regions with their
#' annotations, as the CARMA vignette does for several loci: each region's
#' annotation matrix is passed as \code{cbind(1, A)} in \code{w.list}.
#' CARMA then links the prior probability of causality to the annotations by
#' a logistic model, \eqn{\mathrm{logit}\,\Pr(\gamma_i = 1) = w_i^\top\theta},
#' fitted by its EM algorithm with a penalised (by default ridge) logistic
#' regression whose data are pooled over all the regions in the call
#' (Yang et al. 2023; CARMA's M-step concatenates the loci). All other
#' settings are those of \code{\link{run_carma}}, including
#' \code{lambda = 1/sqrt(p)} per region, and CARMA's own defaults
#' (\code{EM.dist = "Logistic"}, \code{prior.prob.computation = "Logistic"}).
#'
#' @param z_list List of numeric vectors. Marginal z-scores, one per region.
#' @param ld_list List of LD matrices, one per region.
#' @param annot_list List of numeric annotation matrices (\code{p_r x d}), one
#'   per region, all with the same \code{d} columns in the same order. CARMA
#'   standardises each column within a region.
#' @param rho.index,num.causal,tau,effect.size.prior,outlier.switch,all.iter
#'   As in \code{\link{run_carma}}.
#' @param input.alpha Numeric. Elastic-net mixing parameter of CARMA's
#'   annotation regression. Default: 0 (ridge), CARMA's default.
#'
#' @return A list with one standardised result per region, in input order.
#'   \code{additional} holds \code{outliers}, \code{annotation_mode}
#'   (\code{"annotated"}), \code{annotation_coef} (CARMA's last fitted
#'   intercept and annotation coefficients, when it wrote them) and
#'   \code{n_regions_pooled}. On failure every element carries the same
#'   \code{error}. The runtime of the pooled call is divided equally between
#'   the regions.
#' @export
run_carma_annotated <- function(z_list,
                                ld_list,
                                annot_list,
                                rho.index         = 0.95,
                                num.causal        = 10,
                                tau               = 0.04,
                                effect.size.prior = "Spike-slab",
                                outlier.switch    = TRUE,
                                all.iter          = 3,
                                input.alpha       = 0) {

  R     <- length(z_list)
  p_vec <- lengths(z_list)

  stopifnot(
    "z_list must contain at least one region" = R >= 1L,
    "ld_list and annot_list must have one entry per region" =
      length(ld_list) == R && length(annot_list) == R,
    "rho.index must be in (0, 1)" =
      is.numeric(rho.index) && rho.index > 0 && rho.index < 1,
    "input.alpha must be in [0, 1]" =
      is.numeric(input.alpha) && length(input.alpha) == 1 &&
      input.alpha >= 0 && input.alpha <= 1
  )
  for (r in seq_len(R)) {
    if (!is.matrix(ld_list[[r]]) || nrow(ld_list[[r]]) != p_vec[r] ||
        ncol(ld_list[[r]]) != p_vec[r]) {
      stop("LD matrix ", r, " must be ", p_vec[r], " x ", p_vec[r], call. = FALSE)
    }
  }
  problem <- .carma_annotation_problem(annot_list, p_vec)
  if (!is.null(problem)) stop(problem, call. = FALSE)
  d <- ncol(annot_list[[1]])

  params <- c(.carma_params(rho.index, num.causal, tau, effect.size.prior,
                            outlier.switch, all.iter),
              list(input.alpha = input.alpha, n_annotations = d))

  fail <- function(msg, elapsed = 0) {
    lapply(seq_len(R), function(r) {
      .carma_error_result(p_vec[r], params, elapsed / R, msg)
    })
  }

  if (!.carma_available()) {
    return(fail("CARMA is not installed. Run setup_carma() to install it."))
  }

  # CARMA replaces the first column of each w.list matrix by an intercept and
  # standardises the rest, so the vignette passes cbind(1, annotations).
  w_list <- lapply(annot_list, function(A) {
    A <- as.matrix(A)
    storage.mode(A) <- "double"
    cbind(1, A)
  })

  out_dir <- tempfile(pattern = "carma_out_")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  labels <- as.list(paste0("region", seq_len(R)))

  start_time <- proc.time()
  fits <- tryCatch(
    .carma_call(
      z.list            = z_list,
      ld.list           = ld_list,
      w.list            = w_list,
      lambda.list       = as.list(1 / sqrt(p_vec)),
      label.list        = labels,
      rho.index         = rho.index,
      num.causal        = as.integer(num.causal),
      tau               = tau,
      effect.size.prior = effect.size.prior,
      outlier.switch    = outlier.switch,
      all.iter          = as.integer(all.iter),
      input.alpha       = input.alpha,
      output.labels     = out_dir
    ),
    error = function(e) structure(conditionMessage(e), class = "carma_error")
  )
  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  if (inherits(fits, "carma_error")) return(fail(as.character(fits), elapsed))
  if (!is.list(fits) || length(fits) != R) {
    return(fail("CARMA did not return one result per region.", elapsed))
  }

  # CARMA writes the coefficients of its last annotation fit, the same for
  # every locus, to post_<label>_theta.txt (intercept first).
  theta_path <- file.path(out_dir, paste0("post_", labels[[1]], "_theta.txt"))
  annotation_coef <- if (file.exists(theta_path)) {
    v <- tryCatch(as.numeric(utils::read.table(theta_path)[[1]]),
                  error = function(e) NULL)
    if (!is.null(v) && length(v) == d + 1L) {
      stats::setNames(v, c("(Intercept)",
                           colnames(annot_list[[1]]) %||% paste0("A", seq_len(d))))
    } else v
  } else NULL

  lapply(seq_len(R), function(r) {
    out <- .carma_extract(fits[[r]], p_vec[r])
    list(
      pip             = out$pip,
      credible_sets   = out$credible_sets,
      method          = "carma",
      input_type      = "summary",
      params          = params,
      runtime_seconds = elapsed / R,
      additional      = list(
        outliers         = out$outliers,
        annotation_mode  = "annotated",
        annotation_coef  = annotation_coef,
        n_regions_pooled = R
      )
    )
  })
}


# =============================================================================
# Scenario-level setup: fit CARMA with annotations on all regions together
# =============================================================================

#' Run CARMA with annotations once per scenario
#'
#' Called by \code{\link{run_methods}} once per scenario. When the regions
#' carry annotations (and \code{use_annotations} is not \code{FALSE}), it runs
#' \code{\link{run_carma_annotated}} on all regions of the scenario together,
#' so that CARMA's annotation model is fitted to the pooled regions, and
#' returns the per-region results keyed by a fingerprint of each region's
#' z-scores. Without annotations it returns an empty list and every region is
#' fine-mapped on its own by \code{\link{run_carma}}, as before. A failed
#' annotation run is passed on as \code{.carma_error}, so those regions are
#' reported as failures rather than silently fine-mapped without annotations.
#'
#' @param genotypes List. \code{simulation$genotypes}.
#' @param regions List. One scenario's \code{regions}.
#' @param user_args List. The user's method_args for carma.
#'
#' @return A list with \code{.carma_cache} or \code{.carma_error}, or an empty
#'   list.
#' @export
run_carma_scenario_setup <- function(genotypes, regions, user_args) {

  if (isFALSE(user_args$use_annotations)) return(list())

  R <- length(regions)
  A_list <- lapply(seq_len(R), function(i) {
    .carma_extract_annotations(genotypes[[i]], regions[[i]])
  })
  if (all(vapply(A_list, is.null, logical(1)))) return(list())

  p_vec   <- vapply(regions, function(x) length(x$z), integer(1))
  problem <- .carma_annotation_problem(A_list, p_vec)
  if (!is.null(problem)) {
    return(list(.carma_error = paste("CARMA annotation run not possible:", problem)))
  }

  keep <- intersect(names(user_args),
                    setdiff(names(formals(run_carma_annotated)),
                            c("z_list", "ld_list", "annot_list")))
  fits <- tryCatch(
    do.call(run_carma_annotated, c(
      list(z_list     = lapply(regions, `[[`, "z"),
           ld_list    = lapply(genotypes[seq_len(R)], `[[`, "LD"),
           annot_list = A_list),
      user_args[keep]
    )),
    error = function(e) conditionMessage(e)
  )
  if (is.character(fits)) {
    return(list(.carma_error = paste("CARMA annotation run failed:", fits)))
  }
  if (!is.null(fits[[1]]$error)) {
    return(list(.carma_error = paste("CARMA annotation run failed:",
                                     fits[[1]]$error)))
  }

  keys <- vapply(regions, function(x) .fb_fingerprint(x$z), character(1))
  list(.carma_cache = stats::setNames(fits, keys))
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run CARMA on a single region from simulation data structures
#'
#' Adapter called by \code{\link{run_methods}}. When the scenario hook
#' \code{\link{run_carma_scenario_setup}} has fitted CARMA with annotations,
#' the region's result is taken from it. When the region carries annotations
#' but no hook result is available (a direct call), CARMA is fitted with the
#' annotations of this region alone. Otherwise the region is fine-mapped
#' without annotations by \code{\link{run_carma}}.
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing \code{LD}, \code{n} and optionally \code{annotations_matrix}.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z}.
#' @param .carma_cache,.carma_error Set by
#'   \code{\link{run_carma_scenario_setup}}; not for direct use.
#' @param use_annotations Logical. Use the regions' annotations when present.
#'   Default: TRUE. Set to FALSE to reproduce runs made before the annotation
#'   model was added (Iterations 004 to 007).
#' @param input.alpha Passed to \code{\link{run_carma_annotated}}.
#' @param ... Additional arguments passed to \code{\link{run_carma}}
#'   (e.g. \code{rho.index}, \code{num.causal}, \code{tau}).
#'
#' @return A standardised result list; \code{additional$annotation_mode}
#'   records whether annotations were used.
#' @export
run_carma_region <- function(region_geno, region_pheno,
                             .carma_cache    = NULL,
                             .carma_error    = NULL,
                             use_annotations = TRUE,
                             input.alpha     = 0,
                             ...) {
  p <- length(region_pheno$z)
  if (!is.null(.carma_error)) {
    return(.carma_error_result(p, .carma_dots_params(...), 0, .carma_error))
  }
  if (!is.null(.carma_cache)) {
    hit <- .carma_cache[[.fb_fingerprint(region_pheno$z)]]
    if (is.null(hit)) {
      return(.carma_error_result(
        p, .carma_dots_params(...), 0,
        "CARMA annotation run returned no result for this region."))
    }
    return(hit)
  }

  A <- .carma_extract_annotations(region_geno, region_pheno)
  if (isTRUE(use_annotations) && !is.null(A)) {
    return(run_carma_annotated(
      z_list      = list(region_pheno$z),
      ld_list     = list(region_geno$LD),
      annot_list  = list(A),
      input.alpha = input.alpha,
      ...
    )[[1]])
  }

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

# params for an error result raised before run_carma() is reached.
.carma_dots_params <- function(...) {
  a <- list(...)
  .carma_params(a$rho.index %||% 0.95, a$num.causal %||% 10, a$tau %||% 0.04,
                a$effect.size.prior %||% "Spike-slab",
                a$outlier.switch %||% TRUE, a$all.iter %||% 3)
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

# Separate functions so the tests can stand in for CARMA.
.carma_available <- function() requireNamespace("CARMA", quietly = TRUE)

# CRITICAL: CARMA's output.labels defaults to '.', i.e. the CURRENT WORKING
# DIRECTORY, and CARMA writes four files per locus into it
# (post_<label>.txt, post_<label>_poi_likeli.txt,
#  post_<label>_poi_gamma.mtx, post_<label>_outliers.rds), plus
# post_<label>_theta.txt when annotations are used. Under the PBS array the
# cwd is the project root and ~500 tasks run concurrently with the same
# default label, so leaving this unset both pollutes the repo and races on
# identical filenames. Callers pass a per-fit temp directory that they delete.
# capture.output suppresses CARMA's per-locus timing print statements.
.carma_call <- function(...) {
  res <- NULL
  invisible(utils::capture.output(
    res <- CARMA::CARMA(..., printing.log = FALSE)
  ))
  res
}

# Standardised pieces of one locus of CARMA's output.
.carma_extract <- function(fit, p) {
  pip <- pmax(0, pmin(1, as.numeric(fit$PIPs)))

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

  outliers <- fit$Outliers
  if (is.null(outliers)) outliers <- data.frame()

  list(pip = pip, credible_sets = credible_sets, outliers = outliers)
}

# Prefer region_geno$annotations_matrix, falling back to region_pheno's copy
# (same convention as the PAINTOR and Funmap wrappers).
.carma_extract_annotations <- function(region_geno, region_pheno) {
  A <- region_geno$annotations_matrix
  if (is.null(A)) A <- region_pheno$annotations_matrix
  A
}

# NULL when the annotation matrices can be passed to CARMA, otherwise why not.
.carma_annotation_problem <- function(annot_list, p_vec) {
  if (any(vapply(annot_list, is.null, logical(1)))) {
    return("some regions have no annotation matrix.")
  }
  ok <- vapply(seq_along(annot_list), function(r) {
    A <- annot_list[[r]]
    (is.matrix(A) || is.data.frame(A)) && nrow(A) == p_vec[r] &&
      ncol(A) >= 1L && is.numeric(as.matrix(A)) && all(is.finite(as.matrix(A)))
  }, logical(1))
  if (!all(ok)) {
    return("each annotation matrix must be numeric and finite, with one row per variant and at least one column.")
  }
  if (length(unique(vapply(annot_list, ncol, integer(1)))) != 1L) {
    return("the regions' annotation matrices have different numbers of columns.")
  }
  NULL
}
