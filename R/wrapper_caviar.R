# =============================================================================
# wrapper_caviar.R
#
# Wrapper for CAVIAR (Hormozdiari et al., Genetics 2014), a summary-statistic
# fine-mapping method that enumerates causal configurations under a multivariate
# normal model of the z-scores given LD, and returns per-variant posterior
# probabilities plus a rho-level causal set.
#
# CAVIAR is an external C++ binary. Verified CLI (upstream README):
#   CAVIAR -z <zfile> -l <ldfile> -c <max_causal> -o <out_prefix> [-r <rho>]
#     -z  M lines, TAB separated: <snp_name> <z_score>          (no header)
#     -l  M lines of M SPACE separated LD values                (no header)
#     -c  maximum number of causal variants to enumerate        (default 2)
#     -o  output prefix; writes <prefix>_post and <prefix>_set
#     -r  rho, the causal-set coverage level                    (default 0.95)
#
#   <prefix>_post : per-variant posterior probabilities (header + one row/SNP)
#   <prefix>_set  : the variant names in the rho-level causal set, one per line
#
# STATUS: BUILT AND WORKING, BUT NOT RUN IN THE BENCHMARK (decision 2026-07-29).
# CAVIAR enumerates causal configurations, so its cost is combinatorial in the
# region size: C(p, c) configurations each carrying matrix operations. On the
# verification locus (p = 150, c = 2, ~11k configurations) it was already too
# slow to finish comfortably; the benchmark's regions run to p = 1000, where the
# configuration space is ~45x larger. Across 337,500 fits that would dominate the
# entire allocation for one method.
#
# It is also the least informative addition: CAVIAR (2014) is the direct ancestor
# of FINEMAP and PAINTOR, both already in the benchmark, both using smarter search
# over the same model class.
#
# The wrapper is kept because it is written and structurally correct - if CAVIAR
# is ever wanted on the small (100-200 SNP) regions only, it is ready. It simply
# is not included in FMB_METHODS for the full grid.
#
# This file provides:
#   - setup_caviar()      : verifies the CAVIAR binary is reachable
#   - run_caviar()        : runs CAVIAR on a single region (explicit inputs)
#   - run_caviar_region() : adapter called by run_methods()
#
# Reference:
#   Hormozdiari F, Kostem E, Kang EY, Pasaniuc B, Eskin E (2014). Identifying
#   causal variants at loci with multiple signals of association.
#   Genetics, 198(2), 497-508.
# =============================================================================


#' Set up the CAVIAR binary
#'
#' @param caviar_path Character. Path to the CAVIAR binary or its name on PATH.
#' @return Invisibly, the resolved path.
#' @export
setup_caviar <- function(caviar_path = "CAVIAR") {
  if (file.exists(caviar_path)) {
    message("CAVIAR binary found: ", normalizePath(caviar_path))
    return(invisible(normalizePath(caviar_path)))
  }
  resolved <- Sys.which(caviar_path)
  if (nchar(resolved) > 0) {
    message("CAVIAR binary found on PATH: ", resolved)
    return(invisible(resolved))
  }
  stop(
    "CAVIAR binary not found (searched for: '", caviar_path, "').\n\n",
    "Build it from source:\n",
    "  git clone https://github.com/fhormoz/caviar.git\n",
    "  cd caviar/CAVIAR-C++ && make\n",
    "  # then pass the path via\n",
    "  #   method_args = list(caviar = list(caviar_path = '/path/to/CAVIAR'))\n",
    "CAVIAR needs a C++ toolchain and GSL headers (module load GSL on the HPC).",
    call. = FALSE
  )
}


#' Run CAVIAR fine-mapping on a single region
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param variant_ids Character vector or NULL. Variant names. Embedded
#'   whitespace is collapsed to underscores, because the .z file is
#'   whitespace-delimited and VCF-style ids ("1 12345 . A T") would otherwise
#'   shift the columns (the same trap fixed for FINEMAP and Functional BEATRICE).
#' @param caviar_path Character. Path to the CAVIAR binary. Default "CAVIAR".
#' @param max_causal Integer. Maximum number of causal variants (-c). Default 2.
#' @param rho Numeric. Causal-set coverage level (-r). Default 0.95.
#'
#' @return A list in the benchmark's standard fine-mapping output format.
#' @export
run_caviar <- function(z,
                       LD,
                       variant_ids = NULL,
                       caviar_path = "CAVIAR",
                       max_causal  = 2L,
                       rho         = 0.95) {

  p <- length(z)
  params <- list(max_causal = max_causal, rho = rho, caviar_path = caviar_path)

  stopifnot(
    "LD must be a p x p matrix" = is.matrix(LD) && nrow(LD) == p && ncol(LD) == p,
    "max_causal must be a positive integer" =
      is.numeric(max_causal) && length(max_causal) == 1 && max_causal >= 1
  )

  if (is.null(variant_ids)) variant_ids <- paste0("SNP_", seq_len(p))
  variant_ids <- gsub("\\s+", "_", variant_ids)
  # CAVIAR matches the _set file back to names, so they must be unique.
  if (anyDuplicated(variant_ids)) variant_ids <- paste0(variant_ids, "_", seq_len(p))

  work_dir <- tempfile(pattern = "caviar_run_")
  dir.create(work_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  z_path  <- file.path(work_dir, "region.z")
  ld_path <- file.path(work_dir, "region.ld")
  prefix  <- file.path(work_dir, "out")

  # -z : TAB separated <name> <z>, no header
  utils::write.table(data.frame(variant_ids, as.numeric(z)), z_path,
                     quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  # -l : space separated LD, no header
  utils::write.table(round(LD, 8), ld_path,
                     quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

  args <- c("-z", z_path, "-l", ld_path,
            "-c", as.character(as.integer(max_causal)),
            "-o", prefix, "-r", as.character(rho))

  start_time <- proc.time()
  run_output <- tryCatch(system2(caviar_path, args = args, stdout = TRUE, stderr = TRUE),
                         error = function(e) structure(conditionMessage(e), class = "caviar_error"))
  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  if (inherits(run_output, "caviar_error")) {
    return(.caviar_error_result(p, params, elapsed, as.character(run_output)))
  }

  post_path <- paste0(prefix, "_post")
  set_path  <- paste0(prefix, "_set")
  if (!file.exists(post_path)) {
    return(.caviar_error_result(p, params, elapsed,
      paste(c("CAVIAR produced no _post file.", run_output), collapse = "\n")))
  }

  post <- tryCatch(utils::read.table(post_path, header = TRUE,
                                     stringsAsFactors = FALSE, sep = ""),
                   error = function(e) NULL)
  if (is.null(post) || ncol(post) < 2L) {
    return(.caviar_error_result(p, params, elapsed, "Failed to parse CAVIAR _post file."))
  }

  # Column naming has varied across CAVIAR builds. The per-variant PIP is the
  # LAST numeric column ("Causal_Post._Prob"); the middle column, when present,
  # is Prob_in_pCausalSet. Prefer an explicit name match, else fall back to the
  # last numeric column.
  num_cols <- which(vapply(post, is.numeric, logical(1)))
  pip_col <- grep("causal.*post|post.*prob", names(post), ignore.case = TRUE)
  pip_col <- if (length(pip_col)) pip_col[length(pip_col)] else
             if (length(num_cols)) num_cols[length(num_cols)] else NA_integer_
  if (is.na(pip_col)) {
    return(.caviar_error_result(p, params, elapsed,
      paste("Could not identify the posterior column in CAVIAR _post; columns:",
            paste(names(post), collapse = ", "))))
  }

  # Re-order to the input variant order using the name column (column 1).
  ord <- match(variant_ids, as.character(post[[1L]]))
  if (anyNA(ord)) ord <- seq_len(min(p, nrow(post)))
  pip <- rep(NA_real_, p)
  vals <- suppressWarnings(as.numeric(post[[pip_col]][ord]))
  pip[seq_along(vals)] <- vals
  pip <- pmax(0, pmin(1, pip))
  if (all(is.na(pip))) {
    return(.caviar_error_result(p, params, elapsed, "CAVIAR posterior column parsed to all NA."))
  }

  # _set holds variant NAMES (not indices) in the rho-level causal set.
  credible_sets <- list()
  if (file.exists(set_path)) {
    nm <- readLines(set_path, warn = FALSE)
    nm <- trimws(nm); nm <- nm[nzchar(nm)]
    idx <- match(nm, variant_ids); idx <- idx[!is.na(idx)]
    if (length(idx)) credible_sets <- list(sort(as.integer(idx)))
  }

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = "caviar",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(set_size = if (length(credible_sets)) length(credible_sets[[1]]) else 0L,
                           post_column = names(post)[pip_col])
  )
}


#' Run CAVIAR on a single region from simulation data structures
#'
#' CAVIAR does not use annotations or the sample size - only z-scores and LD.
#'
#' @param region_geno One element of \code{simulation$genotypes}.
#' @param region_pheno One element of a scenario's \code{regions}.
#' @param ... Passed to \code{\link{run_caviar}}.
#' @return The output of \code{\link{run_caviar}}.
#' @export
run_caviar_region <- function(region_geno, region_pheno, ...) {
  run_caviar(
    z           = region_pheno$z,
    LD          = region_geno$LD,
    variant_ids = region_geno$variant_ids,
    ...
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

.caviar_error_result <- function(p, params, elapsed, error_msg) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "caviar",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(set_size = 0L, post_column = NA_character_),
    error           = error_msg
  )
}
