# =============================================================================
# paintor.R
#
# Wrapper for PAINTOR v3.0 (Kichaev et al. 2014, 2015) fine-mapping.
#
# PAINTOR is an external C++ binary that uses a probabilistic model to compute
# posterior inclusion probabilities, optionally integrating functional
# annotations. Unlike most methods it explicitly learns annotation enrichment
# weights via an EM algorithm.
#
# Annotations are optional here: if none are supplied PAINTOR runs with an
# intercept term only (equivalent to a uniform prior).
#
# This file provides:
#   - setup_paintor()        : verifies the PAINTOR binary is accessible
#   - run_paintor()          : runs PAINTOR on a single region (explicit inputs)
#   - run_paintor_region()   : adapter called by run_methods()
#
# Standard output format:
#   pip              Numeric vector (length p). Marginal PIPs.
#   credible_sets    List of integer vectors. Constructed greedily from PIPs
#                    (variants sorted by descending PIP, accumulated until
#                    cumulative sum >= coverage). One credible set returned.
#   method           Character. "paintor".
#   input_type       Character. Always "summary".
#   params           List. Hyperparameters used.
#   runtime_seconds  Numeric. Wall-clock time.
#   additional       List. PAINTOR-specific outputs (see run_paintor() docs).
#
# References:
#   Kichaev G et al. (2014). Integrating functional data to prioritize causal
#   variants in statistical fine-mapping studies. PLoS Genet, 10(10), e1004722.
#   Kichaev G & Pasaniuc B (2015). Leveraging functional-annotation data in
#   trans-ethnic fine-mapping studies. Am J Hum Genet, 97(2), 260-271.
# =============================================================================


# =============================================================================
# Setup
# =============================================================================

#' Set up the PAINTOR binary
#'
#' Checks that the PAINTOR binary can be found at the supplied path or on the
#' system PATH. If not found, prints installation instructions.
#'
#' PAINTOR must be compiled from source or installed via conda - no
#' pre-compiled binary is distributed. Typical installation:
#'
#' \preformatted{
#'   # Option 1: conda (easiest)
#'   conda install -c bioconda paintor
#'
#'   # Option 2: compile from source
#'   git clone https://github.com/gkichaev/PAINTOR_V3.0
#'   cd PAINTOR_V3.0
#'   make
#'   # Then either move the binary to your PATH or pass the full path to
#'   # run_methods() via method_args = list(paintor = list(paintor_path = ...))
#' }
#'
#' @param paintor_path Character. Path to the PAINTOR binary, or \code{"PAINTOR"}
#'   to search the system PATH. Default: \code{"PAINTOR"}.
#'
#' @return Invisibly returns the resolved path to the binary.
#' @export
setup_paintor <- function(paintor_path = "PAINTOR") {

  # Full path supplied: check directly
  if (file.exists(paintor_path)) {
    message("PAINTOR binary found: ", normalizePath(paintor_path))
    return(invisible(normalizePath(paintor_path)))
  }

  # Search PATH
  resolved <- Sys.which(paintor_path)
  if (nchar(resolved) > 0) {
    message("PAINTOR binary found on PATH: ", resolved)
    return(invisible(resolved))
  }

  stop(
    "PAINTOR binary not found (searched for: '", paintor_path, "').\n\n",
    "PAINTOR v3.0 must be installed manually. Options:\n\n",
    "  Option 1 - conda (easiest):\n",
    "    conda install -c bioconda paintor\n\n",
    "  Option 2 - compile from source:\n",
    "    git clone https://github.com/gkichaev/PAINTOR_V3.0\n",
    "    cd PAINTOR_V3.0 && make\n",
    "    # Then add to PATH or pass the full path via:\n",
    "    #   method_args = list(paintor = list(paintor_path = '/path/to/PAINTOR'))\n",
    call. = FALSE
  )
}


# =============================================================================
# Run PAINTOR on a single region
# =============================================================================

#' Run PAINTOR fine-mapping on a single region
#'
#' Writes the required PAINTOR input files to a temporary directory, calls the
#' PAINTOR binary, parses the posterior probability output, and returns results
#' in the standardised format. The temporary directory is deleted on exit.
#'
#' Annotations are optional. When \code{annotations} is supplied, PAINTOR
#' learns annotation enrichment weights via its EM algorithm. When omitted,
#' PAINTOR runs with an intercept only (equivalent to a uniform prior over
#' causal variants).
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param annotations Matrix or NULL. Functional annotations matrix (p x K).
#'   Numeric values; binary or continuous annotations are both supported.
#'   If NULL (default), PAINTOR runs without functional annotations.
#' @param variant_ids Character vector or NULL. Variant identifiers (length p).
#'   If NULL, set to \code{"SNP_1", "SNP_2", ...}.
#' @param paintor_path Character. Path to the PAINTOR binary or name on PATH.
#'   Default: \code{"PAINTOR"}.
#' @param max_causal Integer. Maximum number of causal variants to consider
#'   when using enumeration mode. Default: 2.
#' @param mcmc Logical. If TRUE, use MCMC-based search instead of exact
#'   enumeration. Recommended for large regions (p > 500) or when
#'   \code{max_causal > 3}. Default: FALSE.
#' @param coverage Numeric. Coverage level for constructing credible sets from
#'   PIPs. Default: 0.95.
#'
#' @return A list with the standardised fine-mapping output:
#' \describe{
#'   \item{pip}{Numeric vector (length p). Marginal posterior inclusion
#'     probabilities.}
#'   \item{credible_sets}{List containing one integer vector: the indices of
#'     variants in the credible set (ordered by decreasing PIP), constructed
#'     greedily from PIPs until cumulative PIP >= \code{coverage}.}
#'   \item{method}{Character. Always \code{"paintor"}.}
#'   \item{input_type}{Character. Always \code{"summary"}.}
#'   \item{params}{List. Hyperparameters used.}
#'   \item{runtime_seconds}{Numeric. Wall-clock time in seconds.}
#'   \item{additional}{List of PAINTOR-specific outputs:
#'     \describe{
#'       \item{annotations_used}{Integer. Number of annotation columns passed
#'         (0 if none).}
#'       \item{log_bayes_factor}{Numeric vector (length p) or NULL. Log Bayes
#'         factor per variant if present in PAINTOR output, otherwise NULL.}
#'     }
#'   }
#'   \item{error}{Character or NULL. Error message if PAINTOR failed.}
#' }
#'
#' @export
run_paintor <- function(z,
                        LD,
                        annotations  = NULL,
                        variant_ids  = NULL,
                        paintor_path = "PAINTOR",
                        max_causal   = 2,
                        mcmc         = FALSE,
                        coverage     = 0.95) {

  # --- Validate inputs --------------------------------------------------------

  p <- length(z)

  stopifnot(
    "LD must be a p x p matrix" =
      is.matrix(LD) && nrow(LD) == p && ncol(LD) == p,
    "max_causal must be a positive integer" =
      is.numeric(max_causal) && length(max_causal) == 1 &&
      max_causal == floor(max_causal) && max_causal >= 1,
    "coverage must be a single number in (0, 1)" =
      is.numeric(coverage) && length(coverage) == 1 &&
      coverage > 0 && coverage < 1
  )

  if (!is.null(annotations)) {
    stopifnot(
      "annotations must be a numeric matrix with p rows" =
        is.matrix(annotations) && is.numeric(annotations) &&
        nrow(annotations) == p
    )
  }

  if (is.null(variant_ids)) variant_ids <- paste0("SNP_", seq_len(p))

  n_annots <- if (!is.null(annotations)) ncol(annotations) else 0L

  params <- list(
    max_causal   = max_causal,
    mcmc         = mcmc,
    coverage     = coverage,
    paintor_path = paintor_path,
    n_annotations = n_annots
  )

  # --- Set up temp working directories ----------------------------------------

  work_dir <- tempfile(pattern = "paintor_run_")
  in_dir   <- file.path(work_dir, "input")
  out_dir  <- file.path(work_dir, "output")
  dir.create(in_dir,  recursive = TRUE)
  dir.create(out_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  locus_name <- "locus1"

  # --- Write locus file (z-scores only) ---------------------------------------
  # PAINTOR expects a single-column file with header "Zscore" and one value
  # per line. Annotations and rsids are NOT placed here.

  locus_path <- file.path(in_dir, locus_name)
  writeLines(c("Zscore", as.character(z)), locus_path)

  # --- Write LD file ----------------------------------------------------------
  # Named <locus_name>.LD (matches -LDname LD). Space-separated, no header.

  ld_path <- file.path(in_dir, paste0(locus_name, ".LD"))
  write.table(round(LD, 8), ld_path,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

  # --- Write annotations file -------------------------------------------------
  # PAINTOR always requires a <locus_name>.annotations file, even when no
  # functional annotations are used. When annotations are supplied, write them;
  # otherwise write a single intercept column (all 1s) as a placeholder that
  # satisfies the requirement without affecting the model (we do not pass it
  # via -annotations so PAINTOR ignores it).

  annot_names <- NULL
  annot_path  <- file.path(in_dir, paste0(locus_name, ".annotations"))

  if (n_annots > 0L) {
    annot_names <- paste0("ANNOT", seq_len(n_annots))
    annot_df <- as.data.frame(annotations)
    names(annot_df) <- annot_names
  } else {
    annot_df <- data.frame(Enrich = rep(1L, p))
  }

  write.table(annot_df, annot_path,
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = " ")

  # --- Write input file (list of locus names) ---------------------------------

  input_file <- file.path(work_dir, "input.txt")
  writeLines(locus_name, input_file)

  # --- Build PAINTOR arguments ------------------------------------------------

  args <- c(
    "-input",  input_file,
    "-in",     paste0(in_dir, "/"),
    "-out",    paste0(out_dir, "/"),
    "-Zhead",  "Zscore",
    "-LDname", "LD"
  )

  if (n_annots > 0L) {
    args <- c(args, "-annotations", paste(annot_names, collapse = ","))
  }

  if (mcmc) {
    args <- c(args, "-mcmc")
  } else {
    args <- c(args, "-enumerate", as.character(as.integer(max_causal)))
  }

  # --- Run PAINTOR ------------------------------------------------------------

  start_time <- proc.time()

  run_output <- tryCatch({
    system2(paintor_path, args = args, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    structure(conditionMessage(e), class = "paintor_error")
  })

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  if (inherits(run_output, "paintor_error")) {
    return(.paintor_error_result(p, params, elapsed, run_output))
  }

  # PAINTOR exits 0 even on some failures; check for the output file
  result_path <- file.path(out_dir, paste0(locus_name, ".results"))

  if (!file.exists(result_path)) {
    err_msg <- paste(
      c("PAINTOR produced no output file.", run_output),
      collapse = "\n"
    )
    return(.paintor_error_result(p, params, elapsed, err_msg))
  }

  # --- Parse results ----------------------------------------------------------

  res_df <- tryCatch(
    utils::read.table(result_path, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  if (is.null(res_df) || !"Posterior_Prob" %in% names(res_df)) {
    return(.paintor_error_result(
      p, params, elapsed,
      "Failed to parse PAINTOR results file or 'Posterior_Prob' column missing."
    ))
  }

  # Re-order to match the original variant order (PAINTOR preserves order, but
  # be defensive in case it doesn't).
  if ("rsid" %in% names(res_df)) {
    ord <- match(variant_ids, res_df$rsid)
    if (anyNA(ord)) ord <- seq_len(p)   # fallback: assume order preserved
  } else {
    ord <- seq_len(p)
  }

  pip <- as.numeric(res_df$Posterior_Prob[ord])

  # Clamp to [0, 1] - floating-point rounding can push values just outside
  pip <- pmax(0, pmin(1, pip))

  # Log Bayes factor column (present in some PAINTOR versions)
  log_bayes_factor <- if ("log_BF" %in% names(res_df)) {
    as.numeric(res_df$log_BF[ord])
  } else {
    NULL
  }

  # --- Build credible set from PIPs -------------------------------------------
  # Greedy: sort by descending PIP, accumulate until cumsum >= coverage.

  ord_pip    <- order(pip, decreasing = TRUE)
  cumulative <- cumsum(pip[ord_pip])
  n_in_cs    <- which(cumulative >= coverage)[1L]
  if (is.na(n_in_cs)) n_in_cs <- p    # safety: include all if PIPs don't sum to 1
  cs_indices <- sort(ord_pip[seq_len(n_in_cs)])

  credible_sets <- list(cs_indices)

  # --- Return -----------------------------------------------------------------

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = "paintor",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      annotations_used  = n_annots,
      log_bayes_factor  = log_bayes_factor
    )
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run PAINTOR on a single region from simulation data structures
#'
#' Thin adapter that extracts the appropriate inputs from the simulation's
#' \code{region_geno} and \code{region_pheno} objects and calls
#' \code{\link{run_paintor}}. This is the function registered in the method
#' registry and called by \code{\link{run_methods}}.
#'
#' Annotations are passed if present in \code{region_pheno}
#' (\code{annotations_matrix}). If absent, PAINTOR runs without annotations
#' (uniform prior). Unlike Funmap, annotations are not required.
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing \code{LD}, \code{n}, and optionally \code{variant_ids}.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z} and optionally \code{annotations_matrix}.
#' @param ... Additional arguments passed to \code{\link{run_paintor}}
#'   (e.g. \code{max_causal}, \code{mcmc}, \code{paintor_path}).
#'
#' @return The output of \code{\link{run_paintor}}.
#' @export
run_paintor_region <- function(region_geno, region_pheno,
                               .paintor_cache = NULL, ...) {
  # Prefer the cross-loci result computed once per scenario by
  # run_paintor_scenario_setup() - that is how PAINTOR is designed to be run
  # (its EM learns annotation enrichment ACROSS loci). Fall back to a
  # single-locus run if the pooled call was unavailable or failed.
  if (!is.null(.paintor_cache)) {
    hit <- .paintor_cache[[.paintor_fingerprint(region_pheno$z)]]
    if (!is.null(hit)) return(hit)
  }
  out <- run_paintor(
    z           = region_pheno$z,
    LD          = region_geno$LD,
    annotations = .paintor_extract_annotations(region_geno, region_pheno),
    variant_ids = region_geno$variant_ids,
    ...
  )
  out$additional$pooled_across_loci <- FALSE
  out
}


# =============================================================================
# Scenario-level setup: pool annotation enrichment across ALL loci
# =============================================================================

#' Run PAINTOR once across every locus in a scenario
#'
#' PAINTOR's EM algorithm estimates annotation enrichment weights from ALL loci
#' listed in its \code{-input} file; this is how it is run on real data and is
#' the point of the method. Running it one locus at a time - as this wrapper
#' originally did - estimates the enrichment from a single locus, which both
#' wastes the cross-loci information PAINTOR is built to exploit and makes the
#' enrichment estimate extremely noisy.
#'
#' This hook is called once per scenario by \code{\link{run_methods}}. It writes
#' every region as its own locus, lists them all in one \code{input.txt}, makes a
#' single PAINTOR call, and returns the per-region results keyed by a
#' fingerprint of the region's z-vector. \code{run_paintor_region()} then looks
#' up its own region.
#'
#' @param genotypes List. \code{simulation$genotypes}.
#' @param regions List. One scenario's \code{regions}.
#' @param user_args List. The user's method_args for paintor.
#'
#' @return A list with \code{.paintor_cache}, or an empty list if the pooled run
#'   could not be performed (callers then fall back to per-locus runs).
#' @export
run_paintor_scenario_setup <- function(genotypes, regions, user_args) {

  n_regions <- length(regions)
  if (n_regions < 2L) return(list())    # nothing to pool

  gv <- function(k, d) { v <- user_args[[k]]; if (is.null(v)) d else v }
  paintor_path <- gv("paintor_path", "PAINTOR")
  max_causal   <- gv("max_causal", 2)
  mcmc         <- gv("mcmc", FALSE)
  coverage     <- gv("coverage", 0.95)

  # ---- collect the loci ------------------------------------------------------
  specs <- vector("list", n_regions)
  keys  <- character(n_regions)
  n_ann <- NULL
  for (i in seq_len(n_regions)) {
    z  <- regions[[i]]$z
    LD <- genotypes[[i]]$LD
    A  <- .paintor_extract_annotations(genotypes[[i]], regions[[i]])
    if (is.null(z) || is.null(LD) || nrow(LD) != length(z)) return(list())
    # every locus must expose the same annotation columns for a joint fit
    k <- if (is.null(A)) 0L else ncol(A)
    if (is.null(n_ann)) n_ann <- k else if (k != n_ann) return(list())
    ids <- genotypes[[i]]$variant_ids
    if (is.null(ids)) ids <- paste0("SNP_", seq_along(z))
    specs[[i]] <- list(z = z, LD = LD, A = A, ids = ids, p = length(z))
    keys[i] <- .paintor_fingerprint(z)
  }

  work_dir <- tempfile(pattern = "paintor_pool_")
  in_dir   <- file.path(work_dir, "input")
  out_dir  <- file.path(work_dir, "output")
  dir.create(in_dir, recursive = TRUE); dir.create(out_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  annot_names <- if (n_ann > 0L) paste0("ANNOT", seq_len(n_ann)) else NULL
  locus_names <- sprintf("locus%d", seq_len(n_regions))

  for (i in seq_len(n_regions)) {
    s <- specs[[i]]; ln <- locus_names[i]
    writeLines(c("Zscore", as.character(s$z)), file.path(in_dir, ln))
    utils::write.table(round(s$LD, 8), file.path(in_dir, paste0(ln, ".LD")),
                       quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
    adf <- if (n_ann > 0L) {
      x <- as.data.frame(s$A); names(x) <- annot_names; x
    } else data.frame(Enrich = rep(1L, s$p))
    utils::write.table(adf, file.path(in_dir, paste0(ln, ".annotations")),
                       quote = FALSE, row.names = FALSE, col.names = TRUE, sep = " ")
  }
  # THE fix: every locus is listed, so PAINTOR's EM pools across all of them.
  input_file <- file.path(work_dir, "input.txt")
  writeLines(locus_names, input_file)

  args <- c("-input", input_file, "-in", paste0(in_dir, "/"),
            "-out", paste0(out_dir, "/"), "-Zhead", "Zscore", "-LDname", "LD")
  if (n_ann > 0L) args <- c(args, "-annotations", paste(annot_names, collapse = ","))
  args <- if (mcmc) c(args, "-mcmc") else c(args, "-enumerate", as.character(as.integer(max_causal)))

  run_output <- tryCatch(system2(paintor_path, args = args, stdout = TRUE, stderr = TRUE),
                         error = function(e) structure(conditionMessage(e), class = "paintor_error"))
  if (inherits(run_output, "paintor_error")) return(list())

  # ---- parse each locus ------------------------------------------------------
  cache <- list()
  for (i in seq_len(n_regions)) {
    rp <- file.path(out_dir, paste0(locus_names[i], ".results"))
    if (!file.exists(rp)) next
    df <- tryCatch(utils::read.table(rp, header = TRUE, stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (is.null(df) || !"Posterior_Prob" %in% names(df)) next
    s <- specs[[i]]
    ord <- if ("rsid" %in% names(df)) {
      o <- match(s$ids, df$rsid); if (anyNA(o)) seq_len(s$p) else o
    } else seq_len(s$p)
    pip <- pmax(0, pmin(1, as.numeric(df$Posterior_Prob[ord])))
    if (length(pip) != s$p) next
    op <- order(pip, decreasing = TRUE)
    n_cs <- which(cumsum(pip[op]) >= coverage)[1L]
    if (is.na(n_cs)) n_cs <- s$p
    cache[[keys[i]]] <- list(
      pip = pip, credible_sets = list(sort(op[seq_len(n_cs)])),
      method = "paintor", input_type = "summary",
      params = list(max_causal = max_causal, mcmc = mcmc, coverage = coverage,
                    paintor_path = paintor_path, n_annotations = n_ann),
      runtime_seconds = NA_real_,
      additional = list(annotations_used = n_ann,
                        log_bayes_factor = if ("log_BF" %in% names(df))
                          as.numeric(df$log_BF[ord]) else NULL,
                        pooled_across_loci = TRUE,
                        n_loci_pooled = n_regions))
  }
  if (length(cache) == 0L) return(list())
  list(.paintor_cache = cache)
}


# =============================================================================
# Internal helpers
# =============================================================================

# Prefer region_geno$annotations_matrix, falling back to region_pheno's copy.
# Both simulators populate the geno-side matrix; the locus pipeline additionally
# copies it onto the pheno object, but simulate_gwfm_data() does NOT - so reading
# only from the pheno side silently drops annotations under the genome-wide
# simulator, which is exactly the bug fixed for Functional BEATRICE in PR #19.
# PAINTOR had the same latent bug.
.paintor_extract_annotations <- function(region_geno, region_pheno) {
  A <- region_geno$annotations_matrix
  if (is.null(A)) A <- region_pheno$annotations_matrix
  A
}

# Deterministic per-region key from the z vector, used to hand a pooled result
# back to the right region. Same construction as R/wrapper_fb_joint.R.
.paintor_fingerprint <- function(z) {
  z <- as.numeric(z)
  paste(length(z), format(sum(z), digits = 15), format(sum(z * z), digits = 15),
        format(z[1], digits = 15), format(z[length(z)], digits = 15), sep = "|")
}

.paintor_error_result <- function(p, params, elapsed, error_msg) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "paintor",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      annotations_used = params$n_annotations %||% 0L,
      log_bayes_factor = NULL
    ),
    error = error_msg
  )
}

# %||% is defined in R/utils.R.
