# =============================================================================
# functional_beatrice.R
#
# Wrapper for Functional BEATRICE (BEATRICE with annotation-informed prior).
#
# Functional BEATRICE extends BEATRICE with a LassoNet prior network that
# maps per-variant annotation vectors to annotation-informed inclusion
# probabilities p_0(j). When no annotations are supplied it degrades
# gracefully to a uniform prior, behaving like standard BEATRICE.
#
# Statistical model
# -----------------
# Inference network  : Binary Concrete (Gumbel-Softmax) variational
#                      distribution over causal configurations; same as
#                      BEATRICE.
# Prior network      : LassoNetPrior — a shallow MLP with skip connections
#                      θ from each annotation directly to the output logits,
#                      plus hidden layers. An L1 proximal penalty on θ drives
#                      annotation sparsity; a hierarchy constraint ensures
#                      that zeroing θ_j also gates out the hidden path for
#                      annotation j.
# Training objective : ELBO = E[log p(z|γ)] − KL(q || p_0)
#                           + λ_reg·||p_0||² + L1(θ)
#
# Input files  (written to a temp directory)
#   .z     space-separated: variant_name z_score, no header
#   .ld    space-separated p×p matrix, no header
#   .annot space-separated: variant_name ann1 ann2 …, no header (optional)
#
# Output files  (read from temp directory)
#   pip.csv                 variant_index, pip, variant_names
#   credible_set.txt        one credible set per line, 0-based indices
#   feature_importance.csv  annotation importance scores (LassoNet |θ|)
#
# This file provides:
#   setup_functional_beatrice()          : verify beatrice_annot.py and deps
#   run_functional_beatrice()            : run on a single region
#   run_functional_beatrice_region()     : adapter called by run_methods()
#
# Standard output format mirrors run_beatrice():
#   pip              Numeric vector (length p). Marginal PIPs.
#   credible_sets    List of integer vectors (1-based).
#   method           "functional_beatrice"
#   input_type       "summary"
#   params           List. Hyperparameters used.
#   runtime_seconds  Numeric.
#   additional       List (see run_functional_beatrice() docs).
# =============================================================================


# =============================================================================
# Setup
# =============================================================================

#' Set up Functional BEATRICE
#'
#' Verifies that \code{beatrice_annot.py} can be found in \code{beatrice_dir}
#' and that the required Python packages are importable.
#'
#' The Functional BEATRICE code lives in the \code{BEATRICE_annot_sparse/}
#' subdirectory of this repository.  Point \code{beatrice_dir} there and
#' supply the Python executable from your BEATRICE conda environment:
#'
#' \preformatted{
#'   setup_functional_beatrice(
#'     beatrice_dir = "BEATRICE_annot_sparse",
#'     python       = "~/anaconda3/envs/beatrice/bin/python"
#'   )
#' }
#'
#' @param beatrice_dir Character. Path to the directory containing
#'   \code{beatrice_annot.py} (i.e. \code{BEATRICE_annot_sparse/}).
#' @param python Character. Path to the Python executable. Default:
#'   \code{"python"} (searches PATH).
#'
#' @return Invisibly returns a named list with \code{beatrice_script} and
#'   \code{python}.
#' @export
setup_functional_beatrice <- function(beatrice_dir, python = "python") {

  beatrice_dir <- normalizePath(path.expand(beatrice_dir), mustWork = FALSE)

  script_path <- file.path(beatrice_dir, "beatrice_annot.py")
  if (!file.exists(script_path)) {
    stop(
      "beatrice_annot.py not found in: ", beatrice_dir, "\n\n",
      "The Functional BEATRICE code should be in BEATRICE_annot_sparse/.\n",
      "Pass that directory:\n",
      "  setup_functional_beatrice(beatrice_dir = 'BEATRICE_annot_sparse')",
      call. = FALSE
    )
  }

  resolved_python <- if (file.exists(path.expand(python))) {
    normalizePath(path.expand(python))
  } else {
    py <- Sys.which(python)
    if (nchar(py) == 0) {
      stop(
        "Python executable not found: '", python, "'\n\n",
        "Pass the full path:\n",
        "  setup_functional_beatrice(\n",
        "    beatrice_dir = 'BEATRICE_annot_sparse',\n",
        "    python       = '~/anaconda3/envs/beatrice/bin/python'\n",
        "  )",
        call. = FALSE
      )
    }
    py
  }

  required_pkgs <- c("torch", "numpy", "pandas", "tqdm", "absl", "imageio",
                     "seaborn", "matplotlib")
  check_pkg <- function(pkg) {
    system2(resolved_python,
            args = c("-c", shQuote(paste0("import ", pkg))),
            stdout = FALSE, stderr = FALSE) == 0L
  }
  missing_pkgs <- Filter(Negate(check_pkg), required_pkgs)

  if (length(missing_pkgs) > 0) {
    stop(
      "The following Python packages are missing from '", resolved_python, "':\n",
      paste0("  ", missing_pkgs, collapse = "\n"), "\n\n",
      "Install the BEATRICE conda environment:\n",
      "  conda env create -f <path>/conda_environment.yml\n",
      "  conda activate beatrice\n",
      call. = FALSE
    )
  }

  message("Functional BEATRICE ready.")
  message("  Script : ", script_path)
  message("  Python : ", resolved_python)

  invisible(list(beatrice_script = script_path, python = resolved_python))
}


# =============================================================================
# Run Functional BEATRICE on a single region
# =============================================================================

#' Run Functional BEATRICE fine-mapping on a single region
#'
#' Writes the required input files to a temporary directory, calls
#' \code{beatrice_annot.py} via the specified Python executable, parses the
#' output, and returns results in the standardised format. The temporary
#' directory is deleted on exit.
#'
#' When \code{annotations} is \code{NULL} the method runs without a functional
#' prior (equivalent to standard BEATRICE with a uniform prior). When
#' annotations are supplied a LassoNet prior network is trained jointly with
#' the inference network; the learned annotation importance scores are returned
#' in \code{additional$feature_importance}.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p × p).
#' @param n Integer. Sample size (GWAS N).
#' @param annotations Numeric matrix or NULL. Annotation matrix (p × m). Each
#'   row is one variant; each column one annotation (binary or continuous).
#'   Pass \code{NULL} to run without annotations (uniform prior).
#' @param variant_ids Character vector or NULL. Variant identifiers (length p).
#'   If NULL, set to \code{"rs0", "rs1", …}.
#' @param beatrice_dir Character. Path to the directory containing
#'   \code{beatrice_annot.py}.
#' @param python Character. Path to the Python executable. Default:
#'   \code{"python"}.
#' @param max_iter Integer. Training iterations. Default: 2000.
#' @param n_caus Integer. Expected maximum number of causal variants.
#'   Default: 5.
#' @param sigma_sq Numeric. Prior variance on effect sizes. Default: 0.05.
#' @param gamma_coverage Numeric. Coverage threshold for credible sets.
#'   Default: 0.95.
#' @param sparse_concrete Integer. Sparsity of Binary Concrete samples per
#'   iteration. Default: 50.
#' @param prior_regularisation Numeric. L2 regularisation weight on the
#'   prior network outputs. Default: 1.0.
#' @param lambda_l1 Numeric. L1 penalty on LassoNet skip connections.
#'   Controls annotation sparsity. Default: 0.01.
#' @param hierarchy_M Numeric. LassoNet hierarchy constraint multiplier
#'   (\code{||W_j^(1)||} ≤ M·|θ_j|). Default: 10.0.
#'
#' @return A list with the standardised fine-mapping output:
#' \describe{
#'   \item{pip}{Numeric vector (length p). Posterior inclusion probabilities.}
#'   \item{credible_sets}{List of integer vectors (1-based).}
#'   \item{method}{Character. Always \code{"functional_beatrice"}.}
#'   \item{input_type}{Character. Always \code{"summary"}.}
#'   \item{params}{List. Hyperparameters used.}
#'   \item{runtime_seconds}{Numeric.}
#'   \item{additional}{List:
#'     \describe{
#'       \item{cs_pip}{List of numeric vectors. Conditional inclusion
#'         probabilities within each credible set. NULL if not produced.}
#'       \item{feature_importance}{Data frame (annotation_index, importance)
#'         from LassoNet |θ|, sorted descending. NULL if no annotations or
#'         file not produced.}
#'     }
#'   }
#'   \item{error}{Character or NULL. Error message if the method failed.}
#' }
#'
#' @export
run_functional_beatrice <- function(z,
                                    LD,
                                    n,
                                    annotations          = NULL,
                                    variant_ids          = NULL,
                                    beatrice_dir,
                                    python               = "python",
                                    max_iter             = 2000,
                                    n_caus               = 5,
                                    sigma_sq             = 0.05,
                                    gamma_coverage       = 0.95,
                                    sparse_concrete      = 50,
                                    prior_regularisation = 1.0,
                                    lambda_l1            = 0.01,
                                    hierarchy_M          = 10.0) {

  # --- Validate inputs --------------------------------------------------------

  p <- length(z)

  stopifnot(
    "LD must be a p x p matrix" =
      is.matrix(LD) && nrow(LD) == p && ncol(LD) == p,
    "n must be a positive integer" =
      is.numeric(n) && length(n) == 1 && n > 0,
    "max_iter must be >= 500" =
      is.numeric(max_iter) && max_iter >= 500
  )

  if (!is.null(annotations)) {
    if (!is.matrix(annotations) || nrow(annotations) != p) {
      stop("annotations must be a numeric matrix with nrow == length(z).",
           call. = FALSE)
    }
  }

  if (is.null(variant_ids)) variant_ids <- paste0("rs", seq(0, p - 1))

  # CRITICAL: both the .z and .annot files are whitespace-separated, and
  # trainer_annot.py parses them positionally:
  #     Z = pd.read_table(z_file, sep=' ', ...).to_numpy()[:,1]
  #     v = pd.read_table(annot_file, sep=None, ...).to_numpy()[:,1:].astype(float)
  # VCF-derived variant_ids look like "1 40023356 . A T" - embedded spaces.
  # Written verbatim they shift the columns, so z becomes the POSITION and
  # the annotation cast hits the "." rsID ("could not convert string to
  # float: '.'"), failing every annotated run. Collapse whitespace so each
  # variant_id is a single token.
  variant_ids <- gsub("\\s+", "_", variant_ids)

  beatrice_dir <- normalizePath(path.expand(beatrice_dir), mustWork = FALSE)
  python       <- if (file.exists(path.expand(python))) {
    normalizePath(path.expand(python))
  } else {
    python
  }
  script_path <- file.path(beatrice_dir, "beatrice_annot.py")

  if (!file.exists(script_path)) {
    return(.fb_error_result(
      p,
      .fb_params(n, max_iter, n_caus, sigma_sq, gamma_coverage,
                 sparse_concrete, prior_regularisation, lambda_l1,
                 hierarchy_M, beatrice_dir, python, !is.null(annotations)),
      0, paste("beatrice_annot.py not found at:", script_path)
    ))
  }

  params <- .fb_params(n, max_iter, n_caus, sigma_sq, gamma_coverage,
                       sparse_concrete, prior_regularisation, lambda_l1,
                       hierarchy_M, beatrice_dir, python, !is.null(annotations))

  # --- Set up temp working directory ------------------------------------------

  work_dir <- tempfile(pattern = "fb_run_")
  out_dir  <- file.path(work_dir, "output")
  dir.create(work_dir, recursive = TRUE)
  dir.create(out_dir,  recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  # --- Write .z file ----------------------------------------------------------

  z_path  <- file.path(work_dir, "region.z")
  writeLines(paste(variant_ids, z, sep = " "), z_path)

  # --- Write .ld file ---------------------------------------------------------

  ld_path <- file.path(work_dir, "region.ld")
  write.table(round(LD, 8), ld_path,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

  # --- Write .annot file (if annotations provided) ----------------------------

  annot_path <- NULL
  if (!is.null(annotations)) {
    annot_path <- file.path(work_dir, "region.annot")
    annot_df   <- cbind(variant_ids, as.data.frame(annotations))
    write.table(annot_df, annot_path,
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
  }

  # --- Build command-line arguments -------------------------------------------

  args <- c(
    script_path,
    "--z",                    z_path,
    "--LD",                   ld_path,
    "--N",                    as.character(as.integer(n)),
    "--target",               paste0(out_dir, "/"),
    "--max_iter",             as.character(as.integer(max_iter)),
    "--n_caus",               as.character(as.integer(n_caus)),
    "--sigma_sq",             as.character(sigma_sq),
    "--gamma_coverage",       as.character(gamma_coverage),
    "--sparse_concrete",      as.character(as.integer(sparse_concrete)),
    "--prior_regularisation", as.character(prior_regularisation),
    "--lambda_l1",            as.character(lambda_l1),
    "--hierarchy_M",          as.character(hierarchy_M),
    "--plot_loss",            "False",
    "--get_cred",             "True"
  )

  if (!is.null(annot_path)) {
    args <- c(args, "--annot", annot_path)
  }

  # --- Run Functional BEATRICE ------------------------------------------------

  start_time <- proc.time()

  run_output <- tryCatch({
    system2(python, args = args, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    structure(conditionMessage(e), class = "fb_error")
  })

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  if (inherits(run_output, "fb_error")) {
    return(.fb_error_result(p, params, elapsed, run_output))
  }

  # --- Check output -----------------------------------------------------------

  pip_path  <- file.path(out_dir, "pip.csv")
  cred_path <- file.path(out_dir, "credible_set.txt")

  if (!file.exists(pip_path)) {
    err_msg <- paste(
      c("Functional BEATRICE produced no pip.csv.", run_output),
      collapse = "\n"
    )
    return(.fb_error_result(p, params, elapsed, err_msg))
  }

  # --- Parse pip.csv ----------------------------------------------------------

  pip_df <- tryCatch(
    utils::read.csv(pip_path, stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  if (is.null(pip_df) || !"pip" %in% names(pip_df)) {
    return(.fb_error_result(
      p, params, elapsed,
      "Failed to parse pip.csv or 'pip' column missing."
    ))
  }

  if ("variant_names" %in% names(pip_df)) {
    ord <- match(variant_ids, pip_df$variant_names)
    if (anyNA(ord)) ord <- seq_len(p)
  } else {
    ord <- seq_len(p)
  }

  pip <- pmax(0, pmin(1, as.numeric(pip_df$pip[ord])))

  # --- Parse credible_set.txt (0-based → 1-based) ----------------------------

  credible_sets <- list()
  cs_pip        <- NULL

  if (file.exists(cred_path)) {
    cred_lines <- readLines(cred_path, warn = FALSE)
    cred_lines <- cred_lines[nchar(trimws(cred_lines)) > 0]
    if (length(cred_lines) > 0) {
      credible_sets <- lapply(cred_lines, function(line) {
        idx <- suppressWarnings(as.integer(strsplit(trimws(line), "\\s+")[[1]]))
        sort(idx[!is.na(idx)] + 1L)
      })
      credible_sets <- credible_sets[lengths(credible_sets) > 0]
    }
  }

  cond_path <- file.path(out_dir, "conditional_credible_variants_probability.txt")
  if (file.exists(cond_path)) {
    cond_lines <- readLines(cond_path, warn = FALSE)
    cond_lines <- cond_lines[nchar(trimws(cond_lines)) > 0]
    cs_pip <- lapply(cond_lines, function(line) {
      suppressWarnings(as.numeric(strsplit(trimws(line), "\\s+")[[1]]))
    })
  }

  # --- Parse feature_importance.csv ------------------------------------------

  feature_importance <- NULL
  fi_path <- file.path(out_dir, "feature_importance.csv")
  if (file.exists(fi_path)) {
    feature_importance <- tryCatch(
      utils::read.csv(fi_path, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  }

  # --- Return -----------------------------------------------------------------

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = "functional_beatrice",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      cs_pip             = cs_pip,
      feature_importance = feature_importance
    )
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run Functional BEATRICE on a single region from simulation data structures
#'
#' Thin adapter that extracts z-scores, LD matrix, sample size, variant IDs,
#' and (optionally) annotation matrix from the simulation's data structures
#' and calls \code{\link{run_functional_beatrice}}.
#'
#' The annotation matrix is read from \code{region_geno$annotations_matrix}
#' first, falling back to \code{region_pheno$annotations_matrix} for
#' backward compatibility with the locus pipeline (which copies annotations
#' onto both objects). If both are NULL the method runs without a
#' functional prior. The geno-first ordering matters for the genome-wide
#' pipeline (\code{simulate_gwfm_data}), which populates only the geno-side
#' matrix.
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing at minimum \code{LD}, \code{n}, \code{variant_ids}, and
#'   optionally \code{annotations_matrix} (p x m).
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z} and, in the locus pipeline, a copy of
#'   \code{annotations_matrix}.
#' @param ... Additional arguments passed to
#'   \code{\link{run_functional_beatrice}} (e.g. \code{beatrice_dir},
#'   \code{python}, \code{max_iter}).
#'
#' @return The output of \code{\link{run_functional_beatrice}}.
#' @export
run_functional_beatrice_region <- function(region_geno, region_pheno, ...) {
  run_functional_beatrice(
    z            = region_pheno$z,
    LD           = region_geno$LD,
    n            = region_geno$n,
    annotations  = .fb_extract_annotations(region_geno, region_pheno),
    variant_ids  = region_geno$variant_ids,
    ...
  )
}


# Prefer region_geno$annotations_matrix, falling back to
# region_pheno$annotations_matrix. Both simulators (run_simulation,
# simulate_gwfm_data) populate the geno-side matrix; the locus pipeline
# additionally copies it onto the pheno object at scenario-build time
# (see run_simulation.R line ~420), but the genome-wide simulator does
# NOT — so reading only from the pheno side silently drops annotations
# under simulate_gwfm_data. This helper is factored out so the
# annotation-selection contract can be regression-tested without
# needing the Python-only run_functional_beatrice() to actually run.
.fb_extract_annotations <- function(region_geno, region_pheno) {
  A <- region_geno$annotations_matrix
  if (is.null(A)) A <- region_pheno$annotations_matrix
  A
}


# =============================================================================
# Internal helpers
# =============================================================================

.fb_params <- function(n, max_iter, n_caus, sigma_sq, gamma_coverage,
                       sparse_concrete, prior_regularisation, lambda_l1,
                       hierarchy_M, beatrice_dir, python, annotations_given) {
  list(
    n                    = n,
    max_iter             = max_iter,
    n_caus               = n_caus,
    sigma_sq             = sigma_sq,
    gamma_coverage       = gamma_coverage,
    sparse_concrete      = sparse_concrete,
    prior_regularisation = prior_regularisation,
    lambda_l1            = lambda_l1,
    hierarchy_M          = hierarchy_M,
    beatrice_dir         = beatrice_dir,
    python               = python,
    annotations_given    = annotations_given
  )
}

.fb_error_result <- function(p, params, elapsed, error_msg) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "functional_beatrice",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(cs_pip = NULL, feature_importance = NULL),
    error           = error_msg
  )
}
