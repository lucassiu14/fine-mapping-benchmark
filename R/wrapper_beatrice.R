# =============================================================================
# beatrice.R
#
# Wrapper for BEATRICE (Ghosal et al.) fine-mapping.
#
# BEATRICE is a Python script that performs variational Bayesian fine-mapping
# using a neural network to learn posterior inclusion probabilities. It uses
# binary concrete distributions (Gumbel-Softmax) for differentiable
# approximation of discrete causal configurations and ABF-based likelihoods.
#
# Unlike the other external-binary methods (FINEMAP, PAINTOR), BEATRICE is
# called as a Python script rather than a compiled executable:
#   python <beatrice_dir>/beatrice.py --z <z_file> --LD <ld_file> --N <n> ...
#
# This file provides:
#   - setup_beatrice()        : verifies beatrice.py and Python deps are found
#   - run_beatrice()          : runs BEATRICE on a single region
#   - run_beatrice_region()   : adapter called by run_methods()
#
# Standard output format:
#   pip              Numeric vector (length p). Marginal PIPs from pip.csv.
#   credible_sets    List of integer vectors (1-based). Parsed from
#                    credible_set.txt (BEATRICE uses 0-based indices;
#                    converted here).
#   method           Character. "beatrice".
#   input_type       Character. Always "summary".
#   params           List. Hyperparameters used.
#   runtime_seconds  Numeric. Wall-clock time.
#   additional       List. BEATRICE-specific outputs (see run_beatrice() docs).
#
# Reference:
#   Ghosal S et al. BEATRICE: Bayesian Fine-mapping from Summary Data using
#   Deep Variational Inference. https://github.com/sayangsep/Beatrice-Finemapping
# =============================================================================


# =============================================================================
# Setup
# =============================================================================

#' Set up BEATRICE
#'
#' Verifies that the BEATRICE Python script can be found at \code{beatrice_dir}
#' and that the required Python packages (torch, numpy, scipy, pandas) are
#' importable by the specified Python executable.
#'
#' Install BEATRICE by cloning the repository and creating its conda
#' environment:
#'
#' \preformatted{
#'   git clone https://github.com/sayangsep/Beatrice-Finemapping
#'   cd Beatrice-Finemapping
#'   conda env create -f conda_environment.yml
#'   conda activate beatrice
#' }
#'
#' Then pass the repo directory and the conda-env Python to this function:
#'
#' \preformatted{
#'   setup_beatrice(
#'     beatrice_dir = "~/Beatrice-Finemapping",
#'     python       = "~/anaconda3/envs/beatrice/bin/python"
#'   )
#' }
#'
#' @param beatrice_dir Character. Path to the cloned Beatrice-Finemapping
#'   repository root (must contain \code{beatrice.py}).
#' @param python Character. Path to the Python executable to use.
#'   Default: \code{"python"} (searches PATH).
#'
#' @return Invisibly returns a named list with \code{beatrice_script} (full
#'   path to \code{beatrice.py}) and \code{python} (resolved Python path).
#' @export
setup_beatrice <- function(beatrice_dir, python = "python") {

  beatrice_dir <- path.expand(beatrice_dir)

  # Check beatrice.py exists
  script_path <- file.path(beatrice_dir, "beatrice.py")
  if (!file.exists(script_path)) {
    stop(
      "beatrice.py not found in: ", beatrice_dir, "\n\n",
      "Please clone the repository:\n",
      "  git clone https://github.com/sayangsep/Beatrice-Finemapping\n\n",
      "Then pass the path:\n",
      "  setup_beatrice(beatrice_dir = '/path/to/Beatrice-Finemapping')",
      call. = FALSE
    )
  }

  # Check Python executable
  resolved_python <- if (file.exists(path.expand(python))) {
    normalizePath(path.expand(python))
  } else {
    py <- Sys.which(python)
    if (nchar(py) == 0) {
      stop(
        "Python executable not found: '", python, "'\n\n",
        "Pass the full path to your Python:\n",
        "  setup_beatrice(\n",
        "    beatrice_dir = '/path/to/Beatrice-Finemapping',\n",
        "    python       = '~/anaconda3/envs/beatrice/bin/python'\n",
        "  )",
        call. = FALSE
      )
    }
    py
  }

  # Check required Python packages
  check_pkg <- function(pkg) {
    out <- system2(resolved_python,
                   args = c("-c", shQuote(paste0("import ", pkg))),
                   stdout = FALSE, stderr = FALSE)
    out == 0L
  }

  missing_pkgs <- Filter(Negate(check_pkg), c("torch", "numpy", "scipy", "pandas"))

  if (length(missing_pkgs) > 0) {
    stop(
      "The following Python packages are missing from '", resolved_python, "':\n",
      paste0("  ", missing_pkgs, collapse = "\n"), "\n\n",
      "Install the BEATRICE conda environment:\n",
      "  cd /path/to/Beatrice-Finemapping\n",
      "  conda env create -f conda_environment.yml\n",
      "  conda activate beatrice\n\n",
      "Then pass the conda-env Python:\n",
      "  setup_beatrice(\n",
      "    beatrice_dir = '/path/to/Beatrice-Finemapping',\n",
      "    python       = '~/anaconda3/envs/beatrice/bin/python'\n",
      "  )",
      call. = FALSE
    )
  }

  message("BEATRICE ready.")
  message("  Script : ", script_path)
  message("  Python : ", resolved_python)

  invisible(list(beatrice_script = script_path, python = resolved_python))
}


# =============================================================================
# Run BEATRICE on a single region
# =============================================================================

#' Run BEATRICE fine-mapping on a single region
#'
#' Writes the required BEATRICE input files to a temporary directory, calls
#' \code{beatrice.py} via the specified Python executable, parses the output,
#' and returns results in the standardised format. The temporary directory is
#' deleted on exit.
#'
#' BEATRICE trains a neural network for each region; runtime scales with
#' \code{max_iter}. For large-scale benchmarking, reduce \code{max_iter}
#' (minimum 500) or run on a machine with GPU support.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param n Integer. Sample size (GWAS N).
#' @param variant_ids Character vector or NULL. Variant identifiers (length p).
#'   If NULL, set to \code{"rs0", "rs1", ...} (matching BEATRICE convention).
#' @param beatrice_dir Character. Path to the cloned Beatrice-Finemapping
#'   repository root.
#' @param python Character. Path to the Python executable. Default:
#'   \code{"python"}.
#' @param max_iter Integer. Training iterations. Default: 2000. Minimum: 500.
#'   Reduce for faster (but noisier) results.
#' @param n_caus Integer. Expected number of causal variants. Default: 5.
#' @param sigma_sq Numeric. Prior variance on effect sizes. Default: 0.05.
#' @param gamma_coverage Numeric. Coverage threshold for credible sets.
#'   Default: 0.95.
#' @param sparse_concrete Integer. Number of non-zero locations sampled per
#'   iteration (sparsity). Default: 50.
#'
#' @return A list with the standardised fine-mapping output:
#' \describe{
#'   \item{pip}{Numeric vector (length p). Posterior inclusion probabilities.}
#'   \item{credible_sets}{List of integer vectors (1-based). One credible set
#'     per line in \code{credible_set.txt}. Empty list if no credible sets
#'     were produced.}
#'   \item{method}{Character. Always \code{"beatrice"}.}
#'   \item{input_type}{Character. Always \code{"summary"}.}
#'   \item{params}{List. Hyperparameters used.}
#'   \item{runtime_seconds}{Numeric. Wall-clock time in seconds.}
#'   \item{additional}{List of BEATRICE-specific outputs:
#'     \describe{
#'       \item{cs_pip}{List of numeric vectors. Conditional inclusion
#'         probabilities within each credible set, as reported by BEATRICE.
#'         NULL if the file is not produced.}
#'     }
#'   }
#'   \item{error}{Character or NULL. Error message if BEATRICE failed.}
#' }
#'
#' @export
run_beatrice <- function(z,
                         LD,
                         n,
                         variant_ids    = NULL,
                         beatrice_dir,
                         python         = "python",
                         max_iter       = 2000,
                         n_caus         = 5,
                         sigma_sq       = 0.05,
                         gamma_coverage = 0.95,
                         sparse_concrete = 50) {

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

  if (is.null(variant_ids)) variant_ids <- paste0("rs", seq(0, p - 1))

  beatrice_dir   <- path.expand(beatrice_dir)
  python         <- if (file.exists(path.expand(python))) normalizePath(path.expand(python)) else python

  # Prefer beatrice_annot.py (BEATRICE_annot_sparse fork) when present -
  # upstream beatrice.py has a numpy-2.x bug in trainer.py:calculate_pip that
  # crashes late in training. The fork fixes it and accepts the same flags
  # when --annot is omitted, giving vanilla BEATRICE semantics.
  script_path <- NULL
  for (cand in c("beatrice_annot.py", "beatrice.py")) {
    p_cand <- file.path(beatrice_dir, cand)
    if (file.exists(p_cand)) { script_path <- p_cand; break }
  }

  if (is.null(script_path)) {
    return(.beatrice_error_result(
      p, .beatrice_params(n, max_iter, n_caus, sigma_sq, gamma_coverage,
                          sparse_concrete, beatrice_dir, python),
      0,
      paste("No BEATRICE script (beatrice_annot.py or beatrice.py) found in:",
            beatrice_dir)
    ))
  }

  params <- .beatrice_params(n, max_iter, n_caus, sigma_sq, gamma_coverage,
                             sparse_concrete, beatrice_dir, python)

  # --- Set up temp working directory ------------------------------------------

  work_dir <- tempfile(pattern = "beatrice_run_")
  out_dir  <- file.path(work_dir, "output")
  dir.create(work_dir, recursive = TRUE)
  dir.create(out_dir,  recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  # --- Write .z file ----------------------------------------------------------
  # Format: two space-separated columns, no header.
  # Column 1: variant name, Column 2: z-score.

  z_path <- file.path(work_dir, "region.z")
  z_lines <- paste(variant_ids, z, sep = " ")
  writeLines(z_lines, z_path)

  # --- Write .ld file ---------------------------------------------------------
  # Format: space-separated NxN matrix, no header.

  ld_path <- file.path(work_dir, "region.ld")
  write.table(round(LD, 8), ld_path,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

  # --- Build arguments --------------------------------------------------------

  args <- c(
    script_path,
    "--z",               z_path,
    "--LD",              ld_path,
    "--N",               as.character(as.integer(n)),
    "--target",          paste0(out_dir, "/"),
    "--max_iter",        as.character(as.integer(max_iter)),
    "--n_caus",          as.character(as.integer(n_caus)),
    "--sigma_sq",        as.character(sigma_sq),
    "--gamma_coverage",  as.character(gamma_coverage),
    "--sparse_concrete", as.character(as.integer(sparse_concrete)),
    "--plot_loss",       "False",   # suppress PDF plots in batch mode
    "--get_cred",        "True"
  )

  # --- Run BEATRICE -----------------------------------------------------------

  start_time <- proc.time()

  run_output <- tryCatch({
    system2(python, args = args, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    structure(conditionMessage(e), class = "beatrice_error")
  })

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  if (inherits(run_output, "beatrice_error")) {
    return(.beatrice_error_result(p, params, elapsed, run_output))
  }

  # --- Check output -----------------------------------------------------------

  pip_path  <- file.path(out_dir, "pip.csv")
  cred_path <- file.path(out_dir, "credible_set.txt")

  if (!file.exists(pip_path)) {
    err_msg <- paste(
      c("BEATRICE produced no pip.csv output.", run_output),
      collapse = "\n"
    )
    return(.beatrice_error_result(p, params, elapsed, err_msg))
  }

  # --- Parse pip.csv ----------------------------------------------------------

  pip_df <- tryCatch(
    utils::read.csv(pip_path, stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  if (is.null(pip_df) || !"pip" %in% names(pip_df)) {
    return(.beatrice_error_result(
      p, params, elapsed,
      "Failed to parse BEATRICE pip.csv or 'pip' column missing."
    ))
  }

  # Re-order to match the original variant order.
  # pip.csv has a variant_names column; use it if available.
  if ("variant_names" %in% names(pip_df)) {
    ord <- match(variant_ids, pip_df$variant_names)
    if (anyNA(ord)) ord <- seq_len(p)   # fallback: assume order preserved
  } else {
    ord <- seq_len(p)
  }

  pip <- as.numeric(pip_df$pip[ord])
  pip <- pmax(0, pmin(1, pip))   # clamp for safety

  # --- Parse credible_set.txt -------------------------------------------------
  # Format: one credible set per line, space-separated 0-based variant indices.
  # Convert to 1-based.

  credible_sets <- list()
  cs_pip        <- NULL

  if (file.exists(cred_path)) {
    cred_lines <- readLines(cred_path, warn = FALSE)
    cred_lines <- cred_lines[nchar(trimws(cred_lines)) > 0]

    if (length(cred_lines) > 0) {
      credible_sets <- lapply(cred_lines, function(line) {
        idx_0based <- suppressWarnings(as.integer(strsplit(trimws(line), "\\s+")[[1]]))
        idx_0based <- idx_0based[!is.na(idx_0based)]
        sort(idx_0based + 1L)   # 0-based → 1-based
      })
      credible_sets <- credible_sets[lengths(credible_sets) > 0]
    }
  }

  # Parse conditional inclusion probabilities (optional output)
  cond_path <- file.path(out_dir, "conditional_credible_variants_probability.txt")
  if (file.exists(cond_path)) {
    cond_lines <- readLines(cond_path, warn = FALSE)
    cond_lines <- cond_lines[nchar(trimws(cond_lines)) > 0]
    cs_pip <- lapply(cond_lines, function(line) {
      suppressWarnings(as.numeric(strsplit(trimws(line), "\\s+")[[1]]))
    })
  }

  # --- Return -----------------------------------------------------------------

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = "beatrice",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      cs_pip = cs_pip
    )
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run BEATRICE on a single region from simulation data structures
#'
#' Thin adapter that extracts the appropriate inputs from the simulation's
#' \code{region_geno} and \code{region_pheno} objects and calls
#' \code{\link{run_beatrice}}. Registered in the method registry and called
#' by \code{\link{run_methods}}.
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing \code{LD}, \code{n}, and optionally \code{variant_ids}.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z}.
#' @param ... Additional arguments passed to \code{\link{run_beatrice}}
#'   (e.g. \code{beatrice_dir}, \code{python}, \code{max_iter}).
#'
#' @return The output of \code{\link{run_beatrice}}.
#' @export
run_beatrice_region <- function(region_geno, region_pheno, ...) {
  run_beatrice(
    z           = region_pheno$z,
    LD          = region_geno$LD,
    n           = region_geno$n,
    variant_ids = region_geno$variant_ids,
    ...
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

.beatrice_params <- function(n, max_iter, n_caus, sigma_sq, gamma_coverage,
                             sparse_concrete, beatrice_dir, python) {
  list(
    n               = n,
    max_iter        = max_iter,
    n_caus          = n_caus,
    sigma_sq        = sigma_sq,
    gamma_coverage  = gamma_coverage,
    sparse_concrete = sparse_concrete,
    beatrice_dir    = beatrice_dir,
    python          = python
  )
}

.beatrice_error_result <- function(p, params, elapsed, error_msg) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "beatrice",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(cs_pip = NULL),
    error           = error_msg
  )
}
