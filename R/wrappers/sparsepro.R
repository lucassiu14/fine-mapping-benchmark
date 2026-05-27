# =============================================================================
# sparsepro.R
#
# Wrapper for SparsePro (Zhang et al. 2023) fine-mapping.
#
# SparsePro is a Python CLI tool that performs variational fine-mapping using
# a sparse-projection prior on effect groups. Unlike Funmap (which exposes a
# Python module importable via reticulate), SparsePro is distributed as the
# script `sparsepro_zld.py` in the upstream GitHub repo, so this wrapper
# follows the BEATRICE / FINEMAP / PAINTOR pattern: write input files to a
# temp directory, invoke the script via system2(), parse the output files,
# and return results in the standardised format.
#
# The upstream repo is at https://github.com/zhwm/SparsePro. The user clones
# it once and passes the path via `sparsepro_dir`. See setup_sparsepro().
#
# This file provides:
#   - setup_sparsepro()       : verifies sparsepro_zld.py is reachable and that
#                                the Python env has numpy / scipy / pandas
#   - run_sparsepro()         : runs SparsePro on a single region (explicit
#                                inputs)
#   - run_sparsepro_region()  : adapter called by run_methods()
#
# Standard output format:
#   pip              Numeric vector (length p). Marginal PIPs from <prefix>.pip.
#   credible_sets    List of integer vectors (1-based). One per effect group
#                    in <prefix>.cs. Empty list if no credible sets reported.
#   method           Character. "sparsepro".
#   input_type       Character. Always "summary".
#   params           List. Hyperparameters used.
#   runtime_seconds  Numeric. Wall-clock time.
#   additional       List. SparsePro-specific outputs (see run_sparsepro() docs).
#
# Reference:
#   Zhang W, Najafabadi H, Li Y (2023). SparsePro: an efficient fine-mapping
#   method integrating summary statistics and functional annotations.
#   PLoS Genetics, 19(12), e1011104.
# =============================================================================


# =============================================================================
# Setup
# =============================================================================

#' Set up SparsePro
#'
#' Verifies that the SparsePro script (\code{sparsepro_zld.py}) can be found
#' at \code{sparsepro_dir} and that the required Python packages
#' (\code{numpy}, \code{scipy}, \code{pandas}) are importable by the specified
#' Python executable.
#'
#' Install SparsePro by cloning the upstream repository and installing its
#' requirements into your Python environment:
#'
#' \preformatted{
#'   git clone https://github.com/zhwm/SparsePro
#'   cd SparsePro
#'   pip install -r requirements.txt
#' }
#'
#' SparsePro shares its Python dependencies (\code{numpy}, \code{scipy},
#' \code{pandas}) with BEATRICE and Funmap, so the existing benchmark conda
#' environment is sufficient — \code{pandas} may not be present in older
#' BEATRICE envs and can be added with \code{pip install pandas}.
#'
#' Then pass the repo directory and (if not on PATH) the Python executable:
#'
#' \preformatted{
#'   setup_sparsepro(
#'     sparsepro_dir = "~/SparsePro",
#'     python        = "~/anaconda3/envs/beatrice/bin/python"
#'   )
#' }
#'
#' @param sparsepro_dir Character. Path to the cloned SparsePro repository
#'   root (must contain \code{sparsepro_zld.py}).
#' @param python Character. Path to the Python executable to use.
#'   Default: \code{"python"} (searches PATH).
#'
#' @return Invisibly returns a named list with \code{sparsepro_script} (full
#'   path to \code{sparsepro_zld.py}) and \code{python} (resolved Python path).
#' @export
setup_sparsepro <- function(sparsepro_dir, python = "python") {

  sparsepro_dir <- path.expand(sparsepro_dir)

  script_path <- file.path(sparsepro_dir, "sparsepro_zld.py")
  if (!file.exists(script_path)) {
    stop(
      "sparsepro_zld.py not found in: ", sparsepro_dir, "\n\n",
      "Please clone the upstream repository:\n",
      "  git clone https://github.com/zhwm/SparsePro\n",
      "  cd SparsePro && pip install -r requirements.txt\n\n",
      "Then pass the path:\n",
      "  setup_sparsepro(sparsepro_dir = '/path/to/SparsePro')",
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
        "Pass the full path to your Python:\n",
        "  setup_sparsepro(\n",
        "    sparsepro_dir = '/path/to/SparsePro',\n",
        "    python        = '~/anaconda3/envs/beatrice/bin/python'\n",
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

  missing_pkgs <- Filter(Negate(check_pkg), c("numpy", "scipy", "pandas"))

  if (length(missing_pkgs) > 0) {
    stop(
      "The following Python packages are missing from '", resolved_python, "':\n",
      paste0("  ", missing_pkgs, collapse = "\n"), "\n\n",
      "Install SparsePro's requirements:\n",
      "  pip install -r ", file.path(sparsepro_dir, "requirements.txt"), "\n\n",
      "If you are using the BEATRICE conda env, you may need to add pandas:\n",
      "  conda run -n beatrice pip install pandas",
      call. = FALSE
    )
  }

  message("SparsePro ready.")
  message("  Script : ", script_path)
  message("  Python : ", resolved_python)

  invisible(list(sparsepro_script = script_path, python = resolved_python))
}


# =============================================================================
# Run SparsePro on a single region
# =============================================================================

#' Run SparsePro fine-mapping on a single region
#'
#' Writes the required input files to a temporary directory (z-score table,
#' LD matrix, and the --zld summary file SparsePro expects), invokes
#' \code{sparsepro_zld.py} via the specified Python executable, parses the
#' output \code{.pip} and \code{.cs} files, and returns results in the
#' standardised format. The temporary directory is deleted on exit.
#'
#' SparsePro fits one model per region (one entry in the --zld summary file).
#' The benchmark calls this wrapper one region at a time, so each call
#' produces a single-locus --zld file. This is less efficient than batching
#' many loci together, but keeps the wrapper interface symmetric with the
#' other per-region methods.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param n Integer. Sample size (GWAS N).
#' @param variant_ids Character vector or NULL. Variant identifiers (length p).
#'   Required for SparsePro's output (the \code{.pip} and \code{.cs} files
#'   key on rsids). If NULL, synthetic IDs of the form \code{"snp1", "snp2",
#'   ...} are generated.
#' @param sparsepro_dir Character. Path to the cloned SparsePro repository
#'   (must contain \code{sparsepro_zld.py}).
#' @param python Character. Path to the Python executable. Default:
#'   \code{"python"}.
#' @param K Integer. Maximum number of effect groups (causal signals)
#'   SparsePro will consider. Default: 5.
#' @param cthres Numeric. Coverage level for credible sets. Default: 0.95.
#'
#' @return A list with the standardised fine-mapping output:
#' \describe{
#'   \item{pip}{Numeric vector (length p). Posterior inclusion probabilities,
#'     re-ordered to match the input \code{variant_ids}.}
#'   \item{credible_sets}{List of integer vectors (1-based). One per effect
#'     group in the \code{.cs} file. Empty list if no credible sets were
#'     reported.}
#'   \item{method}{Character. Always \code{"sparsepro"}.}
#'   \item{input_type}{Character. Always \code{"summary"}.}
#'   \item{params}{List. Hyperparameters used.}
#'   \item{runtime_seconds}{Numeric. Wall-clock time in seconds.}
#'   \item{additional}{List of SparsePro-specific outputs:
#'     \describe{
#'       \item{cs_pip}{List of numeric vectors. Per-variant inclusion
#'         probabilities within each credible set (the \code{pip} column of
#'         the \code{.cs} file), in the same order as \code{credible_sets}.
#'         NULL if not parseable.}
#'       \item{cs_effect_size}{List of numeric vectors. Posterior effect
#'         sizes per variant per credible set (the \code{effect_size} column
#'         of the \code{.cs} file). NULL if not parseable.}
#'       \item{run_log}{Character. Captured stdout/stderr from the SparsePro
#'         invocation; useful for debugging unexpected failures.}
#'     }
#'   }
#'   \item{error}{Character or NULL. Error message if SparsePro failed.}
#' }
#'
#' @export
run_sparsepro <- function(z,
                          LD,
                          n,
                          variant_ids   = NULL,
                          sparsepro_dir,
                          python        = "python",
                          K             = 5,
                          cthres        = 0.95) {

  # --- Validate inputs --------------------------------------------------------

  p <- length(z)

  stopifnot(
    "LD must be a p x p matrix" =
      is.matrix(LD) && nrow(LD) == p && ncol(LD) == p,
    "n must be a positive integer" =
      is.numeric(n) && length(n) == 1 && n > 0,
    "K must be a positive integer" =
      is.numeric(K) && length(K) == 1 && K >= 1,
    "cthres must be in (0, 1)" =
      is.numeric(cthres) && length(cthres) == 1 &&
      cthres > 0 && cthres < 1
  )

  if (is.null(variant_ids)) variant_ids <- paste0("snp", seq_len(p))

  sparsepro_dir <- path.expand(sparsepro_dir)
  python        <- if (file.exists(path.expand(python))) normalizePath(path.expand(python)) else python
  script_path   <- file.path(sparsepro_dir, "sparsepro_zld.py")

  params <- .sparsepro_params(n, K, cthres, sparsepro_dir, python)

  if (!file.exists(script_path)) {
    return(.sparsepro_error_result(
      p, params, 0,
      paste("sparsepro_zld.py not found at:", script_path,
            "\nRun setup_sparsepro() for instructions.")
    ))
  }

  # --- Set up temp working directory ------------------------------------------

  work_dir <- tempfile(pattern = "sparsepro_run_")
  out_dir  <- file.path(work_dir, "output")
  dir.create(work_dir, recursive = TRUE)
  dir.create(out_dir,  recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  region_label <- "region"
  zfile_name   <- paste0(region_label, ".z")
  ldfile_name  <- paste0(region_label, ".ld")

  # --- Write z-score file -----------------------------------------------------
  # Format: 2-column TSV, no header. Column 1 = variant ID, column 2 = z-score.

  z_path  <- file.path(work_dir, zfile_name)
  z_df    <- data.frame(rsid = variant_ids, z = as.numeric(z),
                        stringsAsFactors = FALSE)
  write.table(z_df, z_path, sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  # --- Write LD matrix file ---------------------------------------------------
  # SparsePro expects whitespace-separated Pearson correlations, no header.

  ld_path <- file.path(work_dir, ldfile_name)
  write.table(round(LD, 8), ld_path,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

  # --- Write --zld summary file -----------------------------------------------
  # One row per locus: <zscore_file>\t<ld_file> (plus optional annotation
  # column, which we leave out). File paths are relative to --zdir, so we
  # use the bare filenames here.

  zld_path <- file.path(work_dir, "zld.txt")
  writeLines(paste(zfile_name, ldfile_name, sep = "\t"), zld_path)

  # --- Build arguments --------------------------------------------------------

  args <- c(
    script_path,
    "--zld",    zld_path,
    "--zdir",   work_dir,
    "--N",      as.character(as.integer(n)),
    "--save",   out_dir,
    "--prefix", region_label,
    "--K",      as.character(as.integer(K)),
    "--cthres", as.character(cthres),
    "--verbose"
  )

  # --- Run SparsePro ----------------------------------------------------------

  start_time <- proc.time()

  run_output <- tryCatch({
    system2(python, args = args, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    structure(conditionMessage(e), class = "sparsepro_error")
  })

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  if (inherits(run_output, "sparsepro_error")) {
    return(.sparsepro_error_result(p, params, elapsed, as.character(run_output),
                                    run_log = character(0)))
  }

  run_log <- if (is.character(run_output)) run_output else character(0)

  # --- Locate output files ----------------------------------------------------
  # SparsePro writes <prefix>.pip and <prefix>.cs into --save.

  pip_path <- file.path(out_dir, paste0(region_label, ".pip"))
  cs_path  <- file.path(out_dir, paste0(region_label, ".cs"))

  if (!file.exists(pip_path)) {
    err_msg <- paste(
      c("SparsePro produced no .pip output.", run_log),
      collapse = "\n"
    )
    return(.sparsepro_error_result(p, params, elapsed, err_msg, run_log = run_log))
  }

  # --- Parse .pip file --------------------------------------------------------
  # Three columns: variant_id, z-score, pip. Whitespace-separated. The order
  # may differ from the input variant order; we re-index against variant_ids.

  pip_df <- tryCatch(
    utils::read.table(pip_path, header = FALSE, sep = "",
                      stringsAsFactors = FALSE,
                      col.names = c("rsid", "z", "pip")),
    error = function(e) NULL
  )

  if (is.null(pip_df) || !"pip" %in% names(pip_df)) {
    return(.sparsepro_error_result(
      p, params, elapsed,
      "Failed to parse SparsePro .pip output.",
      run_log = run_log
    ))
  }

  # Re-index PIPs to match the input variant_ids order. Variants present in
  # the input but missing from .pip get NA.
  pip <- numeric(p)
  ord <- match(variant_ids, pip_df$rsid)
  pip <- pip_df$pip[ord]
  pip <- pmax(0, pmin(1, suppressWarnings(as.numeric(pip))))

  # --- Parse .cs file ---------------------------------------------------------
  # SparsePro emits one row per effect group (credible set). Columns:
  #   cs            : whitespace- or comma-separated rsids in the group
  #   pip           : whitespace- or comma-separated within-group PIPs
  #   effect_size   : whitespace- or comma-separated effect sizes
  #
  # The exact intra-cell separator varies between SparsePro versions. We
  # parse defensively: try comma-split first, fall back to whitespace.

  credible_sets  <- list()
  cs_pip         <- NULL
  cs_effect_size <- NULL

  if (file.exists(cs_path) && file.info(cs_path)$size > 0L) {
    cs_lines <- readLines(cs_path, warn = FALSE)
    cs_lines <- cs_lines[nchar(trimws(cs_lines)) > 0L]

    # Skip a header row if present (column names commonly start with "cs").
    if (length(cs_lines) > 0L && grepl("^(cs|set|group)\\b", cs_lines[1L],
                                       ignore.case = TRUE)) {
      cs_lines <- cs_lines[-1L]
    }

    if (length(cs_lines) > 0L) {
      parsed <- lapply(cs_lines, function(line) {
        # Split top-level on tabs (between the three columns).
        cells <- strsplit(line, "\t", fixed = TRUE)[[1L]]
        if (length(cells) < 1L) return(NULL)

        # Robust intra-cell split.
        split_cell <- function(s) {
          s <- trimws(s)
          if (!nzchar(s)) return(character(0))
          if (grepl(",", s, fixed = TRUE)) strsplit(s, ",\\s*")[[1L]]
          else strsplit(s, "\\s+")[[1L]]
        }

        rsids  <- split_cell(cells[1L])
        pips   <- if (length(cells) >= 2L) suppressWarnings(as.numeric(split_cell(cells[2L]))) else NA_real_
        effs   <- if (length(cells) >= 3L) suppressWarnings(as.numeric(split_cell(cells[3L]))) else NA_real_

        idx <- match(rsids, variant_ids)
        idx <- idx[!is.na(idx)]
        list(idx = sort(unique(idx)), pip = pips, effect = effs)
      })

      parsed <- Filter(function(p) !is.null(p) && length(p$idx) > 0L, parsed)

      if (length(parsed) > 0L) {
        credible_sets  <- lapply(parsed, `[[`, "idx")
        cs_pip         <- lapply(parsed, `[[`, "pip")
        cs_effect_size <- lapply(parsed, `[[`, "effect")
      }
    }
  }

  # --- Return -----------------------------------------------------------------

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = "sparsepro",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      cs_pip         = cs_pip,
      cs_effect_size = cs_effect_size,
      run_log        = run_log
    )
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run SparsePro on a single region from simulation data structures
#'
#' Thin adapter that extracts the appropriate inputs from the simulation's
#' \code{region_geno} and \code{region_pheno} objects and calls
#' \code{\link{run_sparsepro}}.
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing \code{LD}, \code{n}, and optionally \code{variant_ids}.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z}.
#' @param ... Additional arguments passed to \code{\link{run_sparsepro}}
#'   (e.g. \code{sparsepro_dir}, \code{python}, \code{K}, \code{cthres}).
#'
#' @return The output of \code{\link{run_sparsepro}}.
#' @export
run_sparsepro_region <- function(region_geno, region_pheno, ...) {
  run_sparsepro(
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

.sparsepro_params <- function(n, K, cthres, sparsepro_dir, python) {
  list(
    n             = n,
    K             = K,
    cthres        = cthres,
    sparsepro_dir = sparsepro_dir,
    python        = python
  )
}

.sparsepro_error_result <- function(p, params, elapsed, error_msg,
                                     run_log = character(0)) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "sparsepro",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      cs_pip         = NULL,
      cs_effect_size = NULL,
      run_log        = run_log
    ),
    error           = error_msg
  )
}
