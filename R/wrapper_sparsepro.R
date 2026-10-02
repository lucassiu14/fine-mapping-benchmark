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
#   - setup_sparsepro()              : verifies sparsepro_zld.py is reachable
#                                       and that the Python env has numpy /
#                                       scipy / pandas
#   - run_sparsepro()                : runs SparsePro on a single region
#                                       without annotations (explicit inputs)
#   - run_sparsepro_annotated()      : runs SparsePro's annotation workflow
#                                       over several regions (explicit inputs)
#   - run_sparsepro_scenario_setup() : per-scenario hook for run_methods()
#   - run_sparsepro_region()         : adapter called by run_methods()
#
# Annotation mode
# ---------------
# When the simulated regions carry annotations, run_methods() uses SparsePro's
# own two-run workflow, exactly as the upstream README runs it:
#   1. "Fine-mapping with GWAS summary statistics and estimating annotation
#      enrichment": every region of the scenario is listed in one --zld file
#      with an annotation column, and the script is called with
#      --anno anno --pthres <pthres>. It fine-maps each region under a uniform
#      prior, pools the PIPs over all listed regions, runs a G-test per
#      annotation and writes joint enrichment weights for the annotations with
#      p < pthres (<prefix>.W<pthres>) and for all annotations (<prefix>.W1.0).
#   2. "Fine-mapping with both GWAS summary statistics and functional
#      annotations selected by G-test": the script is called again with
#      --zld <prefix>.h2 --anno anno --aW <prefix>.W<pthres>, which fine-maps
#      each region under the prior softmax(A W).
# If no annotation passes the G-test, the W<pthres> file is not written and,
# as in the paper ("SparsePro+ was not implemented"), the run-1 result is kept.
# prior = "all" uses W1.0 instead (the paper's "SparsePro+1.0").
# The script's enrichment estimates are relative risks of binary annotation
# status, so the annotation workflow is used only for 0/1 annotations; other
# annotations are fine-mapped without them and the reason is recorded.
# The script is only called, never copied: SparsePro is licensed
# CC BY-NC-ND 4.0.
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
#'       \item{annotation_mode}{Character. Always \code{"none"} here; see
#'         \code{\link{run_sparsepro_annotated}} for the annotation workflow.}
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

  variant_ids <- .sparsepro_ids(variant_ids, p)

  sparsepro_dir <- path.expand(sparsepro_dir)
  python        <- .sparsepro_resolve_python(python)
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

  # --- Write inputs -----------------------------------------------------------

  files <- .sparsepro_write_locus(work_dir, "region", z, LD, variant_ids)

  # SparsePro's sparsepro_zld.py reads this via
  #   pd.read_csv(args.zld, sep='\s+')
  # which treats the first row as a HEADER. Without a header line the
  # loop over ldlists iterates zero times, no .pip is written, and the
  # wrapper reports "produced no .pip output". Include the header explicitly.
  # File paths are relative to --zdir, so the bare filenames are used.
  zld_path <- file.path(work_dir, "zld.txt")
  writeLines(c("z\tld", paste(files$z, files$ld, sep = "\t")), zld_path)

  # --- Run SparsePro ----------------------------------------------------------

  args <- .sparsepro_args(script_path, zld_path, work_dir, n, out_dir,
                          prefix = "region", K = K, cthres = cthres)

  start_time <- proc.time()
  run        <- .sparsepro_exec(python, args)
  elapsed    <- as.numeric((proc.time() - start_time)["elapsed"])

  if (!is.null(run$error)) {
    return(.sparsepro_error_result(p, params, elapsed, run$error,
                                    run_log = character(0)))
  }

  parsed <- .sparsepro_parse_locus(out_dir, files$z, variant_ids)
  if (!is.null(parsed$error)) {
    return(.sparsepro_error_result(
      p, params, elapsed,
      paste(c(parsed$error, run$log), collapse = "\n"),
      run_log = run$log
    ))
  }

  # --- Return -----------------------------------------------------------------

  list(
    pip             = parsed$pip,
    credible_sets   = parsed$credible_sets,
    method          = "sparsepro",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      cs_pip          = parsed$cs_pip,
      cs_effect_size  = parsed$cs_effect_size,
      run_log         = run$log,
      annotation_mode = "none"
    )
  )
}


# =============================================================================
# Run SparsePro's annotation workflow over several regions
# =============================================================================

#' Run SparsePro with functional annotations over several regions
#'
#' Runs the two-step annotation workflow of the upstream SparsePro README on
#' all supplied regions at once, so that enrichment is estimated from the
#' regions pooled, as \code{sparsepro_zld.py} does for every locus listed in
#' its \code{--zld} file.
#'
#' Step 1 calls the script with \code{--anno} and \code{--pthres}: each region
#' is fine-mapped under a uniform prior, a G-test is run for every annotation
#' on the pooled posterior inclusion probabilities, and joint enrichment
#' weights are written for the annotations with p-value below
#' \code{pthres}. Step 2 calls the script with \code{--aW} on that weight
#' file, which fine-maps each region again under the prior
#' \code{softmax(A W)}. If no annotation passes the G-test, no weight file is
#' written and the step-1 result is returned, as in Zhang et al. (2023).
#' With \code{prior = "all"}, step 2 uses the weights of all annotations
#' (the \code{W1.0} file) instead.
#'
#' SparsePro's enrichment estimates are relative risks of binary annotation
#' status, so every annotation must be 0/1.
#'
#' @param z_list List of numeric vectors. Marginal z-scores, one per region.
#' @param ld_list List of LD matrices, one per region.
#' @param annot_list List of annotation matrices (\code{p_r x d}, entries 0 or
#'   1), one per region, all with the same \code{d} columns in the same order.
#' @param n Integer. GWAS sample size, shared by all regions.
#' @param variant_ids_list List of character vectors, or NULL. Variant
#'   identifiers per region (must be unique within a region). Synthetic
#'   identifiers are used when NULL.
#' @param sparsepro_dir Character. Path to the cloned SparsePro repository.
#' @param python Character. Python executable. Default: \code{"python"}.
#' @param K Integer. Maximum number of effect groups. Default: 5.
#' @param cthres Numeric. Credible-set coverage. Default: 0.95.
#' @param pthres Numeric. G-test p-value threshold for selecting annotations.
#'   Default: \code{1e-5}, the value used in the paper and README.
#' @param prior Character. \code{"significant"} (default) uses the annotations
#'   selected by the G-test; \code{"all"} uses all annotations.
#'
#' @return A list with one standardised result per region, in input order.
#'   \code{additional} holds \code{annotation_mode} (\code{"significant"},
#'   \code{"all"}, or \code{"no_significant_annotation"} when step 2 was not
#'   run), \code{gtest} (the step-1 G-test table), \code{enrichment} (the
#'   weights used in step 2, or NULL), \code{pip_without_annotations} (the
#'   step-1 PIPs), \code{n_regions_pooled} and \code{run_log}. On failure
#'   every element carries the same \code{error}. The runtime of the pooled
#'   call is divided equally between the regions.
#' @export
run_sparsepro_annotated <- function(z_list,
                                    ld_list,
                                    annot_list,
                                    n,
                                    variant_ids_list = NULL,
                                    sparsepro_dir,
                                    python   = "python",
                                    K        = 5,
                                    cthres   = 0.95,
                                    pthres   = 1e-5,
                                    prior    = c("significant", "all")) {

  prior <- match.arg(prior)
  R     <- length(z_list)
  p_vec <- lengths(z_list)

  stopifnot(
    "z_list must contain at least one region" = R >= 1L,
    "ld_list and annot_list must have one entry per region" =
      length(ld_list) == R && length(annot_list) == R,
    "n must be a positive number" =
      is.numeric(n) && length(n) == 1 && n > 0,
    "K must be a positive integer" =
      is.numeric(K) && length(K) == 1 && K >= 1,
    "cthres must be in (0, 1)" =
      is.numeric(cthres) && length(cthres) == 1 &&
      cthres > 0 && cthres < 1,
    "pthres must be in (0, 1]" =
      is.numeric(pthres) && length(pthres) == 1 &&
      pthres > 0 && pthres <= 1
  )
  for (r in seq_len(R)) {
    if (!is.matrix(ld_list[[r]]) || nrow(ld_list[[r]]) != p_vec[r] ||
        ncol(ld_list[[r]]) != p_vec[r]) {
      stop("LD matrix ", r, " must be ", p_vec[r], " x ", p_vec[r], call. = FALSE)
    }
  }
  problem <- .sparsepro_annotation_problem(annot_list, p_vec)
  if (!is.null(problem)) stop(problem, call. = FALSE)

  d <- ncol(annot_list[[1]])
  annot_names <- paste0("ANNOT", seq_len(d))
  original_names <- colnames(annot_list[[1]])
  if (is.null(original_names)) original_names <- annot_names

  ids <- lapply(seq_len(R), function(r) {
    .sparsepro_ids(if (is.null(variant_ids_list)) NULL else variant_ids_list[[r]],
                   p_vec[r])
  })
  for (r in seq_len(R)) {
    if (anyDuplicated(ids[[r]])) {
      stop("variant identifiers must be unique within region ", r, call. = FALSE)
    }
  }

  sparsepro_dir <- path.expand(sparsepro_dir)
  python        <- .sparsepro_resolve_python(python)
  script_path   <- file.path(sparsepro_dir, "sparsepro_zld.py")

  params <- c(.sparsepro_params(n, K, cthres, sparsepro_dir, python),
              list(pthres = pthres, prior = prior, n_annotations = d))

  fail <- function(msg, elapsed = 0, run_log = character(0)) {
    lapply(seq_len(R), function(r) {
      .sparsepro_error_result(p_vec[r], params, elapsed / R, msg,
                              run_log = run_log)
    })
  }

  if (!file.exists(script_path)) {
    return(fail(paste("sparsepro_zld.py not found at:", script_path,
                      "\nRun setup_sparsepro() for instructions.")))
  }

  # --- Write every region as one locus of a single --zld list -----------------

  work_dir <- tempfile(pattern = "sparsepro_annot_")
  out1     <- file.path(work_dir, "step1")
  out2     <- file.path(work_dir, "step2")
  dir.create(out1, recursive = TRUE)
  dir.create(out2, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  files <- lapply(seq_len(R), function(r) {
    A <- annot_list[[r]]
    colnames(A) <- annot_names
    .sparsepro_write_locus(work_dir, paste0("region", r), z_list[[r]],
                           ld_list[[r]], ids[[r]], A = A)
  })
  zld_path <- file.path(work_dir, "zld.txt")
  writeLines(c("z\tld\tanno",
               vapply(files, function(f) paste(f$z, f$ld, f$anno, sep = "\t"),
                      character(1))),
             zld_path)

  # The weight file is named with Python's str() of --pthres
  # ('{}'.format(args.pthres)), so the same string is passed and formatted.
  pthres_arg <- format(pthres, digits = 15)

  start_time <- proc.time()
  elapsed_now <- function() as.numeric((proc.time() - start_time)["elapsed"])

  # --- Step 1: fine-map without a prior and estimate enrichment ---------------

  args1 <- c(.sparsepro_args(script_path, zld_path, work_dir, n, out1,
                             prefix = "step1", K = K, cthres = cthres),
             "--anno", "anno", "--pthres", pthres_arg)
  run1 <- .sparsepro_exec(python, args1)
  if (!is.null(run1$error)) return(fail(run1$error, elapsed_now()))

  step1 <- lapply(seq_len(R), function(r) {
    .sparsepro_parse_locus(out1, files[[r]]$z, ids[[r]])
  })
  bad <- Filter(function(x) !is.null(x$error), step1)
  if (length(bad) > 0L) {
    return(fail(paste(c(paste("Step 1:", bad[[1]]$error), run1$log),
                      collapse = "\n"),
                elapsed_now(), run1$log))
  }

  gtest <- .sparsepro_read_table(file.path(out1, "step1.wsep"))
  if (!is.null(gtest)) gtest$annotation <- .sparsepro_rename(gtest$index,
                                                             annot_names,
                                                             original_names)

  w_name <- if (prior == "all") {
    "step1.W1.0"
  } else {
    fmt <- .sparsepro_exec(python, c("-c", shQuote(sprintf(
      "print('{}'.format(float('%s')))", pthres_arg))))
    if (!is.null(fmt$error) || length(fmt$log) < 1L) {
      return(fail(paste("Could not format pthres with", python),
                  elapsed_now(), run1$log))
    }
    paste0("step1.W", trimws(fmt$log[length(fmt$log)]))
  }
  w_path <- file.path(out1, w_name)

  # --- Step 2: fine-map under the prior softmax(A W) ---------------------------

  if (!file.exists(w_path)) {
    final      <- step1
    mode       <- "no_significant_annotation"
    enrichment <- NULL
    run_log    <- run1$log
  } else {
    h2_path <- file.path(out1, "step1.h2")
    args2 <- c(.sparsepro_args(script_path, h2_path, work_dir, n, out2,
                               prefix = "step2", K = K, cthres = cthres),
               "--anno", "anno", "--aW", shQuote(w_path))
    run2 <- .sparsepro_exec(python, args2)
    run_log <- c(run1$log, run2$log)
    if (!is.null(run2$error)) return(fail(run2$error, elapsed_now(), run_log))

    final <- lapply(seq_len(R), function(r) {
      .sparsepro_parse_locus(out2, files[[r]]$z, ids[[r]])
    })
    bad <- Filter(function(x) !is.null(x$error), final)
    if (length(bad) > 0L) {
      return(fail(paste(c(paste("Step 2:", bad[[1]]$error), run_log),
                        collapse = "\n"),
                  elapsed_now(), run_log))
    }
    mode       <- prior
    enrichment <- .sparsepro_read_table(w_path)
    if (!is.null(enrichment)) {
      enrichment$annotation <- .sparsepro_rename(enrichment$index,
                                                 annot_names, original_names)
    }
  }

  elapsed <- elapsed_now()

  lapply(seq_len(R), function(r) {
    list(
      pip             = final[[r]]$pip,
      credible_sets   = final[[r]]$credible_sets,
      method          = "sparsepro",
      input_type      = "summary",
      params          = params,
      runtime_seconds = elapsed / R,
      additional      = list(
        cs_pip                  = final[[r]]$cs_pip,
        cs_effect_size          = final[[r]]$cs_effect_size,
        run_log                 = run_log,
        annotation_mode         = mode,
        gtest                   = gtest,
        enrichment              = enrichment,
        pip_without_annotations = step1[[r]]$pip,
        n_regions_pooled        = R
      )
    )
  })
}


# =============================================================================
# Scenario-level setup: run the annotation workflow on all regions together
# =============================================================================

#' Run SparsePro's annotation workflow once per scenario
#'
#' Called by \code{\link{run_methods}} once per scenario. When the regions
#' carry annotations (and \code{use_annotations} is not \code{FALSE}), it runs
#' \code{\link{run_sparsepro_annotated}} on all regions of the scenario
#' together, so that enrichment is estimated from the pooled regions, and
#' returns the per-region results keyed by a fingerprint of each region's
#' z-scores. Without annotations it returns an empty list and every region is
#' fine-mapped on its own by \code{\link{run_sparsepro}}, as before.
#'
#' Annotations that are not 0/1 cannot be used by SparsePro's enrichment
#' estimator; the regions are then fine-mapped without them and the reason is
#' passed on as \code{.sparsepro_note}. A failed annotation run is passed on
#' as \code{.sparsepro_error}, so those regions are reported as failures
#' rather than silently fine-mapped without annotations.
#'
#' @param genotypes List. \code{simulation$genotypes}.
#' @param regions List. One scenario's \code{regions}.
#' @param user_args List. The user's method_args for sparsepro.
#'
#' @return A list with one of \code{.sparsepro_cache}, \code{.sparsepro_note}
#'   or \code{.sparsepro_error}, or an empty list.
#' @export
run_sparsepro_scenario_setup <- function(genotypes, regions, user_args) {

  if (isFALSE(user_args$use_annotations)) return(list())

  R <- length(regions)
  A_list <- lapply(seq_len(R), function(i) {
    .sparsepro_extract_annotations(genotypes[[i]], regions[[i]])
  })
  if (all(vapply(A_list, is.null, logical(1)))) return(list())

  p_vec   <- vapply(regions, function(x) length(x$z), integer(1))
  problem <- .sparsepro_annotation_problem(A_list, p_vec)
  if (!is.null(problem)) return(list(.sparsepro_note = problem))

  n_vals <- unique(vapply(genotypes[seq_len(R)], function(g) as.numeric(g$n),
                          numeric(1)))
  if (length(n_vals) != 1L) {
    return(list(.sparsepro_error = paste(
      "SparsePro takes one sample size for all regions, but the regions",
      "have different n; the annotation workflow cannot pool them.")))
  }
  if (is.null(user_args$sparsepro_dir)) {
    return(list(.sparsepro_error = "sparsepro_dir must be supplied."))
  }

  fits <- tryCatch(
    run_sparsepro_annotated(
      z_list           = lapply(regions, `[[`, "z"),
      ld_list          = lapply(genotypes[seq_len(R)], `[[`, "LD"),
      annot_list       = A_list,
      n                = n_vals,
      variant_ids_list = lapply(genotypes[seq_len(R)], `[[`, "variant_ids"),
      sparsepro_dir    = user_args$sparsepro_dir,
      python           = user_args$python %||% "python",
      K                = user_args$K %||% 5,
      cthres           = user_args$cthres %||% 0.95,
      pthres           = user_args$pthres %||% 1e-5,
      prior            = user_args$prior %||% "significant"
    ),
    error = function(e) conditionMessage(e)
  )
  if (is.character(fits)) {
    return(list(.sparsepro_error = paste("SparsePro annotation run failed:", fits)))
  }
  if (!is.null(fits[[1]]$error)) {
    return(list(.sparsepro_error = paste("SparsePro annotation run failed:",
                                         fits[[1]]$error)))
  }

  keys <- vapply(regions, function(x) .fb_fingerprint(x$z), character(1))
  list(.sparsepro_cache = stats::setNames(fits, keys))
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run SparsePro on a single region from simulation data structures
#'
#' Adapter called by \code{\link{run_methods}}. When the scenario hook
#' \code{\link{run_sparsepro_scenario_setup}} has run SparsePro's annotation
#' workflow, the region's result is taken from it. When the region carries
#' 0/1 annotations but no hook result is available (a direct call), the
#' annotation workflow is run on this region alone. Otherwise the region is
#' fine-mapped without annotations by \code{\link{run_sparsepro}}.
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing \code{LD}, \code{n}, optionally \code{variant_ids} and
#'   \code{annotations_matrix}.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z}.
#' @param .sparsepro_cache,.sparsepro_note,.sparsepro_error Set by
#'   \code{\link{run_sparsepro_scenario_setup}}; not for direct use.
#' @param use_annotations Logical. Use the regions' annotations when present.
#'   Default: TRUE. Set to FALSE to reproduce runs made before the annotation
#'   workflow was added (Iterations 004 to 007).
#' @param pthres,prior Passed to \code{\link{run_sparsepro_annotated}}.
#' @param ... Additional arguments passed to \code{\link{run_sparsepro}}
#'   (e.g. \code{sparsepro_dir}, \code{python}, \code{K}, \code{cthres}).
#'
#' @return A standardised result list; \code{additional$annotation_mode}
#'   records whether and how annotations were used.
#' @export
run_sparsepro_region <- function(region_geno, region_pheno,
                                 .sparsepro_cache = NULL,
                                 .sparsepro_note  = NULL,
                                 .sparsepro_error = NULL,
                                 use_annotations  = TRUE,
                                 pthres           = 1e-5,
                                 prior            = "significant",
                                 ...) {
  p <- length(region_pheno$z)
  if (!is.null(.sparsepro_error)) {
    return(.sparsepro_error_result(p, .sparsepro_dots_params(region_geno, ...),
                                   0, .sparsepro_error))
  }
  if (!is.null(.sparsepro_cache)) {
    hit <- .sparsepro_cache[[.fb_fingerprint(region_pheno$z)]]
    if (is.null(hit)) {
      return(.sparsepro_error_result(
        p, .sparsepro_dots_params(region_geno, ...), 0,
        "SparsePro annotation run returned no result for this region."))
    }
    return(hit)
  }

  A <- .sparsepro_extract_annotations(region_geno, region_pheno)
  if (isTRUE(use_annotations) && !is.null(A) && is.null(.sparsepro_note)) {
    .sparsepro_note <- .sparsepro_annotation_problem(list(A), p)
    if (is.null(.sparsepro_note)) {
      return(run_sparsepro_annotated(
        z_list           = list(region_pheno$z),
        ld_list          = list(region_geno$LD),
        annot_list       = list(A),
        n                = region_geno$n,
        variant_ids_list = list(region_geno$variant_ids),
        pthres           = pthres,
        prior            = prior,
        ...
      )[[1]])
    }
  }

  fit <- run_sparsepro(
    z           = region_pheno$z,
    LD          = region_geno$LD,
    n           = region_geno$n,
    variant_ids = region_geno$variant_ids,
    ...
  )
  if (!is.null(.sparsepro_note) && isTRUE(use_annotations)) {
    fit$additional$annotation_note <- .sparsepro_note
  }
  fit
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

# params for an error result raised before run_sparsepro() is reached.
.sparsepro_dots_params <- function(region_geno, ...) {
  a <- list(...)
  .sparsepro_params(region_geno$n, a$K %||% 5, a$cthres %||% 0.95,
                    a$sparsepro_dir %||% NA_character_, a$python %||% "python")
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

# Prefer region_geno$annotations_matrix, falling back to region_pheno's copy
# (same convention as the PAINTOR and Funmap wrappers).
.sparsepro_extract_annotations <- function(region_geno, region_pheno) {
  A <- region_geno$annotations_matrix
  if (is.null(A)) A <- region_pheno$annotations_matrix
  A
}

# NULL when the annotation matrices can be used by SparsePro's workflow,
# otherwise the reason they cannot.
.sparsepro_annotation_problem <- function(annot_list, p_vec) {
  if (any(vapply(annot_list, is.null, logical(1)))) {
    return("Some regions have no annotation matrix; SparsePro's annotation workflow needs one for every region.")
  }
  ok_shape <- vapply(seq_along(annot_list), function(r) {
    A <- annot_list[[r]]
    (is.matrix(A) || is.data.frame(A)) && nrow(A) == p_vec[r] && ncol(A) >= 1L
  }, logical(1))
  if (!all(ok_shape)) {
    return("Each annotation matrix must have one row per variant and at least one column.")
  }
  if (length(unique(vapply(annot_list, ncol, integer(1)))) != 1L) {
    return("The regions' annotation matrices have different numbers of columns.")
  }
  binary <- vapply(annot_list, function(A) {
    v <- as.matrix(A)
    is.numeric(v) && !anyNA(v) && all(v == 0 | v == 1)
  }, logical(1))
  if (!all(binary)) {
    return("SparsePro's enrichment estimator (relative risk and G-test of annotation status) needs 0/1 annotations, so these annotations were not used.")
  }
  NULL
}

# VCF-derived variant_ids are "1 40023356 . A T" (embedded spaces).
# SparsePro echoes the id back as the INDEX of its .pip output, and we
# parse that file with read.table(sep = "") - i.e. split on ANY
# whitespace - expecting exactly 3 fields (rsid, z, pip). A spaced id
# yields 7 fields and the parse fails. Collapsing whitespace keeps each
# id a single token on both the way out and the way back, and keeps
# match(variant_ids, pip_df$rsid) consistent since both sides are
# sanitised the same way.
.sparsepro_ids <- function(variant_ids, p) {
  if (is.null(variant_ids)) variant_ids <- paste0("snp", seq_len(p))
  gsub("\\s+", "_", variant_ids)
}

.sparsepro_resolve_python <- function(python) {
  if (file.exists(path.expand(python))) normalizePath(path.expand(python)) else python
}

# Writes <label>.z (id, z; no header), <label>.ld (whitespace-separated
# correlations; no header) and, when A is given, <label>.anno (tab-separated,
# header "SNP" then one column per annotation, as in the upstream
# dat/anno.txt). Returns the file names relative to `dir`.
.sparsepro_write_locus <- function(dir, label, z, LD, variant_ids, A = NULL) {
  zfile  <- paste0(label, ".z")
  ldfile <- paste0(label, ".ld")
  utils::write.table(
    data.frame(rsid = variant_ids, z = as.numeric(z), stringsAsFactors = FALSE),
    file.path(dir, zfile), sep = "\t",
    quote = FALSE, row.names = FALSE, col.names = FALSE)
  utils::write.table(round(LD, 8), file.path(dir, ldfile),
                     quote = FALSE, row.names = FALSE, col.names = FALSE,
                     sep = " ")
  annofile <- NULL
  if (!is.null(A)) {
    annofile <- paste0(label, ".anno")
    A   <- as.matrix(A)
    adf <- data.frame(SNP = variant_ids, stringsAsFactors = FALSE)
    for (k in seq_len(ncol(A))) adf[[colnames(A)[k]]] <- as.integer(A[, k])
    utils::write.table(adf, file.path(dir, annofile), sep = "\t",
                       quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
  list(z = zfile, ld = ldfile, anno = annofile)
}

.sparsepro_args <- function(script_path, zld_path, zdir, n, save, prefix,
                            K, cthres) {
  c(
    shQuote(script_path),
    "--zld",    shQuote(zld_path),
    "--zdir",   shQuote(zdir),
    "--N",      as.character(as.integer(n)),
    "--save",   shQuote(save),
    "--prefix", prefix,
    "--K",      as.character(as.integer(K)),
    "--cthres", as.character(cthres),
    "--verbose"
  )
}

# Runs the script; returns list(error = NULL or message, log = output lines).
.sparsepro_exec <- function(python, args) {
  out <- tryCatch(
    system2(python, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) structure(conditionMessage(e), class = "sparsepro_error")
  )
  if (inherits(out, "sparsepro_error")) {
    return(list(error = as.character(out), log = character(0)))
  }
  list(error = NULL, log = if (is.character(out)) as.character(out) else character(0))
}

.sparsepro_read_table <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(utils::read.delim(path, stringsAsFactors = FALSE),
           error = function(e) NULL)
}

.sparsepro_rename <- function(x, internal, original) {
  out <- original[match(x, internal)]
  ifelse(is.na(out), x, out)
}

# Parses <zfile>.pip and <zfile>.cs in `out_dir`. Returns list(pip,
# credible_sets, cs_pip, cs_effect_size), or list(error = message).
.sparsepro_parse_locus <- function(out_dir, zfile, variant_ids) {
  p <- length(variant_ids)

  # SparsePro names .pip/.cs after the Z FILENAME from the --zld list, not
  # after --prefix (sparsepro_zld.py):
  #     z.to_csv(os.path.join(args.save, "{}.pip".format(zfile)), ...)
  # so with zfile = "region.z" it writes "region.z.pip" / "region.z.cs".
  # (--prefix only names the .h2/.wsep/.W* files.)
  pip_path <- file.path(out_dir, paste0(zfile, ".pip"))
  cs_path  <- file.path(out_dir, paste0(zfile, ".cs"))

  if (!file.exists(pip_path)) {
    return(list(error = "SparsePro produced no .pip output."))
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
    return(list(error = "Failed to parse SparsePro .pip output."))
  }

  # Variants present in the input but missing from .pip get NA.
  pip <- pip_df$pip[match(variant_ids, pip_df$rsid)]
  pip <- pmax(0, pmin(1, suppressWarnings(as.numeric(pip))))

  # --- Parse .cs file ---------------------------------------------------------
  # One row per effect group after a "cs  pip  beta" header. Within a cell the
  # upstream script joins values with "/" (sparsepro_zld.py:
  # '/'.join(...)); commas and whitespace are also accepted.

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
      split_cell <- function(s) {
        s <- trimws(s)
        if (!nzchar(s)) return(character(0))
        if (grepl("/", s, fixed = TRUE)) strsplit(s, "/", fixed = TRUE)[[1L]]
        else if (grepl(",", s, fixed = TRUE)) strsplit(s, ",\\s*")[[1L]]
        else strsplit(s, "\\s+")[[1L]]
      }

      parsed <- lapply(cs_lines, function(line) {
        # Split top-level on tabs (between the three columns).
        cells <- strsplit(line, "\t", fixed = TRUE)[[1L]]
        if (length(cells) < 1L) return(NULL)

        rsids <- split_cell(cells[1L])
        pips  <- if (length(cells) >= 2L) suppressWarnings(as.numeric(split_cell(cells[2L]))) else NA_real_
        effs  <- if (length(cells) >= 3L) suppressWarnings(as.numeric(split_cell(cells[3L]))) else NA_real_

        idx <- match(rsids, variant_ids)
        idx <- idx[!is.na(idx)]
        list(idx = sort(unique(idx)), pip = pips, effect = effs)
      })

      parsed <- Filter(function(x) !is.null(x) && length(x$idx) > 0L, parsed)

      if (length(parsed) > 0L) {
        credible_sets  <- lapply(parsed, `[[`, "idx")
        cs_pip         <- lapply(parsed, `[[`, "pip")
        cs_effect_size <- lapply(parsed, `[[`, "effect")
      }
    }
  }

  list(error = NULL, pip = pip, credible_sets = credible_sets,
       cs_pip = cs_pip, cs_effect_size = cs_effect_size)
}
