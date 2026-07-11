# =============================================================================
# finemap.R
#
# Wrapper for FINEMAP (Benner et al. 2016) fine-mapping.
#
# FINEMAP is an external C++ binary. This wrapper handles writing the required
# input files, calling the binary, and parsing the output back into the
# standardised result format.
#
# This file provides:
#   - setup_finemap()        : verifies the FINEMAP binary is accessible
#   - run_finemap()          : runs FINEMAP on a single region (explicit inputs)
#   - run_finemap_region()   : adapter called by run_methods(); extracts inputs
#                              from simulation data structures
#
# Standard output format:
#   pip              Numeric vector (length p). Marginal PIPs.
#   credible_sets    List of integer vectors (variant indices). Credible sets
#                    from the most probable number of causal variants k*.
#                    Empty list if k* = 0 or no credible set file found.
#   method           Character. "finemap".
#   input_type       Character. Always "summary".
#   params           List. Hyperparameters used.
#   runtime_seconds  Numeric. Wall-clock time.
#   additional       List. FINEMAP-specific outputs (see run_finemap() docs).
# =============================================================================


# =============================================================================
# Setup
# =============================================================================

#' Set up the FINEMAP binary
#'
#' Checks that the FINEMAP binary can be found at the supplied path or on the
#' system PATH. If not found and \code{download = TRUE} (the default), the
#' binary is downloaded automatically for your OS and cached locally — no
#' manual installation required.
#'
#' The downloaded binary is saved to your R user cache directory
#' (see \code{tools::R_user_dir("fmbenchmark", "cache")}). The resolved path
#' is returned invisibly, so you can capture it for use in \code{run_methods()}:
#'
#' \preformatted{
#'   fp  <- setup_finemap()
#'   out <- run_methods(sim, methods = "finemap",
#'                      method_args = list(finemap = list(finemap_path = fp)))
#' }
#'
#' @param finemap_path Character. Path to the FINEMAP binary, or \code{"finemap"}
#'   to search the system PATH. Default: \code{"finemap"}.
#' @param download Logical. If \code{TRUE} (default) and the binary is not
#'   found, attempt to download it automatically. Set to \code{FALSE} to
#'   disable downloading and receive manual installation instructions instead.
#'
#' @return Invisibly returns the resolved path to the binary.
#' @export
setup_finemap <- function(finemap_path = "finemap", download = TRUE) {

  # Full path supplied: check directly
  if (file.exists(finemap_path)) {
    message("FINEMAP binary found: ", normalizePath(finemap_path))
    return(invisible(normalizePath(finemap_path)))
  }

  # Search PATH
  resolved <- Sys.which(finemap_path)
  if (nchar(resolved) > 0) {
    message("FINEMAP binary found on PATH: ", resolved)
    return(invisible(resolved))
  }

  # Not found — try to download, or explain how to install manually
  if (download) {
    path <- .download_finemap()
    return(invisible(path))
  }

  sysname <- Sys.info()["sysname"]
  os_hint <- if (.Platform$OS.type == "windows") {
    "  FINEMAP does not provide a Windows binary. Consider using WSL (Linux)."
  } else if (sysname == "Darwin") {
    "  chmod +x finemap_v1.4.1_MacOSX && mv finemap_v1.4.1_MacOSX /usr/local/bin/finemap"
  } else {
    "  chmod +x finemap_v1.4.1_x86_64 && mv finemap_v1.4.1_x86_64 /usr/local/bin/finemap"
  }

  stop(
    "FINEMAP binary not found (searched for: '", finemap_path, "').\n\n",
    "To download automatically, run:\n",
    "  setup_finemap()   # download = TRUE is the default\n\n",
    "To install manually:\n",
    "  1. Go to: http://www.christianbenner.com\n",
    "  2. Download the binary for your OS.\n",
    "  3. ", os_hint, "\n",
    call. = FALSE
  )
}


#' Download the FINEMAP binary for the current OS
#'
#' Downloads a pre-compiled FINEMAP v1.4.1 binary from the echofinemap GitHub
#' releases page, saves it to the R user cache directory, makes it executable,
#' and returns the path.
#'
#' Called automatically by \code{\link{setup_finemap}} when the binary is not
#' found. Can also be called directly to force a fresh download.
#'
#' @return Character. Path to the downloaded binary.
#' @export
download_finemap <- function() {
  .download_finemap()
}


# Internal implementation
.download_finemap <- function() {

  sysname <- Sys.info()["sysname"]

  if (.Platform$OS.type == "windows") {
    stop(
      "Automatic FINEMAP download is not supported on Windows as no Windows\n",
      "binary is available. Consider using WSL (Linux subsystem).",
      call. = FALSE
    )
  }

  fname <- if (sysname == "Darwin") {
    "finemap_v1.4.1_MacOSX.tgz"
  } else {
    "finemap_v1.4.1_x86_64.tgz"
  }

  url <- paste0(
    "https://github.com/RajLabMSSM/echofinemap/releases/download/latest/",
    fname
  )

  # Cache directory
  cache_dir <- tools::R_user_dir("fmbenchmark", which = "cache")
  save_dir  <- file.path(cache_dir, "finemap")
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  tgz_path <- file.path(save_dir, fname)

  message("Downloading FINEMAP binary from:\n  ", url)
  tryCatch(
    utils::download.file(url, destfile = tgz_path, mode = "wb", quiet = FALSE),
    error = function(e) {
      stop(
        "Download failed: ", conditionMessage(e), "\n",
        "Check your internet connection or download manually from:\n",
        "  http://www.christianbenner.com",
        call. = FALSE
      )
    }
  )

  # Extract
  utils::untar(tgz_path, exdir = save_dir)

  # Locate the extracted binary (strip .tgz suffix)
  binary_name <- sub("\\.tgz$", "", fname)
  binary_path <- file.path(save_dir, binary_name, binary_name)

  # Fallback: search for any file named like the binary in save_dir
  if (!file.exists(binary_path)) {
    candidates <- list.files(save_dir, pattern = binary_name,
                             full.names = TRUE, recursive = TRUE)
    candidates <- candidates[!grepl("\\.tgz$", candidates)]
    if (length(candidates) == 0) {
      stop(
        "Downloaded and extracted the archive but could not locate the binary.\n",
        "Archive saved at: ", tgz_path,
        call. = FALSE
      )
    }
    binary_path <- candidates[[1]]
  }

  # Make executable
  Sys.chmod(binary_path, "0755")

  if (sysname == "Darwin") {
    arch <- system("uname -m", intern = TRUE)
    if (identical(trimws(arch), "arm64")) {
      message(
        "Note: you are on Apple Silicon (arm64). The downloaded binary is\n",
        "x86_64 and will run via Rosetta 2. Two prerequisites:\n",
        "  1. Rosetta must be installed: softwareupdate --install-rosetta\n",
        "  2. libzstd must be at /usr/local/lib/libzstd.1.dylib (x86_64 build).\n",
        "     If missing, install x86_64 Homebrew and run:\n",
        "       arch -x86_64 /bin/bash -c \"$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\"\n",
        "       arch -x86_64 /usr/local/bin/brew install zstd"
      )
    }
  }

  message("FINEMAP binary ready at:\n  ", binary_path)
  message(
    "\nTo use it, pass the path to run_methods():\n",
    "  run_methods(sim, methods = \"finemap\",\n",
    "              method_args = list(finemap = list(finemap_path = \"",
    binary_path, "\")))"
  )

  binary_path
}


# =============================================================================
# Run FINEMAP on a single region
# =============================================================================

#' Run FINEMAP fine-mapping on a single region
#'
#' Writes the required FINEMAP input files to a temporary directory, calls the
#' FINEMAP binary, parses the output, and returns results in the standardised
#' format. The temporary directory is deleted on exit.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param n Integer. Sample size.
#' @param beta Numeric vector or NULL. Marginal effect size estimates (length
#'   p). If NULL, derived as \code{z * se}. When both \code{beta} and
#'   \code{se} are NULL, \code{se} is set to 1 and \code{beta = z}.
#' @param se Numeric vector or NULL. Standard errors (length p). See above.
#' @param maf Numeric vector or NULL. Minor allele frequencies (length p).
#'   Used in the FINEMAP input file. If NULL, set to 0.3 for all variants.
#' @param variant_ids Character vector or NULL. Variant identifiers (length p).
#'   Used as rsids in FINEMAP input/output. If NULL, set to
#'   \code{"SNP_1", "SNP_2", ...}.
#' @param finemap_path Character. Path to the FINEMAP binary or name on PATH.
#'   Default: \code{"finemap"}.
#' @param n_causal Integer. Maximum number of causal variants considered.
#'   Default: 5.
#' @param n_iter Integer. Number of stochastic shotgun search iterations.
#'   Default: 100000.
#' @param prior_std Numeric. Prior standard deviation for effect sizes.
#'   Default: 0.05.
#' @param coverage Numeric. Coverage level for credible sets. Default: 0.95.
#'
#' @return A list with the standardised fine-mapping output:
#' \describe{
#'   \item{pip}{Numeric vector (length p). Marginal posterior inclusion
#'     probabilities.}
#'   \item{credible_sets}{List of integer vectors. Credible sets derived from
#'     the most probable number of causal variants k*. Each element is the
#'     set of variant indices for one causal signal. Empty list if k* = 0.}
#'   \item{method}{Character. Always \code{"finemap"}.}
#'   \item{input_type}{Character. Always \code{"summary"}.}
#'   \item{params}{List. Hyperparameters used.}
#'   \item{runtime_seconds}{Numeric. Wall-clock time in seconds.}
#'   \item{additional}{List of FINEMAP-specific outputs:
#'     \describe{
#'       \item{log10bf}{Numeric vector (length p). Log10 Bayes factor per
#'         variant.}
#'       \item{posterior_mean}{Numeric vector (length p). Posterior mean
#'         effect size.}
#'       \item{posterior_sd}{Numeric vector (length p). Posterior SD of
#'         effect size.}
#'       \item{k_posterior}{Named numeric vector. Posterior probability that
#'         exactly k variants are causal (names are k = 0, 1, 2, ...).}
#'       \item{best_k}{Integer. Most probable number of causal variants.}
#'       \item{configs}{Data frame. Top causal configurations with their
#'         posterior probabilities, as returned in FINEMAP's .config file.}
#'     }
#'   }
#'   \item{error}{Character or NULL. Error message if FINEMAP failed.}
#' }
#'
#' @export
run_finemap <- function(z,
                        LD,
                        n,
                        beta         = NULL,
                        se           = NULL,
                        maf          = NULL,
                        variant_ids  = NULL,
                        finemap_path = "finemap",
                        n_causal     = 5,
                        n_iter       = 100000,
                        prior_std    = 0.05,
                        coverage     = 0.95) {

  # Expand ~ so system2 can exec the binary. Users often pass "~/tools/..."
  # and R's default tilde-expansion doesn't apply in system2's first arg.
  finemap_path <- path.expand(finemap_path)

  # --- Validate inputs --------------------------------------------------------

  p <- length(z)

  stopifnot(
    "LD must be a p x p matrix" = is.matrix(LD) && nrow(LD) == p && ncol(LD) == p,
    "n must be a positive integer" = is.numeric(n) && length(n) == 1 && n > 0
  )

  if (is.null(variant_ids)) variant_ids <- paste0("SNP_", seq_len(p))
  if (is.null(maf))         maf         <- rep(0.3, p)

  # Derive beta/se if not supplied
  if (is.null(se) && is.null(beta)) {
    se   <- rep(1.0, p)
    beta <- z
  } else if (is.null(beta)) {
    beta <- z * se
  } else if (is.null(se)) {
    # se = beta / z; guard against z = 0
    se <- ifelse(abs(z) > 1e-10, abs(beta / z), 1.0)
  }

  params <- list(
    n_causal     = n_causal,
    n_iter       = n_iter,
    prior_std    = prior_std,
    coverage     = coverage,
    finemap_path = finemap_path
  )

  # --- Set up temp working directory ------------------------------------------

  work_dir <- tempfile(pattern = "finemap_run_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  prefix <- file.path(work_dir, "fm")

  # --- Write .z file ----------------------------------------------------------

  z_file <- paste0(prefix, ".z")
  z_df <- data.frame(
    rsid       = variant_ids,
    chromosome = 1L,
    position   = seq_len(p),
    allele1    = "A",
    allele2    = "T",
    maf        = round(maf, 6),
    beta       = beta,
    se         = se,
    stringsAsFactors = FALSE
  )
  write.table(z_df, z_file,
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = " ")

  # --- Write .ld file ---------------------------------------------------------

  ld_file <- paste0(prefix, ".ld")
  write.table(round(LD, 6), ld_file,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

  # --- Write .master file -----------------------------------------------------

  snp_file    <- paste0(prefix, ".snp")
  config_file <- paste0(prefix, ".config")
  cred_base   <- paste0(prefix, ".cred")   # FINEMAP appends the k (e.g. .cred1)
  # File suffix must equal the master-file header column name ("log"),
  # else FINEMAP v1.4.2 errors with 'Extension X of file Y ... cannot be
  # found in the header of the master file'. Previously ".log_sss" was
  # used and no longer accepted.
  log_file    <- paste0(prefix, ".log")

  master_file <- file.path(work_dir, "master")
  writeLines(
    c(
      "z;ld;snp;config;cred;log;n_samples",
      paste(z_file, ld_file, snp_file, config_file, cred_base, log_file, n,
            sep = ";")
    ),
    master_file
  )

  # --- Run FINEMAP ------------------------------------------------------------

  # FINEMAP v1.4.2 CLI: --n-iterations was renamed to --n-iter, and
  # --credible-config-value was replaced by --prob-cred-set. --log directs
  # FINEMAP to write to the log path specified in the master's 'log' column.
  args <- c(
    "--sss",
    "--in-files",       master_file,
    "--n-causal-snps",  as.character(n_causal),
    "--n-iter",         as.character(n_iter),
    "--prior-std",      as.character(prior_std),
    "--prob-cred-set",  as.character(coverage),
    "--log"
  )

  start_time <- proc.time()

  run_output <- tryCatch({
    system2(finemap_path, args = args, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    structure(conditionMessage(e), class = "finemap_error")
  })

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  # Check for binary-level failure
  if (inherits(run_output, "finemap_error")) {
    return(.finemap_error_result(p, params, elapsed, run_output))
  }

  # Check that the primary output file was produced
  if (!file.exists(snp_file)) {
    err_msg <- paste(
      c("FINEMAP produced no output.", run_output),
      collapse = "\n"
    )
    return(.finemap_error_result(p, params, elapsed, err_msg))
  }

  # --- Parse .snp file --------------------------------------------------------

  snp_df <- tryCatch(
    utils::read.table(snp_file, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  if (is.null(snp_df)) {
    return(.finemap_error_result(p, params, elapsed,
                                 "Failed to parse FINEMAP .snp output file."))
  }

  # Re-order rows to match the original variant order
  ord <- match(variant_ids, snp_df$rsid)

  pip            <- snp_df$prob[ord]
  posterior_mean <- snp_df$mean[ord]
  posterior_sd   <- snp_df$sd[ord]
  log10bf        <- snp_df$log10bf[ord]

  # --- Parse .config file to get posterior over k -----------------------------

  config_df  <- tryCatch(
    utils::read.table(config_file, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  k_posterior <- NULL
  best_k      <- 0L   # default when config_df is missing/unparseable

  if (!is.null(config_df) && "config" %in% names(config_df) &&
      "prob" %in% names(config_df)) {

    # Prefer n_snps column (present in FINEMAP v1.4.1) — more reliable than
    # parsing the config string, which uses "0" for the null model (0 causal
    # variants) and would be misclassified as k=1 by a comma-split approach.
    if ("n_snps" %in% names(config_df)) {
      config_df$k <- as.integer(config_df$n_snps)
    } else {
      config_df$k <- vapply(config_df$config, function(cfg) {
        cfg <- trimws(cfg)
        if (is.na(cfg) || nchar(cfg) == 0 || cfg == "0") 0L
        else length(strsplit(cfg, ",")[[1L]])
      }, integer(1L))
    }

    k_probs     <- tapply(config_df$prob, config_df$k, sum)
    k_posterior <- setNames(as.numeric(k_probs), names(k_probs))
    best_k      <- as.integer(names(which.max(k_probs)))
  }

  # --- Parse .cred<k*> file for credible sets ---------------------------------

  credible_sets <- list()

  if (best_k > 0L) {
    cred_k_file <- paste0(cred_base, best_k)

    if (file.exists(cred_k_file)) {
      cred_df <- tryCatch(
        utils::read.table(cred_k_file, header = TRUE, stringsAsFactors = FALSE,
                          fill = TRUE, na.strings = c("", "NA")),
        error = function(e) NULL
      )

      if (!is.null(cred_df) && ncol(cred_df) > 0L) {
        credible_sets <- lapply(seq_len(ncol(cred_df)), function(j) {
          rsids <- cred_df[[j]]
          rsids <- rsids[!is.na(rsids) & nchar(trimws(rsids)) > 0L]
          idx   <- match(rsids, variant_ids)
          sort(idx[!is.na(idx)])
        })
        credible_sets <- credible_sets[lengths(credible_sets) > 0L]
      }
    }
  }

  # --- Return -----------------------------------------------------------------

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = "finemap",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      log10bf        = log10bf,
      posterior_mean = posterior_mean,
      posterior_sd   = posterior_sd,
      k_posterior    = k_posterior,
      best_k         = best_k,
      configs        = config_df
    )
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run FINEMAP on a single region from simulation data structures
#'
#' Thin adapter that extracts the appropriate inputs from the simulation's
#' \code{region_geno} and \code{region_pheno} objects and calls
#' \code{\link{run_finemap}}. This is the function registered in the method
#' registry and called by \code{\link{run_methods}}.
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing \code{LD}, \code{n}, and optionally \code{maf} and
#'   \code{variant_ids}.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z}, \code{beta_hat}, and \code{se}.
#' @param ... Additional arguments passed to \code{\link{run_finemap}}
#'   (e.g. \code{n_causal}, \code{prior_std}, \code{finemap_path}).
#'
#' @return The output of \code{\link{run_finemap}}.
#' @export
run_finemap_region <- function(region_geno, region_pheno, ...) {
  run_finemap(
    z           = region_pheno$z,
    LD          = region_geno$LD,
    n           = region_geno$n,
    beta        = region_pheno$beta_hat,
    se          = region_pheno$se,
    maf         = region_geno$maf,
    variant_ids = region_geno$variant_ids,
    ...
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

# Construct a standardised error result when FINEMAP fails, so the batch
# runner can continue rather than abort.
.finemap_error_result <- function(p, params, elapsed, error_msg) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "finemap",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      log10bf        = rep(NA_real_, p),
      posterior_mean = rep(NA_real_, p),
      posterior_sd   = rep(NA_real_, p),
      k_posterior    = NULL,
      best_k         = NA_integer_,
      configs        = NULL
    ),
    error = error_msg
  )
}
