# =============================================================================
# wrapper_finemap_inf.R
#
# Wrapper for FINEMAP-inf (Cui et al., Nature Genetics 2023), which extends
# FINEMAP with an infinitesimal (polygenic) effect component alongside a small
# number of large sparse effects.
#
# WHY THIS METHOD IS HERE. Iteration 003 found that susie_inf trails plain susie
# in EVERY stratum, including the sparse+infinitesimal model it was built for -
# but susie_inf inherits SuSiE's LD fragility (LD explains 16% of its AUPRC
# variance, 56% of its reliability variance). FINEMAP-inf puts the SAME
# infinitesimal extension on the FINEMAP backbone, which is uniquely LD-robust
# in this benchmark (LD explains <1% of finemap's AUPRC variance). Running both
# separates "modelling the infinitesimal term is unhelpful" from "the SuSiE
# backbone is the problem" - a question the current results cannot answer.
#
# FINEMAP-inf is a Python program (FinucaneLab/fine-mapping-inf). CLI verified
# against run_fine_mapping.py's argparse definitions:
#   python run_fine_mapping.py \
#     --sumstats <file> --z-col-name <col> --ld-file <file> --n <int> \
#     --method susieinf,finemapinf --num-sparse-effects <int> --coverage <num> \
#     --purity <num> --save-tsv --output-prefix <prefix>
#
# NOTE on --method (verified in run_fine_mapping.py, 2026-07-30). The script's
# own default is 'susieinf,finemapinf' and it branches on
#   if 'susieinf' not in methods or args.est_finemapinf_tausq:
# so running FINEMAP-inf ALONE cold-starts it from sigmasq=1, tausq=0 and makes
# it estimate its own variance components, whereas the default path has it
# INHERIT SuSiE-inf's tau^2 / sigma^2. Both are supported, but only the latter is
# the configuration the method is distributed to run, so we use it and parse the
# .finemapinf output. This costs a SuSiE-inf fit per region.
#
# Output: <prefix>.finemapinf.bgz - the input sumstats plus columns
#   prob, post_mean_cond, post_sd_cond, alpha, post_mean, tausq, sigmasq
# where `prob` is the per-variant PIP. With --save-tsv a plain .tsv is written
# alongside; we read whichever exists (gzfile() handles the .bgz transparently
# for reading, since bgzip output is gzip-compatible).
#
# This file provides:
#   - setup_finemap_inf()      : checks the script and Python are reachable
#   - run_finemap_inf()        : runs FINEMAP-inf on a single region
#   - run_finemap_inf_region() : adapter called by run_methods()
#
# Reference:
#   Cui R, Elzur RA, Kanai M, et al. (2024). Improving fine-mapping by modeling
#   infinitesimal effects. Nature Genetics, 56, 162-169.
# =============================================================================


#' Check that the FINEMAP-inf script and its Python are available
#'
#' @param finemap_inf_dir Character. Directory containing run_fine_mapping.py.
#' @param python Character. Python executable (or a wrapper script).
#' @return Invisible TRUE if the script is found.
#' @export
setup_finemap_inf <- function(finemap_inf_dir, python = "python") {
  script <- file.path(finemap_inf_dir, "run_fine_mapping.py")
  if (!file.exists(script)) {
    stop(
      "run_fine_mapping.py not found at: ", script, "\n\n",
      "Install FINEMAP-inf:\n",
      "  git clone https://github.com/FinucaneLab/fine-mapping-inf.git\n",
      "  cd fine-mapping-inf && pip install -r requirements.txt\n\n",
      "Then pass the directory via\n",
      "  method_args = list(finemap_inf = list(finemap_inf_dir = '/path/to/fine-mapping-inf'))",
      call. = FALSE
    )
  }
  message("FINEMAP-inf script found: ", script)
  invisible(TRUE)
}


#' Run FINEMAP-inf on a single region
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param n Integer. GWAS sample size.
#' @param variant_ids Character vector or NULL.
#' @param finemap_inf_dir Character. Directory holding run_fine_mapping.py.
#' @param python Character. Python executable / wrapper script.
#' @param num_sparse_effects Integer. Maximum number of sparse large effects
#'   (\code{--num-sparse-effects}). Default 10, matching upstream.
#' @param coverage Numeric. Credible-set coverage. Default 0.95.
#' @param purity Numeric. Credible-set purity threshold. Default 0.5.
#' @param method Character. Which method's OUTPUT to return: "finemapinf"
#'   (default) or "susieinf". This is the file that gets parsed; it is not what
#'   is passed to \code{--method} (see \code{run_methods_str}).
#' @param run_methods_str Character. The value passed to \code{--method}.
#'   Defaults to \code{"susieinf,finemapinf"}, which is the SCRIPT'S OWN
#'   DEFAULT and the authors' intended pipeline: run_fine_mapping.py branches on
#'   \code{if 'susieinf' not in methods or args.est_finemapinf_tausq}, so with
#'   FINEMAP-inf alone it cold-starts from \code{sigmasq=1, tausq=0} and
#'   estimates its own variance components, whereas the default path has
#'   FINEMAP-inf INHERIT SuSiE-inf's tau^2 and sigma^2. Running FINEMAP-inf on
#'   its own is supported but is a different configuration from the one the
#'   method is distributed to run.
#' @param est_tausq Logical. Force FINEMAP-inf to estimate tau^2 itself even
#'   when SuSiE-inf is in the run (\code{--est-finemapinf-tausq}). Default FALSE.
#'
#' @return A list in the benchmark's standard fine-mapping output format.
#' @export
run_finemap_inf <- function(z,
                            LD,
                            n,
                            variant_ids        = NULL,
                            finemap_inf_dir,
                            python             = "python",
                            num_sparse_effects = 10L,
                            coverage           = 0.95,
                            purity             = 0.5,
                            method             = "finemapinf",
                            run_methods_str    = "susieinf,finemapinf",
                            est_tausq          = FALSE) {

  p <- length(z)
  params <- list(num_sparse_effects = num_sparse_effects, coverage = coverage,
                 purity = purity, method = method, run_methods_str = run_methods_str,
                 est_tausq = est_tausq,
                 n = n, finemap_inf_dir = finemap_inf_dir, python = python)

  stopifnot("LD must be a p x p matrix" = is.matrix(LD) && nrow(LD) == p && ncol(LD) == p,
            "n must be a positive integer" = is.numeric(n) && length(n) == 1 && n > 0)

  script <- file.path(finemap_inf_dir, "run_fine_mapping.py")
  if (!file.exists(script)) {
    return(.fminf_error_result(p, params, 0,
      paste("run_fine_mapping.py not found at:", script)))
  }

  if (is.null(variant_ids)) variant_ids <- paste0("SNP_", seq_len(p))
  variant_ids <- gsub("\\s+", "_", variant_ids)

  work_dir <- tempfile(pattern = "fminf_run_")
  dir.create(work_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  ss_path <- file.path(work_dir, "sumstats.tsv")
  # VERIFIED 2026-07-30: run_fine_mapping.py's read_large_file() dispatches on the
  # FILE EXTENSION and accepts only .npy / .npz / .gz / .bgz - a plain-text matrix
  # (e.g. "region.ld") is rejected with
  #   ValueError: File extension .ld of file region.ld currently not supported
  # .gz is the one we can write directly from R (numpy reads it via np.loadtxt).
  ld_path <- file.path(work_dir, "region.ld.gz")
  prefix  <- file.path(work_dir, "out")

  # --sumstats is a table with a named z column (--z-col-name).
  utils::write.table(data.frame(SNP = variant_ids, Z = as.numeric(z)),
                     ss_path, quote = FALSE, row.names = FALSE,
                     col.names = TRUE, sep = "\t")
  # --ld-file: full numeric matrix, gzip-compressed (see the note above).
  ld_con <- gzfile(ld_path, "w")
  utils::write.table(round(LD, 8), ld_con, quote = FALSE,
                     row.names = FALSE, col.names = FALSE, sep = " ")
  close(ld_con)

  args <- c(script,
            "--sumstats",           ss_path,
            "--z-col-name",         "Z",
            "--ld-file",            ld_path,
            "--n",                  as.character(as.integer(n)),
            "--method",             run_methods_str,
            "--num-sparse-effects", as.character(as.integer(num_sparse_effects)),
            "--coverage",           as.character(coverage),
            "--purity",             as.character(purity),
            "--save-tsv",
            "--output-prefix",      prefix)
  if (isTRUE(est_tausq)) args <- c(args, "--est-finemapinf-tausq")

  start_time <- proc.time()
  run_output <- tryCatch(system2(python, args = args, stdout = TRUE, stderr = TRUE),
                         error = function(e) structure(conditionMessage(e), class = "fminf_error"))
  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])
  if (inherits(run_output, "fminf_error")) {
    return(.fminf_error_result(p, params, elapsed, as.character(run_output)))
  }

  # Upstream writes <prefix>.<method>.bgz; --save-tsv also yields a plain .tsv.
  cands <- c(paste0(prefix, ".", method, ".tsv"),
             paste0(prefix, ".", method, ".bgz"),
             paste0(prefix, ".", method, ".tsv.gz"))
  out_file <- cands[file.exists(cands)][1]
  if (is.na(out_file)) {
    return(.fminf_error_result(p, params, elapsed,
      paste(c("FINEMAP-inf produced no output table. Looked for:",
              cands, run_output), collapse = "\n")))
  }

  df <- tryCatch({
    con <- if (grepl("\\.(bgz|gz)$", out_file)) gzfile(out_file) else out_file
    utils::read.table(con, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, comment.char = "")
  }, error = function(e) NULL)

  if (is.null(df) || !"prob" %in% names(df)) {
    return(.fminf_error_result(p, params, elapsed,
      paste("FINEMAP-inf output lacks a 'prob' column; columns:",
            if (is.null(df)) "<unreadable>" else paste(names(df), collapse = ", "))))
  }

  # Re-order to the input variant order when the id column survives.
  ord <- if ("SNP" %in% names(df)) {
    o <- match(variant_ids, as.character(df$SNP)); if (anyNA(o)) seq_len(p) else o
  } else seq_len(p)
  pip <- pmax(0, pmin(1, suppressWarnings(as.numeric(df$prob[ord]))))
  if (length(pip) != p) {
    return(.fminf_error_result(p, params, elapsed,
      sprintf("FINEMAP-inf returned %d probabilities, expected %d.", length(pip), p)))
  }

  # SuSiE-inf output carries a `cs` column; FINEMAP-inf's table does not, so a
  # coverage set is built greedily from the PIPs (as for PAINTOR).
  credible_sets <- list()
  if ("cs" %in% names(df)) {
    csv <- df$cs[ord]
    for (lab in setdiff(unique(csv), c(NA, -1, "", "-1", "NA"))) {
      idx <- which(csv == lab)
      if (length(idx)) credible_sets[[length(credible_sets) + 1L]] <- sort(as.integer(idx))
    }
  } else if (any(is.finite(pip))) {
    op <- order(pip, decreasing = TRUE)
    n_cs <- which(cumsum(pip[op]) >= coverage)[1L]
    if (is.na(n_cs)) n_cs <- p
    credible_sets <- list(sort(op[seq_len(n_cs)]))
  }

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = if (identical(method, "susieinf")) "susie_inf_py" else "finemap_inf",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      tausq    = if ("tausq"    %in% names(df)) df$tausq[1]    else NA_real_,
      sigmasq  = if ("sigmasq"  %in% names(df)) df$sigmasq[1]  else NA_real_,
      post_mean = if ("post_mean" %in% names(df)) as.numeric(df$post_mean[ord]) else NULL
    )
  )
}


#' Run FINEMAP-inf on a single region from simulation data structures
#'
#' @param region_geno One element of \code{simulation$genotypes}.
#' @param region_pheno One element of a scenario's \code{regions}.
#' @param ... Passed to \code{\link{run_finemap_inf}}.
#' @return The output of \code{\link{run_finemap_inf}}.
#' @export
run_finemap_inf_region <- function(region_geno, region_pheno, ...) {
  run_finemap_inf(
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

.fminf_error_result <- function(p, params, elapsed, error_msg) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "finemap_inf",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(tausq = NA_real_, sigmasq = NA_real_, post_mean = NULL),
    error           = error_msg
  )
}
