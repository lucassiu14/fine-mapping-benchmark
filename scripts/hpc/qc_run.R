#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/qc_run.R
#
# Quality control for a benchmark run, to be used WHILE it is still going or
# once it has finished. Answers three questions in order of importance:
#
#   1. Is any method silently producing nothing?
#      This is the failure mode this project has actually suffered. CARMA was
#      registered from Iteration 001 but never installed, so it returned 100%
#      NA for three iterations without anyone noticing, and PAINTOR ran
#      mis-configured for the same three. A method that errors loudly is not
#      the danger; a method that returns NA fits forever is.
#
#   2. How far through is the run, and is progress even across the grid?
#      A row that is far behind usually means its tasks are dying, not that
#      they are slow.
#
#   3. Did any scenario write a results.rds but no evaluation.rds?
#      That means evaluate_methods() failed after the expensive part
#      succeeded - recoverable with collect_results.R, but it needs knowing.
#
# Usage (from the project root):
#   Rscript scripts/hpc/qc_run.R $EPHEMERAL/fmbench_iter004/results/benchmark
#   Rscript scripts/hpc/qc_run.R <root> --sample 300     # fast partial check
#
# Reading every results.rds is I/O heavy, so --sample reads a random subset for
# the per-method statistics. Completion counting always covers everything.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
if (is.na(root) || !nzchar(root) || !dir.exists(root)) {
  stop("usage: Rscript scripts/hpc/qc_run.R <output_root> [--sample N]", call. = FALSE)
}
sample_n <- NA_integer_
if ("--sample" %in% args) {
  sample_n <- suppressWarnings(as.integer(args[which(args == "--sample") + 1L]))
}
`%||%` <- function(x, y) if (is.null(x)) y else x

# --- 1. completion ----------------------------------------------------------
job_dirs <- list.dirs(root, recursive = FALSE)
job_dirs <- job_dirs[grepl("/job_[0-9]+", job_dirs)]
if (!length(job_dirs)) stop("no job_* directories under ", root, call. = FALSE)

eval_files <- list.files(root, pattern = "^evaluation\\.rds$",
                         recursive = TRUE, full.names = TRUE)
res_files  <- list.files(root, pattern = "^results\\.rds$",
                         recursive = TRUE, full.names = TRUE)

cat("=============================================================\n")
cat(sprintf("QC: %s\n", root))
cat("=============================================================\n\n")
cat(sprintf("rows (job_*) present : %d\n", length(job_dirs)))
cat(sprintf("scenarios evaluated  : %d\n", length(eval_files)))
cat(sprintf("scenarios with results.rds : %d\n", length(res_files)))

orphan <- setdiff(dirname(res_files), dirname(eval_files))
if (length(orphan)) {
  cat(sprintf("\n!! %d scenario(s) have results.rds but NO evaluation.rds.\n",
              length(orphan)))
  cat("   evaluate_methods() failed after the expensive fitting succeeded.\n")
  cat("   collect_results.R can re-evaluate these; the compute is not lost.\n")
  cat(sprintf("   e.g. %s\n", orphan[1]))
}

# --- 2. progress per row ----------------------------------------------------
per_row <- vapply(job_dirs, function(d)
  length(list.files(d, pattern = "^evaluation\\.rds$", recursive = TRUE)),
  integer(1))
names(per_row) <- basename(job_dirs)
cat(sprintf("\nscenarios per row: min %d  median %d  max %d\n",
            min(per_row), as.integer(median(per_row)), max(per_row)))
lag <- sort(per_row)[seq_len(min(5L, length(per_row)))]
if (length(unique(per_row)) > 1L) {
  cat("slowest rows (a row far behind is usually dying, not slow):\n")
  for (i in seq_along(lag)) cat(sprintf("  %-46s %d\n", names(lag)[i], lag[i]))
}

# --- 3. per-method fit success ---------------------------------------------
if (!length(res_files)) {
  cat("\nNo results.rds yet - nothing to check per method.\n")
  quit(save = "no", status = 0)
}
pick <- res_files
if (!is.na(sample_n) && sample_n > 0 && sample_n < length(res_files)) {
  set.seed(1); pick <- sample(res_files, sample_n)
  cat(sprintf("\nper-method check on a random sample of %d/%d scenarios\n",
              length(pick), length(res_files)))
} else {
  cat(sprintf("\nper-method check across all %d scenarios\n", length(pick)))
}

acc     <- new.env(parent = emptyenv())
skipped <- new.env(parent = emptyenv())
n_read  <- 0L
for (f in pick) {
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(r)) {
    acc[["__unreadable__"]] <- (acc[["__unreadable__"]] %||% 0L) + 1L
    next
  }
  n_read <- n_read + 1L
  for (m in names(r)) {
    v <- r[[m]]
    # results.rds can carry top-level entries that are not per-method lists
    # (metadata, timings). Indexing those with $ errors on atomic vectors, so
    # only count entries that actually look like a method result.
    if (!is.list(v) || (is.null(v$n_total) && is.null(v$n_failed))) {
      skipped[[m]] <- (skipped[[m]] %||% 0L) + 1L
      next
    }
    cur <- acc[[m]] %||% c(scen = 0, fits = 0, failed = 0)
    acc[[m]] <- c(scen   = cur[["scen"]]   + 1,
                  fits   = cur[["fits"]]   + (v$n_total  %||% 0L),
                  failed = cur[["failed"]] + (v$n_failed %||% 0L))
  }
}
if (!is.null(acc[["__unreadable__"]])) {
  cat(sprintf("\n!! %d results.rds file(s) could not be read (truncated write?)\n",
              acc[["__unreadable__"]]))
}

sk <- setdiff(ls(skipped), character(0))
if (length(sk)) {
  cat(sprintf("\nnon-method top-level entries in results.rds (ignored): %s\n",
              paste(sk, collapse = ", ")))
}

meths <- setdiff(ls(acc), "__unreadable__")
tab <- do.call(rbind, lapply(meths, function(m) {
  a <- acc[[m]]
  data.frame(method = m, scenarios = a[["scen"]], fits = a[["fits"]],
             failed = a[["failed"]],
             pct_failed = if (a[["fits"]] > 0) 100 * a[["failed"]] / a[["fits"]] else NA_real_,
             stringsAsFactors = FALSE)
}))
tab <- tab[order(-tab$pct_failed, tab$method), ]

cat("\n")
cat(sprintf("%-22s %10s %10s %8s %9s\n",
            "method", "scenarios", "fits", "failed", "%failed"))
cat(strrep("-", 63), "\n")
for (i in seq_len(nrow(tab))) with(tab[i, ],
  cat(sprintf("%-22s %10d %10d %8d %8.1f%%\n",
              method, scenarios, fits, failed, pct_failed)))

# --- verdict ---------------------------------------------------------------
cat("\n--- verdict ---------------------------------------------------\n")
expected <- c("susie","susie_inf","finemap","finemap_inf","abf","marginal_z",
              "carma","finimom","sparsepro","funmap","paintor","polyfun_oracle",
              "polyfun_est","polyfun_ldsc","sbayesrc","beatrice",
              "functional_beatrice","fb_pooled","fb_xregion")
missing <- setdiff(expected, tab$method)
if (length(missing)) {
  cat(sprintf("MISSING ENTIRELY (%d): %s\n", length(missing),
              paste(missing, collapse = ", ")))
  cat("  These produced no entry at all - not registered, or the wrapper\n")
  cat("  errored before returning. Investigate before trusting the run.\n")
}
dead <- tab$method[!is.na(tab$pct_failed) & tab$pct_failed >= 99.9]
if (length(dead)) {
  cat(sprintf("100%% FAILED (%d): %s\n", length(dead), paste(dead, collapse = ", ")))
  cat("  This is the CARMA failure mode: registered, running, returning\n")
  cat("  nothing usable. Fix the wrapper and re-run with FMB_METHODS.\n")
}
degraded <- tab$method[!is.na(tab$pct_failed) & tab$pct_failed > 5 & tab$pct_failed < 99.9]
if (length(degraded)) {
  cat(sprintf("PARTIAL FAILURES >5%% (%d): %s\n", length(degraded),
              paste(degraded, collapse = ", ")))
  cat("  Usually size-dependent - check whether the failures concentrate at\n")
  cat("  the largest regions (p = 2000) before treating the numbers as valid.\n")
}
if (!length(missing) && !length(dead) && !length(degraded)) {
  cat("All expected methods present, none above a 5% failure rate.\n")
}
cat(sprintf("\nCompletion: %d scenarios evaluated.\n", length(eval_files)))
