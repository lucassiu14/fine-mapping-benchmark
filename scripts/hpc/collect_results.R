# =============================================================================
# scripts/hpc/collect_results.R
#
# Collect and combine results from all completed benchmark jobs into a single
# analysis-ready data frame.
#
# Usage (from project root, after all jobs complete):
#   Rscript scripts/hpc/collect_results.R
#
# Output:
#   results/benchmark/combined_evaluation.rds  — long data frame, one row per
#                                                 (job, method, stratum, value)
#   results/benchmark/combined_evaluation.csv  — same, CSV form
#   results/benchmark/run_summary.csv          — per-job completion status
# =============================================================================

PROJECT_ROOT <- "."
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "results", "benchmark")
GRID_PATH    <- file.path(PROJECT_ROOT, "scripts", "hpc", "params_grid.csv")

# =============================================================================
# Internal: flatten one evaluate_methods() output to a long data frame.
#
# evaluate_methods() returns a nested list:
#   eval_out$<method>$global       (named list of scalars + se's)
#   eval_out$<method>$by_S         (named list, one entry per S value)
#   eval_out$<method>$by_phi       (named list, one entry per phi value)
#   eval_out$<method>$by_p_causal  (named list, sparse_inf only) or NULL
# Each stratum contains scalar metrics + matching <metric>_se fields.
# =============================================================================

.SCALAR_METRICS <- c("auprc", "cs_coverage", "cs_power",
                     "cs_size_median", "cs_size_mean", "runtime_mean")

.flatten_stratum <- function(stratum_obj, method, stratum, stratum_val) {
  if (is.null(stratum_obj)) return(NULL)
  rows <- lapply(.SCALAR_METRICS, function(metric) {
    val <- stratum_obj[[metric]]
    se  <- stratum_obj[[paste0(metric, "_se")]]
    data.frame(
      method      = method,
      stratum     = stratum,
      stratum_val = stratum_val,
      metric      = metric,
      value       = if (is.null(val) || length(val) == 0L) NA_real_ else as.numeric(val[[1L]]),
      se          = if (is.null(se)  || length(se)  == 0L) NA_real_ else as.numeric(se[[1L]]),
      n_fits      = stratum_obj$n_fits   %||% NA_integer_,
      n_failed    = stratum_obj$n_failed %||% NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.flatten_eval <- function(eval_obj) {
  methods <- eval_obj$methods_evaluated
  if (is.null(methods)) methods <- setdiff(names(eval_obj),
                                            c("methods_evaluated",
                                              "simulation_params",
                                              "pip_thresholds_used"))
  per_method <- lapply(methods, function(m) {
    rows <- list()
    rows[[length(rows) + 1L]] <- .flatten_stratum(eval_obj[[m]]$global,
                                                    m, "global", NA_real_)
    for (strat in c("by_S", "by_phi", "by_p_causal")) {
      strat_obj <- eval_obj[[m]][[strat]]
      if (is.null(strat_obj)) next
      for (key in names(strat_obj)) {
        rows[[length(rows) + 1L]] <- .flatten_stratum(
          strat_obj[[key]], m, strat, suppressWarnings(as.numeric(key))
        )
      }
    }
    do.call(rbind, rows)
  })
  do.call(rbind, per_method)
}

`%||%` <- function(x, y) if (!is.null(x)) x else y

# =============================================================================
# Scan completed jobs
# =============================================================================

grid <- read.csv(GRID_PATH, stringsAsFactors = FALSE)
cat(sprintf("Scanning %d jobs...\n", nrow(grid)))

summary_rows <- vector("list", nrow(grid))
eval_list    <- list()

for (i in seq_len(nrow(grid))) {
  params    <- grid[i, ]
  job_label <- sprintf("model=%s_p=%d_annot=%s",
                       params$model, params$p, params$annot_name)
  job_dir   <- file.path(OUTPUT_DIR, job_label)

  meta_path <- file.path(job_dir, "metadata.rds")
  eval_path <- file.path(job_dir, "evaluation.rds")

  completed <- file.exists(meta_path) && file.exists(eval_path)

  summary_rows[[i]] <- data.frame(
    job_id    = i,
    job_label = job_label,
    model     = params$model,
    p         = params$p,
    annot     = params$annot_name,
    completed = completed,
    stringsAsFactors = FALSE
  )

  if (!completed) {
    cat(sprintf("  [%2d] MISSING : %s\n", i, job_label))
    next
  }

  eval_obj <- tryCatch(readRDS(eval_path), error = function(e) NULL)
  if (is.null(eval_obj)) {
    cat(sprintf("  [%2d] READ ERROR : %s\n", i, job_label))
    next
  }

  flat <- tryCatch(.flatten_eval(eval_obj), error = function(e) NULL)
  if (is.null(flat) || nrow(flat) == 0L) {
    cat(sprintf("  [%2d] EMPTY    : %s\n", i, job_label))
    next
  }

  flat$job_id    <- i
  flat$job_label <- job_label
  flat$model     <- params$model
  flat$p_snps    <- params$p
  flat$annot     <- params$annot_name
  eval_list[[length(eval_list) + 1L]] <- flat

  cat(sprintf("  [%2d] OK : %s  (%d rows)\n", i, job_label, nrow(flat)))
}

# =============================================================================
# Combine and save
# =============================================================================

run_summary <- do.call(rbind, summary_rows)
n_done <- sum(run_summary$completed)

cat(sprintf("\n%d / %d jobs completed.\n", n_done, nrow(grid)))

summary_path <- file.path(OUTPUT_DIR, "run_summary.csv")
write.csv(run_summary, summary_path, row.names = FALSE)
cat(sprintf("Run summary written to: %s\n", summary_path))

if (length(eval_list) > 0L) {
  combined      <- do.call(rbind, eval_list)
  combined_path <- file.path(OUTPUT_DIR, "combined_evaluation.rds")
  csv_path      <- file.path(OUTPUT_DIR, "combined_evaluation.csv")
  saveRDS(combined, combined_path)
  write.csv(combined, csv_path, row.names = FALSE)
  cat(sprintf("Combined evaluation (rds): %s\n", combined_path))
  cat(sprintf("Combined evaluation (csv): %s\n", csv_path))
  cat(sprintf("  Rows: %d, Columns: %s\n", nrow(combined),
              paste(names(combined), collapse = ", ")))
} else {
  cat("No completed evaluations found — nothing to combine.\n")
}

# =============================================================================
# Identify incomplete jobs (for resubmission)
# =============================================================================

incomplete <- run_summary[!run_summary$completed, "job_id"]
if (length(incomplete) > 0L) {
  cat(sprintf("\nIncomplete job IDs (resubmit with --array):\n"))
  cat(paste(incomplete, collapse = ","), "\n")
  cat(sprintf("\nResubmit command:\n"))
  cat(sprintf("  sbatch --array=%s scripts/hpc/submit_benchmark.sh\n",
              paste(incomplete, collapse = ",")))
}
