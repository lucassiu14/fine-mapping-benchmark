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
#   results/benchmark/combined_evaluation.rds  — combined eval data frame
#   results/benchmark/run_summary.csv          — per-job completion status
# =============================================================================

PROJECT_ROOT <- "."
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "results", "benchmark")
GRID_PATH    <- file.path(PROJECT_ROOT, "scripts", "hpc", "params_grid.csv")

# =============================================================================
# Load project
# =============================================================================

r_files <- list.files(file.path(PROJECT_ROOT, "R"),
                       pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
invisible(lapply(r_files, source))

# =============================================================================
# Scan completed jobs
# =============================================================================

grid <- read.csv(GRID_PATH, stringsAsFactors = FALSE)
cat(sprintf("Scanning %d jobs...\n", nrow(grid)))

job_dirs <- list.dirs(OUTPUT_DIR, recursive = FALSE, full.names = TRUE)

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

  # Attach job-level parameters to the evaluation data frame
  if (is.data.frame(eval_obj)) {
    eval_obj$job_id    <- i
    eval_obj$job_label <- job_label
    eval_obj$p_snps    <- params$p
    eval_obj$annot     <- params$annot_name
    eval_list[[length(eval_list) + 1]] <- eval_obj
  } else if (is.list(eval_obj) && !is.null(eval_obj$metrics)) {
    # If evaluate_methods returns a list, pull out the metrics data frame
    df <- eval_obj$metrics
    df$job_id    <- i
    df$job_label <- job_label
    df$p_snps    <- params$p
    df$annot     <- params$annot_name
    eval_list[[length(eval_list) + 1]] <- df
  }

  cat(sprintf("  [%2d] OK : %s\n", i, job_label))
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

if (length(eval_list) > 0) {
  combined <- do.call(rbind, eval_list)
  combined_path <- file.path(OUTPUT_DIR, "combined_evaluation.rds")
  saveRDS(combined, combined_path)
  cat(sprintf("Combined evaluation written to: %s\n", combined_path))
  cat(sprintf("  Rows: %d, Columns: %s\n", nrow(combined),
              paste(names(combined), collapse = ", ")))
} else {
  cat("No completed evaluations found — nothing to combine.\n")
}

# =============================================================================
# Identify incomplete jobs (for resubmission)
# =============================================================================

incomplete <- run_summary[!run_summary$completed, "job_id"]
if (length(incomplete) > 0) {
  cat(sprintf("\nIncomplete job IDs (resubmit with --array):\n"))
  cat(paste(incomplete, collapse = ","), "\n")
  cat(sprintf("\nResubmit command:\n"))
  cat(sprintf("  sbatch --array=%s scripts/hpc/submit_benchmark.sh\n",
              paste(incomplete, collapse = ",")))
}
