#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/collect_results.R
#
# Walk results/benchmark/job_*/ and stitch the per-job evaluation.rds files
# into:
#   results/benchmark/combined_evaluation.rds  - long data frame (one row per
#                                                method x scenario x stratum)
#   results/benchmark/combined_evaluation.csv  - same, CSV
#   results/benchmark/run_summary.csv          - per-job completion status
#
# Usage (from project root, after the array finishes):
#   Rscript scripts/hpc/collect_results.R
# =============================================================================

OUTPUT_ROOT <- "results/benchmark"
PARAMS_CSV  <- "scripts/hpc/params_grid.csv"

if (!dir.exists(OUTPUT_ROOT)) {
  stop("No output dir at ", OUTPUT_ROOT, ". Has the array run?")
}

job_dirs <- list.dirs(OUTPUT_ROOT, recursive = FALSE)
job_dirs <- job_dirs[grepl("^job_", basename(job_dirs))]
if (length(job_dirs) == 0L) {
  stop("No job_*/ subdirectories found under ", OUTPUT_ROOT)
}

grid <- if (file.exists(PARAMS_CSV)) {
  read.csv(PARAMS_CSV, stringsAsFactors = FALSE)
} else {
  NULL
}

cat(sprintf("Collecting %d job directories from %s ...\n",
            length(job_dirs), OUTPUT_ROOT))

# -----------------------------------------------------------------------------
# Per-job status
# -----------------------------------------------------------------------------
status_rows <- list()
combined    <- list()

extract_global <- function(eval_method, method) {
  g <- eval_method$global
  if (is.null(g)) return(NULL)
  data.frame(
    method        = method,
    stratum       = "global",
    stratum_value = NA_character_,
    auprc         = g$auprc       %||% NA_real_,
    auroc         = g$auroc       %||% NA_real_,
    cs_size_mean  = g$cs_size_mean %||% NA_real_,
    cs_coverage   = g$cs_coverage  %||% NA_real_,
    cs_power      = g$cs_power     %||% NA_real_,
    n_fits        = g$n_fits       %||% NA_integer_,
    stringsAsFactors = FALSE
  )
}

extract_stratum <- function(eval_method, method, stratum_name) {
  s <- eval_method[[stratum_name]]
  if (is.null(s) || length(s) == 0L) return(NULL)
  do.call(rbind, lapply(names(s), function(key) {
    v <- s[[key]]
    if (is.null(v)) return(NULL)
    data.frame(
      method        = method,
      stratum       = stratum_name,
      stratum_value = key,
      auprc         = v$auprc       %||% NA_real_,
      auroc         = v$auroc       %||% NA_real_,
      cs_size_mean  = v$cs_size_mean %||% NA_real_,
      cs_coverage   = v$cs_coverage  %||% NA_real_,
      cs_power      = v$cs_power     %||% NA_real_,
      n_fits        = v$n_fits       %||% NA_integer_,
      stringsAsFactors = FALSE
    )
  }))
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

for (jd in job_dirs) {
  job_label <- basename(jd)
  eval_path <- file.path(jd, "evaluation.rds")
  meta_path <- file.path(jd, "params.json")

  if (!file.exists(eval_path)) {
    status_rows[[length(status_rows) + 1L]] <- data.frame(
      job_dir = job_label, status = "missing_evaluation.rds",
      runtime_sec = NA_real_, stringsAsFactors = FALSE)
    next
  }

  evaluation <- tryCatch(readRDS(eval_path),
                         error = function(e) NULL)
  if (is.null(evaluation)) {
    status_rows[[length(status_rows) + 1L]] <- data.frame(
      job_dir = job_label, status = "unreadable_evaluation.rds",
      runtime_sec = NA_real_, stringsAsFactors = FALSE)
    next
  }

  meta <- if (file.exists(meta_path)) {
    tryCatch(jsonlite::read_json(meta_path), error = function(e) NULL)
  } else NULL
  runtime_sec <- if (!is.null(meta) && !is.null(meta$runtime_sec)) {
    as.numeric(meta$runtime_sec)
  } else NA_real_

  status_rows[[length(status_rows) + 1L]] <- data.frame(
    job_dir     = job_label,
    status      = "ok",
    runtime_sec = runtime_sec,
    stringsAsFactors = FALSE)

  for (method in names(evaluation)) {
    em <- evaluation[[method]]
    # Skip the non-method metadata entries (methods_evaluated, simulation_params,
    # pip_thresholds_used). Method entries are lists with a `global` sub-list.
    if (is.null(em) || !is.list(em) || is.null(em[["global"]])) next
    rows <- list(
      extract_global(em, method),
      extract_stratum(em, method, "by_S"),
      extract_stratum(em, method, "by_phi"),
      extract_stratum(em, method, "by_p_causal"),
      extract_stratum(em, method, "by_causal_maf"),
      extract_stratum(em, method, "by_true_annotation_type")
    )
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0L) next
    df <- do.call(rbind, rows)
    df$job_dir <- job_label
    combined[[length(combined) + 1L]] <- df
  }
}

status_df <- do.call(rbind, status_rows)
write.csv(status_df, file.path(OUTPUT_ROOT, "run_summary.csv"), row.names = FALSE)

if (length(combined) == 0L) {
  cat("No readable evaluation files. Wrote run_summary.csv only.\n")
  quit(save = "no")
}

combined_df <- do.call(rbind, combined)

# Join the params_grid info if available so each row carries its job context.
if (!is.null(grid)) {
  grid$job_dir <- sprintf("job_%03d_%s", grid$job_id, grid$label)
  keep_cols <- c("job_dir", "model", "p", "annotation_type",
                 "n_annotations", "enrichment", "n_ref", "n", "n_iter")
  combined_df <- merge(combined_df, grid[, intersect(keep_cols, names(grid))],
                       by = "job_dir", all.x = TRUE, sort = FALSE)
}

saveRDS(combined_df, file.path(OUTPUT_ROOT, "combined_evaluation.rds"))
write.csv(combined_df, file.path(OUTPUT_ROOT, "combined_evaluation.csv"),
          row.names = FALSE)

cat(sprintf("\nWrote:\n  %s/combined_evaluation.rds (%d rows)\n",
            OUTPUT_ROOT, nrow(combined_df)))
cat(sprintf("  %s/combined_evaluation.csv\n", OUTPUT_ROOT))
cat(sprintf("  %s/run_summary.csv\n", OUTPUT_ROOT))

cat("\nPer-job status:\n")
print(status_df)
