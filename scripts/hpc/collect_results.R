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

# Same override the worker uses (submit_benchmark_pbs.sh points this at
# scratch). Set FMB_OUTPUT_ROOT before running collect if the array wrote
# to a non-default location.
OUTPUT_ROOT <- Sys.getenv("FMB_OUTPUT_ROOT", unset = "results/benchmark")
PARAMS_CSV  <- "scripts/hpc/params_grid.csv"

if (!dir.exists(OUTPUT_ROOT)) {
  stop("No output dir at ", OUTPUT_ROOT, ". Has the array run?")
}

# Evaluations live at job_*/scenario_*/evaluation.rds (scenario-split array).
# Fall back to the older flat job_*/evaluation.rds layout if present.
eval_paths <- Sys.glob(file.path(OUTPUT_ROOT, "job_*", "scenario_*", "evaluation.rds"))
if (length(eval_paths) == 0L) {
  eval_paths <- Sys.glob(file.path(OUTPUT_ROOT, "job_*", "evaluation.rds"))
}
if (length(eval_paths) == 0L) {
  stop("No evaluation.rds found under ", OUTPUT_ROOT,
       " (looked in job_*/scenario_*/ and job_*/). Has the array run?")
}

# Supplemental re-runs (FMB_METHODS=...) write evaluation_supp.rds next to the
# original. Those methods were re-run after a wrapper fix, so they SUPERSEDE
# the stale entries in evaluation.rds for the same scenario.
supp_for <- function(eval_path) {
  sp <- file.path(dirname(eval_path), "evaluation_supp.rds")
  if (file.exists(sp)) tryCatch(readRDS(sp), error = function(e) NULL) else NULL
}
n_supp <- length(Sys.glob(file.path(OUTPUT_ROOT, "job_*", "scenario_*",
                                    "evaluation_supp.rds")))
if (n_supp > 0L) {
  cat(sprintf("Found %d supplemental evaluation(s); these override the originals.\n",
              n_supp))
}

grid <- if (file.exists(PARAMS_CSV)) {
  read.csv(PARAMS_CSV, stringsAsFactors = FALSE)
} else {
  NULL
}

cat(sprintf("Collecting %d evaluation file(s) from %s ...\n",
            length(eval_paths), OUTPUT_ROOT))

# -----------------------------------------------------------------------------
# Per-job status
# -----------------------------------------------------------------------------
status_rows <- list()
combined    <- list()

# Curve accumulators. We POOL COUNTS by (job, S, phi, method) - summing
# tp/fp/fn per threshold and n/n_causal/sum_pip per calibration bin - rather
# than averaging per-scenario RATES, which would be badly biased when many
# scenarios have tiny denominators. Counts are additive, so downstream
# analysis can pool these further to ANY coarser stratum (e.g. by
# annotation_correlation) just by summing again. The 5 iterations within a
# (job, S, phi) cell are replicates, so pooling them here is the natural unit
# and keeps the artefact ~1.9M rows instead of ~9.4M.
fdr_acc <- new.env(hash = TRUE, parent = emptyenv())
cal_acc <- new.env(hash = TRUE, parent = emptyenv())
THRESH  <- NULL

.acc <- function(env, key, m) {
  if (!is.null(env[[key]])) env[[key]] <- env[[key]] + m else env[[key]] <- m
}

# Average precision from a stored fdr_power_curve. The trapezoidal `auprc`
# field is an invalid PR-space interpolation (Davis & Goadrich 2006) and
# understates by ~29% at weak signal / ~0% at strong - a SIGNAL-DEPENDENT
# distortion, so it cannot be cancelled by comparison. AP is the step-function
# estimator; from this 201-threshold grid it recovers the exact ranking-based
# AP to within ~1-5%. Recomputed here so existing runs are salvaged without
# re-evaluating (evaluate.R now also stores `ap` natively for future runs).
.ap_from_curve <- function(v) {
  ap <- v$ap %||% NULL
  if (!is.null(ap) && !is.na(ap)) return(as.numeric(ap))
  cur <- v$fdr_power_curve
  if (is.null(cur) || nrow(cur) < 2L) return(NA_real_)
  o <- order(cur$recall, -cur$precision)
  sum(diff(c(0, cur$recall[o])) * cur$precision[o])
}

.metric_row <- function(v, method, stratum, stratum_value) {
  data.frame(
    method        = method,
    stratum       = stratum,
    stratum_value = stratum_value,
    ap            = .ap_from_curve(v),
    auprc         = v$auprc        %||% NA_real_,
    auroc         = v$auroc        %||% NA_real_,
    cs_size_mean  = v$cs_size_mean %||% NA_real_,
    cs_coverage   = v$cs_coverage  %||% NA_real_,
    cs_power      = v$cs_power     %||% NA_real_,
    n_fits        = v$n_fits       %||% NA_integer_,
    n_failed      = v$n_failed     %||% NA_integer_,
    stringsAsFactors = FALSE
  )
}

extract_global <- function(eval_method, method) {
  g <- eval_method$global
  if (is.null(g)) return(NULL)
  .metric_row(g, method, "global", NA_character_)
}

extract_stratum <- function(eval_method, method, stratum_name) {
  s <- eval_method[[stratum_name]]
  if (is.null(s) || length(s) == 0L) return(NULL)
  do.call(rbind, lapply(names(s), function(key) {
    v <- s[[key]]
    if (is.null(v)) return(NULL)
    .metric_row(v, method, stratum_name, key)
  }))
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

for (eval_path in eval_paths) {
  unit_dir_path <- dirname(eval_path)
  is_scenario_layout <- grepl("^scenario_", basename(unit_dir_path))
  job_label      <- if (is_scenario_layout) basename(dirname(unit_dir_path))
                    else                     basename(unit_dir_path)
  scenario_label <- if (is_scenario_layout) basename(unit_dir_path) else NA_character_

  evaluation <- tryCatch(readRDS(eval_path), error = function(e) NULL)
  if (is.null(evaluation)) {
    status_rows[[length(status_rows) + 1L]] <- data.frame(
      job_dir = job_label, scenario = scenario_label,
      status = "unreadable_evaluation.rds",
      runtime_sec = NA_real_, stringsAsFactors = FALSE)
    next
  }

  # Overlay supplemental (re-run) methods on top of the originals.
  supp <- supp_for(eval_path)
  if (!is.null(supp)) {
    for (sm in names(supp)) {
      em <- supp[[sm]]
      if (is.list(em) && !is.null(em[["global"]])) evaluation[[sm]] <- em
    }
  }

  status_rows[[length(status_rows) + 1L]] <- data.frame(
    job_dir = job_label, scenario = scenario_label,
    status = "ok", runtime_sec = NA_real_, stringsAsFactors = FALSE)

  for (method in names(evaluation)) {
    em <- evaluation[[method]]
    # Skip non-method metadata entries (methods_evaluated, simulation_params,
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
    df$job_dir  <- job_label
    df$scenario <- scenario_label
    combined[[length(combined) + 1L]] <- df

    # --- accumulate curves, keyed by (job, S, phi, method) ------------------
    # A scenario has exactly one S and one phi, so read them off the strata.
    S_val   <- names(em$by_S   %||% list())[1] %||% NA_character_
    phi_val <- names(em$by_phi %||% list())[1] %||% NA_character_
    key <- paste(job_label, S_val, phi_val, method, sep = "|")

    cur <- em$global$fdr_power_curve
    if (!is.null(cur) && nrow(cur) > 0L) {
      if (is.null(THRESH)) THRESH <<- cur$threshold
      .acc(fdr_acc, key, cbind(tp = cur$tp, fp = cur$fp, fn = cur$fn))
    }
    pc <- em$global$pip_calibration
    if (!is.null(pc) && nrow(pc) > 0L) {
      # store sum_pip (= n * mean_pip) so mean_pip is recoverable after pooling
      sp <- ifelse(is.na(pc$mean_pip), 0, pc$mean_pip) * pc$n
      .acc(cal_acc, key, cbind(n = pc$n, n_causal = pc$n_causal, sum_pip = sp))
    }
  }
}

# --- Materialise pooled curves ----------------------------------------------
.unkey <- function(keys) {
  parts <- strsplit(keys, "|", fixed = TRUE)
  data.frame(
    job_dir = vapply(parts, `[`, character(1), 1L),
    S       = vapply(parts, `[`, character(1), 2L),
    phi     = vapply(parts, `[`, character(1), 3L),
    method  = vapply(parts, `[`, character(1), 4L),
    stringsAsFactors = FALSE)
}

fdr_keys <- ls(fdr_acc)
if (length(fdr_keys) > 0L) {
  meta <- .unkey(fdr_keys)
  fdr_long <- do.call(rbind, lapply(seq_along(fdr_keys), function(i) {
    m <- fdr_acc[[fdr_keys[i]]]
    data.frame(meta[i, , drop = FALSE], row.names = NULL,
               threshold = THRESH, tp = m[, "tp"], fp = m[, "fp"], fn = m[, "fn"],
               stringsAsFactors = FALSE)
  }))
  saveRDS(fdr_long, file.path(OUTPUT_ROOT, "combined_fdr_curves.rds"))
  cat(sprintf("  combined_fdr_curves.rds       (%d rows, pooled counts)\n",
              nrow(fdr_long)))
}

cal_keys <- ls(cal_acc)
if (length(cal_keys) > 0L) {
  meta <- .unkey(cal_keys)
  cal_long <- do.call(rbind, lapply(seq_along(cal_keys), function(i) {
    m <- cal_acc[[cal_keys[i]]]
    data.frame(meta[i, , drop = FALSE], row.names = NULL,
               bin = seq_len(nrow(m)), n = m[, "n"], n_causal = m[, "n_causal"],
               sum_pip = m[, "sum_pip"], stringsAsFactors = FALSE)
  }))
  cal_long$mean_pip    <- ifelse(cal_long$n > 0, cal_long$sum_pip / cal_long$n, NA_real_)
  cal_long$frac_causal <- ifelse(cal_long$n > 0, cal_long$n_causal / cal_long$n, NA_real_)
  saveRDS(cal_long, file.path(OUTPUT_ROOT, "combined_pip_calibration.rds"))
  cat(sprintf("  combined_pip_calibration.rds  (%d rows, pooled counts)\n",
              nrow(cal_long)))
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
  # NB: "p"/"enrichment"/"n_ref" are from the OLD grid schema and silently
  # vanished via intersect(), taking p_causal + annotation_correlation with
  # them - i.e. two of the plan's headline swept axes (Section 1.3) were
  # missing from the aggregate. Use the current schema.
  keep_cols <- c("job_dir", "model", "p_causal", "annotation_type",
                 "annotation_correlation", "n_annotations", "n_regions",
                 "n", "n_iter")
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
