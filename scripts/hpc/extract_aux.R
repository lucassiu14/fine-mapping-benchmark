#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/extract_aux.R
#
# Rescue RUNTIME and ANNOTATION-IMPORTANCE from a benchmark tree before
# $EPHEMERAL purges it.
#
# WHY. iter004_collect.R stores only the metrics needed for the accuracy
# analysis. Two quantities the paper needs are therefore absent from
# combined_fit_metrics.rds and survive only inside results.rds on ephemeral:
#
#   runtime_seconds        per fit, and total_runtime_seconds per method per
#                          scenario. Needed for the computational-cost section.
#   feature_importance     the LassoNet annotation importances returned by the
#                          Functional BEATRICE family. Needed for the
#                          annotation-selection section.
#
# Neither can be recomputed without re-running the benchmark. This reads
# results.rds only - it does not touch sim.rds, which extract_truth.R handles.
#
#   Rscript scripts/hpc/extract_aux.R <array_idx> <bench_root> <out_dir>
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

args       <- commandArgs(trailingOnly = TRUE)
idx        <- suppressWarnings(as.integer(args[1]))
bench_root <- args[2]
out_dir    <- args[3]
stopifnot(!is.na(idx), dir.exists(bench_root))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

job_dirs <- sort(list.dirs(bench_root, recursive = FALSE))
job_dirs <- job_dirs[grepl("/job_[0-9]+", job_dirs)]
if (idx > length(job_dirs)) {
  message("array index ", idx, " > ", length(job_dirs), " rows; nothing to do")
  quit(save = "no", status = 0L)
}
job_dir <- job_dirs[idx]
label   <- basename(job_dir)
message("row ", idx, ": ", label)

scen_dirs <- sort(list.dirs(job_dir, recursive = FALSE))
scen_dirs <- scen_dirs[grepl("/scenario_[0-9]+", scen_dirs)]

rt_rows  <- vector("list", length(scen_dirs) * 400L)   # per fit
fi_rows  <- vector("list", length(scen_dirs) * 400L)   # per fit, FB family only
kr <- 0L; kf <- 0L; n_read <- 0L

for (sd in scen_dirs) {
  f <- file.path(sd, "results.rds")
  if (!file.exists(f)) next
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(r)) next
  n_read <- n_read + 1L
  sc_idx <- suppressWarnings(as.integer(sub("^scenario_", "", basename(sd))))

  for (m in names(r)) {
    v <- r[[m]]
    # Skip the non-method entries the worker also writes (methods_run,
    # run_timestamp, simulation_params).
    if (!is.list(v) || is.null(v$results)) next

    kr <- kr + 1L
    rt_rows[[kr]] <- list(
      job_dir = label, scenario_id = sc_idx, method = m,
      region_id = NA_integer_,
      runtime_seconds = as.numeric(v$total_runtime_seconds %||% NA_real_),
      scope = "scenario_total",
      n_total  = as.integer(v$n_total  %||% NA_integer_),
      n_failed = as.integer(v$n_failed %||% NA_integer_))

    for (rg in seq_along(v$results)) {
      fit <- v$results[[rg]]
      if (!is.list(fit)) next

      # Per-fit runtime: finer than the scenario total, and lets cost be
      # stratified by region size.
      if (!is.null(fit$runtime_seconds)) {
        kr <- kr + 1L
        rt_rows[[kr]] <- list(
          job_dir = label, scenario_id = sc_idx, method = m,
          region_id = as.integer(fit$region_id %||% rg),
          runtime_seconds = as.numeric(fit$runtime_seconds),
          scope = "fit", n_total = NA_integer_, n_failed = NA_integer_)
      }

      # Annotation importance, where the method returns it.
      fi <- fit$additional$feature_importance
      if (!is.null(fi)) {
        kf <- kf + 1L
        fi_rows[[kf]] <- list(
          job_dir = label, scenario_id = sc_idx, method = m,
          region_id = as.integer(fit$region_id %||% rg),
          # The cross-region wrappers never return importances from the joint
          # run itself, so an importance record labelled fb_xregion / fb_pooled
          # can only have come from the per-region fallback in
          # wrapper_fb_joint.R. Carry the flag or the row is a single-locus
          # result wearing a cross-region label.
          joint_fallback = isTRUE(fit$additional$joint_fallback),
          importance = fi)
      }
    }
  }
}

rt_rows <- rt_rows[seq_len(kr)]
fi_rows <- fi_rows[seq_len(kf)]
out <- file.path(out_dir, sprintf("aux_%s.rds", label))
saveRDS(list(job_dir = label, scenarios_read = n_read,
             runtime = rt_rows, importance = fi_rows), out)
message(sprintf("  wrote %s  (%d scenarios, %d runtime records, %d importance records, %.1f MB)",
                basename(out), n_read, kr, kf, file.size(out) / 1024^2))
