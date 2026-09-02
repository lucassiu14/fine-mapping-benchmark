#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/extract_pip_tail.R
#
# Rescue enough of the PIP vectors to recompute threshold-based metrics at ANY
# threshold, before $EPHEMERAL purges results.rds.
#
# WHY. iter004_collect.R stores counts at five fixed thresholds only
# (0.50, 0.80, 0.90, 0.95, 0.99). A power-FDR curve, or a false discovery rate
# at any other threshold, cannot be computed from those five points. Doing so
# requires the PIP vectors, which exist nowhere but results.rds on ephemeral.
#
# Storing the full vectors would be ~17 GB. But PIPs are extremely sparse -
# in Iteration 004 fewer than 0.5% of variants exceed 0.1 - so retaining only
# the tail above a floor reconstructs every threshold-based quantity above that
# floor exactly, at a small fraction of the size. Below the floor only the
# COUNTS are kept, which is all that any threshold above it needs.
#
#   Rscript scripts/hpc/extract_pip_tail.R <array_idx> <bench_root> <out_dir> [floor]
#
# Default floor 0.01. Thresholds below it cannot be recovered; nothing in
# fine-mapping practice uses them.
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

args       <- commandArgs(trailingOnly = TRUE)
idx        <- suppressWarnings(as.integer(args[1]))
bench_root <- args[2]
out_dir    <- args[3]
floor_pip  <- if (length(args) >= 4) as.numeric(args[4]) else 0.01
stopifnot(!is.na(idx), dir.exists(bench_root), is.finite(floor_pip))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

job_dirs <- sort(list.dirs(bench_root, recursive = FALSE))
job_dirs <- job_dirs[grepl("/job_[0-9]+", job_dirs)]
if (idx > length(job_dirs)) {
  message("array index ", idx, " > ", length(job_dirs), " rows; nothing to do")
  quit(save = "no", status = 0L)
}
job_dir <- job_dirs[idx]; label <- basename(job_dir)
message("row ", idx, ": ", label, "   floor = ", floor_pip)

# Truth comes from sim.rds; without it the tail cannot be scored.
sim_file <- file.path(job_dir, "sim.rds")
sim <- if (file.exists(sim_file)) readRDS(sim_file) else NULL
if (is.null(sim)) message("  WARNING: no sim.rds; storing tails without causal flags")

scen_dirs <- sort(list.dirs(job_dir, recursive = FALSE))
scen_dirs <- scen_dirs[grepl("/scenario_[0-9]+", scen_dirs)]

rows <- vector("list", length(scen_dirs) * 400L)
k <- 0L; n_kept <- 0; n_tot <- 0

for (sd in scen_dirs) {
  f <- file.path(sd, "results.rds")
  if (!file.exists(f)) next
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(r)) next
  sc_idx <- suppressWarnings(as.integer(sub("^scenario_", "", basename(sd))))
  sc <- if (!is.null(sim) && sc_idx <= length(sim$scenarios)) sim$scenarios[[sc_idx]] else NULL

  for (m in names(r)) {
    v <- r[[m]]
    if (!is.list(v) || is.null(v$results)) next
    for (rg in seq_along(v$results)) {
      fit <- v$results[[rg]]
      pip <- fit$pip
      if (is.null(pip) || !length(pip)) next
      region_id <- as.integer(fit$region_id %||% rg)
      causal <- if (!is.null(sc)) as.integer(sc$regions[[region_id]]$truth$causal_indices) else integer(0)

      keep <- which(is.finite(pip) & pip >= floor_pip)
      n_tot <- n_tot + length(pip); n_kept <- n_kept + length(keep)
      k <- k + 1L
      rows[[k]] <- list(
        job_dir = label, scenario_id = sc_idx, region_id = region_id, method = m,
        n_variants = length(pip),
        n_causal   = length(causal),
        # Causal variants are ALWAYS retained regardless of their PIP: a causal
        # variant the method assigned low probability is exactly the case a
        # power calculation must count, and dropping it would silently inflate
        # recall at every threshold.
        idx        = sort(unique(c(keep, causal))),
        pip        = pip[sort(unique(c(keep, causal)))],
        is_causal  = as.integer(sort(unique(c(keep, causal))) %in% causal),
        # Everything below the floor, summarised: enough to normalise any rate.
        n_below    = length(pip) - length(keep),
        sum_pip_below = sum(pip[setdiff(seq_along(pip), keep)], na.rm = TRUE))
    }
  }
}
rows <- rows[seq_len(k)]
out <- file.path(out_dir, sprintf("piptail_%s.rds", label))
saveRDS(list(job_dir = label, floor = floor_pip, fits = k, tail = rows), out)
message(sprintf("  wrote %s  (%d fits, %.2f%% of PIPs retained, %.1f MB)",
                basename(out), k, 100 * n_kept / max(n_tot, 1), file.size(out) / 1024^2))
