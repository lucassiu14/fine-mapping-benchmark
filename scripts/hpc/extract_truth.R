#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/extract_truth.R
#
# Rescue the GROUND TRUTH from a benchmark output tree before ephemeral purges.
#
# WHY. $EPHEMERAL is deleted 30 days after write and is not backed up. The
# per-row sim.rds files hold the only copy of the causal variant indices; every
# metric in the analysis is scored against them, so once they go the surviving
# results.rds files cannot be scored at all. sim.rds is also the OLDEST object
# in the tree - it is written once per row at the start and never rewritten,
# even when results are recomputed - so it expires first, potentially leaving
# PIPs with no truth to score them against.
#
# This extracts only the truth. It does not touch results.rds, so it is small
# and fast: one integer vector per (scenario, region) plus a few scalars.
#
#   Rscript scripts/hpc/extract_truth.R <array_idx> <bench_root> <out_dir>
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
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

sim_file <- file.path(job_dir, "sim.rds")
if (!file.exists(sim_file)) {
  stop("no sim.rds in ", job_dir, " - truth already lost for this row.", call. = FALSE)
}
message("row ", idx, ": ", label)
sim <- readRDS(sim_file)

p_nominal <- sim$params$p
n_reg     <- length(p_nominal)
# Same convention as iter004_collect.R: two same-size regions per class, in
# order, so region_id 1,2 are the first class, 3,4 the second, and so on.
region_idx_of <- function(region_id) ((region_id - 1L) %% 2L) + 1L

rows <- vector("list", length(sim$scenarios) * n_reg)
k <- 0L
for (sc_idx in seq_along(sim$scenarios)) {
  sc <- sim$scenarios[[sc_idx]]
  for (rg in seq_len(min(n_reg, length(sc$regions)))) {
    tr <- sc$regions[[rg]]$truth
    if (is.null(tr)) next
    k <- k + 1L
    rows[[k]] <- list(
      job_dir         = label,
      scenario_id     = sc_idx,
      region_id       = rg,
      region_idx      = region_idx_of(rg),
      n_variants      = p_nominal[rg],
      S               = tr$S               %||% sc$S,
      phi             = tr$phi             %||% sc$phi,
      model           = tr$model           %||% sc$model,
      p_causal        = tr$p_causal        %||% NA_real_,
      annotation_type = tr$annotation_type %||% NA_character_,
      causal_indices  = as.integer(tr$causal_indices),
      # Present only where the simulator stored it (Iteration 005 onward);
      # polyfun_oracle's true prior cannot be reconstructed without it.
      causal_probs    = tr$causal_probs
    )
  }
}
rows <- rows[seq_len(k)]
out <- file.path(out_dir, sprintf("truth_%s.rds", label))
saveRDS(list(job_dir = label, n_regions = n_reg, p = p_nominal,
             n_scenarios = length(sim$scenarios), truth = rows),
        out)
message(sprintf("  wrote %s  (%d scenario-region records, %.1f KB)",
                basename(out), k, file.size(out) / 1024))
