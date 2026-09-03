#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/extract_credible_sets.R
#
# Rescue the methods' OWN credible sets before $EPHEMERAL purges results.rds.
#
#   Rscript scripts/hpc/extract_credible_sets.R <array_idx> <bench_root> <out_dir>
#
# WHY. iter004_collect.R stores only the uniform reconstruction - set_size_*,
# set_hit_*, set_prec_*, set_reached_* from .rank_set_metrics(), which cuts the
# ranked PIP vector at a fixed cumulative mass. The sets each method actually
# reported are dropped at collection: there is no cs_* column anywhere in
# combined_fit_metrics.rds, and none of the other three rescue jobs touches
# them. They exist only inside results.rds.
#
# They are needed to compare set CONSTRUCTIONS - SuSiE's purity-filtered
# per-component sets against FINEMAP's configuration-posterior sets against a
# uniform rule - which is a question about the methods as shipped, and cannot be
# answered from a reconstruction that deliberately discards their conventions.
#
# Cheap: a few sets of a few variants per fit, against ~50 stored PIPs per fit
# in the tail extractor. Also records n_sets and each set's size so a run can be
# summarised without unpacking every index vector.
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
job_dir <- job_dirs[idx]; label <- basename(job_dir)
message("row ", idx, ": ", label)

scen_dirs <- sort(list.dirs(job_dir, recursive = FALSE))
scen_dirs <- scen_dirs[grepl("/scenario_[0-9]+", scen_dirs)]

rows <- vector("list", length(scen_dirs) * 400L)
k <- 0L; n_fits <- 0L; n_with_sets <- 0L
for (sd in scen_dirs) {
  f <- file.path(sd, "results.rds")
  if (!file.exists(f)) next
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(r)) next
  sc_idx <- suppressWarnings(as.integer(sub("^scenario_", "", basename(sd))))

  for (m in names(r)) {
    v <- r[[m]]
    if (!is.list(v) || is.null(v$results)) next
    for (rg in seq_along(v$results)) {
      fit <- v$results[[rg]]
      if (!is.list(fit)) next
      n_fits <- n_fits + 1L
      cs <- fit$credible_sets
      # An empty list is a RESULT - the method reported no set - and is stored
      # as such. NULL means the wrapper never populated the field at all; the
      # two are distinguished by n_sets = 0 versus reported = FALSE.
      reported <- !is.null(cs)
      if (reported && length(cs)) n_with_sets <- n_with_sets + 1L
      k <- k + 1L
      rows[[k]] <- list(
        job_dir = label, scenario_id = sc_idx,
        region_id = as.integer(fit$region_id %||% rg), method = m,
        n_variants = length(fit$pip %||% integer(0)),
        reported = reported,
        n_sets = if (reported) length(cs) else NA_integer_,
        set_sizes = if (reported) vapply(cs, length, 1L) else integer(0),
        sets = if (reported) lapply(cs, as.integer) else list())
    }
  }
}
rows <- rows[seq_len(k)]
out <- file.path(out_dir, sprintf("cs_%s.rds", label))
saveRDS(list(job_dir = label, fits = k, rows = rows), out)
message(sprintf("  wrote %s  (%d fits, %d with at least one set, %.1f MB)",
                basename(out), n_fits, n_with_sets, file.size(out) / 1024^2))
