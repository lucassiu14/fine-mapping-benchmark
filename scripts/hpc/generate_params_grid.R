#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/generate_params_grid.R
#
# Build params_grid.csv for the auto-research Phase 1 simulation (see
# docs/autoresearch/iteration-001-phase1-grid.md and §1.1 / §1.3 of the plan).
#
# Encoding:
#   - Each row of the resulting CSV is ONE SLURM array task.
#   - Within a task, run_benchmark_job.R sweeps the WITHIN_JOB axes (S, phi,
#     n_iter) so a single run_simulation() call produces all 125 scenarios for
#     that task's (model, annotation_regime, p_causal) combination.
#   - The per-region p vector, enrichment vector, and n_annotations are FIXED
#     across all tasks per §1.1; they are serialised into the row for
#     reproducibility but never vary between rows.
#
# Total tasks: 25
#   - sparse         x 5 annotation regimes                =  5 tasks
#   - sparse_inf     x 5 annotation regimes x 4 p_causal   = 20 tasks
# Total scenarios: 25 * 125 = 3125  (matches §1.3)
#
# Usage:
#   Rscript scripts/hpc/generate_params_grid.R [output_path]
#   Default output: scripts/hpc/params_grid.csv
# =============================================================================

out_path <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_path) || !nzchar(out_path)) {
  out_path <- "scripts/hpc/params_grid.csv"
}

# -----------------------------------------------------------------------------
# Fixed simulation constants (§1.1)
# -----------------------------------------------------------------------------

# Per-region p vector: 4 regions per size class {100, 200, 400, 500, 1000},
# length 20 total. NOT a swept axis - the same vector applies to every task.
P_VECTOR <- c(rep(100L, 4L), rep(200L, 4L), rep(400L, 4L),
              rep(500L, 4L), rep(1000L, 4L))

# Fixed annotation setup: 20 annotations. 4 truly enriched (2 at fold 7.4,
# 2 at fold 2.7), 16 null decoys at fold 1. Same across every task.
N_ANNOTATIONS      <- 20L
ENRICHMENT_VECTOR  <- c(7.4, 7.4, 2.7, 2.7, rep(1, 16L))
stopifnot(length(ENRICHMENT_VECTOR) == N_ANNOTATIONS)

# Within-job sweep: run_simulation() expands the full grid internally when
# handed vector S / phi values. p_causal is a scalar per task (see below).
WITHIN_JOB <- list(
  n_regions = length(P_VECTOR),     # 20
  n         = 1000L,
  n_iter    = 5L,
  S         = c(1L, 2L, 3L, 5L, 10L),
  phi       = c(0.0075, 0.05, 0.1, 0.2, 0.4)
)

# p_causal only applies to the sparse_inf model. Four values from §1.1.
P_CAUSAL_SPARSE_INF <- c(0.5, 0.7, 0.9, 1.0)

# Annotation regimes: none + binary at four correlation levels (§1.3).
ANNOTATION_REGIMES <- data.frame(
  annotation_type        = c("none",
                             rep("binary", 4L)),
  annotation_correlation = c(NA_real_,
                             0, 0.25, 0.5, 0.75),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Build the grid: 25 rows total
# -----------------------------------------------------------------------------

rows <- list()
job_id <- 0L

for (m in c("sparse", "sparse_inf")) {

  # p_causal expansion: sparse has one row per annotation regime; sparse_inf
  # has one row per (annotation regime, p_causal) pair.
  p_causal_values <- if (m == "sparse") NA_real_ else P_CAUSAL_SPARSE_INF

  for (pc in p_causal_values) {
    for (ar in seq_len(nrow(ANNOTATION_REGIMES))) {
      job_id <- job_id + 1L

      annot_type <- ANNOTATION_REGIMES$annotation_type[ar]
      annot_corr <- ANNOTATION_REGIMES$annotation_correlation[ar]

      label <- sprintf(
        "%s_%s%s%s",
        m,
        c(none = "anNone", binary = "anBinary")[annot_type],
        if (annot_type == "binary") sprintf("_ac%s", format(annot_corr, nsmall = 2)) else "",
        if (m == "sparse_inf") sprintf("_pc%g", pc) else ""
      )

      rows[[job_id]] <- data.frame(
        job_id                 = job_id,
        label                  = label,
        model                  = m,
        p_causal               = pc,
        annotation_type        = annot_type,
        annotation_correlation = annot_corr,
        n_annotations          = N_ANNOTATIONS,
        n_regions              = WITHIN_JOB$n_regions,
        n                      = WITHIN_JOB$n,
        n_iter                 = WITHIN_JOB$n_iter,
        S_values               = paste(WITHIN_JOB$S,   collapse = "|"),
        phi_values             = paste(WITHIN_JOB$phi, collapse = "|"),
        p_values               = paste(P_VECTOR,       collapse = "|"),
        enrichment_values      = paste(ENRICHMENT_VECTOR, collapse = "|"),
        stringsAsFactors       = FALSE
      )
    }
  }
}

grid <- do.call(rbind, rows)

# -----------------------------------------------------------------------------
# Write out
# -----------------------------------------------------------------------------

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(grid, out_path, row.names = FALSE)

cat(sprintf("Wrote %d tasks to %s\n", nrow(grid), out_path))
n_scenarios_per <- length(WITHIN_JOB$S) * length(WITHIN_JOB$phi) * WITHIN_JOB$n_iter
cat(sprintf("Scenarios per task: %d (S=%d, phi=%d, n_iter=%d)\n",
            n_scenarios_per, length(WITHIN_JOB$S), length(WITHIN_JOB$phi),
            WITHIN_JOB$n_iter))
cat(sprintf("Total scenarios: %d\n", nrow(grid) * n_scenarios_per))
cat(sprintf("Total region-scenario-method fits (with 9 methods): %d\n",
            nrow(grid) * n_scenarios_per * WITHIN_JOB$n_regions * 9L))

cat("\nFirst few rows (compact):\n")
print(grid[seq_len(min(6, nrow(grid))),
           c("job_id", "label", "model", "p_causal",
             "annotation_type", "annotation_correlation")],
      row.names = FALSE)
