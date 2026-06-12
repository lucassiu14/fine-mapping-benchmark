#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/generate_params_grid.R
#
# Build params_grid.csv: one row per *experimental condition* (model x region
# size x annotation regime). Each row becomes one SLURM array task; within a
# task we sweep S x phi x n_iter to produce many scenarios.
#
# Usage (from project root):
#   Rscript scripts/hpc/generate_params_grid.R [output_path]
#
# Default output: scripts/hpc/params_grid.csv
# =============================================================================

out_path <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_path) || !nzchar(out_path)) {
  out_path <- "scripts/hpc/params_grid.csv"
}

# -----------------------------------------------------------------------------
# Job axes - edit these to change what the array covers.
# Each row of the resulting CSV is one job; you'll get prod(lengths) jobs total.
# -----------------------------------------------------------------------------
AXES <- list(
  model            = c("sparse", "sparse_inf"),
  p                = c(200L, 500L, 1000L),
  annotation_type  = c("none", "binary")
)

# Within-job sweep (recorded for traceability; the worker script also reads
# these so a single source of truth controls every job's inner grid).
WITHIN_JOB <- list(
  n_regions = 50L,
  S         = c(1L, 2L, 3L, 5L),
  phi       = c(0.1, 0.2, 0.4),
  n_iter    = 20L,
  n         = 500L
)

# Annotation parameters (only used when annotation_type != "none")
N_ANNOTATIONS <- 3L
ENRICHMENT    <- 5

# LD-mismatch experiment: set TRUE to add a second row per condition that
# fits methods with an *independent* reference panel of n_ref individuals.
ADD_LD_MISMATCH <- FALSE
N_REF           <- 500L

# -----------------------------------------------------------------------------
grid <- expand.grid(AXES, stringsAsFactors = FALSE)
grid$n_annotations <- ifelse(grid$annotation_type == "none", NA_integer_, N_ANNOTATIONS)
grid$enrichment    <- ifelse(grid$annotation_type == "none", NA_real_,    ENRICHMENT)
grid$n_ref         <- NA_integer_

if (ADD_LD_MISMATCH) {
  mm <- grid
  mm$n_ref <- N_REF
  grid <- rbind(grid, mm)
}

# Record the within-job sweep as fixed columns (one value per row, repeated).
grid$S_values      <- paste(WITHIN_JOB$S,   collapse = "|")
grid$phi_values    <- paste(WITHIN_JOB$phi, collapse = "|")
grid$n_iter        <- WITHIN_JOB$n_iter
grid$n             <- WITHIN_JOB$n
grid$n_regions     <- WITHIN_JOB$n_regions

# Human-readable label for each job
grid$label <- sprintf(
  "%s_p%d_an%s%s",
  grid$model,
  grid$p,
  c(none = "None", binary = "Binary", continuous = "Continuous")[grid$annotation_type],
  ifelse(is.na(grid$n_ref), "", sprintf("_ref%d", grid$n_ref))
)

grid$job_id <- seq_len(nrow(grid))

# Final column order
grid <- grid[, c("job_id", "label",
                 "model", "p", "annotation_type",
                 "n_annotations", "enrichment", "n_ref",
                 "n_regions", "n", "n_iter", "S_values", "phi_values")]

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(grid, out_path, row.names = FALSE)

cat(sprintf("Wrote %d jobs to %s\n", nrow(grid), out_path))
cat("\nFirst few rows:\n")
print(head(grid))
