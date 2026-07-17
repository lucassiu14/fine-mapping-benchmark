#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/generate_params_grid.R
#
# Build params_grid.csv for the auto-research simulation grid.
#
# ***** ITERATION 002 (2026-07-17) *****
# Changes from Iteration 001 (documented so the two grids are diffable):
#   - annotation arms: none / binary / CONTINUOUS  (was none / binary)
#   - annotation_correlation: 0 only               (was 0, 0.25, 0.5, 0.75)
#   - n_annotations: 10                            (was 20)
#   - ENRICHMENT is now a SWEPT axis: {2.7, 5.4, 8.1, 10.8}, applied to HALF
#     the annotations (first 5 at fold E, last 5 inert at fold 1)  (was fixed)
#   - n_regions: 10  (5 length classes x 2 each)   (was 20, 5 x 4)
#   - n_iter: 10                                   (was 5)
#     [per-scenario replication 10 regions x 10 iter = 100, same as 20 x 5]
#   - NEW AXIS n_ref (LD reference-panel size = "LD with noise"):
#     NA = in-sample LD (perfect), 500, 200. Smaller panel = noisier LD.
#   Unchanged: model, p_causal, S, phi, region-length classes, n = 1000.
#
# Encoding:
#   - Each CSV row is ONE array task's SIMULATION (a run_simulation() call).
#   - Within a task, run_benchmark_job.R sweeps the WITHIN_JOB axes (S, phi,
#     n_iter) so one run_simulation() produces all 250 scenarios for that
#     task's (model, annotation regime, enrichment, p_causal, n_ref) combo.
#
# Row count:
#   model x p_causal  : sparse (1) + sparse_inf x 4 p_causal (4)   =  5
#   annotation regime : none (1) + binary x 4 E (4) + cont x 4 E (4) = 9
#   n_ref             : {in-sample, 500, 200}                       = 3
#   -> 5 * 9 * 3 = 135 rows
# Scenarios: 135 * (5 S * 5 phi * 10 iter = 250) = 33,750
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
# Fixed simulation constants
# -----------------------------------------------------------------------------

# Per-region p vector: 2 regions per size class {100, 200, 400, 500, 1000},
# length 10 total. Region length is therefore a within-sim axis with 5 levels
# and 2 regions each. NOT swept between rows.
P_VECTOR <- c(rep(100L, 2L), rep(200L, 2L), rep(400L, 2L),
              rep(500L, 2L), rep(1000L, 2L))

N_ANNOTATIONS <- 10L

# Enrichment is swept. For a given fold E, half the annotations (the first
# N_ANNOTATIONS/2) carry fold E and the rest are inert (fold 1). This builds
# the length-N_ANNOTATIONS vector the simulator expects.
ENRICHMENT_FOLDS <- c(2.7, 5.4, 8.1, 10.8)
N_ENRICHED       <- N_ANNOTATIONS %/% 2L          # 5 of 10
enrichment_vector_for <- function(E) {
  c(rep(E, N_ENRICHED), rep(1, N_ANNOTATIONS - N_ENRICHED))
}
# Placeholder vector for the annotation-free arm (parsed but never used by the
# worker, which passes enrichment = NULL when annotation_type == "none").
ENRICHMENT_NONE <- rep(1, N_ANNOTATIONS)

# Within-job sweep: run_simulation() expands the full grid internally when
# handed vector S / phi values. p_causal is a scalar per task.
WITHIN_JOB <- list(
  n_regions = length(P_VECTOR),     # 10
  n         = 1000L,
  n_iter    = 10L,
  S         = c(1L, 2L, 3L, 5L, 10L),
  phi       = c(0.0075, 0.05, 0.1, 0.2, 0.4)
)

# p_causal only applies to the sparse_inf model.
P_CAUSAL_SPARSE_INF <- c(0.5, 0.7, 0.9, 1.0)

# n_ref: LD reference-panel size. NA = in-sample LD (methods get cor(X) of the
# GWAS sample itself); 500 / 200 = independent panel of that size, so methods
# get cor(X_ref) -- a noisier LD estimate. Smaller panel = more noise.
N_REF_LEVELS <- c(NA_integer_, 500L, 200L)

# Annotation regimes: (type, enrichment fold). correlation is 0 throughout;
# none carries NA fold / NA correlation.
ANNOTATION_REGIMES <- rbind(
  data.frame(annotation_type = "none",       enrichment_fold = NA_real_,
             annotation_correlation = NA_real_, stringsAsFactors = FALSE),
  data.frame(annotation_type = "binary",     enrichment_fold = ENRICHMENT_FOLDS,
             annotation_correlation = 0,        stringsAsFactors = FALSE),
  data.frame(annotation_type = "continuous", enrichment_fold = ENRICHMENT_FOLDS,
             annotation_correlation = 0,        stringsAsFactors = FALSE)
)

# -----------------------------------------------------------------------------
# Build the grid
# -----------------------------------------------------------------------------

nref_token <- function(nr) if (is.na(nr)) "refInsample" else sprintf("ref%d", nr)
arm_token  <- c(none = "anNone", binary = "anBinary", continuous = "anCont")

rows <- list()
job_id <- 0L

for (m in c("sparse", "sparse_inf")) {
  p_causal_values <- if (m == "sparse") NA_real_ else P_CAUSAL_SPARSE_INF

  for (pc in p_causal_values) {
    for (ar in seq_len(nrow(ANNOTATION_REGIMES))) {
      annot_type <- ANNOTATION_REGIMES$annotation_type[ar]
      annot_corr <- ANNOTATION_REGIMES$annotation_correlation[ar]
      E          <- ANNOTATION_REGIMES$enrichment_fold[ar]

      enrich_vec <- if (annot_type == "none") ENRICHMENT_NONE else enrichment_vector_for(E)
      stopifnot(length(enrich_vec) == N_ANNOTATIONS)

      for (nr in N_REF_LEVELS) {
        job_id <- job_id + 1L

        label <- sprintf(
          "%s_%s%s%s_%s",
          m,
          arm_token[[annot_type]],
          if (annot_type == "none") "" else sprintf("_e%g", E),
          if (m == "sparse_inf") sprintf("_pc%g", pc) else "",
          nref_token(nr)
        )

        rows[[job_id]] <- data.frame(
          job_id                 = job_id,
          label                  = label,
          model                  = m,
          p_causal               = pc,
          annotation_type        = annot_type,
          annotation_correlation = annot_corr,
          enrichment_fold        = E,               # scalar fold for readability
          n_ref                  = nr,              # NA = in-sample LD
          n_annotations          = N_ANNOTATIONS,
          n_regions              = WITHIN_JOB$n_regions,
          n                      = WITHIN_JOB$n,
          n_iter                 = WITHIN_JOB$n_iter,
          S_values               = paste(WITHIN_JOB$S,   collapse = "|"),
          phi_values             = paste(WITHIN_JOB$phi, collapse = "|"),
          p_values               = paste(P_VECTOR,       collapse = "|"),
          enrichment_values      = paste(enrich_vec,     collapse = "|"),
          stringsAsFactors       = FALSE
        )
      }
    }
  }
}

grid <- do.call(rbind, rows)

# -----------------------------------------------------------------------------
# Write out
# -----------------------------------------------------------------------------

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(grid, out_path, row.names = FALSE)

n_scenarios_per <- length(WITHIN_JOB$S) * length(WITHIN_JOB$phi) * WITHIN_JOB$n_iter
cat(sprintf("Wrote %d tasks to %s\n", nrow(grid), out_path))
cat(sprintf("  model x p_causal combos : %d\n",
            length(unique(paste(grid$model, grid$p_causal)))))
cat(sprintf("  annotation regimes      : %d (none + binary x4E + continuous x4E)\n",
            nrow(ANNOTATION_REGIMES)))
cat(sprintf("  n_ref levels            : %d {in-sample, 500, 200}\n", length(N_REF_LEVELS)))
cat(sprintf("Scenarios per task: %d (S=%d, phi=%d, n_iter=%d)\n",
            n_scenarios_per, length(WITHIN_JOB$S), length(WITHIN_JOB$phi),
            WITHIN_JOB$n_iter))
cat(sprintf("Total scenarios: %d\n", nrow(grid) * n_scenarios_per))
cat(sprintf("Region-scenario-method fits (14 methods): %s\n",
            format(nrow(grid) * n_scenarios_per * WITHIN_JOB$n_regions * 14L,
                   big.mark = ",")))

cat("\nSanity: annotation-type x n_ref coverage\n")
print(table(grid$annotation_type, grid$n_ref, useNA = "ifany"))
cat("\nFirst rows:\n")
print(grid[seq_len(min(6, nrow(grid))),
           c("job_id", "label", "model", "p_causal",
             "annotation_type", "enrichment_fold", "n_ref")],
      row.names = FALSE)
