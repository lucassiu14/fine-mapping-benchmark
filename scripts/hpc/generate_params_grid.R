# =============================================================================
# scripts/hpc/generate_params_grid.R
#
# Generate the parameter grid CSV used by the SLURM job array.
#
# Each row = one SLURM array task (one HPC job).
# Structure: 2 models × 4 p values × 5 annotation settings = 40 jobs.
# Within each job, run_simulation() sweeps all S × phi (× p_causal) values.
#
# Usage (from project root):
#   Rscript scripts/hpc/generate_params_grid.R
#
# Output: scripts/hpc/params_grid.csv
# =============================================================================

# =============================================================================
# Benchmark design
# =============================================================================

MODELS       <- c("sparse", "sparse_inf")
P_VALUES     <- c(100L, 200L, 400L, 800L)

# S, phi, p_causal are swept inside each job via run_simulation()
S_VALUES     <- c(1L, 2L, 3L, 5L, 8L)
PHI_VALUES   <- c(0.05, 0.1, 0.2, 0.4, 0.6)
P_CAUSAL     <- c(0.1, 0.2, 0.4, 0.8)   # sparse_inf only

N            <- 1000L
N_REGIONS    <- 4L
N_ITER       <- 10L
N_ANNOTATIONS <- 3L

# Five annotation settings:
#   none            : no functional annotations (methods run without annotations)
#   binary_low      : 3 binary annotations, low enrichment  (~2-fold)
#   binary_high     : 3 binary annotations, high enrichment (~10-fold)
#   continuous_low  : 3 continuous annotations, low enrichment (~2-fold)
#   continuous_high : 3 continuous annotations, high enrichment (~10-fold)

ANNOT_SETTINGS <- data.frame(
  annot_name        = c("none", "binary_low", "binary_high",
                        "continuous_low", "continuous_high"),
  annot_type        = c("none",   "binary",    "binary",
                        "continuous", "continuous"),
  annot_enrichment  = c(NA_real_, 2.0,         10.0,
                        2.0,         10.0),
  stringsAsFactors  = FALSE
)

# =============================================================================
# Build grid
# =============================================================================

grid_rows <- vector("list", nrow(ANNOT_SETTINGS) * length(MODELS) * length(P_VALUES))
k <- 0L

for (model in MODELS) {
  for (p in P_VALUES) {
    for (ai in seq_len(nrow(ANNOT_SETTINGS))) {
      k <- k + 1L
      grid_rows[[k]] <- data.frame(
        job_id            = k,
        model             = model,
        n                 = N,
        p                 = p,
        n_regions         = N_REGIONS,
        n_iter            = N_ITER,
        n_annotations     = N_ANNOTATIONS,
        annot_name        = ANNOT_SETTINGS$annot_name[ai],
        annot_type        = ANNOT_SETTINGS$annot_type[ai],
        annot_enrichment  = ANNOT_SETTINGS$annot_enrichment[ai],
        # Embed S/phi/p_causal as semicolon-separated strings for the job to parse
        S_values          = paste(S_VALUES, collapse = ";"),
        phi_values        = paste(PHI_VALUES, collapse = ";"),
        p_causal_values   = paste(P_CAUSAL, collapse = ";"),
        stringsAsFactors  = FALSE
      )
    }
  }
}

grid <- do.call(rbind, grid_rows)

# Sanity check: should be 2 × 4 × 5 = 40 rows
stopifnot(nrow(grid) == 40L)

out_path <- file.path("scripts", "hpc", "params_grid.csv")
write.csv(grid, out_path, row.names = FALSE, quote = FALSE)

cat(sprintf("Parameter grid written to: %s\n", out_path))
cat(sprintf("  Total jobs : %d\n", nrow(grid)))
cat(sprintf("  Models     : %s\n", paste(unique(grid$model), collapse = ", ")))
cat(sprintf("  p values   : %s\n", paste(unique(grid$p), collapse = ", ")))
cat(sprintf("  Annot sets : %s\n", paste(unique(grid$annot_name), collapse = ", ")))
cat(sprintf("\nS values swept per job  : %s\n", paste(S_VALUES, collapse = ", ")))
cat(sprintf("phi values swept per job: %s\n", paste(PHI_VALUES, collapse = ", ")))
cat(sprintf("p_causal (sparse_inf)   : %s\n", paste(P_CAUSAL, collapse = ", ")))
cat(sprintf("\nEstimated method calls per job:\n"))
cat(sprintf("  sparse     : %d S × %d phi × %d iter × %d regions × 6 methods = %d\n",
            length(S_VALUES), length(PHI_VALUES), N_ITER, N_REGIONS,
            length(S_VALUES) * length(PHI_VALUES) * N_ITER * N_REGIONS * 6))
cat(sprintf("  sparse_inf : %d S × %d phi × %d p_causal × %d iter × %d regions × 6 = %d\n",
            length(S_VALUES), length(PHI_VALUES), length(P_CAUSAL), N_ITER, N_REGIONS,
            length(S_VALUES) * length(PHI_VALUES) * length(P_CAUSAL) * N_ITER * N_REGIONS * 6))
