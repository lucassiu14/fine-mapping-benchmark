#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/make_test_fixture.R
#
# Build a SYNTHETIC benchmark output tree in exactly the format the real run
# produces, so the analysis pipeline can be exercised end to end without waiting
# for the HPC array and without touching its output.
#
# Mirrors, in particular, the two structural details that are easy to get wrong
# and that would silently corrupt the collect:
#
#   1. Every saved fit carries scenario_id = 1, NOT the global scenario index.
#      run_benchmark_job.R subsets sim to one scenario and resets it, because
#      evaluate_methods() uses scenario_id as a LIST INDEX into the subset. The
#      true index survives only in the scenario_<sc> directory name.
#   2. results.rds carries non-method top-level entries (methods_run,
#      run_timestamp, simulation_params) alongside the per-method lists.
#
# A fixture that got either of those wrong would let a broken collect pass.
#
#   Rscript scripts/analysis/make_test_fixture.R <out_root>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
root <- args[1] %||% "/tmp/fmb_fixture"
`%||%` <- function(x, y) if (is.null(x)) y else x
set.seed(404)

bench <- file.path(root, "results", "benchmark")
dir.create(bench, recursive = TRUE, showWarnings = FALSE)

# Small but structurally complete: 2 model x annotation combinations, a real
# S x phi x iter sweep, and 2 regions per size class so region_idx and hence
# sigma^2_u are identifiable.
S_LV    <- c(1L, 3L)
PHI_LV  <- c(0.1, 0.4)
N_ITER  <- 10L
P_VEC   <- c(60L, 60L, 120L, 120L)          # 2 regions x 2 size classes
METHODS <- c("susie", "polyfun_oracle", "beatrice", "functional_beatrice",
             "fb_xregion", "funmap")

ROWS <- list(
  list(job_id = 1L, label = "sparse_anNone_refInsample",
       model = "sparse", p_causal = NA_real_, annotation_type = "none",
       enrichment_fold = NA_real_),
  list(job_id = 2L, label = "sparse_anBinary_e5.4_refInsample",
       model = "sparse", p_causal = NA_real_, annotation_type = "binary",
       enrichment_fold = 5.4),
  list(job_id = 3L, label = "sparse_anBinary_e10.8_refInsample",
       model = "sparse", p_causal = NA_real_, annotation_type = "binary",
       enrichment_fold = 10.8)
)

grid <- do.call(rbind, lapply(ROWS, function(r) data.frame(
  job_id = r$job_id, label = r$label, model = r$model, p_causal = r$p_causal,
  annotation_type = r$annotation_type, enrichment_fold = r$enrichment_fold,
  annotation_correlation = 0, n_ref = NA_integer_, n_annotations = 10L,
  n_regions = length(P_VEC), n = 5000L, n_iter = N_ITER,
  S_values = paste(S_LV, collapse = "|"),
  phi_values = paste(PHI_LV, collapse = "|"),
  p_values = paste(P_VEC, collapse = "|"),
  stringsAsFactors = FALSE)))
grid_csv <- file.path(root, "params_grid.csv")
write.csv(grid, grid_csv, row.names = FALSE)

# Method "skill": a known ordering so the analysis has real structure to find.
# susie is the annotation-blind baseline, polyfun_oracle the ceiling, and the FB
# family sits between - so kappa is well defined and bounded.
SKILL <- c(susie = 0.00, polyfun_oracle = 0.45, beatrice = 0.05,
           functional_beatrice = 0.20, fb_xregion = 0.28, funmap = 0.18)

for (r in ROWS) {
  jd <- file.path(bench, sprintf("job_%03d_%s", r$job_id, r$label))
  dir.create(jd, recursive = TRUE, showWarnings = FALSE)

  # expand.grid(S, phi, iter): S varies FASTEST, matching run_simulation().
  pg <- expand.grid(S = S_LV, phi = PHI_LV, iter = seq_len(N_ITER))
  scenarios <- vector("list", nrow(pg))
  # One region draw per region, shared across ALL scenarios of the job - this is
  # what creates sigma^2_u and makes the split-plot structure real.
  region_effect <- rnorm(length(P_VEC), 0, 0.05)

  for (sc in seq_len(nrow(pg))) {
    regs <- lapply(seq_along(P_VEC), function(rg) {
      p <- P_VEC[rg]
      list(truth = list(causal_indices = sort(sample.int(p, pg$S[sc])),
                        S = pg$S[sc], phi = pg$phi[sc], model = r$model))
    })
    scenarios[[sc]] <- list(scenario_id = sc, S = pg$S[sc], phi = pg$phi[sc],
                            p_causal = r$p_causal, iter = pg$iter[sc],
                            model = r$model, regions = regs)
  }
  saveRDS(list(scenarios = scenarios,
               genotypes = lapply(P_VEC, function(p) list(maf = runif(p, 0.05, 0.5))),
               params = list(p = P_VEC, n = 5000L)),
          file.path(jd, "sim.rds"))

  for (sc in seq_len(nrow(pg))) {
    sd <- file.path(jd, sprintf("scenario_%03d", sc))
    dir.create(sd, recursive = TRUE, showWarnings = FALSE)
    res <- list()

    for (m in METHODS) {
      # funmap CANNOT run without annotations - 100% NA on the none arm,
      # structurally. The fixture must reproduce this or the pipeline's handling
      # of a legitimately absent method goes untested.
      dead <- (m == "funmap" && r$annotation_type == "none")

      fits <- lapply(seq_along(P_VEC), function(rg) {
        p  <- P_VEC[rg]
        ci <- scenarios[[sc]]$regions[[rg]]$truth$causal_indices
        if (dead) {
          return(list(pip = rep(NA_real_, p), scenario_id = 1L, region_id = rg,
                      method = m, error = "funmap requires annotations"))
        }
        # Signal grows with phi and skill, shrinks with S and region size, and
        # carries the region's shared draw.
        sig <- (0.35 + SKILL[[m]]) * pg$phi[sc] * 6 /
               (1 + 0.35 * pg$S[sc]) * (60 / p)^0.3 + region_effect[rg]
        z <- rnorm(p, 0, 1); z[ci] <- z[ci] + max(sig, 0.05) * 6
        pip <- exp(z) / sum(exp(z))
        list(pip = pip, scenario_id = 1L, region_id = rg, method = m)
      })
      res[[m]] <- list(results = fits, n_total = length(P_VEC),
                       n_failed = if (dead) length(P_VEC) else 0L,
                       method_args = list(), total_runtime_seconds = 1)
    }
    # Non-method top-level entries, exactly as the real worker writes them.
    res$methods_run       <- METHODS
    res$run_timestamp     <- "2026-08-06 00:00:00"
    res$simulation_params <- list(n = 5000L, p = P_VEC)

    saveRDS(res, file.path(sd, "results.rds"))
    saveRDS(list(), file.path(sd, "evaluation.rds"))
  }
}

cat(sprintf("fixture: %d rows x %d scenarios x %d regions x %d methods\n",
            length(ROWS), N_ITER * length(S_LV) * length(PHI_LV),
            length(P_VEC), length(METHODS)))
cat("bench root: ", bench, "\n", sep = "")
cat("grid csv  : ", grid_csv, "\n", sep = "")
