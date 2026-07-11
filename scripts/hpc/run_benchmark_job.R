#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/run_benchmark_job.R
#
# SLURM array worker. Run with:
#   Rscript scripts/hpc/run_benchmark_job.R <job_id>
#
# <job_id> is the 1-based row of params_grid.csv (= $SLURM_ARRAY_TASK_ID).
#
# Each invocation:
#   1. Reads its parameter row from params_grid.csv
#   2. Runs simulate -> run_methods -> evaluate_methods
#   3. Saves results to results/benchmark/job_<id>_<label>/
# =============================================================================

# -----------------------------------------------------------------------------
# Configuration - edit to match your environment
# -----------------------------------------------------------------------------
PARAMS_CSV  <- "scripts/hpc/params_grid.csv"
OUTPUT_ROOT <- "results/benchmark"
VCF_DIR     <- "data/vcf"          # per-locus benchmark; set to NULL to use
                                    # the bundled chr4 sim1000G example VCF
GENETIC_MAP_DIR <- "data/genetic_maps"

# Methods. Full 15-method benchmark set: 9 Tier-1 (R-only) + 6 Tier-3
# (external binaries via FINEMAP/PAINTOR, and Python via a venv).
# Any missing binary degrades to NA fits without breaking the array.
METHODS <- c("susie", "susie_inf", "abf", "carma",
             "marginal_z", "polyfun_oracle", "polyfun_est",
             "polyfun_ldsc", "sbayesrc",
             "finemap", "paintor", "beatrice",
             "functional_beatrice", "sparsepro", "funmap")

# Tool locations - edit if you move ~/tools to project space.
TOOLS_ROOT <- normalizePath("~/tools", mustWork = FALSE)
PY_VENV    <- file.path(TOOLS_ROOT, "py-venv-runner.sh")

# BEATRICE and Functional BEATRICE both point at the in-repo
# BEATRICE_annot_sparse/ fork. It ships a numpy-2.x-safe calculate_pip and
# accepts --annot None, so it serves as vanilla BEATRICE too. The upstream
# sayangsep/Beatrice-Finemapping repo has a late-training crash under
# numpy>=2 and is intentionally not used.
FB_DIR <- normalizePath("BEATRICE_annot_sparse", mustWork = FALSE)

METHOD_ARGS <- list(
  susie               = list(L = 10, coverage = 0.95),
  susie_inf           = list(L = 10),
  abf                 = list(prior_variance = 0.04),
  carma               = list(num.causal = 5),
  marginal_z          = list(coverage = 0.95),
  polyfun_oracle      = list(L = 10),
  polyfun_est         = list(L = 10),
  polyfun_ldsc        = list(L = 10),
  sbayesrc            = list(n_iter = 300L, burn_in = 150L,
                             gamma_update_every = 10L),
  finemap             = list(finemap_path = file.path(TOOLS_ROOT,
                                "finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64"),
                             n_causal = 5, n_iter = 100000,
                             prior_std = 0.05, coverage = 0.95),
  paintor             = list(paintor_path = file.path(TOOLS_ROOT,
                                "PAINTOR_V3.0/PAINTOR"),
                             max_causal = 2, mcmc = FALSE, coverage = 0.95),
  # BEATRICE / FB: max_iter = 500 (the wrapper's minimum) - critical for
  # Phase 1 wall-clock. At p = 1000 each fit is O(minutes); the grid has
  # ~2500 BEATRICE fits per task, so max_iter = 2000 (the previous value)
  # would multiply per-task runtime by ~4x. 500 iters is what BEATRICE's
  # own smoke tests use.
  beatrice            = list(beatrice_dir = FB_DIR,
                             python = PY_VENV, max_iter = 500, n_caus = 5,
                             sigma_sq = 0.05, gamma_coverage = 0.95,
                             sparse_concrete = 50),
  functional_beatrice = list(beatrice_dir = FB_DIR,
                             python = PY_VENV, max_iter = 500, n_caus = 5,
                             sigma_sq = 0.05, gamma_coverage = 0.95,
                             sparse_concrete = 50,
                             prior_regularisation = 1.0),
  sparsepro           = list(sparsepro_dir = file.path(TOOLS_ROOT, "SparsePro"),
                             python = PY_VENV, cthres = 0.95),
  funmap              = list(python = PY_VENV, max_iter = 100, tol = 5e-5)
)

# -----------------------------------------------------------------------------
# Parse job_id
# -----------------------------------------------------------------------------
args   <- commandArgs(trailingOnly = TRUE)
job_id <- suppressWarnings(as.integer(args[1]))
if (is.na(job_id) || job_id < 1L) {
  stop("Usage: Rscript scripts/hpc/run_benchmark_job.R <job_id>")
}

# -----------------------------------------------------------------------------
# Load the package
# -----------------------------------------------------------------------------
ok <- requireNamespace("fmbenchmark", quietly = TRUE)
if (!ok) {
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(quiet = TRUE)
  } else {
    stop("fmbenchmark is not installed. Run from the project root after:\n",
         "  R -e 'renv::restore(); install.packages(\".\", repos = NULL, type = \"source\")'")
  }
} else {
  library(fmbenchmark)
}

# -----------------------------------------------------------------------------
# Read this job's params row
# -----------------------------------------------------------------------------
if (!file.exists(PARAMS_CSV)) {
  stop("params_grid.csv not found at ", PARAMS_CSV,
       ". Run: Rscript scripts/hpc/generate_params_grid.R")
}
grid <- read.csv(PARAMS_CSV, stringsAsFactors = FALSE)
if (job_id > nrow(grid)) {
  stop(sprintf("job_id %d is past end of grid (nrow = %d)", job_id, nrow(grid)))
}
job <- as.list(grid[job_id, ])

S_vec         <- as.integer(strsplit(job$S_values,          "\\|")[[1]])
phi_vec       <- as.numeric(strsplit(job$phi_values,        "\\|")[[1]])
p_vec         <- as.integer(strsplit(job$p_values,          "\\|")[[1]])
enrichment_vec <- as.numeric(strsplit(job$enrichment_values, "\\|")[[1]])

job_dir <- file.path(OUTPUT_ROOT, sprintf("job_%03d_%s", job$job_id, job$label))
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

# Logs go to SLURM's --output / --error files (see submit_benchmark.sh).
# No per-job-dir log file: avoids on.exit-at-top-level pitfalls and keeps
# log management in one place.
cat("=====================================================================\n")
cat(sprintf("Job %d / %d: %s\n", job$job_id, nrow(grid), job$label))
cat("=====================================================================\n")
cat(sprintf("Started:        %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Output dir:     %s\n", job_dir))
cat(sprintf("Model:          %s\n", job$model))
if (identical(job$model, "sparse_inf") && !is.null(job$p_causal) && !is.na(job$p_causal)) {
  cat(sprintf("p_causal:       %g\n", job$p_causal))
}
cat(sprintf("n_regions:      %d\n", job$n_regions))
cat(sprintf("p (per region): %s\n", paste(p_vec, collapse = ",")))
cat(sprintf("n:              %d\n", job$n))
cat(sprintf("n_iter:         %d\n", job$n_iter))
cat(sprintf("S values:       %s\n", paste(S_vec,   collapse = ", ")))
cat(sprintf("phi values:     %s\n", paste(phi_vec, collapse = ", ")))
cat(sprintf("Annotation:     %s\n", job$annotation_type))
if (job$annotation_type != "none") {
  cat(sprintf("  tracks:       %d\n", job$n_annotations))
  cat(sprintf("  enrichment:   %s\n", paste(enrichment_vec, collapse = ",")))
  cat(sprintf("  correlation:  %s\n",
              if (is.null(job$annotation_correlation) ||
                  is.na(job$annotation_correlation)) "NA"
              else format(job$annotation_correlation, nsmall = 2)))
}
# LD mismatch panel (n_ref) is not swept in Phase 1; kept for backward compat
# with older grids that still carry the column.
if (!is.null(job$n_ref) && !is.na(job$n_ref)) {
  cat(sprintf("LD mismatch:    n_ref = %d (independent reference panel)\n", job$n_ref))
}
cat(sprintf("Methods:        %s\n", paste(METHODS, collapse = ", ")))
cat("\n")

# Reproducible seed: derived from job_id so jobs differ but each is deterministic.
seed <- 1000L + job$job_id

# -----------------------------------------------------------------------------
# Pipeline
# -----------------------------------------------------------------------------
t0 <- Sys.time()

cat("[1/3] Simulating ... ")

# p_causal is a scalar in the grid; run_simulation() ignores it for the
# "sparse" model. NA sentinels in the CSV become NA in R - only forward a
# real value.
p_causal_arg <- if (identical(job$model, "sparse_inf") &&
                    !is.null(job$p_causal) && !is.na(job$p_causal))
  as.numeric(job$p_causal) else NULL

annotation_corr_arg <- if (!is.null(job$annotation_correlation) &&
                           !is.na(job$annotation_correlation))
  as.numeric(job$annotation_correlation) else 0

sim_args <- list(
  n_regions              = as.integer(job$n_regions),
  n                      = as.integer(job$n),
  p                      = p_vec,
  n_iter                 = as.integer(job$n_iter),
  S                      = S_vec,
  phi                    = phi_vec,
  model                  = job$model,
  annotations            = job$annotation_type,
  n_annotations          = if (is.null(job$n_annotations) ||
                               is.na(job$n_annotations)) 0L
                           else as.integer(job$n_annotations),
  enrichment             = if (identical(job$annotation_type, "none")) NULL
                           else enrichment_vec,
  annotation_correlation = annotation_corr_arg,
  vcf_dir                = if (is.null(VCF_DIR) || !nzchar(VCF_DIR)) NULL else VCF_DIR,
  genetic_map_dir        = GENETIC_MAP_DIR,
  n_ref                  = if (is.null(job$n_ref) || is.na(job$n_ref)) NULL
                           else as.integer(job$n_ref),
  seed                   = seed,
  save                   = FALSE,
  verbose                = FALSE
)
if (!is.null(p_causal_arg)) sim_args$p_causal <- p_causal_arg
sim <- do.call(run_simulation, sim_args)
t1 <- Sys.time()
cat(sprintf("done (%.1f min)\n", as.numeric(t1 - t0, units = "mins")))
saveRDS(sim, file.path(job_dir, "sim.rds"))

cat(sprintf("[2/3] Running %d methods ... \n", length(METHODS)))
t1 <- Sys.time()
results <- run_methods(
  sim,
  methods     = METHODS,
  method_args = METHOD_ARGS[intersect(names(METHOD_ARGS), METHODS)],
  save        = FALSE,
  verbose     = FALSE
)
t2 <- Sys.time()
for (m in METHODS) {
  cat(sprintf("       %-20s n_fits=%d  failed=%d\n",
              m, results[[m]]$n_total, results[[m]]$n_failed))
}
cat(sprintf("       Methods total: %.1f min\n",
            as.numeric(t2 - t1, units = "mins")))
saveRDS(results, file.path(job_dir, "results.rds"))

cat("[3/3] Evaluating ... ")
t2 <- Sys.time()
evaluation <- evaluate_methods(sim, results, save = FALSE, verbose = FALSE)
t3 <- Sys.time()
cat(sprintf("done (%.1f s)\n", as.numeric(t3 - t2, units = "secs")))
saveRDS(evaluation, file.path(job_dir, "evaluation.rds"))

# Headline AUPRC per method
cat("\nGlobal AUPRC by method:\n")
for (m in METHODS) {
  a <- evaluation[[m]]$global$auprc
  cat(sprintf("  %-20s %s\n", m,
              if (is.null(a) || is.na(a)) "  NA" else sprintf("%.3f", a)))
}

# Reproducibility metadata
meta <- list(
  job_id      = job$job_id,
  label       = job$label,
  params_row  = job,
  methods     = METHODS,
  method_args = METHOD_ARGS[intersect(names(METHOD_ARGS), METHODS)],
  seed        = seed,
  vcf_dir     = VCF_DIR,
  genetic_map_dir = GENETIC_MAP_DIR,
  started_at  = format(t0, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  runtime_sec = as.numeric(Sys.time() - t0, units = "secs"),
  R_version   = R.version.string,
  fmbenchmark_version =
    tryCatch(as.character(utils::packageVersion("fmbenchmark")),
             error = function(e) "source-checkout"),
  sessionInfo = capture.output(sessionInfo())
)
jsonlite::write_json(meta, file.path(job_dir, "params.json"),
                     auto_unbox = TRUE, pretty = TRUE, null = "null")

cat(sprintf("\nDone in %.1f min. Outputs at:\n  %s\n",
            as.numeric(Sys.time() - t0, units = "mins"), job_dir))
