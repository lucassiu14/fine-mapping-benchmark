# =============================================================================
# scripts/hpc/run_benchmark_job.R
#
# SLURM array worker script.
# Run by submit_benchmark.sh as:
#   Rscript scripts/hpc/run_benchmark_job.R <job_id>
#
# <job_id> is the 1-based row index into params_grid.csv (= $SLURM_ARRAY_TASK_ID).
#
# Each job:
#   1. Reads its parameter row from params_grid.csv
#   2. Calls run_simulation() sweeping all S × phi (× p_causal) values
#   3. Calls run_methods() with all 6 methods
#   4. Saves results to results/benchmark/<job_label>/
#
# =============================================================================

# =============================================================================
# Configuration — edit these paths to match your HPC environment
# =============================================================================

PROJECT_ROOT  <- "."        # path to the project root (. if running from there)

# Python executable with BEATRICE dependencies installed
PYTHON        <- "python"   # e.g. "~/miniconda3/envs/beatrice/bin/python"

# Path to BEATRICE_annot_sparse/ (Functional BEATRICE)
FB_DIR        <- file.path(PROJECT_ROOT, "BEATRICE_annot_sparse")

# Path to PAINTOR binary
PAINTOR_PATH  <- "PAINTOR"  # or full path, e.g. "/usr/local/bin/PAINTOR"

# Path to FUNMAP binary / Python
FUNMAP_DIR    <- file.path(PROJECT_ROOT, "alt_methods", "Funmap_main")

# VCF directory (pre-downloaded by prepare_vcfs.R)
VCF_DIR       <- file.path(PROJECT_ROOT, "data", "gwfm_vcf")

# Output root
OUTPUT_DIR    <- file.path(PROJECT_ROOT, "results", "benchmark")

# Methods to run
METHODS <- c("susie", "susie_inf", "beatrice", "funmap",
             "functional_beatrice", "paintor")

# =============================================================================
# Parse job ID from command line
# =============================================================================

args   <- commandArgs(trailingOnly = TRUE)
job_id <- as.integer(args[1])

if (is.na(job_id) || job_id < 1L) {
  stop("Usage: Rscript run_benchmark_job.R <job_id>  (job_id >= 1)", call. = FALSE)
}

# =============================================================================
# Load project (source all R files)
# =============================================================================

r_files <- list.files(file.path(PROJECT_ROOT, "R"),
                       pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
invisible(lapply(r_files, source))

# =============================================================================
# Read parameter grid
# =============================================================================

grid_path <- file.path(PROJECT_ROOT, "scripts", "hpc", "params_grid.csv")
if (!file.exists(grid_path)) {
  stop("params_grid.csv not found. Run scripts/hpc/generate_params_grid.R first.",
       call. = FALSE)
}

grid <- read.csv(grid_path, stringsAsFactors = FALSE)

if (job_id > nrow(grid)) {
  stop(sprintf("job_id %d exceeds grid size %d.", job_id, nrow(grid)), call. = FALSE)
}

params <- grid[job_id, ]

# Parse semicolon-delimited sweep vectors
S_values      <- as.integer(strsplit(params$S_values,     ";")[[1]])
phi_values    <- as.numeric(strsplit(params$phi_values,   ";")[[1]])
p_causal_vals <- as.numeric(strsplit(params$p_causal_values, ";")[[1]])

cat(sprintf("=== Job %d / %d ===\n", job_id, nrow(grid)))
cat(sprintf("  model       : %s\n",  params$model))
cat(sprintf("  p           : %d\n",  params$p))
cat(sprintf("  annot       : %s\n",  params$annot_name))
cat(sprintf("  S values    : %s\n",  paste(S_values, collapse = ", ")))
cat(sprintf("  phi values  : %s\n",  paste(phi_values, collapse = ", ")))
if (params$model == "sparse_inf") {
  cat(sprintf("  p_causal    : %s\n", paste(p_causal_vals, collapse = ", ")))
}
cat("\n")

# =============================================================================
# Derive simulation parameters
# =============================================================================

annot_type       <- params$annot_type        # "none", "binary", "continuous"
annot_enrichment <- params$annot_enrichment  # NA for "none"

enrichment_arg   <- if (is.na(annot_enrichment)) NULL else annot_enrichment
p_causal_arg     <- if (params$model == "sparse_inf") p_causal_vals else NULL

# Job-specific seed for reproducibility
job_seed <- job_id * 1000L

# =============================================================================
# Construct output label and directory
# =============================================================================

job_label <- sprintf("model=%s_p=%d_annot=%s",
                     params$model, params$p, params$annot_name)
job_out   <- file.path(OUTPUT_DIR, job_label)
if (!dir.exists(job_out)) dir.create(job_out, recursive = TRUE)

cat(sprintf("Output directory: %s\n\n", job_out))

# =============================================================================
# Run simulation
# =============================================================================

cat("--- Running simulation ---\n")
sim_start <- proc.time()

sim <- run_simulation(
  n_regions            = params$n_regions,
  n                    = params$n,
  p                    = params$p,
  n_iter               = params$n_iter,
  S                    = S_values,
  phi                  = phi_values,
  model                = params$model,
  p_causal             = p_causal_arg,
  effect_distribution  = "normal",
  effect_variance      = 0.36,
  annotations          = annot_type,
  n_annotations        = params$n_annotations,
  enrichment           = enrichment_arg,
  vcf_dir              = if (dir.exists(VCF_DIR)) VCF_DIR else NULL,
  seed                 = job_seed,
  save                 = FALSE,
  verbose              = TRUE
)

sim_elapsed <- as.numeric((proc.time() - sim_start)["elapsed"])
cat(sprintf("\nSimulation complete in %.1f s.\n", sim_elapsed))
cat(sprintf("  Scenarios: %d\n\n", length(sim$scenarios)))

# =============================================================================
# Build per-method argument lists
# =============================================================================

method_args <- list(
  susie = list(
    L        = max(S_values) + 2L,   # generous upper bound on causal number
    coverage = 0.95
  ),
  susie_inf = list(
    L        = max(S_values) + 2L,
    coverage = 0.95
  ),
  beatrice = list(
    beatrice_dir    = file.path(PROJECT_ROOT, "alt_methods",
                                "Beatrice-Finemapping"),
    python          = PYTHON,
    max_iter        = 2000L,
    n_caus          = max(S_values) + 2L,
    sparse_concrete = min(50L, params$p)
  ),
  funmap = list(
    funmap_dir = FUNMAP_DIR,
    python     = PYTHON
  ),
  functional_beatrice = list(
    beatrice_dir         = FB_DIR,
    python               = PYTHON,
    max_iter             = 2000L,
    n_caus               = max(S_values) + 2L,
    sparse_concrete      = min(50L, params$p),
    prior_regularisation = 1.0,
    lambda_l1            = 0.01,
    hierarchy_M          = 10.0
  ),
  paintor = list(
    paintor_path = PAINTOR_PATH,
    max_causal   = max(S_values),
    coverage     = 0.95
  )
)

# =============================================================================
# Run methods
# =============================================================================

cat("--- Running fine-mapping methods ---\n")
methods_start <- proc.time()

results <- run_methods(
  simulation  = sim,
  methods     = METHODS,
  method_args = method_args,
  save        = FALSE,
  verbose     = TRUE
)

methods_elapsed <- as.numeric((proc.time() - methods_start)["elapsed"])
cat(sprintf("\nMethods complete in %.1f s.\n\n", methods_elapsed))

# =============================================================================
# Evaluate
# =============================================================================

cat("--- Evaluating ---\n")
eval_out <- evaluate_methods(results, sim)

# =============================================================================
# Save outputs
# =============================================================================

# Save full results (can be large — includes all PIPs and credible sets)
saveRDS(results, file.path(job_out, "results.rds"))

# Save evaluation summary (compact)
saveRDS(eval_out, file.path(job_out, "evaluation.rds"))

# Save job metadata
metadata <- list(
  job_id          = job_id,
  job_label       = job_label,
  params          = params,
  S_values        = S_values,
  phi_values      = phi_values,
  p_causal_vals   = p_causal_vals,
  methods_run     = METHODS,
  sim_elapsed_s   = sim_elapsed,
  method_elapsed_s = methods_elapsed,
  n_scenarios     = length(sim$scenarios),
  timestamp       = Sys.time()
)
saveRDS(metadata, file.path(job_out, "metadata.rds"))

cat(sprintf("Results saved to: %s\n", job_out))
cat(sprintf("Total job time: %.1f s\n",
            as.numeric((proc.time() - sim_start)["elapsed"])))
