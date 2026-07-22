#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/run_benchmark_job.R
#
# Array worker. Run with:
#   Rscript scripts/hpc/run_benchmark_job.R <array_index>
#
# The array is subdivided BY SCENARIO, keeping all regions together in each
# task. A "scenario" is one (S, phi, iter) draw; there are SCENARIOS_PER_ROW
# (= |S| * |phi| * n_iter) of them per grid row, and every scenario contains
# ALL n_regions regions. Splitting by scenario (never by region) is essential:
# methods like SBayesRC, polyfun_ldsc and the cross-fitted-prior candidates
# (plan Idea A / Idea C) POOL annotation + summary-stat information across all
# regions inside run_<method>_scenario_setup(). That pooling only works if a
# task holds the whole region panel. Scenarios, by contrast, are statistically
# independent, so they are the correct unit to parallelise and checkpoint over
# (plan lines: "the shared prior is a single artifact written once per
# scenario"; "you cannot remove specific regions of an iteration").
#
# Each task processes a CONTIGUOUS BLOCK of FMB_SCENARIOS_PER_TASK scenarios
# from ONE row (all from the same sim, so it is simulated/loaded once and
# amortised over the block). With TASKS_PER_ROW = ceil(SCENARIOS_PER_ROW /
# FMB_SCENARIOS_PER_TASK), array_index i (1-based) maps to:
#   row         = ((i - 1) %/% TASKS_PER_ROW) + 1
#   task_in_row = ((i - 1) %%  TASKS_PER_ROW)
#   scenarios   = [task_in_row*chunk + 1 .. min(+chunk, SCENARIOS_PER_ROW)]
# giving nrow(grid) * TASKS_PER_ROW tasks. Iteration 002 (135 rows x 250
# scenarios) at chunk=25 -> 10 tasks/row -> 1350 tasks, each ~10-25 scenarios
# well inside a 72h walltime. chunk=1 recovers one-task-per-scenario (3125 in
# Iteration 001). Every scenario is still checkpointed and skipped on resume, so
# a walltime kill loses only the scenario in flight - which is exactly what
# makes a fat, few-task array safe.
#
# Each invocation:
#   1. Decodes (row, scenario block) from the array index.
#   2. Simulates the row's full region panel once (cached per row in sim.rds).
#   3. For each scenario in the block: runs run_methods on that ONE scenario
#      (all regions -> cross-region pooling intact), evaluates, and saves to
#      <OUTPUT_ROOT>/job_<row>_<label>/scenario_<sc>/ before starting the next.
# =============================================================================

# -----------------------------------------------------------------------------
# Configuration - edit to match your environment
# -----------------------------------------------------------------------------
PARAMS_CSV  <- "scripts/hpc/params_grid.csv"
# OUTPUT_ROOT holds per-task sim/results/evaluation RDS files (~100 MB+
# each). It MUST live on a filesystem with room - on a 1 TB home quota the
# full array overflows and saveRDS() fails with "error writing to
# connection". submit_benchmark_pbs.sh sets FMB_OUTPUT_ROOT to a scratch
# path; the "results/benchmark" fallback is for laptop / small runs only.
OUTPUT_ROOT <- Sys.getenv("FMB_OUTPUT_ROOT", unset = "results/benchmark")
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

# Supplemental mode: FMB_METHODS="beatrice,sparsepro" re-runs ONLY those
# methods and writes results_supp.rds / evaluation_supp.rds alongside the
# originals, leaving the existing (good) results untouched. Used to re-run
# a method after a wrapper fix without redoing the whole 15-method grid -
# the per-row sim.rds cache means no re-simulation either. collect_results.R
# prefers _supp results for any method present in both.
SUPP_METHODS <- Sys.getenv("FMB_METHODS", unset = "")
IS_SUPP      <- nzchar(SUPP_METHODS)
if (IS_SUPP) {
  METHODS <- trimws(strsplit(SUPP_METHODS, ",")[[1]])
  METHODS <- METHODS[nzchar(METHODS)]
}
OUT_SUFFIX <- if (IS_SUPP) "_supp" else ""

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
  # BEATRICE / FB: max_iter = 1500 balances convergence against per-task
  # wall-clock. 500 iters is not enough for the variational objective to
  # stabilise; 1500 matches BEATRICE's paper experiments. A scenario runs
  # 20 regions x 2 BEATRICE-family fits; these dominate the per-task time.
  beatrice            = list(beatrice_dir = FB_DIR,
                             python = PY_VENV, max_iter = 1500, n_caus = 5,
                             sigma_sq = 0.05, gamma_coverage = 0.95,
                             sparse_concrete = 50),
  functional_beatrice = list(beatrice_dir = FB_DIR,
                             python = PY_VENV, max_iter = 1500, n_caus = 5,
                             sigma_sq = 0.05, gamma_coverage = 0.95,
                             sparse_concrete = 50,
                             prior_regularisation = 1.0),
  sparsepro           = list(sparsepro_dir = file.path(TOOLS_ROOT, "SparsePro"),
                             python = PY_VENV, cthres = 0.95),
  funmap              = list(python = PY_VENV, max_iter = 100, tol = 5e-5)
)

# --- Iteration 003: FB / BEATRICE hyperparameter variants --------------------
# Reuse the tested FB / BEATRICE wrappers with one knob changed each. Run
# supplementally (FMB_METHODS="fb_l1hi,fb_prreg5,...") against the cached row
# sims so the 14 baselines (already collected) are not re-run.
.FBv <- METHOD_ARGS$functional_beatrice
.BEv <- METHOD_ARGS$beatrice
METHOD_ARGS$fb_l1hi           <- modifyList(.FBv, list(lambda_l1 = 0.1))          # feature sparsity x10
METHOD_ARGS$fb_l1vhi          <- modifyList(.FBv, list(lambda_l1 = 0.5))          # feature sparsity x50
METHOD_ARGS$fb_prreg5         <- modifyList(.FBv, list(prior_regularisation = 5.0))   # shrink toward BEATRICE
METHOD_ARGS$fb_prreg20        <- modifyList(.FBv, list(prior_regularisation = 20.0))  # strong shrink
METHOD_ARGS$fb_ncaus2         <- modifyList(.FBv, list(n_caus = 2))               # less mass inflation
METHOD_ARGS$fb_concrete       <- modifyList(.FBv, list(sparse_concrete = 200))    # softer posteriors
METHOD_ARGS$fb_sigma_hi       <- modifyList(.FBv, list(sigma_sq = 0.2))           # larger effect prior
METHOD_ARGS$fb_reg_combo      <- modifyList(.FBv, list(lambda_l1 = 0.1,
                                            prior_regularisation = 5.0, n_caus = 3))  # combined
METHOD_ARGS$beatrice_ncaus2   <- modifyList(.BEv, list(n_caus = 2))
METHOD_ARGS$beatrice_sigma_hi <- modifyList(.BEv, list(sigma_sq = 0.2))

# --- Iteration 003 Track B: cross-region JOINT prior models ------------------
# Same FB base args (beatrice_dir, python, max_iter, n_caus, ...). The joint
# training is driven by each method's scenario_setup hook (run_fb_*_scenario_setup
# in R/wrapper_fb_joint.R), which invokes BEATRICE_annot_sparse/beatrice_joint.py
# once per scenario. fb_pooled = shared logistic head (idea #1); fb_xregion =
# shared LassoNet head (idea #2). On the `none` arm they fall back to plain FB.
METHOD_ARGS$fb_pooled         <- .FBv
METHOD_ARGS$fb_xregion        <- .FBv

# -----------------------------------------------------------------------------
# Parse array index
# -----------------------------------------------------------------------------
args      <- commandArgs(trailingOnly = TRUE)
array_idx <- suppressWarnings(as.integer(args[1]))
if (is.na(array_idx) || array_idx < 1L) {
  stop("Usage: Rscript scripts/hpc/run_benchmark_job.R <array_index>")
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

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# -----------------------------------------------------------------------------
# Read grid + decode this task's (row, scenario)
# -----------------------------------------------------------------------------
if (!file.exists(PARAMS_CSV)) {
  stop("params_grid.csv not found at ", PARAMS_CSV,
       ". Run: Rscript scripts/hpc/generate_params_grid.R")
}
grid <- read.csv(PARAMS_CSV, stringsAsFactors = FALSE)

# Scenarios per row = |S| * |phi| * n_iter, identical for every row. Derive it
# from the grid so the decode and the submit-side array size always agree.
scenarios_for_row <- function(r) {
  length(strsplit(r$S_values,   "\\|")[[1]]) *
  length(strsplit(r$phi_values, "\\|")[[1]]) *
  as.integer(r$n_iter)
}
per_row <- vapply(seq_len(nrow(grid)),
                  function(i) scenarios_for_row(grid[i, ]), integer(1))
if (length(unique(per_row)) != 1L) {
  stop("All grid rows must yield the same scenario count for the scenario-split array.")
}
SCENARIOS_PER_ROW <- per_row[1]

# Scenario chunking. Each array task processes a CONTIGUOUS BLOCK of
# SCENARIOS_PER_TASK scenarios from ONE row, so it simulates/loads that row's
# sim.rds just once and amortises it over the whole block. Every scenario is
# still checkpointed individually (its own scenario_<sc>/evaluation.rds) and
# skipped on resume, so a walltime kill costs at most the single scenario in
# flight - the guarantee that lets us pack many scenarios per task and raise
# walltime rather than run 33,750 one-scenario tasks past the array cap. The
# worker and submit script MUST agree on this value: the submit script exports
# FMB_SCENARIOS_PER_TASK into the job. Set to 1 to recover one-task-per-scenario.
SCENARIOS_PER_TASK <- suppressWarnings(as.integer(
  Sys.getenv("FMB_SCENARIOS_PER_TASK", unset = "25")))
if (is.na(SCENARIOS_PER_TASK) || SCENARIOS_PER_TASK < 1L) SCENARIOS_PER_TASK <- 1L
SCENARIOS_PER_TASK <- min(SCENARIOS_PER_TASK, SCENARIOS_PER_ROW)
TASKS_PER_ROW <- as.integer(ceiling(SCENARIOS_PER_ROW / SCENARIOS_PER_TASK))
n_tasks <- nrow(grid) * TASKS_PER_ROW
if (array_idx > n_tasks) {
  stop(sprintf("array_index %d past end (%d rows x %d tasks/row = %d tasks)",
               array_idx, nrow(grid), TASKS_PER_ROW, n_tasks))
}
row_idx        <- ((array_idx - 1L) %/% TASKS_PER_ROW) + 1L
task_in_row    <-  (array_idx - 1L) %%  TASKS_PER_ROW           # 0-based
scenario_start <- task_in_row * SCENARIOS_PER_TASK + 1L
scenario_end   <- min(scenario_start + SCENARIOS_PER_TASK - 1L, SCENARIOS_PER_ROW)
scenario_ids   <- scenario_start:scenario_end

job <- as.list(grid[row_idx, ])

S_vec          <- as.integer(strsplit(job$S_values,          "\\|")[[1]])
phi_vec        <- as.numeric(strsplit(job$phi_values,        "\\|")[[1]])
p_vec          <- as.integer(strsplit(job$p_values,          "\\|")[[1]])
enrichment_vec <- as.numeric(strsplit(job$enrichment_values, "\\|")[[1]])

# job_<row>/ holds the shared per-row sim.rds; each scenario_<sc>/ holds one
# scenario's results. This task owns scenarios scenario_start..scenario_end.
job_dir <- file.path(OUTPUT_ROOT, sprintf("job_%03d_%s", job$job_id, job$label))
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

cat("=====================================================================\n")
cat(sprintf("Task %d  ->  row %d (%s), scenarios %d-%d of %d  (chunk=%d)\n",
            array_idx, row_idx, job$label, scenario_start, scenario_end,
            SCENARIOS_PER_ROW, SCENARIOS_PER_TASK))
cat("=====================================================================\n")
cat(sprintf("Started:  %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Out dir:  %s\n", job_dir))
cat(sprintf("Model:    %s%s\n", job$model,
            if (identical(job$model, "sparse_inf") &&
                !is.na(job$p_causal)) sprintf("  p_causal=%g", job$p_causal) else ""))
cat(sprintf("Regions:  %d  (p = %s)\n", job$n_regions, paste(p_vec, collapse = ",")))
cat(sprintf("Annot:    %s%s\n", job$annotation_type,
            if (job$annotation_type != "none")
              sprintf("  (%s tracks, corr=%s, enrich_fold=%s)", job$n_annotations,
                      if (is.na(job$annotation_correlation)) "NA"
                      else format(job$annotation_correlation, nsmall = 2),
                      if (is.null(job$enrichment_fold) || is.na(job$enrichment_fold)) "NA"
                      else format(job$enrichment_fold)) else ""))
cat(sprintf("LD:       %s\n",
            if (is.null(job$n_ref) || is.na(job$n_ref)) "in-sample (perfect)"
            else sprintf("reference panel n_ref=%d", as.integer(job$n_ref))))
cat(sprintf("Methods:  %s\n", paste(METHODS, collapse = ", ")))
if (IS_SUPP) cat(sprintf("SUPPLEMENTAL run: methods = %s -> evaluation%s.rds\n",
                         paste(METHODS, collapse = ", "), OUT_SUFFIX))
cat("\n")

# Scenario-level resume. Determine which of this task's scenarios still need
# doing BEFORE touching the sim: if the whole block is already checkpointed we
# can exit without paying the (1-3h) simulation/load cost.
scenario_eval_file <- function(sc)
  file.path(job_dir, sprintf("scenario_%03d", sc),
            paste0("evaluation", OUT_SUFFIX, ".rds"))
todo <- scenario_ids[!vapply(scenario_ids,
                             function(sc) file.exists(scenario_eval_file(sc)),
                             logical(1))]
if (length(todo) == 0L) {
  cat(sprintf("All %d scenarios (%d-%d) already complete. Nothing to do.\n",
              length(scenario_ids), scenario_start, scenario_end))
  quit(save = "no", status = 0)
}
cat(sprintf("Scenarios to run this task: %d of %d (rest already checkpointed)\n\n",
            length(todo), length(scenario_ids)))

# -----------------------------------------------------------------------------
# Pipeline
# -----------------------------------------------------------------------------
t0 <- Sys.time()

# --- [1/3] Simulate the row's FULL region panel (cached per row) ------------
# All scenarios of a row share the same genotypes + the same deterministic
# per-scenario phenotypes (seed = 1000 + row). Simulating the whole row and
# then subsetting to one scenario keeps the region panel intact for pooling.
# The sim.rds is shared across all of the row's tasks; the atomic write
# (tmp + rename) makes concurrent first-run writers safe, and resubmits skip
# re-simulation entirely.
sim_file <- file.path(job_dir, "sim.rds")
if (file.exists(sim_file)) {
  cat("[1/3] loading cached row sim.rds\n")
  sim <- readRDS(sim_file)
} else {
  cat("[1/3] simulating row region panel ... ")

  p_causal_arg <- if (identical(job$model, "sparse_inf") &&
                      !is.null(job$p_causal) && !is.na(job$p_causal))
    as.numeric(job$p_causal) else NULL
  annotation_corr_arg <- if (!is.null(job$annotation_correlation) &&
                             !is.na(job$annotation_correlation))
    as.numeric(job$annotation_correlation) else 0
  # n_ref: LD reference-panel size. Absent column or NA -> NULL (in-sample LD),
  # preserving Iteration 001 behaviour. Integer -> methods receive cor(X_ref).
  n_ref_arg <- if (!is.null(job$n_ref) && !is.na(job$n_ref))
    as.integer(job$n_ref) else NULL

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
    seed                   = 1000L + row_idx,
    save                   = FALSE,
    verbose                = FALSE
  )
  if (!is.null(p_causal_arg)) sim_args$p_causal <- p_causal_arg
  if (!is.null(n_ref_arg))    sim_args$n_ref    <- n_ref_arg
  sim <- do.call(run_simulation, sim_args)

  tmp <- paste0(sim_file, ".tmp.", Sys.getpid())
  saveRDS(sim, tmp)
  file.rename(tmp, sim_file)   # atomic on the same filesystem
  cat(sprintf("done (%.1f min)\n", as.numeric(Sys.time() - t0, units = "mins")))
}

n_sc_total <- length(sim$scenarios)
if (scenario_end > n_sc_total) {
  stop(sprintf("scenario_end %d exceeds simulated scenarios %d",
               scenario_end, n_sc_total))
}

# --- [2/3 + 3/3] Loop over this task's scenarios, checkpointing each ---------
# For each scenario: keep every region (sim$genotypes) and just that one
# scenario, so scenario_setup pools across the full region panel as designed.
# evaluate_methods() maps fits back via scenario_id used as a LIST INDEX
# (simulation$scenarios[[f$scenario_id]]); after subsetting the single scenario
# sits at position 1, so its scenario_id must be reset to 1. The true scenario
# index is preserved in the scenario_<sc> dir name and params.json. Each
# scenario is saved before the next begins, so a walltime kill only loses the
# scenario in flight and a resubmit resumes at `todo`.
n_done <- 0L
for (scenario_idx in todo) {
  scenario_dir <- file.path(job_dir, sprintf("scenario_%03d", scenario_idx))
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  evaluation_file <- file.path(scenario_dir, paste0("evaluation", OUT_SUFFIX, ".rds"))
  results_file    <- file.path(scenario_dir, paste0("results",    OUT_SUFFIX, ".rds"))
  # Re-check just before the work: a resubmit that overlaps could have filled it.
  if (file.exists(evaluation_file)) {
    cat(sprintf("  scenario %d already complete; skip\n", scenario_idx)); next
  }

  ts   <- Sys.time()
  mini <- sim
  mini$scenarios <- sim$scenarios[scenario_idx]   # length-1 list
  mini$scenarios[[1]]$scenario_id <- 1L
  sc_meta <- sim$scenarios[[scenario_idx]]
  cat(sprintf("[scenario %d/%d] S=%s phi=%s iter=%s  x %d regions x %d methods\n",
              scenario_idx, SCENARIOS_PER_ROW, sc_meta$S %||% "?", sc_meta$phi %||% "?",
              sc_meta$iter %||% "?", job$n_regions, length(METHODS)))

  results <- run_methods(mini, methods = METHODS,
                         method_args = METHOD_ARGS[intersect(names(METHOD_ARGS), METHODS)],
                         save = FALSE, verbose = FALSE)
  saveRDS(results, results_file)
  for (m in METHODS) {
    cat(sprintf("       %-20s n_fits=%d  failed=%d\n",
                m, results[[m]]$n_total %||% 0L, results[[m]]$n_failed %||% 0L))
  }

  # Wrapped so a structural hiccup can't discard the expensive method compute -
  # results.rds is already on disk and collect_results.R can re-evaluate.
  evaluation <- tryCatch(
    evaluate_methods(mini, results, save = FALSE, verbose = FALSE),
    error = function(e) {
      cat(sprintf("  evaluate_methods failed: %s\n", conditionMessage(e)))
      cat("  results.rds is saved; collect_results.R can evaluate later.\n")
      NULL
    }
  )
  if (!is.null(evaluation)) {
    saveRDS(evaluation, evaluation_file)
    meta <- list(
      array_index = array_idx,
      job_id      = job$job_id,
      label       = job$label,
      scenario    = scenario_idx,
      methods     = METHODS,
      seed        = 1000L + row_idx,
      vcf_dir     = VCF_DIR,
      runtime_sec = as.numeric(Sys.time() - ts, units = "secs"),
      R_version   = R.version.string,
      finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    )
    jsonlite::write_json(meta, file.path(scenario_dir, "params.json"),
                         auto_unbox = TRUE, pretty = TRUE, null = "null")
    n_done <- n_done + 1L
    cat(sprintf("  [done] scenario %d in %.1f min  (%d/%d this task)\n",
                scenario_idx, as.numeric(Sys.time() - ts, units = "mins"),
                n_done, length(todo)))
  }
}
cat(sprintf("\n[task complete] %d scenario(s) evaluated in %.1f min total. Output: %s\n",
            n_done, as.numeric(Sys.time() - t0, units = "mins"), job_dir))
