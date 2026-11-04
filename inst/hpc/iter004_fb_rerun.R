#!/usr/bin/env Rscript
# Re-run fb_xregion ONLY over the Iteration 004 grid, with the four fixes that
# landed on Functional BEATRICE main at f3bd29a (the version that produced the
# UK Biobank cholesterol application results).
#
# One PBS array task = one of the 45 grid rows, read from iter004_grid.csv,
# which was generated from results/iter004/combined_scenario_metrics.rds and is
# therefore the grid as actually run, not a reconstruction from memory.
#
# Output goes to its own directory. results/iter004 is NEVER written to: the
# LSR was submitted from those numbers and they stay reproducible from the
# directory they came from, exactly as Iteration 008 handled the comparator
# re-run (see results/iter004/SUPERSEDED_COMPARATORS.md).
#
# IMPORTANT - the simulated data must be IDENTICAL to Iteration 004's, or the
# new fb_xregion numbers cannot be compared against the other 18 methods'
# existing numbers; they would only be internally comparable. The seed is the
# one parameter that cannot be recovered from the stored results, so it must be
# supplied and is checked for below.

suppressPackageStartupMessages(library(fmbenchmark))

req <- function(name, cast = identity) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || !nzchar(v))
    stop("Set ", name, ". See the header of iter004_fb_rerun.pbs.", call. = FALSE)
  cast(v)
}

# --- the parameters that must match Iteration 004 exactly --------------------
SEED       <- req("ITER004_SEED",      as.integer)
N          <- req("ITER004_N",         as.integer)
N_ITER     <- req("ITER004_N_ITER",    as.integer)
N_REGIONS  <- req("ITER004_N_REGIONS", as.integer)
VCF_DIR    <- req("ITER004_VCF_DIR")
# --- where the method lives --------------------------------------------------
FB_DIR     <- req("FB_DIR")
PYTHON     <- req("FMB_PYTHON")
OUT_DIR    <- Sys.getenv("OUT_DIR", "results/iter009_fb_xregion")

# The grid levels, read off Iteration 004's own metrics.
S_LEVELS    <- c(1L, 2L, 3L, 5L, 10L)
PHI_LEVELS  <- c(0.05, 0.1, 0.2, 0.4, 0.6)
REGION_SIZE <- c(500L, 750L, 1000L, 1500L, 2000L)

task <- as.integer(Sys.getenv("PBS_ARRAY_INDEX", "1"))
grid <- utils::read.csv(system.file("hpc", "iter004_grid.csv",
                                    package = "fmbenchmark"),
                        stringsAsFactors = FALSE)
if (!nrow(grid)) grid <- utils::read.csv("inst/hpc/iter004_grid.csv",
                                         stringsAsFactors = FALSE)
stopifnot("PBS_ARRAY_INDEX outside 1..45" = task >= 1L && task <= nrow(grid))
row <- grid[grid$job == task, ]
stopifnot("no grid row for this task" = nrow(row) == 1L)

message(sprintf("[iter004-fb] task %d  %s", task, row$job_dir))
message(sprintf("[iter004-fb] seed %d  n %d  n_iter %d  n_regions %d",
                SEED, N, N_ITER, N_REGIONS))

job_out <- file.path(OUT_DIR, row$job_dir)
dir.create(job_out, recursive = TRUE, showWarnings = FALSE)

# fb_xregion shares ONE prior head across the regions of a scenario, so the
# method is given the whole scenario; run_methods handles that through the
# registry's scenario_setup hook.
margs <- list(fb_xregion = list(beatrice_dir = FB_DIR, python = PYTHON,
                                prior_head = "lassonet"))

all_metrics <- list()
for (p in REGION_SIZE) {
  message(sprintf("[iter004-fb]   region_size %d", p))
  sim <- run_simulation(
    n_regions   = N_REGIONS,
    n           = N,
    p           = p,
    n_iter      = N_ITER,
    S           = S_LEVELS,
    phi         = PHI_LEVELS,
    model       = row$model,
    p_causal    = if (is.na(row$p_causal)) 0.1 else row$p_causal,
    annotations = row$annotations,
    enrichment  = if (is.na(row$enrichment)) NULL else row$enrichment,
    vcf_dir     = VCF_DIR,
    seed        = SEED + task,          # per-row, deterministic given SEED
    save        = FALSE,
    verbose     = FALSE
  )
  fits <- run_methods(sim, methods = "fb_xregion", method_args = margs,
                      save = FALSE, verbose = TRUE)
  ev <- evaluate_methods(sim, fits, per_region = TRUE,
                         save = FALSE, verbose = FALSE)
  saveRDS(list(region_size = p, sim_label = row$job_dir, eval = ev, fits = fits),
          file.path(job_out, sprintf("fb_xregion_p%d.rds", p)))
  all_metrics[[as.character(p)]] <- ev
}

saveRDS(all_metrics, file.path(job_out, "metrics.rds"))
message("[iter004-fb] done -> ", job_out)
