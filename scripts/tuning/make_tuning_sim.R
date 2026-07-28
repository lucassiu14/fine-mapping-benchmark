#!/usr/bin/env Rscript
# =============================================================================
# scripts/tuning/make_tuning_sim.R
#
# Build and cache the FIXED tuning simulation that every Optuna trial is scored
# on. Fixing the data matters: trials must differ only in hyperparameters, so
# the same simulation (same seed) is reused for every trial and every worker.
#
# The tuning set is deliberately SMALL - a few scenarios x a few regions - so a
# trial costs minutes, not hours. It is drawn from the same generative model as
# the Iteration 002 benchmark grid, restricted to ONE stratum, because the
# never-pool rule applies here too: an annotation prior tuned on binary
# annotations under in-sample LD has no reason to be optimal for continuous
# annotations under a reference panel. Run one study per stratum.
#
# Usage (defaults = the stratum where Functional BEATRICE actually broke:
#        binary annotations, high enrichment, reference-panel LD):
#   Rscript scripts/tuning/make_tuning_sim.R --out tuning/sim_binary_ref500.rds \
#     --annotations binary --enrichment 10.8 --n_ref 500 --model sparse \
#     --scenarios 4 --regions 4
# =============================================================================
suppressWarnings(suppressMessages({
  if (requireNamespace("fmbenchmark", quietly = TRUE)) {
    library(fmbenchmark)
  } else if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else stop("fmbenchmark not installed and pkgload/DESCRIPTION unavailable")
}))

args <- commandArgs(trailingOnly = TRUE)
opt <- list(); i <- 1L
while (i <= length(args)) {
  if (startsWith(args[i], "--")) {
    k <- sub("^--", "", args[i])
    v <- if (i + 1L <= length(args) && !startsWith(args[i + 1L], "--")) args[i + 1L] else "TRUE"
    opt[[k]] <- v; i <- i + 2L
  } else i <- i + 1L
}
gets <- function(k, d) if (!is.null(opt[[k]])) opt[[k]] else d
getn <- function(k, d) if (!is.null(opt[[k]])) as.numeric(opt[[k]]) else d
geti <- function(k, d) if (!is.null(opt[[k]])) as.integer(as.numeric(opt[[k]])) else d

OUT         <- gets("out", "tuning/tuning_sim.rds")
ANNOT       <- gets("annotations", "binary")          # none | binary | continuous
ENRICH      <- getn("enrichment", 10.8)
N_REF       <- gets("n_ref", "500")                   # "NA" / "" -> in-sample LD
MODEL       <- gets("model", "sparse")                # sparse | sparse_inf
P_CAUSAL    <- gets("p_causal", "")                   # sparse_inf only
N_SCEN_S    <- gets("S", "1,3")                       # scenario grid: S values
N_SCEN_PHI  <- gets("phi", "0.05,0.2")                #                phi values
N_REGIONS   <- geti("regions", 4L)
N_ANNOT     <- geti("n_annotations", 10L)
N_SUB       <- geti("n", 1000L)
SEED        <- geti("seed", 20260801L)
VCF_DIR     <- gets("vcf_dir", Sys.getenv("FMB_VCF_DIR", unset = ""))
MAP_DIR     <- gets("genetic_map_dir", Sys.getenv("FMB_GENETIC_MAP_DIR", unset = ""))

S_vec   <- as.integer(strsplit(N_SCEN_S,   ",")[[1]])
phi_vec <- as.numeric(strsplit(N_SCEN_PHI, ",")[[1]])
# region sizes: spread across the benchmark's length classes, truncated to N_REGIONS
p_vec <- rep(c(100L, 200L, 400L, 500L, 1000L), length.out = N_REGIONS)

# Half the annotations enriched, matching the Iteration 002 grid convention.
enrich_vec <- if (identical(ANNOT, "none")) NULL else
  c(rep(ENRICH, N_ANNOT %/% 2L), rep(1, N_ANNOT - N_ANNOT %/% 2L))

sim_args <- list(
  n_regions              = N_REGIONS,
  n                      = N_SUB,
  p                      = p_vec,
  n_iter                 = 1L,          # one draw per (S, phi): the tuning set is fixed
  S                      = S_vec,
  phi                    = phi_vec,
  model                  = MODEL,
  annotations            = ANNOT,
  n_annotations          = if (identical(ANNOT, "none")) 0L else N_ANNOT,
  enrichment             = enrich_vec,
  annotation_correlation = 0,
  vcf_dir                = if (nzchar(VCF_DIR)) VCF_DIR else NULL,
  genetic_map_dir        = if (nzchar(MAP_DIR)) MAP_DIR else NULL,
  seed                   = SEED,
  save                   = FALSE,
  verbose                = FALSE
)
if (identical(MODEL, "sparse_inf") && nzchar(P_CAUSAL)) sim_args$p_causal <- as.numeric(P_CAUSAL)
if (nzchar(N_REF) && !identical(toupper(N_REF), "NA")) sim_args$n_ref <- as.integer(N_REF)

cat(sprintf("Building tuning sim: model=%s annot=%s enrich=%s LD=%s regions=%d S={%s} phi={%s}\n",
            MODEL, ANNOT, ifelse(is.null(enrich_vec), "-", ENRICH),
            ifelse(nzchar(N_REF) && !identical(toupper(N_REF), "NA"),
                   paste0("n_ref=", N_REF), "in-sample"),
            N_REGIONS, N_SCEN_S, N_SCEN_PHI))

t0 <- Sys.time()
sim <- do.call(run_simulation, sim_args)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
saveRDS(sim, OUT)

cat(sprintf("wrote %s\n  scenarios: %d   regions/scenario: %d   (%.1f min)\n",
            OUT, length(sim$scenarios), length(sim$genotypes),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("  scenario grid:\n")
for (k in seq_along(sim$scenarios)) {
  s <- sim$scenarios[[k]]
  cat(sprintf("    [%d] S=%s phi=%s iter=%s\n", k, s$S, s$phi, s$iter))
}
