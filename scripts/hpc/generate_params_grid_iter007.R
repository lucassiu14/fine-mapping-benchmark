#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/generate_params_grid_iter007.R
#
# >>> ITERATION 007 ONLY - TEMPORARY. See docs/autoresearch/iteration-007-REVERT.md
#
# A SEPARATE generator, like Iterations 005 and 006, so the standing generator
# (scripts/hpc/generate_params_grid.R, Iteration 004's design) is untouched and
# this file can be deleted when the iteration is reported.
#
# QUESTION. How does LD misspecification - a method receiving the LD of an
# independent reference panel rather than the GWAS sample's own - change
# fine-mapping, and how does that depend on the size of the panel?
#
# DESIGN (user-specified 2026-09-14). 20 rows:
#   LD           in-sample; reference panels of 500, 750, 1500, 2000       5
#   model        sparse; sparse_inf at p_causal 0.6                        2
#   annotations  binary, continuous; fold 5.4 on the first 5 of 10 tracks  2
# Within each row: S {1,3,5} x phi {0.1,0.4} x 10 iterations = 60 scenarios.
#
# GWAS n = 1000. Regions are 2 each of {100,150,200,300,400} variants: Iteration
# 004's sizes divided by 5, so p/n <= 0.4 as in Iteration 004. Every region is
# smaller than every sample (GWAS or panel), so no LD matrix loses rank through
# sample size; variants in perfect LD can still make one singular, as in
# Iteration 004.
#
# UNPAIRED. Each row is its own run_simulation() call with seed 1000 + row
# (run_benchmark_job.R), so the LD levels do not share genotypes, causal variants
# or z-scores. The panel is drawn in the same sim1000G session as the GWAS
# sample, independently of it, and standardised by its own column means and SDs
# (R/simulate_genotypes.R); methods then receive cor(X_ref).
#
# Not comparable to Iteration 004: n, region sizes, S and phi levels, enrichment
# and LD all differ.
#
#   Rscript scripts/hpc/generate_params_grid_iter007.R <out.csv>
# =============================================================================

`%||%`  <- function(x, y) if (is.null(x) || is.na(x) || !nzchar(x)) y else x
args    <- commandArgs(trailingOnly = TRUE)
out_csv <- args[1] %||% "scripts/hpc/params_grid_iter007.csv"

P_VECTOR <- c(rep(100L, 2L), rep(150L, 2L), rep(200L, 2L),
              rep(300L, 2L), rep(400L, 2L))
WITHIN_JOB <- list(
  n_regions = length(P_VECTOR),
  n         = 1000L,
  n_iter    = 10L,
  S         = c(1L, 3L, 5L),
  phi       = c(0.1, 0.4)
)
N_REF_LEVELS        <- c(NA_integer_, 500L, 750L, 1500L, 2000L)   # NA = in-sample
P_CAUSAL_SPARSE_INF <- 0.6
ANNOTATION_TYPES    <- c("binary", "continuous")
ENRICHMENT          <- 5.4
N_ANNOTATIONS       <- 10L
N_ENRICHED          <- 5L

# The premise of the region sizes above: no LD matrix is rank-limited by n.
stopifnot(max(P_VECTOR) / WITHIN_JOB$n <= 0.4,
          max(P_VECTOR) < min(N_REF_LEVELS, na.rm = TRUE))

# Labels follow Iteration 004's format, so scripts that read the model and the
# annotation type from the job directory name work unchanged.
nref_token <- function(nr) if (is.na(nr)) "refInsample" else sprintf("ref%d", nr)
arm_token  <- c(binary = "anBinary", continuous = "anCont")
enrich_vec <- c(rep(ENRICHMENT, N_ENRICHED), rep(1, N_ANNOTATIONS - N_ENRICHED))

rows <- list()
job_id <- 0L
for (m in c("sparse", "sparse_inf")) {
  pc <- if (m == "sparse") NA_real_ else P_CAUSAL_SPARSE_INF
  for (at in ANNOTATION_TYPES) {
    for (nr in N_REF_LEVELS) {
      job_id <- job_id + 1L
      rows[[job_id]] <- data.frame(
        job_id                 = job_id,
        label                  = sprintf("%s_%s_e%g%s_%s", m, arm_token[[at]], ENRICHMENT,
                                         if (m == "sparse_inf") sprintf("_pc%g", pc) else "",
                                         nref_token(nr)),
        model                  = m,
        p_causal               = pc,
        annotation_type        = at,
        annotation_correlation = 0,
        enrichment_fold        = ENRICHMENT,
        n_ref                  = nr,
        n_annotations          = N_ANNOTATIONS,
        n_regions              = WITHIN_JOB$n_regions,
        n                      = WITHIN_JOB$n,
        n_iter                 = WITHIN_JOB$n_iter,
        S_values               = paste(WITHIN_JOB$S,   collapse = "|"),
        phi_values             = paste(WITHIN_JOB$phi, collapse = "|"),
        p_values               = paste(P_VECTOR,       collapse = "|"),
        enrichment_values      = paste(enrich_vec,     collapse = "|"),
        stringsAsFactors       = FALSE)
    }
  }
}
grid <- do.call(rbind, rows)
stopifnot(nrow(grid) == 20L, !anyDuplicated(grid$label))

dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(grid, out_csv, row.names = FALSE)

spr <- length(WITHIN_JOB$S) * length(WITHIN_JOB$phi) * WITHIN_JOB$n_iter
cat(sprintf("ITERATION 007 grid -> %s\n", out_csv))
cat(sprintf("  %d rows: LD {%s} x model {sparse, sparse_inf pc=%g} x annotations {%s}\n",
            nrow(grid), paste(vapply(N_REF_LEVELS, nref_token, ""), collapse = ", "),
            P_CAUSAL_SPARSE_INF, paste(ANNOTATION_TYPES, collapse = ", ")))
cat(sprintf("  S={%s}  phi={%s}  iters=%d  ->  %d scenarios/row, %d total\n",
            paste(WITHIN_JOB$S, collapse = ","), paste(WITHIN_JOB$phi, collapse = ","),
            WITHIN_JOB$n_iter, spr, spr * nrow(grid)))
cat(sprintf("  n=%d; regions {%s}; %d annotations, first %d at fold %g\n",
            WITHIN_JOB$n, paste(P_VECTOR, collapse = ","), N_ANNOTATIONS, N_ENRICHED, ENRICHMENT))
