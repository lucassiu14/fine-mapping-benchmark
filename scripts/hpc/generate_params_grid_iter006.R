#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/generate_params_grid_iter006.R
#
# >>> ITERATION 006 ONLY - TEMPORARY. See docs/autoresearch/iteration-006-REVERT.md
#
# A SEPARATE generator, like Iteration 005's, so the standing generator
# (scripts/hpc/generate_params_grid.R, Iteration 004's design) is untouched and
# this file can be deleted when the iteration is reported.
#
# QUESTION. Iteration 005's, with relationships the user specified: how do the
# annotation-aware methods behave when the annotation -> causality relationship
# is not the log-linear form they assume?
#
# DESIGN. 8 relationships, continuous annotations only = 8 rows. Everything else
# is Iteration 005's grid.
#   null             annotations handed to every method, no effect on causality
#   additive         the log-linear control - the form every method assumes
#   wavg2            weighted average of the variant and 2 neighbours each side
#   valthresh        a value counts only above 1
#   valthresh_wavg2  threshold each variant, then the weighted average
#   square           a^2
#   cubic            a^3
#   cosine10         damped-cosine weighted average over +-10 variants
# The forms are defined in R/simulate_phenotypes.R and checked against
# brute force in scripts/analysis/test_iter006_relationships.R.
#
# STRENGTH. Every non-null arm is scaled so the top 10% of variants hold 34.9% of
# the prior probability - Iteration 005's lower concentration level, chosen by
# the user over the 60.5% that Iteration 005 ran at. select_causal_variants()
# looks the target up from the informative tracks' enrichment value, so
# ENRICHMENT = 2.7 below is what selects it. As in Iteration 005,
# enrichment_fold therefore indexes that ladder, not a fold.
#
#   Rscript scripts/hpc/generate_params_grid_iter006.R <out.csv>
# =============================================================================

`%||%`  <- function(x, y) if (is.null(x) || is.na(x) || !nzchar(x)) y else x
args    <- commandArgs(trailingOnly = TRUE)
out_csv <- args[1] %||% "scripts/hpc/params_grid_iter006.csv"

RELATIONSHIPS <- c("null", "additive", "wavg2", "valthresh", "valthresh_wavg2",
                   "square", "cubic", "cosine10")

P_VECTOR <- rep(1000L, 10L)            # Iteration 005: 10 regions, all p = 1000
WITHIN_JOB <- list(
  n_regions = length(P_VECTOR),
  n         = 5000L,
  n_iter    = 25L,
  S         = c(1L, 3L),
  phi       = c(0.1, 0.4)
)
ENRICHMENT    <- 2.7                   # selects top-decile target 0.349 - see STRENGTH
N_ANNOTATIONS <- 10L
N_INFORMATIVE <- 5L

rows <- lapply(seq_along(RELATIONSHIPS), function(k) {
  rel <- RELATIONSHIPS[k]
  data.frame(
    job_id                 = k,
    label                  = sprintf("rel%s_anCont", rel),
    model                  = "sparse",
    p_causal               = NA_real_,
    annotation_type        = "continuous",
    enrichment_fold        = if (rel == "null") NA_real_ else ENRICHMENT,
    # What the worker parses: the first N_INFORMATIVE tracks carry the value, the
    # rest are inert at 1; the null arm is all ones (see Iteration 005's generator).
    enrichment_values      = paste(
                               if (rel == "null") rep(1, N_ANNOTATIONS)
                               else c(rep(ENRICHMENT, N_INFORMATIVE),
                                      rep(1, N_ANNOTATIONS - N_INFORMATIVE)),
                               collapse = "|"),
    relationship           = rel,
    n_informative          = N_INFORMATIVE,
    annotation_correlation = 0,
    n_ref                  = NA_integer_,
    n_annotations          = N_ANNOTATIONS,
    n_regions              = WITHIN_JOB$n_regions,
    n                      = WITHIN_JOB$n,
    n_iter                 = WITHIN_JOB$n_iter,
    S_values               = paste(WITHIN_JOB$S,   collapse = "|"),
    phi_values             = paste(WITHIN_JOB$phi, collapse = "|"),
    p_values               = paste(P_VECTOR,       collapse = "|"),
    stringsAsFactors = FALSE)
})
grid <- do.call(rbind, rows)
write.csv(grid, out_csv, row.names = FALSE)

spr <- length(WITHIN_JOB$S) * length(WITHIN_JOB$phi) * WITHIN_JOB$n_iter
cat(sprintf("ITERATION 006 grid -> %s\n", out_csv))
cat(sprintf("  %d rows: %s, continuous annotations\n", nrow(grid), paste(RELATIONSHIPS, collapse = ", ")))
cat(sprintf("  S={%s}  phi={%s}  iters=%d  ->  %d scenarios/row, %d total\n",
            paste(WITHIN_JOB$S, collapse = ","), paste(WITHIN_JOB$phi, collapse = ","),
            WITHIN_JOB$n_iter, spr, spr * nrow(grid)))
cat(sprintf("  %d regions of p=%d, n=%d, in-sample LD; %d annotations, first %d informative\n",
            length(P_VECTOR), P_VECTOR[1], WITHIN_JOB$n, N_ANNOTATIONS, N_INFORMATIVE))
cat(sprintf("  strength: top-decile target %.3f (enrichment value %.1f)\n", 0.349, ENRICHMENT))
