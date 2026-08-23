#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/generate_params_grid_iter005.R
#
# >>> ITERATION 005 ONLY - TEMPORARY. See docs/autoresearch/iteration-005-REVERT.md
#
# A SEPARATE generator, so scripts/hpc/generate_params_grid.R (Iteration 004's,
# and the benchmark's standing design) is untouched and this file can simply be
# deleted when the iteration is reported.
#
# QUESTION. How do annotation-aware methods behave when the annotation ->
# causality relationship is NOT the log-linear form they all assume? Every
# competitor's prior is linear or log-linear in the annotations - PolyFun-est
# regresses chi^2 on A, PAINTOR uses a logistic-linear model, SBayesRC a
# multinomial logit, fb_pooled a linear head. Only functional_beatrice and
# fb_xregion use a LassoNet, which has hidden layers and CAN represent
# nonlinearity. Relationships a linear model cannot express are therefore the
# sharpest available test of whether that architecture earns its complexity.
#
# DESIGN. 6 relationships x 2 annotation types = 12 rows.
#   additive   the assumed form. Control, and the bridge to Iteration 004.
#   cooccur    a variant needs two marks TOGETHER.
#   nonmono    enrichment peaks at INTERMEDIATE annotation load.
#   mixed      three annotations enrich, two DEPLETE.
#   threshold  a step in annotation count.
#   null       annotations independent of causality; methods still receive them.
#
# Every row carries 10 annotations of which only the first 5 are informative;
# the other 5 are independent noise, present in the data handed to every method
# but never entering selection. Annotation SELECTION is therefore a constant
# challenge across all arms, and only the functional form varies.
#
# ENRICHMENT IS FIXED AT 5.4 AND THIS IS DELIBERATE. The relationships produce
# very different enrichment strengths at the same nominal fold; the selection
# weights are calibrated to a common top-decile concentration so that SHAPE is
# the only free variable (see .calibrate_log_weights in R/simulate_phenotypes.R).
# Fold 5.4 is the value at which all twelve arms reach the same target - at
# higher folds nonmono/binary saturates near 0.69 and cannot match.
#
# GRID. Depth over breadth: Iteration 004 found only 4.7% of cells had a
# statistically decidable winner at 10 iterations, so replication buys
# resolution that additional factor levels do not.
#
#   Rscript scripts/hpc/generate_params_grid_iter005.R <out.csv>
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
out_csv <- args[1] %||% "scripts/hpc/params_grid_iter005.csv"
`%||%`  <- function(x, y) if (is.null(x) || is.na(x) || !nzchar(x)) y else x

RELATIONSHIPS <- c("additive", "cooccur", "nonmono", "mixed", "threshold", "null")
ANNOT_TYPES   <- c("binary", "continuous")

# 10 regions, ALL at p = 1000. A single size class means the Iteration 004
# pair-difference estimator for sigma^2_u does not apply; the generalised
# estimator uses the variance across all ten draws directly, which is better
# conditioned. See iteration-005.md.
P_VECTOR <- rep(1000L, 10L)

WITHIN_JOB <- list(
  n_regions = length(P_VECTOR),
  n         = 5000L,
  n_iter    = 25L,              # 2.5x Iteration 004, for per-cell resolution
  S         = c(1L, 3L),
  phi       = c(0.1, 0.4)
)
ENRICHMENT    <- 5.4            # see the banner above - not a free choice
N_ANNOTATIONS <- 10L
N_INFORMATIVE <- 5L

rows <- list(); k <- 0L
for (rel in RELATIONSHIPS) for (at in ANNOT_TYPES) {
  k <- k + 1L
  rows[[k]] <- data.frame(
    job_id                 = k,
    label                  = sprintf("rel%s_an%s", rel,
                                     if (at == "binary") "Binary" else "Cont"),
    model                  = "sparse",
    p_causal               = NA_real_,
    annotation_type        = at,
    # The null arm still SIMULATES annotations and hands them to every method;
    # they simply do not enter causal selection. That is the point of the arm.
    enrichment_fold        = if (rel == "null") NA_real_ else ENRICHMENT,
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
}
grid <- do.call(rbind, rows)
write.csv(grid, out_csv, row.names = FALSE)

spr <- length(WITHIN_JOB$S) * length(WITHIN_JOB$phi) * WITHIN_JOB$n_iter
cat(sprintf("ITERATION 005 grid -> %s\n", out_csv))
cat(sprintf("  %d rows (%d relationships x %d annotation types)\n",
            nrow(grid), length(RELATIONSHIPS), length(ANNOT_TYPES)))
cat(sprintf("  S={%s}  phi={%s}  iters=%d  ->  %d scenarios/row, %d total\n",
            paste(WITHIN_JOB$S, collapse=","), paste(WITHIN_JOB$phi, collapse=","),
            WITHIN_JOB$n_iter, spr, spr * nrow(grid)))
cat(sprintf("  %d regions of p=%d, n=%d, in-sample LD\n",
            length(P_VECTOR), P_VECTOR[1], WITHIN_JOB$n))
cat(sprintf("  %d annotations, first %d informative, enrichment fold %.1f\n",
            N_ANNOTATIONS, N_INFORMATIVE, ENRICHMENT))
