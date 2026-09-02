#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_annotation_recovery.R
#
# Does the LassoNet prior select the annotations that actually carry signal?
#
#   Rscript scripts/analysis/iter004_annotation_recovery.R <aux_dir> <l3_rds> [out_rds]
#
# The simulator generates ten annotation tracks of which the FIRST FIVE carry
# the fold-enrichment and the last five are inert (generate_params_grid.R:
# n_annotations = 10, n_enriched = 5). Recovery is therefore a clean ranking
# problem with a known answer, and it is the one claim in the paper no
# competitor can make in the same form: PolyFun and Funmap consume annotations
# but do not report which ones they used.
#
# Metrics, per fit, over the importance vector (LassoNet |theta| contrast):
#   precision@5   how many of the top five ranked tracks are truly enriched
#   auc           P(a random enriched track outranks a random inert one),
#                 computed exactly from the rank sum (Mann-Whitney)
#   mean_rank_e   mean rank of the enriched tracks, 1 = best
#
# WHAT THIS CAN AND CANNOT COVER
# ------------------------------
# functional_beatrice returns feature_importance and is fully analysable.
#
# fb_xregion IS NOT. wrapper_fb_joint.R's .fb_joint_parse_output reads only
# pip.csv, credible_set.txt and conditional_..._probability.txt; the shared
# head's importances are written by joint_trainer.py into <target>/res, and
# <target> lives under work_dir <- tempfile(), which does not survive the R
# session. They were computed and discarded.
#
# So any importance record labelled fb_xregion can only have come from the
# per-region FALLBACK path (wrapper_fb_joint.R:247), which calls
# run_functional_beatrice_region and tags additional$joint_fallback = TRUE.
# extract_aux.R does not carry that flag, so such rows are single-locus
# results wearing a cross-region label. This script REFUSES to score them and
# reports the count instead, because reporting them as cross-region recovery
# would be a straightforward misattribution.
# =============================================================================

N_ANNOT    <- 10L
N_ENRICHED <- 5L
JOINT_METHODS <- c("fb_xregion", "fb_pooled")

args    <- commandArgs(trailingOnly = TRUE)
aux_dir <- args[1]
l3_file <- if (length(args) >= 2) args[2] else "results/iter004/combined_scenario_metrics.rds"
out_file<- if (length(args) >= 3) args[3] else "results/iter004/annotation_recovery.rds"
stopifnot(dir.exists(aux_dir), file.exists(l3_file))

files <- sort(list.files(aux_dir, pattern = "^aux_.*\\.rds$", full.names = TRUE))
message("aux files: ", length(files), " of 45 expected")
if (!length(files)) stop("nothing to read", call. = FALSE)

design <- unique(readRDS(l3_file)[, c("job_dir", "model", "annotation_type")])
truth  <- c(rep(TRUE, N_ENRICHED), rep(FALSE, N_ANNOT - N_ENRICHED))

rows <- list(); k <- 0L
n_joint_labelled <- 0L; n_malformed <- 0L
for (f in files) {
  o <- readRDS(f)
  d <- design[match(o$job_dir, design$job_dir), ]
  for (r in o$importance) {
    if (r$method %in% JOINT_METHODS) { n_joint_labelled <- n_joint_labelled + 1L; next }
    fi <- r$importance
    # data frame (annotation_index, importance), sorted descending by the
    # wrapper - reorder by index so position i is annotation i.
    imp <- if (is.data.frame(fi) && all(c("annotation_index", "importance") %in% names(fi))) {
      v <- rep(NA_real_, N_ANNOT)
      ix <- as.integer(fi$annotation_index)
      ix <- if (min(ix, na.rm = TRUE) == 0L) ix + 1L else ix        # 0- or 1-based
      ok <- ix >= 1L & ix <= N_ANNOT
      v[ix[ok]] <- as.numeric(fi$importance)[ok]; v
    } else if (is.numeric(fi) && length(fi) == N_ANNOT) as.numeric(fi) else NULL
    if (is.null(imp) || anyNA(imp)) { n_malformed <- n_malformed + 1L; next }

    rk   <- rank(-imp, ties.method = "average")     # 1 = most important
    top5 <- order(-imp)[seq_len(N_ENRICHED)]
    # Mann-Whitney AUC from the rank sum of the enriched group.
    Re <- sum(rk[truth]); ne <- sum(truth); ni <- N_ANNOT - ne
    auc <- (ne * ni + ne * (ne + 1) / 2 - Re) / (ne * ni)

    k <- k + 1L
    rows[[k]] <- data.frame(
      job_dir = o$job_dir, model = d$model[1], annotation_type = d$annotation_type[1],
      scenario_id = r$scenario_id, region_id = r$region_id, method = r$method,
      precision5 = mean(truth[top5]), auc = auc, mean_rank_e = mean(rk[truth]),
      stringsAsFactors = FALSE)
  }
}
if (!k) stop("no scorable importance records", call. = FALSE)
tab <- do.call(rbind, rows)

message("scorable fits: ", k)
if (n_joint_labelled)
  message("SKIPPED ", n_joint_labelled, " records labelled ",
          paste(JOINT_METHODS, collapse = "/"), " - these are per-region fallbacks, ",
          "not cross-region fits (see the header)")
if (n_malformed) message("skipped ", n_malformed, " malformed importance vectors")

saveRDS(tab, out_file); message("wrote ", out_file)

cat("\n== recovery by stratum (mean over fits) ==\n")
agg <- aggregate(cbind(precision5, auc, mean_rank_e) ~ model + annotation_type + method,
                 tab[tab$annotation_type != "none", ], mean)
agg$n <- aggregate(auc ~ model + annotation_type + method,
                   tab[tab$annotation_type != "none", ], length)$auc
print(agg[order(agg$model, agg$annotation_type), ], row.names = FALSE, digits = 3)
cat("\nchance level: precision@5 = 0.5, auc = 0.5, mean rank of enriched = 5.5\n")
