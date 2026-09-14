#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter006_ap_ties.R
#
# >>> ITERATION 006 ONLY - TEMPORARY. See docs/autoresearch/iteration-006-REVERT.md
#
# How much of Iteration 006's average-precision picture is tied PIPs?
#
# The pipeline's AP (.compute_ap_exact in R/evaluate_extras.R, matching
# sklearn.metrics.average_precision_score) takes precision over whole groups of
# tied PIPs - equivalent to ranking the causal variant LAST among the variants it
# ties with. A method whose PIPs tie often is therefore marked down against one
# whose PIPs never tie, even when neither knows more.
#
# This script recomputes AP from the full PIPs (results/iter006/piptail, floor 0)
# with ties broken at random, averaged over 20 draws (the tie-neutral value). It
# prints that beside the pipeline's value, which it reproduces exactly first.
#
#   Rscript scripts/analysis/iter006_ap_ties.R [iter006_dir] [out_dir]
#
# Unit and conventions as iter006_relationship_results.R: per-draw means over the
# ten regions, then mean ± 2 SE over the 100 draws in a relationship; paired
# contrasts carry one SE.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
I6  <- if (length(args) >= 1L) args[1] else "results/iter006"
OUT <- if (length(args) >= 2L) args[2] else file.path(I6, "analysis")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
options(width = 250, stringsAsFactors = FALSE)
REL <- c("null", "additive", "wavg2", "valthresh", "valthresh_wavg2", "square", "cubic", "cosine10")
METHODS <- c("susie", "beatrice", "polyfun_ldsc", "paintor", "sbayesrc", "funmap",
             "functional_beatrice", "fb_xregion", "polyfun_oracle")
CONTRASTS <- list(c("polyfun_oracle", "susie"), c("polyfun_ldsc", "susie"), c("funmap", "susie"),
                  c("fb_xregion", "polyfun_ldsc"), c("fb_xregion", "beatrice"))

ap_block <- function(pip, y) {                       # the pipeline's definition
  o <- order(-pip); p <- pip[o]; yy <- y[o]; g <- cumsum(c(TRUE, diff(p) != 0))
  tp <- as.numeric(tapply(yy, g, sum)); n <- as.numeric(tapply(yy, g, length))
  sum((cumsum(tp) / cumsum(n)) * tp) / sum(y)
}
ap_rank   <- function(o, y) { k <- which(y[o] == 1); sum(seq_along(k) / k) / sum(y) }
ap_random <- function(pip, y, R = 20L) mean(vapply(seq_len(R), function(i) ap_rank(order(-pip, runif(length(pip))), y), 0))

set.seed(2026)
grid <- read.csv(file.path(I6, "params_grid_iter006.csv"))
job_rel <- setNames(grid$relationship, sprintf("job_%03d_%s", grid$job_id, grid$label))
fs <- sort(list.files(file.path(I6, "piptail"), "^piptail_.*[.]rds$", full.names = TRUE))
message("reading ", length(fs), " full-PIP files ...")
FITS <- do.call(rbind, lapply(fs, function(f) {
  o <- readRDS(f)
  if (!identical(as.numeric(o$floor), 0)) stop(basename(f), " was not extracted at floor 0", call. = FALSE)
  do.call(rbind, lapply(o$tail, function(t) {
    if (!length(t$pip) || !any(t$is_causal == 1L)) return(NULL)
    if (t$n_below != 0) stop("PIPs missing below the floor in ", basename(f), call. = FALSE)
    y <- t$is_causal; pip <- t$pip; cp <- pip[y == 1L]
    data.frame(job_dir = t$job_dir, scenario_id = t$scenario_id, region_id = t$region_id, method = t$method,
               causal_tied = mean(vapply(cp, function(v) any(pip[y == 0L] == v), TRUE)),
               ap_pipeline = ap_block(pip, y), ap_tie_neutral = ap_random(pip, y))
  }))
}))
FITS$relationship <- job_rel[FITS$job_dir]

L1 <- readRDS(file.path(I6, "combined_fit_metrics.rds"))
chk <- merge(FITS, L1[, c("job_dir", "scenario_id", "region_id", "method", "ap")],
             by = c("job_dir", "scenario_id", "region_id", "method"))
dev <- max(abs(chk$ap_pipeline - chk$ap), na.rm = TRUE)
DRAW <- aggregate(cbind(causal_tied, ap_pipeline, ap_tie_neutral) ~ relationship + job_dir + scenario_id + method, FITS, mean)

sink(file.path(OUT, "iter006_ap_ties.txt"), split = TRUE)
cat("Iteration 006 - average precision with ties as the pipeline scores them, and broken at random\n")
cat(sprintf("fits read: %d (expect 72000); pipeline AP reproduced from the full PIPs, max |diff| %.1e over %d fits\n",
            nrow(FITS), dev, nrow(chk)))

cat("\n== Share of fits whose causal variant's PIP is exactly tied with a non-causal variant's\n")
t0 <- data.frame(method = METHODS)
for (r in REL) t0[[r]] <- sprintf("%.1f%%", 100 * vapply(METHODS, function(m) mean(FITS$causal_tied[FITS$relationship == r & FITS$method == m]), 0))
print(t0, row.names = FALSE, right = FALSE)

cat("\n== Average precision: pipeline -> ties broken at random (mean over 100 draws; ± 2 SE is about 0.02 for both)\n")
t1 <- data.frame(method = METHODS)
for (r in REL) t1[[r]] <- vapply(METHODS, function(m) { d <- DRAW[DRAW$relationship == r & DRAW$method == m, ]
  sprintf("%.3f -> %.3f", mean(d$ap_pipeline), mean(d$ap_tie_neutral)) }, "")
print(t1, row.names = FALSE, right = FALSE)

cat("\n== Paired differences by draw: pipeline (SE) | ties broken at random (SE)\n")
KEYD <- c("relationship", "job_dir", "scenario_id")
con <- do.call(rbind, lapply(CONTRASTS, function(ab) do.call(rbind, lapply(REL, function(r) {
  x <- merge(DRAW[DRAW$relationship == r & DRAW$method == ab[1], ], DRAW[DRAW$relationship == r & DRAW$method == ab[2], ], by = KEYD)
  dp <- x$ap_pipeline.x - x$ap_pipeline.y; dt <- x$ap_tie_neutral.x - x$ap_tie_neutral.y
  data.frame(contrast = paste(ab[1], "-", ab[2]), relationship = r,
             pipeline = mean(dp), pipeline_se = sd(dp) / sqrt(length(dp)),
             tie_neutral = mean(dt), tie_neutral_se = sd(dt) / sqrt(length(dt)))
}))))
t2 <- data.frame(contrast = vapply(CONTRASTS, paste, "", collapse = " - "))
for (r in REL) { s <- con[con$relationship == r, ]
  t2[[r]] <- sprintf("%+.3f (%.3f) | %+.3f (%.3f)", s$pipeline, s$pipeline_se, s$tie_neutral, s$tie_neutral_se)[match(t2$contrast, s$contrast)] }
print(t2, row.names = FALSE, right = FALSE)
sink()
saveRDS(list(fits = FITS, draws = DRAW, contrasts = con, reproduction_max_diff = dev), file.path(OUT, "iter006_ap_ties.rds"))
message("wrote ", file.path(OUT, "iter006_ap_ties.txt"))
