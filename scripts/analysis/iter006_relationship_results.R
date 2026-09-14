#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter006_relationship_results.R
#
# >>> ITERATION 006 ONLY - TEMPORARY. See docs/autoresearch/iteration-006-REVERT.md
#
# First results for Iteration 006, one relationship at a time. Annotations are
# continuous throughout, so the relationship IS the stratum; nothing is pooled
# across relationships.
#
#   Rscript scripts/analysis/iter006_relationship_results.R [iter006_dir] [out_dir]
#
# Defaults: results/iter006 and results/iter006/analysis. Inputs are only read.
#
#   1. average precision for every method under each relationship
#   2. how much of the oracle's gain over susie each method recovers
#   3. paired contrasts: fb_xregion against polyfun_ldsc, functional_beatrice and
#      beatrice; functional_beatrice against beatrice; polyfun_ldsc against susie
#   4. annotation selection by the BEATRICE family: tracks zeroed, and tie-aware
#      precision@5 and AUC for the five informative tracks, read against the
#      null arm, where no track is informative
#   5. each method's advantage over susie MINUS the same advantage under null.
#      Under null the annotations carry no information, so an advantage there is
#      not an annotation gain. polyfun_ldsc and funmap both show one on continuous
#      annotations, in Iteration 005 as well.
#
# AP here is the pipeline's: .compute_ap_exact scores tied PIPs as one block, as
# sklearn does. susie's causal variant is tied with a non-causal one in 22-33% of
# fits, so every comparison against susie - including the headroom in section 2 -
# is inflated. iter006_ap_ties.R recomputes AP with ties broken at random.
#
# Conventions, with sources:
#   * The unit is one scenario draw - S x phi x iteration, ten regions - i.e. L2
#     (combined_replicate_metrics.rds), the per-iteration unit the LSR adopted for
#     calibration and FDR. Each relationship has only four S x phi cells, so an SE
#     across L3 cells would rest on four values. The mean over draws equals the
#     mean over cells because the design is balanced.
#   * Means carry +- 2 SE, as the LSR text prints them; paired differences carry
#     one SE in brackets. Pairs match on job_dir x S x phi x region_size x iter.
#   * Every non-null arm is strength-matched: the top 10% of variants hold 34.9%
#     of the prior (generate_params_grid_iter006.R), so the arms differ in shape.
#   * The learnability shares printed above the AP table were measured for this
#     design (docs/autoresearch/iteration-006.md, section 4): the share of the log
#     selection weight a log-linear prior, and any function of a variant's own
#     annotations, can explain.
#   * polyfun_oracle reads the true selection probabilities, so it is the ceiling.
#     susie is annotation-blind, so it is the floor. Recovered share =
#     mean(method - susie) / mean(oracle - susie), with a ratio-estimator SE.
#   * Importance. fb_xregion's shared head writes one vector per scenario into all
#     ten regions; it is scored once per scenario, and fallback fits are excluded.
#     functional_beatrice's is per region. Tracks 1-5 are informative
#     (n_informative = 5). Under null nothing is informative, so about 0.5 is
#     expected. Ties are broken at random in expectation, as in
#     fb_lambda_sweep_analysis.R.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
I6  <- if (length(args) >= 1L) args[1] else "results/iter006"
OUT <- if (length(args) >= 2L) args[2] else file.path(I6, "analysis")
if (!dir.exists(I6)) stop("no such folder: ", I6, call. = FALSE)
if (grepl("(^|/)iter00[45](/|$)", OUT)) stop("refusing to write into another iteration's folder: ", OUT, call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
options(width = 250, stringsAsFactors = FALSE)

REL <- c("null", "additive", "wavg2", "valthresh", "valthresh_wavg2", "square", "cubic", "cosine10")
METHODS <- c("susie", "beatrice", "polyfun_ldsc", "paintor", "sbayesrc", "funmap",
             "functional_beatrice", "fb_xregion", "polyfun_oracle")
LEARN <- data.frame(relationship = REL,                      # iteration-006.md, section 4
                    linear = c(NA, 1.00, 0.48, 0.48, 0.22, 0.01, 0.61, 0.38),
                    own    = c(NA, 1.00, 0.49, 0.94, 0.46, 1.00, 0.99, 0.39))
CONTRASTS <- list(c("fb_xregion", "polyfun_ldsc"), c("fb_xregion", "functional_beatrice"),
                  c("fb_xregion", "beatrice"), c("functional_beatrice", "beatrice"),
                  c("polyfun_ldsc", "susie"))
KEY  <- c("job_dir", "S", "phi", "region_size", "iter")
N_ANNOT <- 10L; TRUTH <- rep(c(TRUE, FALSE), each = 5L); ZERO_TOL <- 1e-8

mse <- function(x) { x <- x[is.finite(x)]; n <- length(x)
  c(mean = if (n) mean(x) else NA_real_, se = if (n > 1L) sd(x) / sqrt(n) else NA_real_, n = n) }
pm  <- function(m, se, d = 3L) ifelse(is.finite(m), sprintf(paste0("%.", d, "f ± %.", d, "f"), m, se), "-")
section <- function(t) cat("\n", strrep("=", 110), "\n", t, "\n", strrep("=", 110), "\n", sep = "")
show <- function(d) print(d, row.names = FALSE, right = FALSE)

# ---- data ----------------------------------------------------------------------
L2 <- readRDS(file.path(I6, "combined_replicate_metrics.rds"))
L3 <- readRDS(file.path(I6, "combined_scenario_metrics.rds"))
need <- c(KEY, "relationship", "annotation_type", "method", "ap")
if (length(setdiff(need, names(L2)))) stop("L2 lacks: ", paste(setdiff(need, names(L2)), collapse = ", "), call. = FALSE)
if (!identical(sort(unique(L2$annotation_type)), "continuous")) stop("expected continuous annotations only", call. = FALSE)
bad <- c(setdiff(METHODS, unique(L2$method)), setdiff(unique(L2$method), METHODS),
         setdiff(REL, unique(L2$relationship)), setdiff(unique(L2$relationship), REL))
if (length(bad)) stop("design mismatch: ", paste(bad, collapse = ", "), call. = FALSE)
units <- function(m, r) L2[L2$method == m & L2$relationship == r, c(KEY, "ap")]
paired <- function(a, b, r) {
  x <- merge(units(a, r), units(b, r), by = KEY, suffixes = c("_a", "_b"))
  mse(x$ap_a - x$ap_b)
}

sink(file.path(OUT, "iter006_report.txt"), split = TRUE)
cat("Iteration 006 - annotation-to-causality relationships, continuous annotations -",
    format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("source:", normalizePath(I6), "\n")

# ---- 0. collected ----------------------------------------------------------------
section("0. What was collected")
per <- with(L2, tapply(ap, list(method, relationship), function(v) sum(is.finite(v))))
cat(sprintf("scenario draws per method per relationship: %s\n", paste(sort(unique(as.vector(per))), collapse = "/")))
cat(sprintf("L3 cells: %d; failed fits: %d of %d\n", nrow(L3), sum(L3$n_failed), sum(L3$n_fits)))
xo <- merge(units("polyfun_oracle", "null"), units("susie", "null"), by = KEY)
cat(sprintf("polyfun_oracle equals susie under null (uniform prior) in %.1f%% of draws - the sanity check\n",
            100 * mean(abs(xo$ap.x - xo$ap.y) < 1e-9)))

# ---- 1. average precision ----------------------------------------------------
section("1. Average precision, mean ± 2 SE over the 100 scenario draws in each relationship")
ap <- do.call(rbind, lapply(REL, function(r) do.call(rbind, lapply(METHODS, function(m) {
  v <- mse(units(m, r)$ap); data.frame(relationship = r, method = m, ap = v[["mean"]], se = v[["se"]], n = v[["n"]])
}))))
tab <- data.frame(row = c("learnable: log-linear prior", "learnable: any own-annotation prior", METHODS))
for (r in REL) {
  l <- LEARN[LEARN$relationship == r, ]
  tab[[r]] <- c(ifelse(is.na(l$linear), "-", sprintf("%.2f", l$linear)), ifelse(is.na(l$own), "-", sprintf("%.2f", l$own)),
                with(ap[ap$relationship == r, ], pm(ap, 2 * se))[match(METHODS, ap$method[ap$relationship == r])])
}
show(tab)

# ---- 2. share of the oracle's gain -------------------------------------------------
section("2. Share of the oracle's gain over susie recovered: mean(method - susie) / mean(oracle - susie), ± 2 SE")
share <- do.call(rbind, lapply(REL, function(r) {
  base <- units("susie", r); orc <- units("polyfun_oracle", r)
  do.call(rbind, lapply(setdiff(METHODS, c("susie", "polyfun_oracle")), function(m) {
    x <- merge(merge(units(m, r), base, by = KEY, suffixes = c("", "_susie")), orc, by = KEY, suffixes = c("", "_oracle"))
    d <- x$ap - x$ap_susie; h <- x$ap_oracle - x$ap_susie; R <- mean(d) / mean(h)
    data.frame(relationship = r, method = m, headroom = mean(h), headroom_se = sd(h) / sqrt(nrow(x)),
               share = R, share_se = sd(d - R * h) / (sqrt(nrow(x)) * abs(mean(h))), n = nrow(x))
  }))
}))
hr <- unique(share[, c("relationship", "headroom", "headroom_se")])
tab2 <- data.frame(row = c("headroom: oracle - susie", setdiff(METHODS, c("susie", "polyfun_oracle"))))
for (r in REL) {
  h <- hr[hr$relationship == r, ]; s <- share[share$relationship == r, ]
  # a share is only meaningful where the headroom itself is clearly positive
  ok <- h$headroom > 2 * h$headroom_se
  tab2[[r]] <- c(sprintf("%+.3f (SE %.3f)", h$headroom, h$headroom_se),
                 if (ok) sprintf("%4.0f%% ± %.0f", 100 * s$share, 200 * s$share_se)[match(tab2$row[-1], s$method)]
                 else rep("-", nrow(tab2) - 1L))
}
show(tab2)
cat("A share is shown only where the headroom exceeds 2 SE; under null there is no gain to recover.\n")
cat("CAUTION: susie, the floor here, ties the causal variant in 22-33% of fits and this AP scores a tie as a block,\n",
    "so the headroom and these shares are inflated. With ties broken at random the headroom roughly halves:\n",
    "see iter006_ap_ties.txt.\n", sep = "")

# ---- 3. paired contrasts --------------------------------------------------------
section("3. Paired differences in average precision, draw by draw: mean (SE)")
con <- do.call(rbind, lapply(CONTRASTS, function(ab) do.call(rbind, lapply(REL, function(r) {
  v <- paired(ab[1], ab[2], r)
  data.frame(contrast = paste(ab[1], "-", ab[2]), relationship = r, diff = v[["mean"]], se = v[["se"]], n = v[["n"]])
}))))
tab3 <- data.frame(contrast = unique(con$contrast))
for (r in REL) { s <- con[con$relationship == r, ]; tab3[[r]] <- sprintf("%+.3f (%.3f)", s$diff, s$se)[match(tab3$contrast, s$contrast)] }
show(tab3)

# ---- 4. annotation selection ----------------------------------------------------
section("4. Annotation selection by the BEATRICE family (tracks 1-5 informative; chance 0.5)")
rel_of <- setNames(unique(L3[, c("job_dir", "relationship")])$relationship, unique(L3[, c("job_dir", "relationship")])$job_dir)
read_imp <- function(method) {
  fs <- sort(list.files(file.path(I6, "aux"), "^aux_.*[.]rds$", full.names = TRUE))
  parts <- lapply(fs, function(f) {
    o <- readRDS(f); im <- Filter(function(x) identical(x$method, method), o$importance)
    if (!length(im)) return(NULL)
    V <- matrix(vapply(im, function(x) { fi <- x$importance; v <- rep(NA_real_, N_ANNOT)
      ix <- as.integer(fi$annotation_index); if (min(ix) == 0L) ix <- ix + 1L
      v[ix] <- as.numeric(fi$importance); v }, numeric(N_ANNOT)), ncol = N_ANNOT, byrow = TRUE)
    list(meta = data.frame(job_dir = o$job_dir, scenario_id = vapply(im, function(x) as.integer(x$scenario_id), 0L),
                           region_id = vapply(im, function(x) as.integer(x$region_id), 0L),
                           fallback = vapply(im, function(x) isTRUE(x$joint_fallback), TRUE)), V = V)
  })
  parts <- Filter(Negate(is.null), parts)
  list(meta = do.call(rbind, lapply(parts, `[[`, "meta")), V = do.call(rbind, lapply(parts, `[[`, "V")))
}
rank_metrics <- function(V) {
  V[V < ZERO_TOL] <- 0
  as.data.frame(t(apply(V, 1, function(v) {
    z <- v == 0; th <- sort(v, decreasing = TRUE)[5]; above <- v > th; tied <- v == th
    p5 <- (sum(TRUTH & above) + (5 - sum(above)) * mean(TRUTH[tied])) / 5
    rk <- rank(-v, ties.method = "average")
    c(zeroed = sum(z), precision5 = p5, auc = (25 + 15 - sum(rk[TRUTH])) / 25)
  })))
}
imp_units <- list()
fbx <- read_imp("fb_xregion")
ok  <- which(!fbx$meta$fallback & rowSums(is.na(fbx$V)) == 0L)
key <- paste(fbx$meta$job_dir, fbx$meta$scenario_id)[ok]
first <- ok[match(key, key)]
n_disagree <- sum(tapply(apply(abs(fbx$V[ok, , drop = FALSE] - fbx$V[first, , drop = FALSE]), 1, max) > 1e-10, key, any))
keep <- ok[!duplicated(key)]
imp_units$fb_xregion <- cbind(fbx$meta[keep, c("job_dir", "scenario_id")], method = "fb_xregion",
                              rank_metrics(fbx$V[keep, , drop = FALSE]))
fbr <- read_imp("functional_beatrice"); ok2 <- which(rowSums(is.na(fbr$V)) == 0L)
imp_units$functional_beatrice <- cbind(fbr$meta[ok2, c("job_dir", "scenario_id")], method = "functional_beatrice",
                                       rank_metrics(fbr$V[ok2, , drop = FALSE]))
IMP <- do.call(rbind, imp_units); IMP$relationship <- rel_of[IMP$job_dir]
cat(sprintf("fb_xregion: %d fallback fits excluded; scenarios whose ten regions disagree on the shared vector: %d\n",
            sum(fbx$meta$fallback), n_disagree))
imp <- do.call(rbind, lapply(c("fb_xregion", "functional_beatrice"), function(m) do.call(rbind, lapply(REL, function(r) {
  d <- IMP[IMP$method == m & IMP$relationship == r, ]; z <- mse(d$zeroed); p <- mse(d$precision5); a <- mse(d$auc)
  data.frame(method = m, relationship = r, n = nrow(d), zeroed = z[["mean"]], precision5 = p[["mean"]],
             precision5_se = p[["se"]], auc = a[["mean"]], auc_se = a[["se"]])
}))))
for (m in c("fb_xregion", "functional_beatrice")) {
  cat(sprintf("\n-- %s (%s)\n", m, if (m == "fb_xregion") "shared head, one vector per scenario" else "per region"))
  s <- imp[imp$method == m, ]
  t4 <- data.frame(row = c("units", "tracks zeroed /10", "precision@5 ± 2 SE", "AUC ± 2 SE"))
  for (r in REL) { x <- s[s$relationship == r, ]
    t4[[r]] <- c(x$n, sprintf("%.2f", x$zeroed), pm(x$precision5, 2 * x$precision5_se), pm(x$auc, 2 * x$auc_se)) }
  show(t4)
}

nb <- IMP[IMP$method == "fb_xregion" & IMP$relationship == "null", ]
zs <- function(v) (mean(v) - 0.5) / (sd(v) / sqrt(length(v)))
cat(sprintf("\nfb_xregion under null, where no track is informative: precision@5 %.3f (%.1f SE above 0.5), AUC %.3f (%.1f SE above 0.5).\n",
            mean(nb$precision5), zs(nb$precision5), mean(nb$auc), zs(nb$auc)))
Vnull <- fbx$V[keep, , drop = FALSE][rel_of[fbx$meta$job_dir[keep]] == "null", , drop = FALSE]
cat("mean |contrast| per track under null:", paste(sprintf("t%d %.4f", 1:10, colMeans(Vnull)), collapse = "  "), "\n")
cat("Every fit starts from the same weights (torch.manual_seed(1) in joint_trainer.py), so read the\n",
    "other relationships against this null row, not against 0.5.\n", sep = "")

# ---- 5. relative to null --------------------------------------------------------
section("5. Advantage over susie in each relationship, minus the same advantage under null: mean (SE)")
cat("Under null the annotations carry no information, so an advantage over susie there is not an annotation gain.\n",
    "The rows are independent simulations, so SE = sqrt(SE_relationship^2 + SE_null^2).\n", sep = "")
cat("With tied PIPs broken at random (iter006_ap_ties.txt), polyfun_ldsc's null-arm gap over susie vanishes (+0.004)\n",
    "and funmap's reverses (-0.019): those gaps are susie's ties, not annotation use.\n", sep = "")
did_rows <- c(setdiff(METHODS, "susie"), "fb_xregion - polyfun_ldsc")
pair_of  <- function(lab) if (grepl(" - ", lab)) strsplit(lab, " - ")[[1]] else c(lab, "susie")
null_gap <- do.call(rbind, lapply(did_rows, function(lab) { ab <- pair_of(lab); b <- paired(ab[1], ab[2], "null")
  data.frame(row = lab, null_gap = b[["mean"]], null_se = b[["se"]]) }))
did <- do.call(rbind, lapply(did_rows, function(lab) do.call(rbind, lapply(setdiff(REL, "null"), function(r) {
  ab <- pair_of(lab); a <- paired(ab[1], ab[2], r); b <- paired(ab[1], ab[2], "null")
  data.frame(row = lab, relationship = r, did = a[["mean"]] - b[["mean"]], se = sqrt(a[["se"]]^2 + b[["se"]]^2))
}))))
tab5 <- data.frame(row = did_rows, "gap under null" = sprintf("%+.3f (%.3f)", null_gap$null_gap, null_gap$null_se), check.names = FALSE)
for (r in setdiff(REL, "null")) { s5 <- did[did$relationship == r, ]; tab5[[r]] <- sprintf("%+.3f (%.3f)", s5$did, s5$se)[match(tab5$row, s5$row)] }
show(tab5)

sink()
saveRDS(list(ap = ap, share = share, contrasts = con, importance = imp, importance_units = IMP,
             learnability = LEARN, shared_head_disagreements = n_disagree,
             null_gaps = null_gap, null_referenced = did, fb_xregion_null_track_means = colMeans(Vnull)),
        file.path(OUT, "iter006_results.rds"))
message("wrote ", file.path(OUT, "iter006_report.txt"), " and iter006_results.rds")
