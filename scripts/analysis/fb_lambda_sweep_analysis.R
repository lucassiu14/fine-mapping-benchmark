#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/fb_lambda_sweep_analysis.R
#
# The lambda_l1 sweep of the cross-region model: eight fb_xregion arms that
# differ only in the prior head's parameterisation and lambda_l1, fitted to
# Iteration 004's own simulations (PBS array 4005367 plus ten reruns) and
# collected with SOURCE=supp into results/fb_lambda_sweep.
#
#   Rscript scripts/analysis/fb_lambda_sweep_analysis.R [sweep_dir] [iter004_dir] [out_dir]
#
# Defaults: results/fb_lambda_sweep, results/iter004, <sweep_dir>/analysis.
# Iteration 004's files are read, never written.
#
# Four questions, each answered per model x annotation type and never pooled
# across them:
#   1. Runtime       seconds per scenario for each arm, against Iteration 004's
#                    fb_xregion on the same scenarios.
#   2. Annotations   how many of the ten tracks each arm zeroes, and how well it
#                    ranks the five enriched tracks above the five inert ones.
#   3. Accuracy      average precision for each arm, and paired against
#                    Iteration 004's fb_xregion, beatrice and polyfun_ldsc.
#   4. Fit by fit    the rerun against Iteration 004's fb_xregion, and _l200
#                    against beatrice.
#
# Conventions, with their sources:
#   * The no-annotation rows are left out. There every arm falls back to
#     per-region Functional BEATRICE with no annotations, which builds no
#     LassoNet (trainer_annot.py, annot_given = False), and the rows were not
#     completed.
#   * The sweep's own fb_xregion is renamed fb_xregion_rerun. fb_xregion alone
#     always means Iteration 004's fit.
#   * Average precision follows iter004_lsr_figures.R: L3 cells with
#     n_failed < n_fits, mean and SE across cells. Paired differences are taken
#     cell by cell on job_dir x S x phi x region_size. Printed as the LSR text
#     prints them: a mean carries +- 2 SE ("0.793 +- 0.013" is two SEs of 0.0065),
#     a paired difference carries one SE in brackets.
#   * Runtime is aux's scenario_total, which includes the cross-region training:
#     run_methods() starts its clock before the scenario-setup hook that runs the
#     joint trainer. Per-region times are NA for joint fits. The two runs were
#     weeks apart and may have landed on different node types.
#   * The first five of ten annotation tracks are enriched: generate_params_grid.R
#     builds c(rep(E, 5), rep(1, 5)), and simulate_gwfm_data.R weights variant j
#     by exp(sum_k A_jk log enrichment_k), so element k belongs to column k.
#   * Importance is |logit contrast| per track from the SHARED head, written
#     identically into all ten regions of a scenario, so it is scored once per
#     scenario (and the ten copies are checked to agree). Fits that fell back to
#     single-region Functional BEATRICE (joint_fallback = TRUE) do not carry a
#     cross-region importance: they are counted and excluded. "Zeroed" is
#     importance < 1e-8, as in the live check during the run.
#   * Ranking metrics break ties at random, in expectation. A large lambda_l1
#     zeroes tracks, and order() would put tied low-index - that is, enriched -
#     tracks first, scoring _l200 (all ten zeroed) a perfect precision@5.
#   * _l200 need not reproduce beatrice exactly. With every skip weight at zero
#     the hierarchy constraint zeroes the first hidden layer, so the prior stops
#     varying between variants - but its level is set by the hidden-layer
#     biases and learned, where beatrice fixes it at 1/p.
# =============================================================================

args  <- commandArgs(trailingOnly = TRUE)
SWEEP <- if (length(args) >= 1L) args[1] else "results/fb_lambda_sweep"
I004  <- if (length(args) >= 2L) args[2] else "results/iter004"
OUT   <- if (length(args) >= 3L) args[3] else file.path(SWEEP, "analysis")
for (d in c(SWEEP, I004)) if (!dir.exists(d)) stop("no such folder: ", d, call. = FALSE)
if (grepl("(^|/)iter004(/|$)", OUT))
  stop("refusing to write into Iteration 004's folder: ", OUT, call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
options(width = 250, stringsAsFactors = FALSE)

# ---- the design --------------------------------------------------------------
LAMBDA_PROD <- 0.1223   # METHOD_ARGS$functional_beatrice$lambda_l1, run_benchmark_job.R
ARMS <- data.frame(
  method    = c("fb_xregion_id_l03", "fb_xregion_id_l06", "fb_xregion_id", "fb_xregion_rerun",
                "fb_xregion_id_l25", "fb_xregion_id_l50", "fb_xregion_id_l100", "fb_xregion_id_l200"),
  head      = c("identifiable", "identifiable", "identifiable", "standard",
                "identifiable", "identifiable", "identifiable", "identifiable"),
  lambda_l1 = c(0.03, 0.06, LAMBDA_PROD, LAMBDA_PROD, 0.25, 0.5, 1.0, 2.0))
REF          <- c("fb_xregion", "functional_beatrice", "beatrice", "polyfun_ldsc", "susie")
PAIR_REF     <- c("fb_xregion", "beatrice", "polyfun_ldsc")
FB_REF_LABEL <- "functional_beatrice (iter004, per region)"
METHOD_ORDER <- c(ARMS$method, REF, FB_REF_LABEL)
N_ANNOT  <- 10L
TRUTH    <- rep(c(TRUE, FALSE), each = N_ANNOT %/% 2L)
ZERO_TOL <- 1e-8
KEY_CELL <- c("job_dir", "S", "phi", "region_size")
KEY_FIT  <- c("job_dir", "scenario_id", "region_id")
STRATA   <- data.frame(model           = c("sparse", "sparse", "sparse_inf", "sparse_inf"),
                       annotation_type = c("binary", "continuous", "binary", "continuous"))
ENRICH   <- c("2.7", "5.4", "8.1", "10.8")

# ---- helpers -----------------------------------------------------------------
rename_sweep <- function(m) replace(m, m == "fb_xregion", "fb_xregion_rerun")

mse <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  c(mean = if (n) mean(x) else NA_real_, se = if (n > 1L) sd(x) / sqrt(n) else NA_real_, n = n)
}
pm <- function(m, se, d = 3L)
  ifelse(is.finite(m), sprintf(paste0("%.", d, "f ± %.", d, "f"), m, se), "-")

# f applied to each model x annotation-type stratum on its own - never pooled.
per_stratum <- function(d, f) {
  do.call(rbind, lapply(seq_len(nrow(STRATA)), function(i) {
    s <- d[d$model == STRATA$model[i] & d$annotation_type == STRATA$annotation_type[i], , drop = FALSE]
    if (!nrow(s)) return(NULL)
    r <- f(s)
    if (is.null(r) || !nrow(r)) return(NULL)
    data.frame(model = STRATA$model[i], annotation_type = STRATA$annotation_type[i], r,
               row.names = NULL, check.names = FALSE)
  }))
}
per_method <- function(s, f) do.call(rbind, lapply(split(s, s$method), f))
ord_methods <- function(d) d[order(match(d$method, METHOD_ORDER)), , drop = FALSE]

section <- function(title)
  cat("\n", strrep("=", 100), "\n", title, "\n", strrep("=", 100), "\n", sep = "")
show_by_stratum <- function(d) {
  if (is.null(d) || !nrow(d)) { cat("  (nothing)\n"); return(invisible(NULL)) }
  for (i in seq_len(nrow(STRATA))) {
    s <- d[d$model == STRATA$model[i] & d$annotation_type == STRATA$annotation_type[i], , drop = FALSE]
    if (!nrow(s)) next
    cat("\n-- ", STRATA$model[i], " / ", STRATA$annotation_type[i], "\n", sep = "")
    s <- ord_methods(s); s$model <- NULL; s$annotation_type <- NULL
    print(s, row.names = FALSE, right = FALSE)
  }
}

# ---- data --------------------------------------------------------------------
message("reading the sweep ...")
L3s <- readRDS(file.path(SWEEP, "combined_scenario_metrics.rds"))
L1s <- readRDS(file.path(SWEEP, "combined_fit_metrics.rds"))
L3s$method <- rename_sweep(L3s$method)
L1s$method <- rename_sweep(L1s$method)
bad <- c(setdiff(ARMS$method, L3s$method), setdiff(unique(L3s$method), ARMS$method))
if (length(bad))
  stop("the sweep's methods are not the eight arms: ", paste(bad, collapse = ", "), call. = FALSE)

design <- unique(L3s[, c("job_dir", "model", "annotation_type", "enrich", "pc")])
JOBS   <- design$job_dir[design$annotation_type != "none"]
L3s <- L3s[L3s$job_dir %in% JOBS, ]
L1s <- L1s[L1s$job_dir %in% JOBS, c(KEY_FIT, "model", "annotation_type", "method", "failed", "ap")]

message("reading Iteration 004 ...")
L3r <- readRDS(file.path(I004, "combined_scenario_metrics.rds"))
L3r <- L3r[L3r$method %in% REF & L3r$job_dir %in% JOBS, ]
L1r <- readRDS(file.path(I004, "combined_fit_metrics.rds"))
L1r <- L1r[L1r$method %in% REF & L1r$job_dir %in% JOBS, c(KEY_FIT, "method", "ap")]

aux_files <- function(dir) {
  fs <- sort(list.files(dir, "^aux_.*[.]rds$", full.names = TRUE))
  if (!length(fs)) stop("no aux_*.rds in ", dir, call. = FALSE)
  fs
}
read_runtime <- function(dir) {
  do.call(rbind, lapply(aux_files(dir), function(f) {
    o <- readRDS(f)
    if (!o$job_dir %in% JOBS) return(NULL)
    rt <- Filter(function(x) identical(x$scope, "scenario_total"), o$runtime)
    if (!length(rt)) return(NULL)
    data.frame(job_dir     = o$job_dir,
               scenario_id = vapply(rt, function(x) as.integer(x$scenario_id), 0L),
               method      = vapply(rt, function(x) as.character(x$method), ""),
               seconds     = vapply(rt, function(x) as.numeric(x$runtime_seconds), 0))
  }))
}
# Importance records -> list(meta = data.frame, V = n x 10 matrix, column k = track k).
read_importance <- function(dir, methods) {
  parts <- lapply(aux_files(dir), function(f) {
    o <- readRDS(f)
    if (!o$job_dir %in% JOBS) return(NULL)
    im <- Filter(function(x) x$method %in% methods, o$importance)
    if (!length(im)) return(NULL)
    V <- matrix(vapply(im, function(x) {
      fi <- x$importance; v <- rep(NA_real_, N_ANNOT)
      if (is.data.frame(fi) && all(c("annotation_index", "importance") %in% names(fi)) && nrow(fi)) {
        ix <- as.integer(fi$annotation_index)
        if (min(ix, na.rm = TRUE) == 0L) ix <- ix + 1L          # 0-based in the CSV
        ok <- !is.na(ix) & ix >= 1L & ix <= N_ANNOT
        v[ix[ok]] <- as.numeric(fi$importance)[ok]
      }
      v
    }, numeric(N_ANNOT)), ncol = N_ANNOT, byrow = TRUE)
    list(meta = data.frame(job_dir     = o$job_dir,
                           scenario_id = vapply(im, function(x) as.integer(x$scenario_id), 0L),
                           region_id   = vapply(im, function(x) as.integer(x$region_id), 0L),
                           method      = vapply(im, function(x) as.character(x$method), ""),
                           fallback    = vapply(im, function(x) isTRUE(x$joint_fallback), TRUE)),
         V = V)
  })
  parts <- Filter(Negate(is.null), parts)
  list(meta = do.call(rbind, lapply(parts, `[[`, "meta")),
       V    = do.call(rbind, lapply(parts, `[[`, "V")))
}
# Per importance vector: zeroed counts and tie-aware ranking of the enriched tracks.
rank_metrics <- function(V) {
  V[V < ZERO_TOL] <- 0
  ne <- sum(TRUTH); ni <- N_ANNOT - ne
  out <- t(apply(V, 1, function(v) {
    z  <- v == 0
    th <- sort(v, decreasing = TRUE)[ne]
    above <- v > th; tied <- v == th
    p5 <- (sum(TRUTH & above) + (ne - sum(above)) * mean(TRUTH[tied])) / ne
    rk <- rank(-v, ties.method = "average")
    auc <- (ne * ni + ne * (ne + 1) / 2 - sum(rk[TRUTH])) / (ne * ni)
    c(zeroed = sum(z), zeroed_enriched = sum(z & TRUTH), zeroed_inert = sum(z & !TRUTH),
      all_zero = as.numeric(all(z)), precision5 = p5, auc = auc)
  }))
  as.data.frame(out)
}

report <- file.path(OUT, "sweep_report.txt")
sink(report, split = TRUE)
cat("fb_xregion lambda_l1 sweep -", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("sweep:", normalizePath(SWEEP), "\niteration 004:", normalizePath(I004), "\n")
cat("Every table is per model x annotation type. The no-annotation rows are left out.\n\n")
print(ARMS, row.names = FALSE, right = FALSE)

# ---- 0. what was collected ---------------------------------------------------
section("0. What was collected (annotated rows)")
collected <- per_stratum(L1s, function(s) per_method(s, function(m) data.frame(
  method = m$method[1], rows = length(unique(m$job_dir)),
  scenarios = nrow(unique(m[c("job_dir", "scenario_id")])),
  fits = nrow(m), failed = sum(m$failed %in% TRUE | !is.finite(m$ap)))))
show_by_stratum(collected)

# ---- 1. runtime --------------------------------------------------------------
section("1. Runtime: seconds per scenario (all ten regions), and x Iteration 004's fb_xregion")
rt_s <- read_runtime(file.path(SWEEP, "aux"))
# aux lays results_supp.rds over results.rds, so its fb_xregion record is the
# rerun only where the scenario has a sweep result. Keep it only where the
# identifiable arm - which exists nowhere else - is present too.
has_sweep <- unique(rt_s[rt_s$method == "fb_xregion_id", c("job_dir", "scenario_id")])
rt_s <- merge(rt_s[rt_s$method %in% c("fb_xregion", ARMS$method), ], has_sweep,
              by = c("job_dir", "scenario_id"))
rt_s$method <- rename_sweep(rt_s$method)
rt_r <- read_runtime(file.path(I004, "aux"))
rt_r <- rt_r[rt_r$method %in% c("fb_xregion", "functional_beatrice", "beatrice"), ]
rt <- merge(rbind(rt_s, rt_r), design[, c("job_dir", "model", "annotation_type")], by = "job_dir")

base <- rt[rt$method == "fb_xregion", c("job_dir", "scenario_id", "seconds")]
names(base)[3] <- "base_s"
rp <- merge(rt[rt$method %in% ARMS$method, ], base, by = c("job_dir", "scenario_id"))
rp$ratio <- rp$seconds / rp$base_s

runtime <- per_stratum(rt, function(s) per_method(s, function(m) {
  q <- quantile(m$seconds, c(.25, .5, .75), na.rm = TRUE, names = FALSE)
  p <- rp[rp$method == m$method[1] & rp$model == m$model[1] &
          rp$annotation_type == m$annotation_type[1] & is.finite(rp$ratio), ]
  r <- if (nrow(p)) quantile(p$ratio, c(.25, .5, .75), names = FALSE) else rep(NA_real_, 3)
  data.frame(method = m$method[1], scenarios = sum(is.finite(m$seconds)),
             median_s = q[2], q25_s = q[1], q75_s = q[3], mean_s = mean(m$seconds, na.rm = TRUE),
             cpu_hours = sum(m$seconds, na.rm = TRUE) / 3600,
             pairs = nrow(p), ratio_median = r[2], ratio_q25 = r[1], ratio_q75 = r[3],
             ratio_of_totals = if (nrow(p)) sum(p$seconds) / sum(p$base_s) else NA_real_)
}))
show_by_stratum(with(runtime, data.frame(
  model, annotation_type, method, scenarios,
  "median s (IQR)" = sprintf("%.0f (%.0f-%.0f)", median_s, q25_s, q75_s),
  "mean s" = sprintf("%.0f", mean_s), "CPU-h" = sprintf("%.0f", cpu_hours),
  "x iter004 fb_xregion, median (IQR)" =
    ifelse(is.finite(ratio_median), sprintf("%.2f (%.2f-%.2f)", ratio_median, ratio_q25, ratio_q75), "-"),
  "total / iter004 total" = ifelse(is.finite(ratio_of_totals), sprintf("%.2f", ratio_of_totals), "-"),
  check.names = FALSE)))

# ---- 2. annotation importance ------------------------------------------------
section("2. Annotation importance from the shared head (tracks 1-5 enriched, 6-10 inert)")
imp <- read_importance(file.path(SWEEP, "aux"), c("fb_xregion", setdiff(ARMS$method, "fb_xregion_rerun")))
meta <- imp$meta; V <- imp$V
meta <- merge(cbind(meta, .i = seq_len(nrow(meta))), has_sweep, by = c("job_dir", "scenario_id"))
V <- V[meta$.i, , drop = FALSE]; meta$.i <- NULL
meta$method <- rename_sweep(meta$method)
meta <- merge(cbind(meta, .i = seq_len(nrow(meta))), design, by = "job_dir")
V <- V[meta$.i, , drop = FALSE]; meta$.i <- NULL

fits <- per_stratum(L1s, function(s) per_method(s, function(m) data.frame(method = m$method[1], fits = nrow(m))))
fallback <- per_stratum(meta, function(s) per_method(s, function(m) data.frame(
  method = m$method[1], fallback_fits = sum(m$fallback), importance_records = nrow(m))))
fallback <- merge(fallback, fits, by = c("model", "annotation_type", "method"), all = TRUE)
fallback$fallback_share <- fallback$fallback_fits / fallback$fits

joint <- which(!meta$fallback & rowSums(is.na(V)) == 0L)
key   <- paste(meta$job_dir, meta$scenario_id, meta$method, sep = "\r")[joint]
first <- joint[match(key, key)]
dev   <- apply(abs(V[joint, , drop = FALSE] - V[first, , drop = FALSE]), 1, max)
n_disagree <- sum(tapply(dev > 1e-10, key, any))
keep  <- joint[!duplicated(key)]
IM    <- cbind(meta[keep, c("job_dir", "scenario_id", "method", "model", "annotation_type", "enrich", "pc")],
               rank_metrics(V[keep, , drop = FALSE]), row.names = NULL)
cat(sprintf("\nscenarios whose ten regions do not carry the same shared-head vector: %d of %d\n",
            n_disagree, length(keep)))

fbr <- read_importance(file.path(I004, "aux"), "functional_beatrice")
fbr_ok <- rowSums(is.na(fbr$V)) == 0L
FBR <- cbind(merge(cbind(fbr$meta[fbr_ok, c("job_dir", "scenario_id", "region_id")], .i = which(fbr_ok)),
                   design, by = "job_dir"))
FBR <- cbind(FBR, rank_metrics(fbr$V[FBR$.i, , drop = FALSE]))
FBR$method <- FB_REF_LABEL; FBR$.i <- NULL

summ_imp <- function(d, unit) per_stratum(d, function(s) per_method(s, function(m) {
  z <- mse(m$zeroed); ze <- mse(m$zeroed_enriched); zi <- mse(m$zeroed_inert)
  p <- mse(m$precision5); a <- mse(m$auc)
  data.frame(method = m$method[1], unit = unit, n = nrow(m),
             zeroed = z[1], zeroed_se = z[2], zeroed_enriched = ze[1], zeroed_inert = zi[1],
             all_zero = mean(m$all_zero), precision5 = p[1], precision5_se = p[2],
             auc = a[1], auc_se = a[2])
}))
importance <- rbind(summ_imp(IM, "scenario"), summ_imp(FBR, "region"))
importance <- merge(importance, fallback[, c("model", "annotation_type", "method", "fallback_share")],
                    by = c("model", "annotation_type", "method"), all.x = TRUE)
show_by_stratum(with(importance, data.frame(
  model, annotation_type, method, unit, n,
  "fell back" = ifelse(is.finite(fallback_share), sprintf("%.1f%%", 100 * fallback_share), "-"),
  "zeroed /10 ± 2 SE" = pm(zeroed, 2 * zeroed_se, 2),
  "enriched zeroed /5" = sprintf("%.2f", zeroed_enriched),
  "inert zeroed /5" = sprintf("%.2f", zeroed_inert),
  "all ten zeroed" = sprintf("%.1f%%", 100 * all_zero),
  "precision@5 ± 2 SE" = pm(precision5, 2 * precision5_se), "AUC ± 2 SE" = pm(auc, 2 * auc_se),
  check.names = FALSE)))
cat("\nchance: precision@5 = 0.5, AUC = 0.5. Ties are broken at random in expectation, so an arm\n",
    "that zeroes all ten tracks scores exactly chance.\n", sep = "")

cat("\n-- precision@5 / AUC by enrichment fold\n")
by_enrich <- per_stratum(rbind(IM[, c("model", "annotation_type", "method", "enrich", "precision5", "auc")],
                               FBR[, c("model", "annotation_type", "method", "enrich", "precision5", "auc")]),
  function(s) do.call(rbind, lapply(split(s, list(s$method, s$enrich), drop = TRUE), function(m)
    data.frame(method = m$method[1], enrich = m$enrich[1], n = nrow(m),
               precision5 = mean(m$precision5), auc = mean(m$auc)))))
show_by_stratum(do.call(rbind, lapply(split(by_enrich, list(by_enrich$model, by_enrich$annotation_type, by_enrich$method), drop = TRUE),
  function(m) {
    row <- data.frame(model = m$model[1], annotation_type = m$annotation_type[1], method = m$method[1])
    for (e in ENRICH) {
      k <- match(e, m$enrich)
      row[[paste0("fold ", e)]] <- if (is.na(k)) "-" else sprintf("%.2f / %.2f", m$precision5[k], m$auc[k])
    }
    row
  })))

# ---- 3. accuracy ---------------------------------------------------------------
section("3. Average precision: mean ± 2 SE across cells, and paired differences (SE) cell by cell")
cols  <- c(KEY_CELL, "model", "annotation_type", "method", "ap", "n_fits", "n_failed")
cells <- rbind(L3s[, cols], L3r[, cols])
cells <- cells[cells$n_failed < cells$n_fits, ]
ap <- per_stratum(cells, function(s) per_method(s, function(m) {
  a <- mse(m$ap); data.frame(method = m$method[1], cells = a[[3]], ap = a[[1]], ap_se = a[[2]])
}))
paired <- do.call(rbind, lapply(PAIR_REF, function(ref) {
  b <- cells[cells$method == ref, c(KEY_CELL, "ap")]; names(b)[names(b) == "ap"] <- "ap_ref"
  x <- merge(cells[cells$method %in% ARMS$method, ], b, by = KEY_CELL)
  x$diff <- x$ap - x$ap_ref
  per_stratum(x, function(s) per_method(s, function(m) {
    d <- mse(m$diff); data.frame(method = m$method[1], versus = ref, cells = d[[3]], diff = d[[1]], diff_se = d[[2]])
  }))
}))
ap_disp <- with(ap, data.frame(model, annotation_type, method, cells,
                               "AP ± 2 SE" = pm(ap, 2 * ap_se), check.names = FALSE))
for (ref in PAIR_REF) {
  p <- paired[paired$versus == ref, ]
  k <- match(paste(ap_disp$model, ap_disp$annotation_type, ap_disp$method),
             paste(p$model, p$annotation_type, p$method))
  ap_disp[[paste("Δ vs", ref)]] <- ifelse(is.na(k), "-", sprintf("%+.3f (SE %.3f)", p$diff[k], p$diff_se[k]))
}
names(ap_disp)[names(ap_disp) == paste("Δ vs", "fb_xregion")] <- "Δ vs fb_xregion (iter004)"
show_by_stratum(ap_disp)

# ---- 4. fit by fit -------------------------------------------------------------
section("4. Fit by fit: the rerun against Iteration 004, and _l200 against beatrice")
fit_compare <- function(a, b) {
  x <- merge(L1s[L1s$method == a, c(KEY_FIT, "model", "annotation_type", "ap")],
             L1r[L1r$method == b, c(KEY_FIT, "ap")], by = KEY_FIT, suffixes = c("", "_ref"))
  x <- x[is.finite(x$ap) & is.finite(x$ap_ref), ]
  per_stratum(x, function(s) data.frame(
    comparison = paste(a, "vs", b), fits = nrow(s),
    identical_ap = mean(abs(s$ap - s$ap_ref) < 1e-9), mean_diff = mean(s$ap - s$ap_ref),
    mean_abs_diff = mean(abs(s$ap - s$ap_ref)),
    correlation = if (nrow(s) > 2L) cor(s$ap, s$ap_ref) else NA_real_))
}
fitwise <- rbind(fit_compare("fb_xregion_rerun", "fb_xregion"),
                 fit_compare("fb_xregion_id_l200", "beatrice"))
print(within(fitwise, {
  identical_ap <- sprintf("%.1f%%", 100 * identical_ap); mean_diff <- sprintf("%+.4f", mean_diff)
  mean_abs_diff <- sprintf("%.4f", mean_abs_diff); correlation <- sprintf("%.4f", correlation)
}), row.names = FALSE, right = FALSE)

sink()
saveRDS(list(arms = ARMS, collected = collected, runtime_per_scenario = rt, runtime = runtime,
             fallback = fallback, importance_per_scenario = IM, importance_fb_iter004 = FBR,
             importance = importance, importance_by_enrichment = by_enrich,
             shared_head_disagreements = n_disagree, ap = ap, ap_paired = paired, fitwise = fitwise),
        file.path(OUT, "sweep_analysis.rds"))
message("wrote ", report, " and ", file.path(OUT, "sweep_analysis.rds"))
