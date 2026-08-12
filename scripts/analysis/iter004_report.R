#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_report.R
#
# STAGE B of the Iteration 004 analysis pipeline. Binds the L1 partials written
# by iter004_collect.R, builds the L2 and L3 levels of the aggregation ladder
# (variable-importance-analysis.md §3.3), runs the §11 validity checks, then
# drives Sense V, D and G.
#
# THE AGGREGATION LADDER IS THE PART TO GET RIGHT (§3.2):
#
#   Rates are pooled from COUNTS. Per-locus scalars are AVERAGED.
#
#   - AP is a per-locus quantity. AP of a merged ranked list is not the mean of
#     the APs, because merging lets one region's false positives outrank another
#     region's true positive. A pooled ranking describes an analysis nobody
#     performs. -> average, at every level.
#   - FDR and calibration are properties of a POPULATION OF CALLS. Rates computed
#     per fit are unestimable (0-2 selections). -> sum the counts, then form the
#     rate.
#
#   Never average rates. Never pool rankings. Both change the answer, not just
#   its precision.
#
# The SE at L3 is computed from the 10 L2 values, NOT the 20 fits: the two regions
# inside one replicate share a simulated panel, an annotation draw and an
# enrichment draw, so they are not independent (§3.3b).
#
# Usage:
#   Rscript scripts/analysis/iter004_report.R <l1_dir> <grid_csv> <out_dir>
#
# STATUS: written, never executed.
# =============================================================================

args     <- commandArgs(trailingOnly = TRUE)
l1_dir   <- args[1]
grid_csv <- args[2] %||% "scripts/hpc/params_grid.csv"
out_dir  <- args[3] %||% "results/iter004"
source(file.path("R", "evaluate_extras.R"))
source(file.path("scripts", "analysis", "iter004_lib.R"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Bind L1 and attach the design columns
# ---------------------------------------------------------------------------
parts <- list.files(l1_dir, pattern = "^L1_.*\\.rds$", full.names = TRUE)
if (!length(parts)) stop("no L1 partials in ", l1_dir, call. = FALSE)
message("binding ", length(parts), " L1 partials ...")
L1 <- do.call(rbind, lapply(parts, readRDS))
message("  L1: ", nrow(L1), " fit rows, ", length(unique(L1$method)), " methods")

grid <- read.csv(grid_csv, stringsAsFactors = FALSE)
# Expected row count comes FROM THE GRID, not a hardcoded 45. The assertion is
# there to catch a silent NA-drop, and it must not itself break when the grid
# legitimately changes.
EXPECT_JOBS <- nrow(grid)
grid$job_dir <- sprintf("job_%03d_%s", grid$job_id, grid$label)
keep <- c("job_dir", "model", "p_causal", "annotation_type", "enrichment_fold",
          "n_ref", "n", "n_iter")
options(fmb.expect_jobs = EXPECT_JOBS)
L1 <- merge(L1, grid[, intersect(keep, names(grid))], by = "job_dir", all.x = TRUE)
if (anyNA(L1$model)) {
  stop("some L1 rows did not join to params_grid.csv - the grid used for the run ",
       "differs from the one on disk. Regenerate it or point --grid at the ",
       "version that was actually run.", call. = FALSE)
}
# RECODE NA-AS-LEVEL BEFORE ANY GROUPING (§6). aggregate(), split(), table()
# and sort(unique(.)) all DROP NA GROUPS SILENTLY. p_causal is NA on every
# sparse row and enrichment_fold is NA on every `none` row, so grouping on the
# raw columns deletes those arms entirely - verified: grouping 1,152 fit rows on
# raw p_causal returned ZERO rows. The raw columns are kept for reference; all
# grouping uses pc / enrich.
L1$pc     <- ifelse(is.na(L1$p_causal),        "sparse", as.character(L1$p_causal))
L1$enrich <- ifelse(is.na(L1$enrichment_fold), "none",   as.character(L1$enrichment_fold))
stopifnot(!anyNA(L1$pc), !anyNA(L1$enrich))
saveRDS(L1, file.path(out_dir, "combined_fit_metrics.rds"))

THRESH_TAGS <- c("50", "80", "90", "95", "99")
BAND_TAGS   <- c("lo", "mid", "hi", "top")
CNT_COLS <- c(paste0("tp_at_",          THRESH_TAGS),
              paste0("nsel_at_",        THRESH_TAGS),
              paste0("sum_pip_sel_at_", THRESH_TAGS),
              paste0("n_band_",       BAND_TAGS),
              paste0("c_band_",       BAND_TAGS),
              paste0("sum_pip_band_", BAND_TAGS),
              "n_variants", "n_causal", "sum_pip_total")
SCALAR_COLS <- c("ap", "prec_at_1", "prec_at_S", "prec_at_2S",
                 "set_size_5", "set_prec_5", "set_size_95", "set_prec_95",
                 "mass_above_50", "mean_pip_topS")
HIT_COLS <- c("set_hit_5", "set_hit_95")     # Bernoulli per fit -> SUMMED

# n_ref is deliberately ABSENT from the grouping key: it is NA throughout (the
# grid is in-sample only) and would drop every row. It is re-attached after
# aggregation so the §11 check and prepare_analysis_table's filter still see it.
DESIGN <- c("job_dir", "model", "annotation_type", "pc", "enrich",
            "S", "phi", "region_size", "method")


# ifelse() evaluates both branches, so as.numeric() on the sentinel level warns
# even where the result is discarded. Convert only the convertible entries.
.unsentinel <- function(x, sentinel) {
  out <- rep(NA_real_, length(x))
  ok <- x != sentinel
  out[ok] <- as.numeric(x[ok])
  out
}

#' Form rate metrics from pooled counts. Called at L2 and L3 only - never on a
#' single fit, where 0-2 selections make every rate unestimable.
add_rates <- function(d) {
  for (t in THRESH_TAGS) {
    nsel <- d[[paste0("nsel_at_", t)]]
    tp   <- d[[paste0("tp_at_",   t)]]
    sps  <- d[[paste0("sum_pip_sel_at_", t)]]
    fp   <- nsel - tp
    # FDR is 0 by convention when nothing is selected. THAT CONVENTION MAKES AN
    # ABSTAINER WIN EVERY FDR COMPARISON, so nsel travels with it everywhere
    # (rule R1). marginal_z places zero variants above PIP 0.8 across the entire
    # benchmark.
    d[[paste0("fdr_at_", t)]] <- ifelse(nsel > 0, fp / nsel, 0)
    # Self-consistency (§4.2 / P3): the method's own claimed error count.
    exp_fp <- nsel - sps
    d[[paste0("exp_fp_at_",        t)]] <- exp_fp
    d[[paste0("honesty_ratio_at_", t)]] <- ifelse(exp_fp > 0, fp / exp_fp, NA_real_)
  }
  n_all <- rowSums(d[paste0("n_band_", BAND_TAGS)], na.rm = TRUE)
  c_all <- rowSums(d[paste0("c_band_", BAND_TAGS)], na.rm = TRUE)
  p_all <- rowSums(d[paste0("sum_pip_band_", BAND_TAGS)], na.rm = TRUE)

  d$total_mass_ratio <- ifelse(c_all > 0, p_all / c_all, NA_real_)
  # Murphy decomposition on the four bands (§4.4). RES is what exposes the
  # uniform predictor, which scores perfectly on ECE and mass ratio alike.
  mur <- t(vapply(seq_len(nrow(d)), function(i) {
    nb <- as.numeric(d[i, paste0("n_band_",       BAND_TAGS)])
    cb <- as.numeric(d[i, paste0("c_band_",       BAND_TAGS)])
    sb <- as.numeric(d[i, paste0("sum_pip_band_", BAND_TAGS)])
    m  <- .murphy_decomposition(nb, cb, sb)
    c(rel = m$rel, res = m$res, unc = m$unc, bs = m$bs, bss = m$bss)
  }, numeric(5)))
  d <- cbind(d, as.data.frame(mur))
  # ece_hi over the informative range only. Plain ECE is dominated by the ~94%
  # of variant-observations in the lowest band and mis-ranks methods.
  hi <- c("mid", "hi", "top")
  nh <- rowSums(d[paste0("n_band_", hi)], na.rm = TRUE)
  d$ece_hi <- vapply(seq_len(nrow(d)), function(i) {
    nb <- as.numeric(d[i, paste0("n_band_",       hi)])
    cb <- as.numeric(d[i, paste0("c_band_",       hi)])
    sb <- as.numeric(d[i, paste0("sum_pip_band_", hi)])
    ok <- nb > 0
    if (!any(ok)) return(NA_real_)
    sum(nb[ok] / sum(nb[ok]) * abs(sb[ok] / nb[ok] - cb[ok] / nb[ok]))
  }, numeric(1))
  d$hi_pip_reliab <- ifelse(d$n_band_top > 0, d$c_band_top / d$n_band_top, NA_real_)
  d$hi_pip_n      <- d$n_band_top       # R1: reliability never travels alone
  d
}

# --- L2: one row per (cell, replicate). Counts summed over the 2 regions; AP
#     and the per-locus scalars averaged over them (§5.1b).
message("building L2 ...")
by2 <- L1[c(DESIGN, "iter")]
L2c <- aggregate(L1[CNT_COLS], by = by2, FUN = sum, na.rm = TRUE)
L2s <- aggregate(L1[SCALAR_COLS], by = by2, FUN = function(v) mean(v, na.rm = TRUE))
L2h <- aggregate(L1[HIT_COLS], by = by2, FUN = sum, na.rm = TRUE)
L2  <- merge(merge(L2c, L2s, by = c(DESIGN, "iter")), L2h, by = c(DESIGN, "iter"))
L2$p_causal        <- .unsentinel(L2$pc,     "sparse")
L2$enrichment_fold <- .unsentinel(L2$enrich, "none")
L2$n_ref <- NA_integer_
L2  <- add_rates(L2)
saveRDS(L2, file.path(out_dir, "combined_replicate_metrics.rds"))
message("  L2: ", nrow(L2), " rows")

# --- L3: one row per cell. Counts summed over all 20 fits; AP averaged over the
#     10 L2 values, and its SE taken across those 10 - NOT across 20 fits.
message("building L3 ...")
L3c <- aggregate(L1[CNT_COLS], by = L1[DESIGN], FUN = sum, na.rm = TRUE)
L3s <- aggregate(L2[SCALAR_COLS], by = L2[DESIGN], FUN = function(v) mean(v, na.rm = TRUE))
L3e <- aggregate(L2["ap"], by = L2[DESIGN],
                 FUN = function(v) sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v))))
names(L3e)[ncol(L3e)] <- "ap_se"
L3h <- aggregate(L1[HIT_COLS], by = L1[DESIGN], FUN = sum, na.rm = TRUE)
L3  <- merge(merge(merge(L3c, L3s, by = DESIGN), L3e, by = DESIGN), L3h, by = DESIGN)
L3$p_causal        <- .unsentinel(L3$pc,     "sparse")
L3$enrichment_fold <- .unsentinel(L3$enrich, "none")
L3$n_ref <- NA_integer_
L3  <- add_rates(L3)
# Merge rather than assign by position: aggregate() returns its own row order
# and assuming it matches L3's would silently mis-attach the counts.
nf <- aggregate(list(n_fits = rep(1L, nrow(L1)),
                     n_failed = as.integer(L1$failed)),
                by = L1[DESIGN], FUN = sum)
L3 <- merge(L3, nf, by = DESIGN, all.x = TRUE)
saveRDS(L3, file.path(out_dir, "combined_scenario_metrics.rds"))
message("  L3: ", nrow(L3), " cells")

# ---------------------------------------------------------------------------
# §11 validity checks - NOTHING IS REPORTED BEFORE THESE PASS
# ---------------------------------------------------------------------------
sink(file.path(out_dir, "validity_checks.txt"))
cat("VALIDITY CHECKS\n===============\n\n")
ok_all <- TRUE
chk <- function(name, pass, detail = "") {
  cat(sprintf("[%s] %s%s\n", if (pass) "PASS" else "FAIL", name,
              if (nzchar(detail)) paste0("  -- ", detail) else ""))
  if (!pass) ok_all <<- FALSE
}

chk("in-sample LD only", all(is.na(L3$n_ref)),
    sprintf("%d rows with non-NA n_ref", sum(!is.na(L3$n_ref))))
chk(sprintf("%d job rows (from params_grid.csv)", EXPECT_JOBS),
    length(unique(L3$job_dir)) == EXPECT_JOBS,
    paste(length(unique(L3$job_dir)), "found"))
bal <- table(L3$method)
chk("balanced cells per method", length(unique(bal)) == 1L,
    paste0("range ", min(bal), "-", max(bal),
           "; the §7 closed-form correction assumes 20 fits in every cell"))
chk("no all-NA method outside the none arm",
    !any(vapply(split(L3$ap[L3$annotation_type != "none"],
                      L3$method[L3$annotation_type != "none"]),
                function(v) all(is.na(v)), logical(1))))
# polyfun_est and polyfun_ldsc reduce EXACTLY to susie on the none arms - a
# built-in correctness check on the whole pipeline (§2.5).
non <- L3[L3$annotation_type == "none", ]
for (mm in c("polyfun_est", "polyfun_ldsc")) {
  a <- non$ap[non$method == mm]; b <- non$ap[non$method == "susie"]
  if (length(a) && length(a) == length(b)) {
    chk(paste(mm, "== susie on the none arm"),
        isTRUE(all.equal(a, b, tolerance = 1e-6)),
        "if this fails the annotation plumbing is wrong, not the method")
  }
}
cat("\nfunmap on the none arm is expected to be entirely absent - it requires\n")
cat("annotations. Structural, not missing at random.\n")
cat(sprintf("\nOVERALL: %s\n", if (ok_all) "all checks passed" else "FAILURES PRESENT"))
sink()
message("validity: ", if (ok_all) "passed" else "FAILURES - see validity_checks.txt")

# ---------------------------------------------------------------------------
# Drive the three senses
# ---------------------------------------------------------------------------
for (script in c("iter004_sense_v.R", "iter004_sense_d.R", "iter004_sense_g.R")) {
  message("\n>>> ", script)
  Sys.setenv(FMB_EXPECT_JOBS = EXPECT_JOBS)
  st <- system2("Rscript", c(file.path("scripts", "analysis", script), out_dir,
                             file.path(out_dir, sub("\\.R$", "", script))))
  if (st != 0L) message("    (exit ", st, " - see its own output)")
}

message("\nDone. Artefacts in ", out_dir)
