#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_sense_g.R   --   SENSE G: diagnostic
#
# "Which variables hurt THIS METHOD, beyond making the task hard for everyone?"
#
# Implements variable-importance-analysis.md §9. This is the section that
# prevents the most common error in the whole analysis:
#
#   A factor affects a method's score in two ways that the raw score confounds -
#   it makes the problem harder for EVERYTHING (a property of the task), and it
#   makes THIS METHOD specifically worse (a property of the method). The factor
#   with the largest omega^2 on raw AP is usually just the difficulty knob.
#   Reporting it as a finding ABOUT THE METHOD is the error this exists to stop.
#
# The payload is the difficulty-vs-sensitivity table (§9.3): the same
# decomposition run on the raw metric and on the headroom-relative kappa, placed
# side by side.
#
#   omega^2 on Y high, on kappa ~0   ->  DIFFICULTY KNOB. Says nothing about the method.
#   omega^2 on Y ~0,  on kappa high  ->  METHOD-SPECIFIC WEAKNESS. The thing to fix.
#   both high                        ->  genuinely hard AND handled badly.
#   both ~0                          ->  irrelevant over the range swept.
#
# Usage:
#   Rscript scripts/analysis/iter004_sense_g.R <collect_dir> [out_dir] [response]
#
# STATUS: written, never executed.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
collect_dir <- args[1] %||% "results/iter004"
out_dir     <- if (length(args) >= 2) args[2] else "results/iter004/sense_g"
RESPONSE    <- if (length(args) >= 3) args[3] else "ap"
source(file.path("scripts", "analysis", "iter004_lib.R"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

CEILING   <- "polyfun_oracle"       # receives the simulator's true prior weights
BASELINE  <- "susie"                # annotation-blind backbone
HEADROOM_EPS <- 0.02                # kappa is unstable when the denominator ~ 0

FOCUS <- c("functional_beatrice", "fb_xregion", "fb_pooled")
# Controlled pairs: annotation use is the ONLY difference within each pair.
PAIRS <- list(
  annot = c(num = "functional_beatrice", den = "beatrice"),
  joint = c(num = "fb_xregion",          den = "functional_beatrice")
)

rep_file <- file.path(collect_dir, "combined_replicate_metrics.rds")
fit_file <- file.path(collect_dir, "combined_fit_metrics.rds")
for (f in c(rep_file, fit_file)) {
  if (!file.exists(f)) stop("missing ", f, call. = FALSE)
}
reps <- prepare_analysis_table(readRDS(rep_file))

#' Reshape wide by method within a (cell, replicate).
#'
#' All methods run on the SAME simulated data, so within a replicate these are
#' exact paired quantities - kappa and Delta are computed at replicate level, not
#' from cell means, which is what makes §3.4b applicable to them.
wide_by_method <- function(d, response) {
  d$.cell <- interaction(d$job_dir, d$S, d$phi, d$region_size, drop = TRUE)
  w <- reshape(d[, c(".cell", "iter", "method", response)],
               idvar = c(".cell", "iter"), timevar = "method", direction = "wide")
  names(w) <- sub(paste0("^", response, "\\."), "", names(w))
  meta <- unique(d[, c(".cell", "job_dir", "S", "phi", "region_size",
                       "pc", "enrich", "model", "annotation_type")])
  merge(w, meta, by = ".cell", all.x = TRUE)
}

all_kappa <- list(); all_delta <- list(); all_excl <- list()

for (st in STRATA) {
  s <- stratum_subset(reps, st)
  if (is.null(s)) next
  message(sprintf("\n=== stratum %s ===", st$key))
  w <- wide_by_method(s$data, RESPONSE)

  have_ceiling  <- CEILING  %in% names(w)
  have_baseline <- BASELINE %in% names(w)

  # --- kappa: fraction of available annotation headroom captured -----------
  if (have_ceiling && have_baseline) {
    w$headroom <- w[[CEILING]] - w[[BASELINE]]
    # GATE. kappa is unstable when the denominator ~ 0, and the cells where the
    # ceiling has COLLAPSED are themselves a finding, not a nuisance - report the
    # excluded fraction rather than silently dropping them.
    keep <- is.finite(w$headroom) & w$headroom > HEADROOM_EPS
    all_excl[[length(all_excl) + 1L]] <- data.frame(
      stratum = st$key, response = RESPONSE,
      n_rows = nrow(w), excluded_frac = mean(!keep),
      median_headroom = median(w$headroom, na.rm = TRUE),
      stringsAsFactors = FALSE)
    message(sprintf("  headroom gate: %.1f%% of replicate-rows excluded (headroom <= %.2f)",
                    100 * mean(!keep), HEADROOM_EPS))

    wk <- w[keep, , drop = FALSE]
    for (m in intersect(FOCUS, names(wk))) {
      k <- (wk[[m]] - wk[[BASELINE]]) / wk$headroom
      kd <- wk[, c(".cell", "job_dir", "S", "phi", "region_size", "pc", "enrich",
                   "model", "annotation_type"), drop = FALSE]
      kd$iter <- wk$iter; kd$kappa <- k; kd$method <- m; kd$stratum <- st$key
      all_kappa[[length(all_kappa) + 1L]] <- kd
    }
  } else {
    message("  ceiling or baseline absent - kappa not computable in this stratum")
  }

  # --- Delta: controlled pairs --------------------------------------------
  for (nm in names(PAIRS)) {
    p <- PAIRS[[nm]]
    if (!all(p %in% names(w))) next
    dd <- w[, c(".cell", "job_dir", "S", "phi", "region_size", "pc", "enrich",
                "model", "annotation_type"), drop = FALSE]
    dd$iter <- w$iter
    dd$delta <- w[[p[["num"]]]] - w[[p[["den"]]]]
    dd$pair <- nm; dd$stratum <- st$key
    all_delta[[length(all_delta) + 1L]] <- dd
  }
}

kappa_t <- do.call(rbind, all_kappa)
delta_t <- do.call(rbind, all_delta)
excl_t  <- do.call(rbind, all_excl)

# --- decompose kappa and Delta exactly as §7 ------------------------------
# NOTE: the §7.2 correction needs THEIR variance components, not the raw
# metric's. kappa and Delta are computed at replicate level above, so
# variance_components() applies directly - but only if the table carries
# region_idx. Without it these decompositions cannot be corrected and must be
# reported as UNCORRECTED, which inflates whole-plot terms by orders of
# magnitude. We refuse rather than mislead.
decompose_derived <- function(tbl, value_col, label) {
  if (is.null(tbl)) return(NULL)
  if (!"region_idx" %in% names(tbl)) {
    message(sprintf(
      "  %s: replicate table has no region_idx, so sigma^2_u is unidentifiable. ",
      label),
      "Skipping the corrected decomposition rather than reporting an ",
      "uncorrected one - whole-plot terms would be inflated by ~m*sigma^2_u/2 ",
      "per df.")
    return(NULL)
  }
  out <- list()
  for (st in STRATA) {
    sub <- tbl[tbl$stratum == st$key, , drop = FALSE]
    if (!nrow(sub)) next
    ss <- stratum_subset(sub, st); if (is.null(ss)) next
    for (m in unique(sub$method %||% sub$pair)) {
      sm <- ss$data[(ss$data$method %||% ss$data$pair) == m, , drop = FALSE]
      if (!nrow(sm)) next
      vc <- tryCatch(variance_components(sm, value_col), error = function(e) NULL)
      if (is.null(vc)) next
      cm <- aggregate(sm[[value_col]],
                      by = sm[c("job_dir", "S", "phi", "region_size",
                                "pc", "enrich")],
                      FUN = mean, na.rm = TRUE)
      names(cm)[ncol(cm)] <- value_col
      dec <- tryCatch(decompose(cm, value_col, ss$free, vc), error = function(e) NULL)
      if (is.null(dec)) next
      tt <- dec$table; tt$stratum <- st$key; tt$method <- m; tt$what <- label
      out[[length(out) + 1L]] <- tt
    }
  }
  do.call(rbind, out)
}

kappa_dec <- decompose_derived(kappa_t, "kappa", "kappa")
delta_dec <- decompose_derived(delta_t, "delta", "delta")

saveRDS(list(kappa = kappa_t, delta = delta_t, excluded = excl_t,
             kappa_decomposition = kappa_dec, delta_decomposition = delta_dec,
             response = RESPONSE),
        file.path(out_dir, paste0("sense_g_", RESPONSE, ".rds")))

sink(file.path(out_dir, paste0("sense_g_", RESPONSE, ".txt")))
cat("SENSE G - diagnostic   (response: ", RESPONSE, ")\n", sep = "")
cat("=================================================\n\n")
cat("kappa = (Y_m - Y_susie) / (Y_polyfun_oracle - Y_susie)\n")
cat("      = the fraction of AVAILABLE ANNOTATION HEADROOM the method captures.\n\n")
cat("Raw metrics cannot answer this question. The factor with the largest\n")
cat("omega^2 on raw AP is usually just the difficulty knob; reporting it as a\n")
cat("finding about the method is the error this section exists to prevent.\n\n")

if (!is.null(excl_t)) {
  cat("HEADROOM GATE - cells where the ceiling has collapsed are a FINDING:\n\n")
  print(excl_t, row.names = FALSE, digits = 3)
}

if (!is.null(kappa_t)) {
  cat("\n\nMEAN kappa BY STRATUM AND METHOD (1.0 = captures all available headroom)\n\n")
  agg <- aggregate(kappa ~ stratum + method, kappa_t, function(v)
    c(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE)))
  print(agg, row.names = FALSE, digits = 3)
}

if (!is.null(delta_t)) {
  cat("\n\nCONTROLLED PAIRS (annotation use is the only difference)\n")
  cat("  annot = functional_beatrice - beatrice\n")
  cat("  joint = fb_xregion - functional_beatrice\n\n")
  agg <- aggregate(delta ~ stratum + pair, delta_t, function(v)
    c(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE)))
  print(agg, row.names = FALSE, digits = 3)
}

cat("\n\nDIFFICULTY vs SENSITIVITY\n")
cat("Place these beside the Sense V table for the SAME stratum and method.\n\n")
cat("  omega^2 on Y high, on kappa ~0  -> DIFFICULTY KNOB - hard for everyone,\n")
cat("                                     says nothing about this method\n")
cat("  omega^2 on Y ~0,  on kappa high -> METHOD-SPECIFIC WEAKNESS - the thing to fix\n")
cat("  both high                       -> genuinely hard AND handled badly\n")
cat("  both ~0                         -> irrelevant over the range swept\n\n")
if (!is.null(kappa_dec)) {
  k <- kappa_dec[kappa_dec$order == 1L, ]
  print(k[order(k$stratum, k$method, -k$omega2),
          c("stratum", "method", "term", "omega2", "span")],
        row.names = FALSE, digits = 3)
} else {
  cat("(kappa decomposition unavailable - see the region_idx note above)\n")
}
sink()
message("\nwrote ", file.path(out_dir, paste0("sense_g_", RESPONSE, ".txt")))
