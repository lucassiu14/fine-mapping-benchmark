#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_sense_d.R   --   SENSE D: decision-relevant
#
# "Which variables change WHICH METHOD WINS?"
#
# Implements variable-importance-analysis.md §8. A variable can explain 25% of
# the variance in AUPRC while never changing which method you would deploy
# (V high, D zero) - which is why this is a separate question with its own
# response, not a re-reading of Sense V.
#
# Two things make or break this section:
#
#   1. GATING. An ungated flip rate counts coin tosses as findings. A change of
#      winner is real only if the margin exceeds the noise.
#   2. CROSS-FITTING. W(c) is an arg-max over noisy means, so the winner's margin
#      over the runner-up is biased upward when measured on the replicates that
#      chose it - the winner's curse. An in-sample gate is anti-conservative in
#      exactly the marginal cells it exists to police.
#
# Usage:
#   Rscript scripts/analysis/iter004_sense_d.R <collect_dir> [out_dir] [response]
#
# STATUS: written, never executed.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
collect_dir <- args[1] %||% "results/iter004"
out_dir     <- if (length(args) >= 2) args[2] else "results/iter004/sense_d"
RESPONSE    <- if (length(args) >= 3) args[3] else "ap"
source(file.path("scripts", "analysis", "iter004_lib.R"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Z_GATE <- 2      # |dbar| > z * SE_paired
set.seed(20260806)

# polyfun_oracle receives the simulator's TRUE prior weights. It is a ceiling,
# not a deployable competitor, and including it would make it win most cells and
# say nothing (§8.1).
EXCLUDE_FROM_DECISION <- c("polyfun_oracle")

rep_file <- file.path(collect_dir, "combined_replicate_metrics.rds")
if (!file.exists(rep_file)) {
  stop("missing ", rep_file, "\n",
       "Sense D needs the REPLICATE-level table. Because all methods run on the\n",
       "same simulated data, the `iter` index is aligned across methods, so\n",
       "d_i = Y_m1,i - Y_m2,i is a valid PAIRED difference and its SE is 3-10x\n",
       "smaller than the unpaired form. Storing only cell means destroys this.",
       call. = FALSE)
}
reps <- prepare_analysis_table(readRDS(rep_file))
reps <- reps[!reps$method %in% EXCLUDE_FROM_DECISION, , drop = FALSE]

cell_key <- function(d) interaction(d$job_dir, d$S, d$phi, d$region_size, drop = TRUE)

# decide_cell() lives in iter004_lib.R so it can be unit-tested in isolation.

all_flip <- list(); all_map <- list()

for (st in STRATA) {
  s <- stratum_subset(reps, st)
  if (is.null(s)) next
  d <- s$data
  message(sprintf("\n=== stratum %s | free: %s ===", st$key,
                  paste(s$free, collapse = ", ")))

  d$.cell <- cell_key(d)
  cells <- split(seq_len(nrow(d)), d$.cell)

  res <- lapply(names(cells), function(k) {
    sub <- d[cells[[k]], , drop = FALSE]
    dc  <- decide_cell(sub, RESPONSE)
    if (is.null(dc)) return(NULL)
    one <- sub[1, , drop = FALSE]
    data.frame(cell = k, winner = dc$winner, runner_up = dc$runner_up,
               decided = dc$decided, dbar = dc$dbar,
               pc = as.character(one$pc), enrich = as.character(one$enrich),
               S = as.character(one$S), phi = as.character(one$phi),
               region_size = as.character(one$region_size),
               stringsAsFactors = FALSE)
  })
  W <- do.call(rbind, Filter(Negate(is.null), res))
  if (is.null(W)) next

  # --- §8.2 flip rate, raw and gated ---------------------------------------
  # Hold everything else fixed, sweep F, ask whether the winner changed, then
  # average over all settings of the others. Report BOTH rates: the gap
  # quantifies how much apparent decision structure is noise.
  for (f in s$free) {
    others <- setdiff(s$free, f)
    grp <- if (length(others)) interaction(W[others], drop = TRUE)
           else factor(rep("all", nrow(W)))
    sp <- split(W, grp)
    raw <- vapply(sp, function(g) as.integer(length(unique(g$winner)) > 1L), integer(1))
    gat <- vapply(sp, function(g) {
      gg <- g[g$decided, , drop = FALSE]
      if (nrow(gg) < 2L) return(NA_integer_)     # not enough decided levels
      as.integer(length(unique(gg$winner)) > 1L)
    }, integer(1))
    all_flip[[length(all_flip) + 1L]] <- data.frame(
      stratum = st$key, factor = f, response = RESPONSE,
      n_settings = length(sp),
      flip_raw = mean(raw),
      flip_gated = mean(gat, na.rm = TRUE),
      undecidable_settings = mean(is.na(gat)),
      pct_cells_decided = mean(W$decided),
      stringsAsFactors = FALSE)
  }

  W$stratum <- st$key
  all_map[[length(all_map) + 1L]] <- W
}

flip <- do.call(rbind, all_flip)
maps <- do.call(rbind, all_map)
saveRDS(list(flip = flip, winners = maps, response = RESPONSE),
        file.path(out_dir, paste0("sense_d_", RESPONSE, ".rds")))

sink(file.path(out_dir, paste0("sense_d_", RESPONSE, ".txt")))
cat("SENSE D - decision-relevant   (response: ", RESPONSE, ")\n", sep = "")
cat("=========================================================\n\n")
cat("flip_raw counts any change of winner. flip_gated counts only changes where\n")
cat("the margin survives a cross-fitted paired test at z = ", Z_GATE, ".\n", sep = "")
cat("THE GAP BETWEEN THEM IS THE POINT: it is how much apparent decision\n")
cat("structure is noise. A large raw rate with a near-zero gated rate means the\n")
cat("winner is effectively arbitrary in that direction.\n\n")
cat("polyfun_oracle is excluded - it reads the truth and is a ceiling, not a\n")
cat("deployable competitor.\n\n")
if (!is.null(flip)) {
  f <- flip[order(flip$stratum, -flip$flip_gated), ]
  cat(sprintf("%-18s %-14s %10s %10s %12s %10s\n", "stratum", "factor",
              "flip_raw", "flip_gated", "undecidable", "decided"))
  cat(strrep("-", 80), "\n")
  for (i in seq_len(nrow(f))) with(f[i, ],
    cat(sprintf("%-18s %-14s %9.1f%% %9.1f%% %11.1f%% %9.1f%%\n",
                stratum, factor, 100*flip_raw, 100*flip_gated,
                100*undecidable_settings, 100*pct_cells_decided)))
}

# --- §8.4 regime map -------------------------------------------------------
cat("\n\nREGIME MAP - modal winner over the two highest-gated-flip factors\n")
cat("A regime map lives WITHIN one stratum; its axes are free factors of that\n")
cat("stratum, never model, annotation_type or n_ref.\n")
if (!is.null(flip) && !is.null(maps)) {
  for (k in unique(flip$stratum)) {
    fk <- flip[flip$stratum == k, ]
    top <- head(fk$factor[order(-fk$flip_gated)], 2)
    if (length(top) < 2L) next
    mk <- maps[maps$stratum == k & maps$decided, , drop = FALSE]
    if (!nrow(mk)) { cat("\n", k, ": no decided cells\n", sep = ""); next }
    cat("\n", k, "  axes: ", paste(top, collapse = " x "), "\n", sep = "")
    tb <- table(mk[[top[1]]], mk[[top[2]]], mk$winner)
    modal <- apply(tb, c(1, 2), function(v)
      if (sum(v) == 0) "-" else dimnames(tb)[[3]][which.max(v)])
    print(modal, quote = FALSE)
  }
}
sink()
message("\nwrote ", file.path(out_dir, paste0("sense_d_", RESPONSE, ".txt")))
