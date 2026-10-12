#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/apply_iter008_comparators.R
#
# Replace polyfun_ldsc, polyfun_oracle and sbayesrc in Iteration 004's tables
# with Iteration 008's, so that one set of files holds the correct numbers for
# every method and the superseded ones are not kept anywhere.
#
# This is safe to do ONLY because the two runs share their simulations exactly:
# a canary compared 5,000 susie fits across rows 1-2 and found identical realised
# region sizes and identical per-fit AP (scripts/hpc/check_sim_reproduction.R).
# The rows being swapped therefore describe the same regions, the same
# phenotypes and the same causal sets as the rows they replace.
#
# IT OVERWRITES IN PLACE AND THE OLD ROWS ARE NOT KEPT. Iteration 004's
# published values for these three methods survive only as the summary tables in
# docs/autoresearch/comparator-fidelity.md. Run with --dry-run first.
#
#   Rscript scripts/analysis/apply_iter008_comparators.R [--dry-run]
#     [--old results/iter004] [--new results/iter008_fidelity]
#
# What it updates
#   combined_fit_metrics.rds                 L1, one row per fit
#   combined_replicate_metrics.rds           L2
#   combined_scenario_metrics.rds            L3, one row per cell
#   combined_scenario_metrics_with_power.rds L3 + credible-set power, rebuilt
#                                            from the merged L1 because power is
#                                            not recoverable from L3
#
# What it CANNOT update, because Iteration 008 has no equivalent extract: the
# per-job stores aux/, credsets/ and piptail/, and everything derived from them
# (calibration_bands9.rds, fdr_thresholds.rds, calib_fdr_periter*.rds,
# native_credible_sets.rds). Those still hold the superseded numbers for these
# three methods. It lists them at the end so the state is never in doubt.
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
DRY     <- "--dry-run" %in% args
argval  <- function(flag, default) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L]
}
OLD <- argval("--old", "results/iter004")
NEW <- argval("--new", "results/iter008_fidelity")
TARGET <- c("polyfun_ldsc", "polyfun_oracle", "sbayesrc")

stopifnot(dir.exists(OLD), dir.exists(NEW))
cat(sprintf("replacing %s\n  in  %s\n  from %s\n%s\n\n",
            paste(TARGET, collapse = ", "), OLD, NEW,
            if (DRY) "  (DRY RUN - nothing is written)" else ""))

# --- the three tables that swap row for row ---------------------------------
TABLES <- c("combined_fit_metrics", "combined_replicate_metrics",
            "combined_scenario_metrics")
# A row's identity, used to prove the swap neither duplicated nor dropped work.
KEYS <- list(
  combined_fit_metrics       = c("job_dir","scenario_id","region_id","method"),
  combined_replicate_metrics = c("job_dir","S","phi","region_size","iter","method"),
  combined_scenario_metrics  = c("job_dir","S","phi","region_size","method"))

merged_L1 <- NULL
for (nm in TABLES) {
  fo <- file.path(OLD, paste0(nm, ".rds")); fn <- file.path(NEW, paste0(nm, ".rds"))
  o <- readRDS(fo); n <- readRDS(fn)

  if (!identical(names(o), names(n)))
    stop(nm, ": the two tables do not share a schema, so the rows are not ",
         "interchangeable. Columns only in the old table: ",
         paste(setdiff(names(o), names(n)), collapse = ", "),
         "; only in the new: ", paste(setdiff(names(n), names(o)), collapse = ", "),
         call. = FALSE)

  n <- n[n$method %in% TARGET, , drop = FALSE]
  drop_n <- sum(o$method %in% TARGET)
  if (drop_n != nrow(n))
    stop(nm, ": would drop ", drop_n, " rows and add ", nrow(n),
         ". The runs do not cover the same grid; refusing to write a table ",
         "that is short or doubled.", call. = FALSE)

  kept <- o[!o$method %in% TARGET, , drop = FALSE]
  out  <- rbind(kept, n)
  # Deterministic order, so the file does not depend on how it was assembled.
  key  <- KEYS[[nm]]
  key  <- key[key %in% names(out)]
  out  <- out[do.call(order, unname(as.list(out[key]))), , drop = FALSE]
  rownames(out) <- NULL

  id <- do.call(paste, c(out[key], sep = "\r"))
  if (anyDuplicated(id))
    stop(nm, ": the merge produced ", sum(duplicated(id)), " duplicate rows on ",
         paste(key, collapse = " x "), call. = FALSE)
  if (nrow(out) != nrow(o))
    stop(nm, ": row count changed, ", nrow(o), " -> ", nrow(out), call. = FALSE)

  cat(sprintf("%-28s %s rows  (%s swapped, %d methods, keys unique)\n",
              nm, format(nrow(out), big.mark = ","),
              format(nrow(n), big.mark = ","), length(unique(out$method))))

  if (nm == "combined_fit_metrics") merged_L1 <- out
  if (!DRY) {
    tmp <- paste0(fo, ".tmp"); saveRDS(out, tmp); file.rename(tmp, fo)
  }
  rm(o, n, kept, out); invisible(gc())
}

# --- rebuild the power column from the merged L1 ----------------------------
# set_prec_95 and set_size_95 are per-fit means at L3, and mean(prec)*mean(k) is
# not mean(prec*k), so power has to be formed at fit level and averaged in.
fo_p <- file.path(OLD, "combined_scenario_metrics_with_power.rds")
if (file.exists(fo_p)) {
  KEY <- c("job_dir","model","annotation_type","S","phi","region_size","method")
  L1  <- if (is.null(merged_L1)) readRDS(file.path(OLD, "combined_fit_metrics.rds")) else merged_L1
  L1  <- L1[!L1$failed & !L1$bad_pip, ]
  L1$power <- ifelse(L1$n_causal > 0,
                     (L1$set_prec_95 * L1$set_size_95) / L1$n_causal, NA_real_)
  pw <- tapply(L1$power, do.call(paste, c(L1[KEY], sep = "\r")), mean, na.rm = TRUE)
  rm(L1); invisible(gc())
  L3 <- readRDS(file.path(OLD, "combined_scenario_metrics.rds"))
  L3$power <- as.numeric(pw[do.call(paste, c(L3[KEY], sep = "\r"))])
  cat(sprintf("%-28s %s rows  (power attached to %.1f%% of cells)\n",
              "…_with_power", format(nrow(L3), big.mark = ","),
              100 * mean(!is.na(L3$power))))
  if (!DRY) { tmp <- paste0(fo_p, ".tmp"); saveRDS(L3, tmp); file.rename(tmp, fo_p) }
  rm(L3); invisible(gc())
}

# --- say plainly what is still superseded -----------------------------------
stale <- c("calibration_bands9.rds", "fdr_thresholds.rds", "calib_fdr_periter.rds",
           "calib_fdr_periter_units.rds", "native_credible_sets.rds",
           "aux/ (45 files)", "credsets/ (45 files)", "piptail/ (45 files)",
           "figures_by_factor/ap_by_factor_data.rds", "iter004_sense_v/sense_v.rds")
present <- stale[file.exists(file.path(OLD, sub(" .*", "", stale)))]
if (length(present)) {
  cat("\nSTILL HOLDING THE SUPERSEDED NUMBERS for these three methods:\n")
  for (p in present) cat("  ", p, "\n")
  cat("Each is derived from the raw per-scenario PIPs, which Iteration 008 wrote\n",
      "as results_supp.rds on scratch. Regenerating them means re-running the\n",
      "extract scripts against that root; they cannot be rebuilt from the tables.\n", sep = "")
}
cat(if (DRY) "\nDRY RUN - nothing was written.\n" else "\nDone.\n")
