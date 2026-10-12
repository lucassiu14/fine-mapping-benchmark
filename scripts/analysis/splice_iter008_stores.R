#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/splice_iter008_stores.R
#
# Second half of folding Iteration 008 into Iteration 004. The first half,
# apply_iter008_comparators.R, swapped rows in the three combined_* data frames.
# This one does the per-job stores, which are not data frames: piptail/,
# credsets/ and aux/ each hold a list of per-fit RECORDS, one row file per
# params_grid row, and each record carries its own `method` field. Splicing
# means dropping the records of the three superseded methods and appending
# Iteration 008's, per file.
#
# Legitimate only because the two runs share their simulations exactly, which
# the canary established on 5,000 susie fits
# (scripts/hpc/check_sim_reproduction.R). The records being swapped describe the
# same regions as the records they replace.
#
# IT OVERWRITES IN PLACE AND THE OLD RECORDS ARE NOT KEPT. Run --dry-run first.
#
#   Rscript scripts/analysis/splice_iter008_stores.R [--dry-run]
#     [--old results/iter004] [--new results/iter008_fidelity]
#     [--only piptail,aux,credsets]
#
# Afterwards these rebuild from the spliced stores, and must be rerun:
#   native_credible_sets.rds   <- credsets/
#   calibration_bands9.rds     <- piptail/ + L1   (iter004_calibration_bands.R)
#   fdr_thresholds.rds         <- piptail/ + L1   (iter004_fdr_thresholds.R)
#   calib_fdr_periter*.rds     <- piptail/
# =============================================================================

args   <- commandArgs(trailingOnly = TRUE)
DRY    <- "--dry-run" %in% args
argval <- function(flag, default) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L]
}
OLD    <- argval("--old", "results/iter004")
NEW    <- argval("--new", "results/iter008_fidelity")
ONLY   <- strsplit(argval("--only", "piptail,credsets,aux"), ",")[[1]]
TARGET <- c("polyfun_ldsc", "polyfun_oracle", "sbayesrc")

# Each store: the subdirectory, the filename stem, and which list components
# hold per-fit records keyed by method.
STORES <- list(
  piptail  = list(dir = "piptail",  stem = "piptail_", parts = "tail",
                  count = "fits"),
  credsets = list(dir = "credsets", stem = "cs_",      parts = "rows",
                  count = "fits"),
  aux      = list(dir = "aux",      stem = "aux_",     parts = c("runtime", "importance"),
                  count = NULL))

methods_of <- function(recs) {
  if (!length(recs)) return(character(0))
  vapply(recs, function(r) if (is.null(r$method)) NA_character_ else
    as.character(r$method)[1], character(1))
}

total <- list(files = 0L, dropped = 0L, added = 0L)
for (s in intersect(ONLY, names(STORES))) {
  cfg <- STORES[[s]]
  od <- file.path(OLD, cfg$dir); nd <- file.path(NEW, cfg$dir)
  if (!dir.exists(od)) { cat(sprintf("%-9s no %s - skipped\n", s, od)); next }
  if (!dir.exists(nd)) {
    cat(sprintf("%-9s no %s - NOT SPLICED, still superseded\n", s, nd)); next
  }
  of <- sort(Sys.glob(file.path(od, paste0(cfg$stem, "*.rds"))))
  nf <- sort(Sys.glob(file.path(nd, paste0(cfg$stem, "*.rds"))))
  if (length(of) != length(nf)) {
    cat(sprintf("%-9s %d old files but %d new - NOT SPLICED\n", s, length(of), length(nf)))
    next
  }
  if (!identical(basename(of), basename(nf))) {
    cat(sprintf("%-9s the row files do not correspond - NOT SPLICED\n", s)); next
  }

  dropped <- 0L; added <- 0L
  for (i in seq_along(of)) {
    o <- readRDS(of[i]); n <- readRDS(nf[i])
    for (p in cfg$parts) {
      if (is.null(o[[p]])) stop(basename(of[i]), ": no component '", p, "'", call. = FALSE)
      mo <- methods_of(o[[p]]); mn <- methods_of(n[[p]])
      keep <- o[[p]][!mo %in% TARGET]
      take <- n[[p]][mn %in% TARGET]
      # A component can legitimately hold nothing for these methods in EITHER
      # run - aux's `importance` is the LassoNet annotation weights, which only
      # the BEATRICE family produces - and then there is nothing to splice.
      # What must never happen is dropping records and putting none back.
      if (sum(mo %in% TARGET) > 0L && !length(take))
        stop(basename(nf[i]), ": component '", p, "' has no records for ",
             paste(TARGET, collapse = ", "), ", but the old file has ",
             sum(mo %in% TARGET), ". Refusing to drop them and put nothing back.",
             call. = FALSE)
      dropped <- dropped + sum(mo %in% TARGET); added <- added + length(take)
      o[[p]] <- c(keep, take)
    }
    if (!is.null(cfg$count)) o[[cfg$count]] <- length(o[[cfg$parts[1]]])
    if (!DRY) { tmp <- paste0(of[i], ".tmp"); saveRDS(o, tmp); file.rename(tmp, of[i]) }
    rm(o, n); if (i %% 10L == 0L) invisible(gc())
  }
  cat(sprintf("%-9s %2d files  %s records dropped, %s added\n", s, length(of),
              format(dropped, big.mark = ","), format(added, big.mark = ",")))
  total$files <- total$files + length(of)
  total$dropped <- total$dropped + dropped; total$added <- total$added + added
}

cat(sprintf("\n%s: %d files, %s records dropped, %s added\n",
            if (DRY) "DRY RUN" else "Done", total$files,
            format(total$dropped, big.mark = ","), format(total$added, big.mark = ",")))
if (!DRY && total$files)
  cat("Now rebuild native_credible_sets.rds, calibration_bands9.rds,\n",
      "fdr_thresholds.rds and calib_fdr_periter*.rds from the spliced stores.\n", sep = "")
