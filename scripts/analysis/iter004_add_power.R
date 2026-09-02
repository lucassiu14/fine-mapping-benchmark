#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_add_power.R
#
# Attach credible-set POWER to the L3 table.
#
#   Rscript scripts/analysis/iter004_add_power.R [results/iter004]
#
# Power - the proportion of causal variants captured by the 95% set - cannot be
# recovered from L3. There, set_prec_95 and set_size_95 are both per-fit MEANS,
# and mean(prec) * mean(k) != mean(prec * k). It has to be formed at fit level:
#
#     captured = set_prec_95 * set_size_95      # exact, per fit
#     power    = captured / n_causal
#
# then averaged into cells on the same key L3 uses.
# =============================================================================

dir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "results/iter004"
KEY <- c("job_dir","model","annotation_type","S","phi","region_size","method")

L1 <- readRDS(file.path(dir, "combined_fit_metrics.rds"))
L1 <- L1[!L1$failed & !L1$bad_pip, ]
L1$power <- ifelse(L1$n_causal > 0,
                   (L1$set_prec_95 * L1$set_size_95) / L1$n_causal, NA_real_)
pw <- tapply(L1$power, do.call(paste, c(L1[KEY], sep = "\r")), mean, na.rm = TRUE)
rm(L1); invisible(gc())

L3 <- readRDS(file.path(dir, "combined_scenario_metrics.rds"))
L3$power <- as.numeric(pw[do.call(paste, c(L3[KEY], sep = "\r"))])
message(sprintf("power attached to %d of %d cells (%.1f%%)",
                sum(!is.na(L3$power)), nrow(L3), 100*mean(!is.na(L3$power))))

out <- file.path(dir, "combined_scenario_metrics_with_power.rds")
saveRDS(L3, out); message("wrote ", out)
