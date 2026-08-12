#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_sense_v.R   --   SENSE V: variance-explaining
#
# "Which simulation variables move this method's metric?"
#
# Implements variable-importance-analysis.md §7. Runs the decomposition per
# stratum x method x response, and reports the variance budget FIRST because it
# bounds every subsequent claim: if the noise floor is 50%, no factor
# decomposition can explain more than half.
#
# Usage:
#   Rscript scripts/analysis/iter004_sense_v.R <collect_dir> [out_dir]
#
# Expects in <collect_dir> (produced by the re-collect, spec §5.5):
#   combined_fit_metrics.rds       L1 - one row per (job,S,phi,iter,region_idx,method)
#   combined_scenario_metrics.rds  L3 - one row per cell
#
# STATUS: written, never executed.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
collect_dir <- args[1] %||% "results/iter004"
out_dir     <- if (length(args) >= 2) args[2] else "results/iter004/sense_v"
source(file.path("scripts", "analysis", "iter004_lib.R"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The three responses. Different factors will top each list, and THAT
# DISAGREEMENT IS THE FINDING - it is the formal version of "this method's
# ranking is fine but its calibration is broken" (§7.6).
RESPONSES <- c("ap", "fdr_at_90", "bss")

fit_file  <- file.path(collect_dir, "combined_fit_metrics.rds")
cell_file <- file.path(collect_dir, "combined_scenario_metrics.rds")
for (f in c(fit_file, cell_file)) {
  if (!file.exists(f)) {
    stop("missing ", f, "\n",
         "Run the re-collect first. Sense V needs BOTH levels: the fit table for\n",
         "the variance components (§3.4b) and the cell table for the SS partition\n",
         "(§7.1). Without region_idx in the fit table, sigma^2_u is unidentifiable\n",
         "and p_causal / enrichment_fold have NO VALID ERROR TERM AT ALL.",
         call. = FALSE)
  }
}

fits  <- prepare_analysis_table(readRDS(fit_file))
cells <- prepare_analysis_table(readRDS(cell_file))

if (!"region_idx" %in% names(fits)) {
  stop("combined_fit_metrics.rds has no `region_idx` column. The two regions in ",
       "a size class cannot be distinguished, so sigma^2_u is unidentifiable ",
       "(§3.4b). Re-collect with region_idx before proceeding.", call. = FALSE)
}

methods <- sort(unique(cells$method))
message(sprintf("methods: %d | strata: %d | responses: %d",
                length(methods), length(STRATA), length(RESPONSES)))

all_budgets <- list(); all_terms <- list(); all_sobol <- list(); all_vc <- list()

for (st in STRATA) {
  sc <- stratum_subset(cells, st); sf <- stratum_subset(fits, st)
  if (is.null(sc) || is.null(sf)) {
    message("  stratum ", st$key, ": absent, skipping"); next
  }
  message(sprintf("\n=== stratum %s | free factors: %s ===",
                  st$key, paste(sc$free, collapse = ", ")))

  for (m in methods) {
    cm <- sc$data[sc$data$method == m, , drop = FALSE]
    fm <- sf$data[sf$data$method == m, , drop = FALSE]
    if (!nrow(cm) || !nrow(fm)) next

    # funmap on the `none` arm is structurally absent, not missing at random.
    # Do not let an NA-dropping aggregation quietly change the method set.
    if (all(is.na(cm$ap))) {
      message(sprintf("  %-22s all-NA in this stratum (structural?) - skipped", m))
      next
    }

    for (resp in RESPONSES) {
      if (!resp %in% names(cm) || all(is.na(cm[[resp]]))) next

      vc <- tryCatch(variance_components(fm, resp), error = function(e) NULL)
      if (is.null(vc) || !is.finite(vc$sigma2_eps)) {
        message(sprintf("  %-22s %-12s variance components failed - skipped", m, resp))
        next
      }

      dec <- tryCatch(decompose(cm, resp, sc$free, vc), error = function(e) {
        message(sprintf("  %-22s %-12s decompose failed: %s", m, resp,
                        conditionMessage(e))); NULL
      })
      if (is.null(dec)) next

      all_vc[[length(all_vc) + 1L]] <- data.frame(
        stratum = st$key, method = m, response = resp,
        sigma2_eps = vc$sigma2_eps, sigma2_u_bar = vc$sigma2_u_bar,
        sigma2_job = vc$sigma2_job, stringsAsFactors = FALSE)

      all_budgets[[length(all_budgets) + 1L]] <- data.frame(
        stratum = st$key, method = m, response = resp,
        ss_tot = dec$ss_tot,
        noise_floor = dec$noise_floor,
        noise_iter = dec$noise_floor_iter,
        noise_region = dec$noise_floor_region,
        structure_mains = sum(dec$table$omega2[dec$table$order == 1L]),
        structure_2way  = sum(dec$table$omega2[dec$table$order == 2L]),
        structure_3way_plus = sum(dec$table$omega2[dec$table$order >= 3L]),
        stringsAsFactors = FALSE)

      tt <- dec$table
      tt$stratum <- st$key; tt$method <- m; tt$response <- resp
      all_terms[[length(all_terms) + 1L]] <- tt

      sb <- sobol_indices(dec, sc$free)
      sb$stratum <- st$key; sb$method <- m; sb$response <- resp
      all_sobol[[length(all_sobol) + 1L]] <- sb
    }
  }
}

budgets <- do.call(rbind, all_budgets)
terms_t <- do.call(rbind, all_terms)
sobol_t <- do.call(rbind, all_sobol)
vc_t    <- do.call(rbind, all_vc)

saveRDS(list(budgets = budgets, terms = terms_t, sobol = sobol_t,
             variance_components = vc_t, caveats = importance_caveats()),
        file.path(out_dir, "sense_v.rds"))

# --- report ----------------------------------------------------------------
sink(file.path(out_dir, "sense_v.txt"))
cat("SENSE V - variance-explaining\n")
cat("=============================\n\n")
cat("THE VARIANCE BUDGET BOUNDS EVERY CLAIM BELOW IT. A noise floor of 50% means\n")
cat("no factor decomposition can explain more than half. If main effects dominate,\n")
cat("marginal plots are honest summaries; if 2-way+ dominates, every claim must be\n")
cat("stated conditionally.\n\n")

if (!is.null(budgets)) {
  b <- budgets[order(budgets$response, budgets$stratum, -budgets$noise_floor), ]
  cat(sprintf("%-18s %-22s %-11s %8s %8s %8s %8s %8s %9s\n",
              "stratum", "method", "response", "noise", "(iter)", "(region)",
              "mains", "2-way", "3-way+"))
  cat(strrep("-", 112), "\n")
  for (i in seq_len(nrow(b))) with(b[i, ],
    cat(sprintf("%-18s %-22s %-11s %7.1f%% %7.1f%% %7.1f%% %7.1f%% %7.1f%% %8.1f%%\n",
                stratum, method, response, 100*noise_floor, 100*noise_iter,
                100*noise_region, 100*structure_mains, 100*structure_2way,
                100*structure_3way_plus)))
}

cat("\n\nTOP TERMS PER (stratum, method, response)\n")
cat("omega^2 is a share of the OBSERVED cell-mean SS; shares plus the noise floor\n")
cat("sum to 1. Negative values are reported as computed and NOT clamped.\n\n")
if (!is.null(terms_t)) {
  sp <- split(terms_t, interaction(terms_t$stratum, terms_t$method,
                                   terms_t$response, drop = TRUE))
  for (k in names(sp)) {
    tt <- head(sp[[k]][order(-sp[[k]]$omega2), ], 6)
    cat("\n", k, "\n", sep = "")
    cat(sprintf("  %-28s %4s %6s %10s %10s\n", "term", "df", "WP", "omega^2", "span"))
    for (i in seq_len(nrow(tt))) with(tt[i, ],
      cat(sprintf("  %-28s %4d %6s %9.3f%% %10s\n", term, df,
                  ifelse(whole_plot, "yes", "-"), 100*omega2,
                  ifelse(is.na(span), "-", sprintf("%.4f", span)))))
  }
}

cat("\n\nINTERACTION INVOLVEMENT (S_Ti - S_i)\n")
cat("~0 means the factor acts additively and its marginal plot tells the whole\n")
cat("story. Large means the effect is conditional and the marginal is misleading.\n\n")
if (!is.null(sobol_t)) print(sobol_t, row.names = FALSE, digits = 3)

cat("\n\nMANDATORY CAVEATS\n")
for (cv in importance_caveats()) cat("\n* ", cv, "\n", sep = "")
sink()

message("\nwrote ", file.path(out_dir, "sense_v.txt"))
