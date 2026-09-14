#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_ap_by_factor_figures.R
#
# Average precision against each simulation factor - S, phi, rho and region
# size - for the LSR's "Main simulation - additional results" section. One
# figure per factor.
#
# POOLED ACROSS STRATA, by the user's explicit one-off decision (2026-09-14).
# Every other LSR figure keeps generative model x annotation type apart; these
# average over them so that each factor gets a single panel. Cells are weighted
# equally, so at each S, phi and region-size level the sparse-plus-infinitesimal
# model supplies 4/5 of the cells and the annotated regimes 8/9 of them.
#
# Otherwise the conventions of fig_auprc (iter004_lsr_figures.R on
# lsr/feedback1-figures): the same L3 table and exclusions (carma, fb_pooled,
# polyfun_est, and cells in which every fit failed), the mean over cells with a
# 95% interval, the seven displayed methods in lsr_palette.R colours and shapes,
# and every other method in grey.
#
# rho is the share of genetic variance carried by the sparse component. It is
# swept (0.2-0.8) only in the sparse-plus-infinitesimal model; the sparse model
# is exactly rho = 1 (generate_params_grid.R leaves 1.0 out of sparse_inf for
# that reason), so it is placed there.
#
#   Rscript scripts/analysis/iter004_ap_by_factor_figures.R [out_dir] [l3_rds]
# =============================================================================

suppressMessages(library(ggplot2))
.a  <- commandArgs(TRUE)
FIG <- if (length(.a) >= 1) .a[1] else "results/iter004/figures_by_factor"
L3F <- if (length(.a) >= 2) .a[2] else "results/iter004/combined_scenario_metrics_with_power.rds"
stopifnot(file.exists(L3F))
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
source("scripts/analysis/lsr_palette.R")

DROP <- c("carma", "fb_pooled", "polyfun_est")
L3 <- readRDS(L3F)
L3 <- L3[!L3$method %in% DROP & L3$n_failed < L3$n_fits, ]
L3 <- L3[is.finite(L3$ap), ]
stopifnot(all(M7 %in% L3$method))
L3$rho <- ifelse(L3$model == "sparse", 1, L3$p_causal)
stopifnot(!anyNA(L3$rho))

# --- the same table as section 6.3, before anything is pooled ---------------
# Results quote these per-stratum means (sparse model, binary annotations).
quoted <- c(polyfun_ldsc = 0.773, susie = 0.723, fb_xregion = 0.748)
got <- vapply(names(quoted), function(m) mean(L3$ap[L3$method == m & L3$model == "sparse" &
                                                     L3$annotation_type == "binary"]), 0)
cat("Section 6.3 check (sparse, binary):\n")
print(data.frame(method = names(quoted), quoted = quoted, recomputed = round(got, 3)), row.names = FALSE)
stopifnot(all(round(got, 3) == quoted))

FACTORS <- list(
  S      = list(col = "S",           lab = "Causal variants per region, S",
                lev = c(1, 2, 3, 5, 10)),
  phi    = list(col = "phi",         lab = "Variance explained, φ",
                lev = c(0.05, 0.1, 0.2, 0.4, 0.6)),
  rho    = list(col = "rho",         lab = "Sparse variance share, ρ (1 = sparse model)",
                lev = c(0.2, 0.4, 0.6, 0.8, 1)),
  region = list(col = "region_size", lab = "Variants per region",
                lev = c(500, 750, 1000, 1500, 2000)))

th <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = .25, colour = "grey92"),
        panel.border = element_rect(colour = "grey70", linewidth = .35),
        axis.title = element_text(size = 7.8),
        axis.text  = element_text(size = 6.9),
        legend.text  = element_text(size = 6.9),
        legend.key.size = unit(9, "pt"),
        legend.position = "right")

OTH <- sprintf("other methods (%d)", length(setdiff(unique(L3$method), M7)))
out <- list()
for (fk in names(FACTORS)) {
  f <- FACTORS[[fk]]
  x <- L3[[f$col]]
  stopifnot(all(x %in% f$lev), all(f$lev %in% x))

  # Balance: cells per level for one method, by stratum. Equal rows across
  # levels mean the pooled weights are the same at every level.
  bal <- with(L3[L3$method == "susie", ], table(level = x[L3$method == "susie"],
                                               stratum = paste(model, annotation_type)))
  cat(sprintf("\n[%s] susie cells per level x stratum:\n", fk)); print(bal)

  d <- do.call(rbind, lapply(split(L3, list(x, L3$method), drop = TRUE), function(r) {
    data.frame(level = r[[f$col]][1], method = r$method[1], m = mean(r$ap),
               se = sd(r$ap) / sqrt(nrow(r)), n = nrow(r))
  }))
  d$pos <- match(d$level, f$lev)
  d$grp <- factor(ifelse(d$method %in% M7, d$method, OTH), c(M7, OTH))
  hi <- d[d$grp != OTH, ]; lo <- d[d$grp == OTH, ]

  p <- ggplot(d, aes(pos, m, group = method, colour = grp, shape = grp)) +
    geom_line(data = lo, linewidth = .35) +
    geom_line(data = hi, linewidth = .55) +
    geom_errorbar(data = hi, aes(ymin = m - 1.96 * se, ymax = m + 1.96 * se),
                  width = .12, linewidth = .4) +
    geom_point(data = hi, size = 1.7, stroke = .6) +
    scale_x_continuous(breaks = seq_along(f$lev), labels = as.character(f$lev),
                       expand = expansion(add = .25)) +
    scale_colour_manual(values = c(PAL7, setNames(GREY, OTH)), limits = c(M7, OTH), name = NULL) +
    scale_shape_manual(values = c(SHP7, setNames(16, OTH)), limits = c(M7, OTH), name = NULL) +
    labs(x = f$lab, y = "Mean average precision over cells (95% CI)") + th
  ggsave(file.path(FIG, sprintf("fig_ap_by_%s.pdf", fk)), p, width = 6.3, height = 3.4, device = cairo_pdf)

  tab <- reshape(hi[, c("level", "method", "m")], idvar = "method", timevar = "level", direction = "wide")
  names(tab) <- sub("^m\\.", "", names(tab))
  cat(sprintf("\n[%s] mean AP of the seven (cells per point: %s):\n", fk,
              paste(range(hi$n), collapse = "-")))
  tab[, -1] <- round(tab[, -1], 3); print(tab[match(M7, tab$method), ], row.names = FALSE)
  out[[fk]] <- d
}
saveRDS(out, file.path(FIG, "ap_by_factor_data.rds"))
cat("\nwrote", length(FACTORS), "figures to", FIG, "\n")
