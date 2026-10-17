#!/usr/bin/env Rscript
# =============================================================================
# lambda_ap.R -- rebuild fig_xregion_lambda_auprc.pdf.
#
# Two changes against fb_lambda_sweep_figures.R:
#   1. all SEVEN lambda arms are shown, not the five chosen for readability;
#   2. the polyfun_ldsc reference row comes from Iteration 008, so it agrees
#      with Figure 2 of the main text. beatrice was not re-run and stays on
#      Iteration 004.
# Everything else -- the estimator, the strata, the layout -- is unchanged.
# =============================================================================
suppressMessages({library(ggplot2); library(grid)})
.a    <- commandArgs(TRUE)
SWEEP <- .a[1]; I004 <- .a[2]; I008 <- .a[3]; OUTPDF <- .a[4]
source("scripts/analysis/lsr_palette.R")

ARMS <- data.frame(
  method = c("fb_xregion_id_l03","fb_xregion_id_l06","fb_xregion_id",
             "fb_xregion_id_l25","fb_xregion_id_l50","fb_xregion_id_l100","fb_xregion_id_l200"),
  lambda = c(0.03, 0.06, 0.1223, 0.25, 0.5, 1.0, 2.0))
SHOW <- ARMS$method                                   # all seven now
LAB  <- c(fb_xregion_id_l03 = "λ = 0.03",  fb_xregion_id_l06  = "λ = 0.06",
          fb_xregion_id     = "λ = 0.1223 (tuned)", fb_xregion_id_l25 = "λ = 0.25",
          fb_xregion_id_l50 = "λ = 0.5",   fb_xregion_id_l100 = "λ = 1",
          fb_xregion_id_l200 = "λ = 2")
COL <- c(fb_xregion_id_l03 = "#56B4E9", fb_xregion_id_l06 = "#009E73",
         fb_xregion_id = PAL7[["fb_xregion"]], fb_xregion_id_l25 = "#0072B2",
         fb_xregion_id_l50 = "#E69F00", fb_xregion_id_l100 = "#000000",
         fb_xregion_id_l200 = "#CC79A7")
SHP <- c(fb_xregion_id_l03 = 17, fb_xregion_id_l06 = 15, fb_xregion_id = SHP7[["fb_xregion"]],
         fb_xregion_id_l25 = 18, fb_xregion_id_l50 = 8, fb_xregion_id_l100 = 6,
         fb_xregion_id_l200 = 4)
stopifnot(!anyDuplicated(COL), !anyDuplicated(SHP))
REFS <- c("polyfun_ldsc", "beatrice")
MODEL_LAB <- c(sparse = "Sparse", sparse_inf = "Sparse + infinitesimal")
ARM_LAB   <- c(binary = "Binary annotations", continuous = "Continuous annotations")
rename_sweep <- function(m) ifelse(m == "fb_xregion", "fb_xregion_rerun", m)

L3s <- readRDS(file.path(SWEEP, "combined_scenario_metrics.rds")); L3s$method <- rename_sweep(L3s$method)
JOBS <- unique(L3s$job_dir[L3s$annotation_type != "none"])
L3s <- L3s[L3s$job_dir %in% JOBS & L3s$n_failed < L3s$n_fits & L3s$method %in% ARMS$method, ]

## references: beatrice from Iteration 004, polyfun_ldsc from Iteration 008
r4 <- readRDS(file.path(I004, "combined_scenario_metrics.rds"))
r4 <- r4[r4$job_dir %in% JOBS & r4$method == "beatrice" & r4$n_failed < r4$n_fits, ]
r8 <- readRDS(file.path(I008, "combined_scenario_metrics.rds"))
r8 <- r8[r8$job_dir %in% JOBS & r8$method == "polyfun_ldsc" & r8$n_failed < r8$n_fits, ]
cat(sprintf("reference cells: beatrice %d (iter004), polyfun_ldsc %d (iter008)\n", nrow(r4), nrow(r8)))

cols <- c("model","annotation_type","method","ap")
apd  <- rbind(L3s[, cols], r4[, cols], r8[, cols])
AP <- do.call(rbind, lapply(split(apd, list(apd$model, apd$annotation_type, apd$method), drop = TRUE),
  function(r) { v <- r$ap[is.finite(r$ap)]
    data.frame(model = r$model[1], annotation_type = r$annotation_type[1], method = r$method[1],
               m = mean(v), se = sd(v)/sqrt(length(v)), cells = length(v)) }))

a1 <- AP[AP$method %in% c(SHOW, REFS), ]
ylab <- c(LAB, polyfun_ldsc = "polyfun_ldsc", beatrice = "beatrice")
a1$row <- factor(ylab[a1$method], ylab[c(REFS, rev(SHOW))])
a1$grp <- factor(ifelse(a1$method %in% SHOW, a1$method, "reference"), c(SHOW, "reference"))
a1$model <- factor(MODEL_LAB[as.character(a1$model)], MODEL_LAB)
a1$arm   <- factor(ARM_LAB[as.character(a1$annotation_type)], ARM_LAB)

th <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linewidth = .25, colour = "grey90"),
        panel.border = element_rect(colour = "grey70", linewidth = .35),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 6.6, face = "bold"),
        axis.title = element_text(size = 6.8), axis.text = element_text(size = 6.2),
        legend.text = element_text(size = 6.2), legend.key.size = unit(8, "pt"),
        legend.position = "right")
p1 <- ggplot(a1, aes(m, row, colour = grp, shape = grp)) +
  geom_errorbar(aes(xmin = m - 1.96*se, xmax = m + 1.96*se), orientation = "y",
                width = .42, linewidth = .4) +
  geom_point(size = 1.6, stroke = .55) +
  facet_grid(model ~ arm) +
  scale_colour_manual(values = c(COL[SHOW], reference = GREY), limits = c(SHOW, "reference"),
                      labels = c(unname(LAB[SHOW]), "main-study method"), name = NULL) +
  scale_shape_manual(values = c(SHP[SHOW], reference = 16), limits = c(SHOW, "reference"),
                     labels = c(unname(LAB[SHOW]), "main-study method"), name = NULL) +
  scale_x_continuous(breaks = seq(.5, .9, .1)) +
  labs(x = "Mean average precision over cells (95% CI)", y = NULL) +
  th + theme(panel.grid.major.y = element_line(linewidth = .2, colour = "grey93"))

if (grepl("\\.png$", OUTPDF)) png(OUTPDF, width = 7.28, height = 4.4, units = "in", res = 200, type = "cairo") else cairo_pdf(OUTPDF, width = 7.28, height = 4.4)
print(p1); invisible(dev.off())
cat("wrote", OUTPDF, "\n\n")
w <- a1[, c("model","arm","row","m")]
print(reshape(w, idvar = c("row"), timevar = "model", direction = "wide")[0, ], row.names = FALSE)
for (s in levels(a1$arm)) { cat("--", s, "--\n")
  d <- a1[a1$arm == s, c("model","row","m")]
  print(reshape(d, idvar = "row", timevar = "model", direction = "wide"), row.names = FALSE, digits = 3) }
