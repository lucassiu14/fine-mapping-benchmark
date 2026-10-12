#!/usr/bin/env Rscript
# =============================================================================
# fig_example3.R  --  one simulated region, end to end, as a 2 x 4 grid.
#
#   Rscript fig_example3.R <out_pdf> <example_rds> <job_label> <scenario> <region>
#
# Top row:    the GWAS z-scores, then the three methods that cannot use
#             annotations.
# Bottom row: the four that can, Functional BEATRICE first.
#
# PIPs come from the stored PIP tails, which keep every PIP at or above 0.01
# and every causal variant whatever its PIP.
# =============================================================================

suppressMessages({library(ggplot2); library(grid); library(gtable)})
.a  <- commandArgs(TRUE)
OUT <- .a[1]; EX <- .a[2]; JOB <- .a[3]; SC <- as.integer(.a[4]); RG <- as.integer(.a[5])
PAL_F <- Sys.getenv("LSR_PALETTE", "scripts/analysis/lsr_palette.R")
if (!file.exists(PAL_F)) stop("cannot find lsr_palette.R at '", PAL_F, "'")
source(PAL_F)

e  <- readRDS(EX)
p  <- length(e$z); cz <- sort(as.integer(e$causal))
BLIND <- c("susie", "finemap_inf", "beatrice")                       # top row, after z
ANN   <- c("fb_xregion", "polyfun_oracle", "polyfun_ldsc", "sbayesrc")  # bottom row
NAMES <- c(susie = "SuSiE", finemap_inf = "FINEMAP-inf", beatrice = "BEATRICE",
           fb_xregion = "Functional BEATRICE", polyfun_oracle = "Oracle PolyFun + SuSiE",
           polyfun_ldsc = "PolyFun + SuSiE", sbayesrc = "SBayesRC")
NEW <- c("polyfun_ldsc", "polyfun_oracle", "sbayesrc")   # these come from iteration 008

get_pip <- function(m) {
  f <- file.path(if (m %in% NEW) "results/iter008_fidelity/piptail" else "results/iter004/piptail",
                 paste0("piptail_", JOB, ".rds"))
  o <- readRDS(f)
  for (z in o$tail)
    if (z$method == m && z$scenario_id == SC && z$region_id == RG) {
      v <- numeric(z$n_variants); v[z$idx] <- z$pip; return(v)
    }
  stop("fit not found: ", m)
}
pips <- lapply(c(BLIND, ANN), get_pip); names(pips) <- c(BLIND, ANN)

base <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = .2, colour = "grey92"),
        panel.border = element_rect(colour = "grey70", linewidth = .35),
        plot.title = element_text(size = 6.8, face = "bold", hjust = .5,
                                  margin = margin(b = 2)),
        axis.title = element_text(size = 6.6), axis.text = element_text(size = 6.0),
        legend.position = "none",
        plot.margin = margin(3, 5, 2, 3))

cmark <- data.frame(x = cz)

## ---- z-scores -------------------------------------------------------------
zd <- data.frame(x = seq_len(p), z = e$z, causal = seq_len(p) %in% cz)
p_z <- ggplot(zd, aes(x, z)) +
  geom_vline(data = cmark, aes(xintercept = x), colour = "grey60", linetype = "22",
             linewidth = .3) +
  geom_point(data = zd[!zd$causal, ], colour = "grey45", size = .5) +
  geom_point(data = zd[zd$causal, ], shape = 21, colour = REF, fill = NA,
             size = 2.4, stroke = .8) +
  scale_x_continuous(limits = c(1, p), breaks = c(1, 500, 1000)) +
  labs(title = "GWAS z-scores", x = "Variant index", y = "z") + base

## ---- PIP panels -----------------------------------------------------------
pip_panel <- function(m, show_y) {
  v  <- pips[[m]]
  d  <- data.frame(x = seq_len(p), y = v, causal = seq_len(p) %in% cz)
  d  <- d[d$y > 0.002 | d$causal, ]
  g <- ggplot(d, aes(x, y)) +
    geom_segment(data = d[!d$causal, ], aes(x = x, xend = x, y = 0, yend = y),
                 colour = PAL7[[m]], linewidth = .3) +
    geom_segment(data = d[d$causal, ], aes(x = x, xend = x, y = 0, yend = y),
                 colour = REF, linewidth = .45) +
    geom_point(data = d[d$causal, ], shape = 21, colour = REF, fill = NA,
               size = 2.2, stroke = .8) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, .5, 1),
                       labels = c("0", "0.5", "1"), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1, p), breaks = c(1, 500, 1000)) +
    labs(title = NAMES[m], x = "Variant index", y = if (show_y) "PIP" else NULL) + base
  if (!show_y) g <- g + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank())
  g
}
# Every PIP panel keeps its own axis. The top row's leftmost panel is the
# z-score panel, so without this the three PIP panels beside it would show no
# PIP scale anywhere in the row; it also keeps all seven panels the same size.
top    <- lapply(BLIND, pip_panel, show_y = TRUE)
bottom <- lapply(ANN,   pip_panel, show_y = TRUE)

## ---- compose: four equal columns, two equal rows ---------------------------
tg <- gtable(widths = unit(rep(1, 4), "null"), heights = unit(c(1, 1), "null"))
tg <- gtable_add_grob(tg, ggplotGrob(p_z), t = 1, l = 1, name = "z")
for (i in seq_along(top))
  tg <- gtable_add_grob(tg, ggplotGrob(top[[i]]), t = 1, l = i + 1, name = paste0("t", i))
for (i in seq_along(bottom))
  tg <- gtable_add_grob(tg, ggplotGrob(bottom[[i]]), t = 2, l = i, name = paste0("b", i))

if (grepl("\\.png$", OUT)) png(OUT, width = 7.28, height = 3.1, units = "in", res = 200, type = "cairo") else PDF_DEVICE(OUT, width = 7.28, height = 3.1)
grid.newpage(); grid.draw(tg)
invisible(dev.off())
cat("wrote", OUT, "\n")
for (m in c(BLIND, ANN)) {
  v <- pips[[m]]
  cat(sprintf("%-22s causal %s | max non-causal %.2f at %d\n", NAMES[m],
      paste(sprintf("%.2f", v[cz]), collapse = " "), max(v[-cz]),
      which.max(replace(v, cz, -1))))
}
