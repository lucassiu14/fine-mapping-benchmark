#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_native_cs_figure.R
#
# What changes when each method's OWN credible sets are used in place of the
# uniform 0.95-mass reconstruction.
#
#   Rscript scripts/analysis/iter004_native_cs_figure.R <out_dir> [native_rds] [l3_rds]
#
# Drawn as a shift rather than as two separate panels: each method is one
# segment running from where the uniform rule puts it to where its own
# convention puts it, with the arrowhead at the latter. The argument of the
# section is about the DIFFERENCE between the two constructions, and a reader
# comparing two scatter plots across a page cannot see a difference that moves
# points in opposite directions for different methods.
#
# Coverage is per reported set under both rules, so a method reporting three
# sets contributes three observations; power is per fit. marginal_z is omitted:
# it reports one set of about 765 variants under either rule, which is a
# statement about the axis rather than the method.
# =============================================================================

.a  <- commandArgs(TRUE)
FIG <- if (length(.a) >= 1) .a[1] else "figures"
NAT <- if (length(.a) >= 2) .a[2] else "results/iter004/native_credible_sets.rds"
L3F <- if (length(.a) >= 3) .a[3] else "results/iter004/combined_scenario_metrics_with_power.rds"
suppressMessages(library(ggplot2))
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

DROP <- c("carma", "fb_pooled", "polyfun_est", "marginal_z")
HI   <- c(susie = "#0072B2", susie_inf = "#56B4E9",
          fb_xregion = "#D55E00", sbayesrc = "#009E73")
GREY <- "#9AA0A6"; REF <- "#C1272D"

nat <- readRDS(NAT)
nat <- nat[!nat$method %in% DROP, ]

L3 <- readRDS(L3F)
L3 <- L3[L3$n_failed < L3$n_fits & !L3$method %in% DROP, ]
uni <- do.call(rbind, lapply(
  split(L3, list(L3$model, L3$annotation_type, L3$method), drop = TRUE),
  function(x) data.frame(model = x$model[1], annotation_type = x$annotation_type[1],
                         method = as.character(x$method[1]),
                         u_cov = mean(x$set_hit_95, na.rm = TRUE) / 20,
                         u_pow = mean(x$power, na.rm = TRUE),
                         u_size = mean(x$set_size_95, na.rm = TRUE))))

d <- merge(nat[, c("model","annotation_type","method",
                   "coverage_per_set","power","mean_size")],
           uni, by = c("model","annotation_type","method"))
d <- d[is.finite(d$mean_size) & is.finite(d$u_size) & d$mean_size > 0, ]
d$model <- factor(d$model, c("sparse","sparse_inf"),
                  c("Sparse","Sparse + infinitesimal"))
d$arm <- factor(d$annotation_type, c("none","binary","continuous"),
                c("No annotations","Binary annotations","Continuous annotations"))
meths <- sort(unique(d$method)); d$id <- match(d$method, meths)
d$grp <- factor(ifelse(d$method %in% names(HI), d$method, "other methods"),
                c(names(HI), "other methods"))
KEY <- paste(sprintf("%d %s", seq_along(meths), meths), collapse = "   ·   ")
.h <- ceiling(length(meths) / 2)
KEY <- paste(paste(sprintf("%d %s", 1:.h, meths[1:.h]), collapse = "   ·   "),
             paste(sprintf("%d %s", (.h+1):length(meths), meths[(.h+1):length(meths)]),
                   collapse = "   ·   "), sep = "\n")

th <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey70", linewidth = .35),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 7.2, face = "bold"),
        axis.title = element_text(size = 7.8), axis.text = element_text(size = 6.8),
        legend.text = element_text(size = 6.9), legend.key.size = unit(8, "pt"),
        legend.position = "right",
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 7.2, colour = "grey35"),
        plot.caption = element_text(size = 6.2, colour = "grey30", hjust = 0,
                                    margin = margin(t = 6)))

shift <- function(y_nat, y_uni, ylab, ttl, sub, hline = NA) {
  g <- ggplot(d, aes(colour = grp))
  if (!is.na(hline))
    g <- g + geom_hline(yintercept = hline, linetype = "22",
                        colour = REF, linewidth = .4)
  g +
    geom_segment(aes(x = u_size, y = .data[[y_uni]],
                     xend = mean_size, yend = .data[[y_nat]]),
                 arrow = arrow(length = unit(3.2, "pt"), type = "closed"),
                 linewidth = .38, alpha = .85) +
    geom_point(aes(u_size, .data[[y_uni]]), size = .7, shape = 1, stroke = .4) +
    geom_text(aes(mean_size, .data[[y_nat]], label = id),
              size = 2.1, fontface = "bold", vjust = -0.8, show.legend = FALSE) +
    facet_grid(model ~ arm) +
    scale_colour_manual(values = c(HI, `other methods` = GREY),
                        limits = c(names(HI), "other methods"), name = NULL) +
    scale_x_log10(breaks = c(1, 2, 5, 10, 50, 200, 1000),
                  labels = c("1","2","5","10","50","200","1000")) +
    labs(x = "Mean credible-set size (variants, log scale)", y = ylab,
         title = ttl, subtitle = sub, caption = KEY) + th
}

p1 <- shift("coverage_per_set", "u_cov", "Coverage: sets containing a causal variant",
  "Coverage under each method's own credible sets, against the uniform rule",
  paste("Each arrow runs from the uniform 0.95-mass reconstruction (open circle)",
        "to the method's own convention (number).\nDashed line is the nominal",
        "0.95 level. marginal_z is omitted: about 765 variants per set under",
        "either rule."), hline = 0.95)
ggsave(file.path(FIG, "fig_cs_native_coverage.pdf"), p1, width = 9.2, height = 5.8)

p2 <- shift("power", "u_pow", "Power: causal variants captured by any set",
  "Power under each method's own credible sets, against the uniform rule",
  paste("Each arrow runs from the uniform 0.95-mass reconstruction (open circle)",
        "to the method's own convention (number).\nA method reporting several",
        "sets can capture several signals; one set cannot."))
ggsave(file.path(FIG, "fig_cs_native_power.pdf"), p2, width = 9.2, height = 5.8)
cat("fig_cs_native_coverage and fig_cs_native_power written\n")
