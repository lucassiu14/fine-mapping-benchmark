#!/usr/bin/env Rscript
# =============================================================================
# cs_figures.R -- the four credible-set figures, rebuilt on the merged tables.
#
#   Rscript cs_figures.R <out_dir> <L3_with_power> <native_credible_sets>
#
# Uniform rule: one 0.95-mass set per fit, from set_hit_95 / set_size_95 and the
# power column. Native rule: each method's own convention, from the native table.
# Seven displayed methods, palette and shapes from lsr_palette.R.
# =============================================================================
suppressMessages({library(ggplot2); library(grid)})
.a <- commandArgs(TRUE); FIG <- .a[1]; L3F <- .a[2]; NATF <- .a[3]
EXT <- if (length(.a) >= 4) .a[4] else "pdf"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
source("scripts/analysis/lsr_palette.R")
FITS <- 20

NAMES <- c(susie = "SuSiE", polyfun_ldsc = "PolyFun + SuSiE", finemap_inf = "FINEMAP-inf",
           polyfun_oracle = "Oracle PolyFun + SuSiE", sbayesrc = "SBayesRC",
           beatrice = "BEATRICE", fb_xregion = "Functional BEATRICE")
MODEL_LAB <- c(sparse="Sparse", sparse_inf="Sparse + infinitesimal")
ARM_LAB   <- c(none="No annotations", binary="Binary annotations",
               continuous="Continuous annotations")
strata <- function(d) { d$model <- factor(MODEL_LAB[as.character(d$model)], MODEL_LAB)
  d$arm <- factor(ARM_LAB[as.character(d$annotation_type)], ARM_LAB); d }

## ---- uniform 0.95-mass rule, aggregated to cells then to strata ------------
L3 <- readRDS(L3F)
L3 <- L3[L3$method %in% M7 & L3$n_failed < L3$n_fits, ]
uni <- do.call(rbind, lapply(split(L3, list(L3$model, L3$annotation_type, L3$method), drop=TRUE),
  function(r) data.frame(model=r$model[1], annotation_type=r$annotation_type[1],
    method=as.character(r$method[1]),
    cov = mean(r$set_hit_95, na.rm=TRUE)/FITS,
    size= mean(r$set_size_95, na.rm=TRUE),
    pow = mean(r$power, na.rm=TRUE))))

th <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="grey70", linewidth=.35),
        strip.background = element_rect(fill="grey94", colour=NA),
        strip.text = element_text(size=6.6, face="bold"),
        axis.title = element_text(size=6.8), axis.text = element_text(size=6.2),
        legend.text = element_text(size=6.0), legend.key.size = unit(7,"pt"),
        legend.position = "bottom", legend.box.margin = margin(t = -4),
        legend.margin = margin(t = 0, b = 0), legend.key.spacing.x = unit(3, "pt"))
pal <- setNames(PAL7[M7], NAMES[M7]); shp <- setNames(SHP7[M7], NAMES[M7])
XS  <- scale_x_log10(breaks=c(1,2,5,10,50,200,1000), labels=c("1","2","5","10","50","200","1000"),
                     expand=expansion(mult=c(.07,.07)))
save <- function(p, f, h=4.3) {
  path <- file.path(FIG, paste0(f, ".", EXT))
  if (EXT == "png") png(path, width=7.28, height=h, units="in", res=200, type="cairo")
  else cairo_pdf(path, width=7.28, height=h)
  print(p); invisible(dev.off()); cat("wrote", path, "\n")
}
u <- strata(uni); u$method <- factor(NAMES[u$method], unname(NAMES[M7]))

mk_uni <- function(yv, ylab, hline) {
  g <- ggplot(u, aes(size, .data[[yv]], colour=method, shape=method))
  if (!is.na(hline)) g <- g + geom_hline(yintercept=hline, linetype="22", colour=REF, linewidth=.4)
  g + geom_point(size=1.7, stroke=.6) + facet_grid(model ~ arm) +
    scale_colour_manual(values=pal, name=NULL, guide=guide_legend(nrow=1)) +
    scale_shape_manual(values=shp, name=NULL, guide=guide_legend(nrow=1)) +
    XS + labs(x="Mean 95% credible-set size (variants, log scale)", y=ylab) + th
}
save(mk_uni("cov","Coverage: sets containing a causal variant", 0.95), "fig_coverage")
save(mk_uni("pow","Power: causal variants captured by a set", NA),     "fig_power")

## ---- each method's own convention, drawn against the uniform rule ----------
nat <- readRDS(NATF); nat <- nat[nat$method %in% M7, ]
k <- c("model","annotation_type","method")
m <- merge(uni, nat[, c(k,"coverage_per_set","power","mean_size")], by=k)
m <- m[is.finite(m$mean_size) & is.finite(m$coverage_per_set), ]
cat("native rows drawn:", nrow(m), " methods:", paste(sort(unique(m$method)), collapse=", "), "\n")
m <- strata(m); m$method <- factor(NAMES[m$method], unname(NAMES[M7]))

mk_nat <- function(y_uni, y_nat, ylab, hline) {
  ggplot(m) +
    { if (!is.na(hline)) geom_hline(yintercept=hline, linetype="22", colour=REF, linewidth=.4) } +
    geom_segment(aes(x=size, y=.data[[y_uni]], xend=mean_size, yend=.data[[y_nat]], colour=method),
                 arrow=arrow(length=unit(3.2,"pt"), type="closed"), linewidth=.45) +
    geom_point(aes(size, .data[[y_uni]], colour=method), shape=1, size=1.5, stroke=.55) +
    geom_point(aes(mean_size, .data[[y_nat]], colour=method, shape=method), size=1.7, stroke=.6) +
    facet_grid(model ~ arm) +
    scale_colour_manual(values=pal, name=NULL, guide=guide_legend(nrow=1)) +
    scale_shape_manual(values=shp, name=NULL, guide=guide_legend(nrow=1)) +
    XS + labs(x="Mean credible-set size (variants, log scale)", y=ylab) + th
}
save(mk_nat("cov","coverage_per_set","Coverage: sets containing a causal variant", 0.95), "fig_cs_native_coverage")
save(mk_nat("pow","power","Power: causal variants captured by a set", NA),               "fig_cs_native_power")

cat("\n--- uniform rule, sparse/binary ---\n")
print(uni[uni$model=="sparse" & uni$annotation_type=="binary", ], row.names=FALSE, digits=3)
