suppressPackageStartupMessages(library(ggplot2))
`%||%` <- function(x,y) if (is.null(x)) y else x
source("scripts/analysis/iter004_lib.R")
L3 <- readRDS("results/iter004/combined_scenario_metrics.rds"); L3<-L3[!L3$method %in% EXCLUDE_METHODS,]
L3$arm <- factor(ifelse(L3$model=="sparse","Sparse","Sparse + infinitesimal"),
                 levels=c("Sparse","Sparse + infinitesimal"))
FAM <- c("polyfun_oracle"="#000000","fb_xregion"="#0072B2","fb_pooled"="#56B4E9",
         "functional_beatrice"="#D55E00","beatrice"="#E69F00","susie"="#009E73","polyfun_est"="#CC79A7")
COL <- function(m) ifelse(m %in% names(FAM), FAM[m], "#BBBBBB")
INK<-"#222222";MUTED<-"#666666"
th<-theme_minimal(base_size=8.5)+theme(panel.grid.minor=element_blank(),
  panel.grid.major=element_line(colour="#E8E8E8",linewidth=.3),
  axis.title=element_text(colour=INK),axis.text=element_text(colour=MUTED),
  strip.text=element_text(colour=INK,face="bold",size=8),
  plot.title=element_text(colour=INK,face="bold",size=10),
  plot.subtitle=element_text(colour=MUTED,size=7.5),legend.position="none")
ag <- function(v) { a <- aggregate(L3[[v]] ~ method + arm, L3, function(z) mean(z,na.rm=TRUE))
  names(a)[3] <- "val"; a$part <- v; a }
mu <- rbind(transform(ag("rel"), part="REL (miscalibration, lower better)"),
            transform(ag("res"), part="RES (resolution, higher better)"))
mu$col <- COL(mu$method)
o <- aggregate(res ~ method, L3, function(z) mean(z,na.rm=TRUE))
mu$method <- factor(mu$method, levels=o$method[order(o$res)])
p <- ggplot(mu, aes(val, method, fill=col)) + geom_col(width=.68) + scale_fill_identity() +
  facet_grid(arm ~ part, scales="free_x") + th +
  labs(title="Murphy decomposition of the Brier score", x="component", y=NULL,
       subtitle="RES exposes calibration achieved by saying nothing: a method emitting S/p everywhere has REL = 0 and RES = 0.")
ggsave("results/iter004_full/report/figs/f17_murphy.pdf", p, width=9.5, height=6.5, device="pdf")
cat("murphy figure written\n")
