`%||%` <- function(x,y) if (is.null(x)) y else x
source("scripts/analysis/iter004_lib.R")
L1 <- readRDS("results/iter004/combined_fit_metrics.rds")
L3 <- readRDS("results/iter004/combined_scenario_metrics.rds")
ok <- function(lab, cond, det="") cat(sprintf("[%s] %-58s %s\n", if(isTRUE(cond)) "PASS" else "FAIL", lab, det))
f <- L1[!is.na(L1$ap), ]

cat("=== C. METRIC INTERNAL CONSISTENCY (per fit) ===\n")
for (t in c("50","90","99")) {
  tp <- f[[paste0("tp_at_",t)]]; ns <- f[[paste0("nsel_at_",t)]]
  ok(sprintf("tp_at_%s <= nsel_at_%s", t, t), all(tp <= ns, na.rm=TRUE),
     sprintf("%d violations", sum(tp > ns, na.rm=TRUE)))
  ok(sprintf("tp_at_%s <= n_causal", t), all(tp <= f$n_causal, na.rm=TRUE),
     sprintf("%d violations", sum(tp > f$n_causal, na.rm=TRUE)))
}
ok("nsel monotone decreasing in threshold (50>=90>=99)",
   all(f$nsel_at_50 >= f$nsel_at_90 & f$nsel_at_90 >= f$nsel_at_99, na.rm=TRUE),
   sprintf("%d violations", sum(!(f$nsel_at_50 >= f$nsel_at_90 & f$nsel_at_90 >= f$nsel_at_99), na.rm=TRUE)))
ok("AP in [0,1]", all(f$ap >= 0 & f$ap <= 1, na.rm=TRUE),
   sprintf("range [%.4f, %.4f]", min(f$ap,na.rm=TRUE), max(f$ap,na.rm=TRUE)))
ok("set_size_95 <= n_variants", all(f$set_size_95 <= f$n_variants, na.rm=TRUE),
   sprintf("%d violations", sum(f$set_size_95 > f$n_variants, na.rm=TRUE)))
ok("band counts sum to n_variants",
   {s <- f$n_band_lo+f$n_band_mid+f$n_band_hi+f$n_band_top; all(abs(s-f$n_variants)<1e-6, na.rm=TRUE)},
   sprintf("%d violations", sum(abs(f$n_band_lo+f$n_band_mid+f$n_band_hi+f$n_band_top - f$n_variants)>1e-6, na.rm=TRUE)))
ok("causal band counts sum to n_causal",
   {s <- f$c_band_lo+f$c_band_mid+f$c_band_hi+f$c_band_top; all(abs(s-f$n_causal)<1e-6, na.rm=TRUE)},
   sprintf("%d violations", sum(abs(f$c_band_lo+f$c_band_mid+f$c_band_hi+f$c_band_top - f$n_causal)>1e-6, na.rm=TRUE)))

cat("\n=== D. SANITY OF THE SCIENCE ===\n")
ap <- tapply(L3$ap, L3$method, mean, na.rm=TRUE)
ok("marginal_z and abf are the two weakest methods",
   all(c("marginal_z","abf") %in% names(sort(ap))[1:2]),
   paste(names(sort(ap))[1:3], collapse=" < "))
d <- L3[L3$model=="sparse" & L3$annotation_type!="none", ]
ok("polyfun_oracle is the best method on annotated sparse",
   names(which.max(tapply(d$ap, d$method, mean, na.rm=TRUE))) == "polyfun_oracle",
   sprintf("best = %s", names(which.max(tapply(d$ap, d$method, mean, na.rm=TRUE)))))
ok("AP decreases monotonically with S (pooled)",
   {m <- tapply(L3$ap, L3$S, mean, na.rm=TRUE); all(diff(m) < 0)},
   paste(sprintf("S=%s:%.3f", names(tapply(L3$ap,L3$S,mean,na.rm=TRUE)), tapply(L3$ap,L3$S,mean,na.rm=TRUE)), collapse=" "))
ok("AP increases monotonically with phi (pooled)",
   {m <- tapply(L3$ap, L3$phi, mean, na.rm=TRUE); all(diff(m) > 0)},
   paste(sprintf("phi=%s:%.3f", names(tapply(L3$ap,L3$phi,mean,na.rm=TRUE)), tapply(L3$ap,L3$phi,mean,na.rm=TRUE)), collapse=" "))
di <- L3[L3$model=="sparse_inf", ]
ok("AP increases with p_causal (more sparse signal = easier)",
   {m <- tapply(di$ap, di$p_causal, mean, na.rm=TRUE); all(diff(m) > 0)},
   paste(sprintf("pc=%s:%.3f", names(tapply(di$ap,di$p_causal,mean,na.rm=TRUE)), tapply(di$ap,di$p_causal,mean,na.rm=TRUE)), collapse=" "))
ok("sparse arm is easier than sparse_inf arm",
   mean(L3$ap[L3$model=="sparse"],na.rm=TRUE) > mean(L3$ap[L3$model=="sparse_inf"],na.rm=TRUE),
   sprintf("sparse %.4f vs sparse_inf %.4f", mean(L3$ap[L3$model=="sparse"],na.rm=TRUE), mean(L3$ap[L3$model=="sparse_inf"],na.rm=TRUE)))
