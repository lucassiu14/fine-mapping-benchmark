# =============================================================================
# scripts/analysis/iter003_mass_location.R
#
# Answers: is a total_mass_ratio of ~2 actually harmful, given AP, FDR and
# high-PIP reliability all look acceptable?
#
# It is not the ratio that matters but WHERE the excess mass sits. SuSiE and
# Functional BEATRICE have near-identical mass ratios in this stratum yet put
# their excess in completely different PIP bands - susie as a thin sub-0.1
# haze, FB concentrated in 0.5-0.9 where it corrupts credible sets.
#
# Second half checks whether the Optuna tuning simulation spans the cells where
# the failure actually lives. It does not.
#
# Run from the project root. Reads results/iter003/combined_pip_calibration.rds.
# NEVER pools across model (sparse/sparse_inf) or annotation type.
# =============================================================================

d <- readRDS("results/iter003/combined_pip_calibration.rds")

# --- strata (NEVER pool across model or annotation type) ---------------------
d$model <- ifelse(grepl("_sparse_inf_", d$job_dir), "sparse_inf", "sparse")
d$annot <- ifelse(grepl("_anNone_",       d$job_dir), "none",
           ifelse(grepl("_anBinary",      d$job_dir), "binary",
           ifelse(grepl("_anCont",  d$job_dir), "continuous", NA)))
d$ld    <- ifelse(grepl("_refInsample$",  d$job_dir), "insample",
           ifelse(grepl("_ref500$",       d$job_dir), "ref500",
           ifelse(grepl("_ref200$",       d$job_dir), "ref200", NA)))
stopifnot(!any(is.na(d$annot)), !any(is.na(d$ld)))

# The stratum the tuning study targets: sparse / binary annotations / n_ref=500.
s <- d[d$model == "sparse" & d$annot == "binary" & d$ld == "ref500", ]
meth <- c("fb_xregion", "fb_pooled", "functional_beatrice", "beatrice",
          "susie", "polyfun_est", "finemap")
s <- s[s$method %in% meth, ]

# PIP band from the bin's own edges (bin 1 = [0,0.1), bin 10 = [0.9,1]).
s$band <- cut(s$bin, breaks = c(0, 1, 5, 9, 10),
              labels = c("<0.1", "0.1-0.5", "0.5-0.9", ">=0.9"))

agg <- aggregate(cbind(sum_pip, n_causal, n) ~ method + band, s, sum)
tot <- aggregate(cbind(sum_pip, n_causal) ~ method,        s, sum)
tot$mass_ratio <- tot$sum_pip / tot$n_causal
tot$excess     <- tot$sum_pip - tot$n_causal

cat("=== stratum: sparse / binary annotations / ref500 (NOT pooled) ===\n\n")
cat(sprintf("%-21s %8s %9s %9s\n", "method", "massrat", "SumPIP", "excess"))
for (i in order(-tot$mass_ratio)) with(tot[i, ],
  cat(sprintf("%-21s %8.2f %9.0f %9.0f\n", method, mass_ratio, sum_pip, excess)))

cat("\n=== where does the EXCESS mass sit? (share of total excess per band) ===\n")
cat(sprintf("%-21s %9s %9s %9s %9s\n", "method", "<0.1", "0.1-0.5", "0.5-0.9", ">=0.9"))
for (m in tot$method[order(-tot$mass_ratio)]) {
  a <- agg[agg$method == m, ]
  ex <- setNames(a$sum_pip - a$n_causal, as.character(a$band))
  ex <- ex[c("<0.1","0.1-0.5","0.5-0.9",">=0.9")]
  ex[is.na(ex)] <- 0
  cat(sprintf("%-21s %8.0f%% %8.0f%% %8.0f%% %8.0f%%\n", m,
              100*ex[1]/sum(ex), 100*ex[2]/sum(ex), 100*ex[3]/sum(ex), 100*ex[4]/sum(ex)))
}

cat("\n=== calibration within band: observed frac causal vs nominal PIP ===\n")
cat(sprintf("%-21s %-9s %9s %9s %10s\n","method","band","meanPIP","obs_frac","n_variants"))
for (m in c("fb_xregion","fb_pooled","susie")) {
  a <- s[s$method == m, ]
  k <- aggregate(cbind(sum_pip, n_causal, n) ~ band, a, sum)
  for (i in seq_len(nrow(k))) with(k[i, ],
    cat(sprintf("%-21s %-9s %9.3f %9.3f %10.0f\n",
                m, as.character(band), sum_pip/n, n_causal/n, n)))
}
d <- readRDS("results/iter003/combined_pip_calibration.rds")
# Match the tuning sim EXACTLY: sparse, binary annotations, enrichment 10.8, ref500.
s <- d[grepl("_sparse_anBinary_e10.8_ref500$", d$job_dir), ]
s$S <- as.numeric(s$S); s$phi <- as.numeric(s$phi); s$p <- as.numeric(s$region_size)
s$corner <- ifelse(s$S %in% c(1,3) & s$phi %in% c(0.05,0.2) & s$p %in% c(100,200,400,500),
                   "tuning-sim cells", "rest of stratum")

f <- function(x, m) {
  a <- x[x$method == m, ]; hi <- a[a$bin == 10, ]
  data.frame(method = m,
             mass_ratio = sum(a$sum_pip)/sum(a$n_causal),
             hi_pip_obs = sum(hi$n_causal)/sum(hi$n),
             n_hi_calls = sum(hi$n))
}
for (grp in c("tuning-sim cells", "rest of stratum")) {
  cat("\n=== ", grp, "  (sparse / binary e10.8 / ref500) ===\n", sep="")
  print(do.call(rbind, lapply(c("fb_xregion","fb_pooled","susie"),
        function(m) f(s[s$corner == grp, ], m))), row.names = FALSE, digits = 3)
}

cat("\n=== fb_xregion hi-PIP reliability, cell by cell ===\n")
a <- s[s$method == "fb_xregion" & s$bin == 10, ]
k <- aggregate(cbind(n, n_causal) ~ S + phi, a, sum)
k$hi_pip_obs <- k$n_causal / k$n
k$in_tuning_sim <- ifelse(k$S %in% c(1,3) & k$phi %in% c(0.05,0.2), "yes", "-")
print(k[order(k$phi, k$S), c("S","phi","n","hi_pip_obs","in_tuning_sim")],
      row.names = FALSE, digits = 3)
