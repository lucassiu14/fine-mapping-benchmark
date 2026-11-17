#!/usr/bin/env Rscript
# Is sigma_sq's ORACLE value known in closed form?
#
# After the RSS z-scale fix (joint_trainer.py rss_rescale_z), FB's ABF is the
# RSS model z ~ N(Rb, R) with b = beta/se and a prior b_gamma ~ N(0, kappa I),
# kappa = n * sigma_sq. Since se_j^2 = sigma^2_resid / (x_j'x_j),
#
#     Var(b_j) = Var(beta_j) * x_j'x_j / sigma^2_resid
#              = n * [ Var(beta_j) * (x_j'x_j / n) / sigma^2_resid ]
#
# so sigma_sq is exactly susie's `scaled_prior_variance`. For this package's
# PVE-calibrated sparse simulator, sigma^2_resid = Var(Xb) * (1 - phi) / phi
# (simulate_phenotypes.R:889) with STANDARDISED columns, so x_j'x_j / n = 1,
# Var(Xb) ~= S * Var(beta_j), and
#
#     sigma_sq* = phi / (S * (1 - phi))
#
# with the effect distribution and effect_variance cancelling, because the
# residual variance is calibrated FROM Var(Xb).
#
# This script checks that algebra against the stored truth of every region of a
# simulation object: it computes Var(beta)/sigma^2_resid from X, y and
# beta_true, with no appeal to the formula, and compares.
#
# Expect per-region scatter: mean(beta^2) over S draws of a normal effect has
# relative sd sqrt(2/S), i.e. 63% at S = 5, so single regions can sit a factor
# of 2 out. The MEDIAN over regions is the testable quantity.
#
# Usage: Rscript scripts/tuning/test_sigma_oracle.R [sim.rds] [tolerance]

args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args) >= 1) args[1] else "tuning/full25.rds"
tol  <- if (length(args) >= 2) as.numeric(args[2]) else 0.10

sim <- readRDS(path)
stopifnot(all(c("genotypes", "scenarios") %in% names(sim)))
gt <- sim$genotypes
sc <- sim$scenarios

rows <- list()
cat(sprintf("%-5s %-3s %-5s %-4s %12s %12s %7s %7s\n",
            "scen", "S", "phi", "reg", "measured", "phi/(S(1-phi))",
            "ratio", "pve"))
cat(strrep("-", 64), "\n")
for (k in seq_along(sc)) {
  s <- sc[[k]]
  if (!identical(s$model, "sparse")) next   # sparse_inf adds a tau^2 component
  for (j in seq_along(s$regions)) {
    r  <- s$regions[[j]]
    bt <- r$truth$beta_true
    if (is.null(bt) || is.null(r$y) || j > length(gt)) next
    g  <- as.numeric(gt[[j]]$X %*% bt)
    s2 <- var(r$y) - var(g)                 # residual variance
    cz <- bt[bt != 0]
    if (!length(cz) || !(s2 > 0)) next
    measured <- mean(cz^2) / s2
    formula  <- s$phi / (length(cz) * (1 - s$phi))
    cat(sprintf("%-5d %-3d %-5.2f %-4d %12.5f %12.5f %7.3f %7.3f\n",
                s$scenario_id, s$S, s$phi, j, measured, formula,
                measured / formula, var(g) / var(r$y)))
    rows[[length(rows) + 1L]] <- c(measured = measured, formula = formula)
  }
}
stopifnot(length(rows) > 0)
d <- as.data.frame(do.call(rbind, rows))
ratio <- d$measured / d$formula
med <- median(ratio)

cat("\n")
cat(sprintf("regions                : %d\n", nrow(d)))
cat(sprintf("median ratio           : %.4f  (target 1, tolerance %.2f)\n", med, tol))
cat(sprintf("IQR of ratio           : [%.3f, %.3f]\n",
            quantile(ratio, 0.25), quantile(ratio, 0.75)))
cat(sprintf("worst single region    : %.3f\n", max(abs(log(ratio)))))
cat(sprintf("oracle span over cells : %.5f to %.5f\n",
            min(d$formula), max(d$formula)))
cat(sprintf("cells above sigma_sq's historical 0.5 ceiling: %d of %d regions\n",
            sum(d$formula > 0.5), nrow(d)))

if (abs(med - 1) > tol) {
  stop(sprintf("FAIL: median ratio %.4f is more than %.2f from 1", med, tol))
}
cat("\nPASS: sigma_sq* = phi/(S(1-phi)) reproduces the measured scaled prior variance.\n")
