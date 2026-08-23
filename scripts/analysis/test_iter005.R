# =============================================================================
# >>> ITERATION 005 - EXPLORATORY FORK. DELETE WITH THE ITERATION. <<<
#
# scripts/analysis/test_iter005.R
#
# The iter004 suite re-pointed at iter005_lib.R, plus the stratification cases
# the fork exists for. Kept as a full copy rather than testing only the delta:
# the fork carries its OWN copy of every shared function, so a change to
# iter005_lib.R must be caught by iter005's tests - test_iter004.R would not see
# it. Iteration 004's suite and library are untouched.
#
#   Rscript scripts/analysis/test_iter005.R
# =============================================================================
#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/test_iter004.R
#
# The §5.8 validation cases, plus a synthetic split-plot test of the §7.2
# lambda_T machinery end to end.
#
# Every expected value below is HAND-COMPUTED from the definitions, not from a
# previous run of this code. A test that asserts what the code currently does is
# not a test.
#
#   Rscript scripts/analysis/test_iter004.R
# =============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x
source(file.path("R", "evaluate_extras.R"))
source(file.path("scripts", "analysis", "iter005_lib.R"))

PASS <- 0L; FAIL <- 0L
ok <- function(name, cond, detail = "") {
  if (isTRUE(cond)) { PASS <<- PASS + 1L; cat(sprintf("  PASS  %s\n", name)) }
  else { FAIL <<- FAIL + 1L
         cat(sprintf("  FAIL  %s%s\n", name,
                     if (nzchar(detail)) paste0("   [", detail, "]") else "")) }
}
eq <- function(a, b, tol = 1e-9) isTRUE(all.equal(a, b, tolerance = tol))

cat("\n=== §5.8 metric cases (expected values hand-computed) ===\n")

# AP = 1 exactly. Sorted 0.9(T) 0.5(F) 0.1(F): tp=1,1,1; prec=1,1/2,1/3;
# rec=1,1,1; AP = (1-0)*1 = 1.
r1 <- .compute_ap_exact(c(0.9, 0.5, 0.1), c(TRUE, FALSE, FALSE))
ok("AP = 1.0 when the causal ranks first", eq(r1, 1), sprintf("got %.6f", r1))

# AP = 1/3. Sorted 0.6(F) 0.4(F) 0.2(T): tp=0,0,1; prec at k=3 is 1/3;
# rec jumps 0 -> 1 there; AP = 1*(1/3).
r2 <- .compute_ap_exact(c(0.6, 0.4, 0.2), c(FALSE, FALSE, TRUE))
ok("AP = 1/3 when the causal ranks last", eq(r2, 1/3), sprintf("got %.6f", r2))

# §5.1b - per-fit AP averaged, NOT the merged ranking.
# Merged: pips .9(T) .6(F) .5(F) .4(F) .2(T) .1(F), S=2
#   tp = 1,1,1,1,2,2 ; prec = 1,.5,1/3,.25,.4,1/3 ; rec = .5,.5,.5,.5,1,1
#   AP_merged = .5*1 + .5*.4 = 0.7
# Correct answer is mean(1, 1/3) = 0.666...  The two MUST differ or the test
# cannot detect the fallacy.
merged <- .compute_ap_exact(c(0.9, 0.6, 0.5, 0.4, 0.2, 0.1),
                            c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE))
ok("merged-ranking AP is 0.7 (the WRONG answer)", eq(merged, 0.7),
   sprintf("got %.6f", merged))
ok("per-fit mean AP is 2/3 and differs from merged", eq(mean(c(r1, r2)), 2/3))

# §5.1 - the grid bias. A causal at 0.0032 with 60 non-causals above it has
# true rank 61, so exact AP = 1/61. The 0.005-grid version returns 1/N.
pip61 <- c(seq(0.9, 0.01, length.out = 60), 0.0032)
y61   <- c(rep(FALSE, 60), TRUE)
r3 <- .compute_ap_exact(pip61, y61)
ok("exact AP = 1/61 for a causal at rank 61", eq(r3, 1/61),
   sprintf("got %.8f, expected %.8f", r3, 1/61))

# Tie safety: all tied at 0.5, one causal among N. One block, k = N, tp = 1,
# so AP = 1/N regardless of input order.
N <- 25L
r4 <- .compute_ap_exact(rep(0.5, N), c(rep(FALSE, 12), TRUE, rep(FALSE, N - 13)))
ok("AP = 1/N when every PIP is tied", eq(r4, 1/N), sprintf("got %.6f", r4))
r4b <- .compute_ap_exact(rep(0.5, N), c(TRUE, rep(FALSE, N - 1)))
ok("tied AP is independent of the causal's position", eq(r4, r4b))

# §4.4 - the uniform predictor scores perfectly on ECE and mass ratio. RES is
# what exposes it.
p <- 1000L; S <- 3L
unif <- rep(S / p, p); yu <- c(rep(TRUE, S), rep(FALSE, p - S))
mu <- .murphy_decomposition(n_b = p, c_b = S, sum_pip_b = sum(unif))
ok("uniform predictor: REL = 0", eq(mu$rel, 0))
ok("uniform predictor: RES = 0  (the gaming hole)", eq(mu$res, 0))
ok("uniform predictor: BSS = 0", eq(mu$bss, 0))
ok("uniform predictor: mass ratio = 1", eq(sum(unif) / S, 1))

# §4.6 R1 - abstention. Nothing above 0.9 must give fdr = 0 AND nsel = 0, so
# the pair is reported together and the abstainer cannot win on FDR alone.
pa <- c(0.4, 0.3, 0.2); ya <- c(TRUE, FALSE, FALSE)
nsel90 <- sum(pa >= 0.9)
ok("abstainer: nsel_at_90 = 0", nsel90 == 0L)
ok("abstainer: fdr_at_90 = 0 by convention (hence R1)", (if (nsel90 > 0) NA else 0) == 0)

# §5.2b - a posterior whose total mass never reaches alpha.
ps <- rep(0.02, 20L)          # total mass 0.4 < 0.95
rs <- .rank_set_metrics(ps, c(TRUE, rep(FALSE, 19L)))
ok("set_reached_95 = FALSE when total mass is 0.4", identical(rs$set_reached_95, FALSE))
ok("set_size_95 = N when alpha is never reached", rs$set_size_95 == 20L)

# §4.2 P3 - honesty. Selected {0.9, 0.8}, one causal => FP = 1.
# exp_fp = (1-0.9)+(1-0.8) = 0.3 ; rho = 1/0.3 = 3.3333...
h <- .honesty(fp = 1, nsel = 2, sum_pip_sel = 0.9 + 0.8)
ok("exp_fp = 0.3 from the method's own claim", eq(h$exp_fp, 0.3))
ok("honesty ratio = 10/3", eq(h$rho, 10/3), sprintf("got %.6f", h$rho))

# §3.2 both directions: identical counts, different rankings.
# Same PIP multiset => identical pooled FDR/calibration counts; APs differ.
pA <- c(0.9, 0.1); pB <- c(0.9, 0.1)
apA <- .compute_ap_exact(pA, c(TRUE, FALSE))
apB <- .compute_ap_exact(pB, c(FALSE, TRUE))
ok("same counts, different ranking: APs differ", !eq(apA, apB))
ok("  and the counts really are identical", eq(sum(pA), sum(pB)))

# Excess mass by band (§4.3 Q2).
em <- .excess_mass_by_band(bin_lower = c(0, 0.1, 0.5, 0.9),
                           c_b = c(0, 0, 1, 2), sum_pip_b = c(5, 3, 2, 2))
ok("excess mass: lo band = 5 - 0", eq(unname(em["lo"]), 5))
ok("excess mass: top band = 2 - 2 = 0", eq(unname(em["top"]), 0))

# ---------------------------------------------------------------------------
cat("\n=== §7.2 split-plot correction: UNBIASEDNESS over replicates ===\n")
# A single draw cannot validate a bias correction. With sigma^2_u = 0.04 the
# whole-plot lambda_T is 0.5045 against SS_tot ~ 14.75, so each df of correction
# removes ~3.4% of the total and one realisation scatters by +/- 0.12 on a null
# factor. What must hold is that the correction is UNBIASED: averaged over
# replicates, a null factor's corrected share goes to zero while its uncorrected
# share stays substantially positive.
#
# Analytic truth for this construction (100 cells per phi level):
#   SS_phi = 100 * sum((0.5+0.05*(1:5) - 0.65)^2) = 2.50
#   SS_u   = 500 * sigma^2_u/2                    = 10.00
#   SS_eps = 500 * sigma^2_eps/20                 = 2.25
#   => true phi share = 2.50 / 14.75 = 0.169 ; enrich and region_size are NULL.
sigma2_u <- 0.04; sigma2_eps <- 0.09
TRUE_PHI_SHARE <- 2.5 / 14.75
NREP <- 40L
set.seed(11)

grid0 <- expand.grid(S = 1:5, phi = 1:5, region_size = 1:5, enrich = 1:4)
mu0   <- 0.5 + 0.05 * grid0$phi
blk   <- interaction(grid0$enrich, grid0$region_size, drop = TRUE)
vc    <- list(sigma2_eps = sigma2_eps, sigma2_u_bar = sigma2_u)

acc <- replicate(NREP, {
  d <- grid0
  u <- rnorm(nlevels(blk), 0, sqrt(sigma2_u))[as.integer(blk)]
  # A cell mean averages the 2 regions of its size class, so its share of the
  # region draw has variance sigma^2_u / 2.
  d$Y <- mu0 + u / sqrt(2) + rnorm(nrow(d), 0, sqrt(sigma2_eps / 20))
  for (v in c("S", "phi", "region_size", "enrich")) d[[v]] <- factor(d[[v]])
  dec <- decompose(d, "Y", c("S", "phi", "region_size", "enrich"), vc)
  tt  <- dec$table
  g   <- function(t) tt$omega2[tt$term == t]
  r   <- function(t) tt$ss[tt$term == t] / dec$ss_tot
  c(phi_c = g("phi"), phi_r = r("phi"),
    en_c  = g("enrich"),      en_r = r("enrich"),
    rs_c  = g("region_size"), rs_r = r("region_size"))
}, simplify = "matrix")
m <- rowMeans(acc)

cat(sprintf("    mean over %d replicates:\n", NREP))
cat(sprintf("      phi          corrected %+.3f   uncorrected %+.3f   (true %.3f)\n",
            m[["phi_c"]], m[["phi_r"]], TRUE_PHI_SHARE))
cat(sprintf("      enrich  NULL corrected %+.3f   uncorrected %+.3f   (true 0)\n",
            m[["en_c"]], m[["en_r"]]))
cat(sprintf("      regsize NULL corrected %+.3f   uncorrected %+.3f   (true 0)\n",
            m[["rs_c"]], m[["rs_r"]]))

ok("phi share recovers its analytic value 0.169",
   abs(m[["phi_c"]] - TRUE_PHI_SHARE) < 0.03,
   sprintf("got %.3f", m[["phi_c"]]))
ok("NULL factor 'enrich': corrected share is unbiased toward 0",
   abs(m[["en_c"]]) < 0.03, sprintf("mean %.4f", m[["en_c"]]))
ok("NULL factor 'region_size': corrected share is unbiased toward 0",
   abs(m[["rs_c"]]) < 0.03, sprintf("mean %.4f", m[["rs_c"]]))
ok("and the UNCORRECTED null shares are substantially inflated",
   m[["en_r"]] > 0.05 && m[["rs_r"]] > 0.05,
   sprintf("enrich %.3f, region_size %.3f - this is what the correction removes",
           m[["en_r"]], m[["rs_r"]]))

# Structural properties, on one draw.
d1 <- grid0
u1 <- rnorm(nlevels(blk), 0, sqrt(sigma2_u))[as.integer(blk)]
d1$Y <- mu0 + u1 / sqrt(2) + rnorm(nrow(d1), 0, sqrt(sigma2_eps / 20))
for (v in c("S", "phi", "region_size", "enrich")) d1[[v]] <- factor(d1[[v]])
dec1 <- decompose(d1, "Y", c("S", "phi", "region_size", "enrich"), vc)
tt1 <- dec1$table
ok("whole-plot terms get the larger lambda",
   unique(tt1$lambda[tt1$whole_plot]) > unique(tt1$lambda[!tt1$whole_plot]))
ok("any term containing S or phi is NOT whole-plot",
   !any(tt1$whole_plot[grepl("(^|:)(S|phi)(:|$)", tt1$term)]))
ok("shares plus noise floor sum to 1", eq(sum(tt1$omega2) + dec1$noise_floor, 1, 1e-8))

cat("\n=== §8.3 cross-fitted decision gate ===\n")
# The gate must do BOTH jobs: decide when a difference is real, and refuse when
# it is a coin toss. A gate that never decides is as useless as one that always
# does, and only the pair of tests distinguishes the two failures.
mk <- function(mA, mB, sdv, R = 10L, seed = 7) {
  set.seed(seed)
  rbind(data.frame(method = "A", iter = 1:R, ap = rnorm(R, mA, sdv)),
        data.frame(method = "B", iter = 1:R, ap = rnorm(R, mB, sdv)))
}
clear <- decide_cell(mk(0.80, 0.50, 0.02), "ap")
ok("gate DECIDES a clear 0.30 gap at sd 0.02",
   isTRUE(clear$decided) && identical(clear$winner, "A"),
   sprintf("decided=%s winner=%s", clear$decided, clear$winner))
toss <- decide_cell(mk(0.80, 0.79, 0.30), "ap")
ok("gate REFUSES a 0.01 gap swamped by sd 0.30",
   isFALSE(toss$decided), sprintf("decided=%s", toss$decided))
ok("gate still names a nominal winner when undecided",
   !is.na(toss$winner))

cat("\n=== variance components recover their inputs ===\n")
# WITH (S,phi) STRUCTURE - the case the POOLED estimator is biased on. The same
# region-draw difference appears in all K = 25 of a job's (S,phi) pairs, so an
# uncorrected pooled variance understates sigma^2_u by a factor
# c = K(J-1)/(JK-1) = 0.758 at J = 4. Again a SINGLE draw cannot test this:
# var() over 4 jobs is wildly variable, so check the MEAN over replicates.
s2u_t <- 0.04; s2e_t <- 0.09; JOBS <- 4L
grid2 <- expand.grid(job = seq_len(JOBS), S = 1:5, phi = 1:5,
                     region_idx = 1:2, iter = 1:10)
one_rep <- function() {
  uj <- matrix(rnorm(JOBS * 2, 0, sqrt(s2u_t)), nrow = JOBS)
  f <- data.frame(
    job_dir     = paste0("job_", grid2$job, "_anBinary_e5.4_refInsample"),
    region_size = 500L, region_idx = grid2$region_idx,
    S = grid2$S, phi = grid2$phi, iter = grid2$iter,
    Y = 0.5 + uj[cbind(grid2$job, grid2$region_idx)] +
        rnorm(nrow(grid2), 0, sqrt(s2e_t)),
    stringsAsFactors = FALSE)
  # Corrected estimate, and the uncorrected pooled one for comparison.
  vc <- variance_components(f, "Y")
  KEY <- c("job_dir", "region_size", "region_idx", "S", "phi")
  rm2 <- aggregate(list(ybar = f$Y), by = f[KEY], FUN = mean)
  w2  <- reshape(rm2, idvar = c("job_dir", "region_size", "S", "phi"),
                 timevar = "region_idx", direction = "wide")
  naive <- max(0, 0.5 * var(w2$ybar.1 - w2$ybar.2) - vc$sigma2_eps / 10)
  c(corrected = vc$sigma2_u_bar, naive = naive)
}
set.seed(21)
reps2 <- replicate(120, one_rep())
mc <- rowMeans(reps2)
cat(sprintf("    mean over 120 replicates: corrected %.4f  naive %.4f  (true %.4f)\n",
            mc[["corrected"]], mc[["naive"]], s2u_t))
ok("debiased sigma^2_u is unbiased with (S,phi) structure",
   abs(mc[["corrected"]] - s2u_t) / s2u_t < 0.12,
   sprintf("ratio %.3f", mc[["corrected"]] / s2u_t))
ok("the NAIVE pooled estimator is biased LOW, as derived (c = 0.758)",
   mc[["naive"]] / s2u_t < 0.90,
   sprintf("ratio %.3f - under-corrects every whole-plot term", mc[["naive"]] / s2u_t))
ok("job labels containing dots do not break the key",
   all(is.finite(reps2["corrected", ])))

# Fit-level table: 2 regions per (block), 10 iterations each.
set.seed(12)
fits <- do.call(rbind, lapply(1:12, function(j)
  do.call(rbind, lapply(1:2, function(ri) {
    uj <- rnorm(1, 0, sqrt(sigma2_u))
    data.frame(job_dir = paste0("job", j), region_size = 500L, region_idx = ri,
               S = 1L, phi = 0.1, iter = 1:10,
               Y = 0.5 + uj + rnorm(10, 0, sqrt(sigma2_eps)),
               stringsAsFactors = FALSE)
  }))))
vc2 <- variance_components(fits, "Y")
ok("sigma^2_eps recovered within 30%",
   abs(vc2$sigma2_eps - sigma2_eps) / sigma2_eps < 0.30,
   sprintf("got %.4f, true %.4f", vc2$sigma2_eps, sigma2_eps))
ok("sigma^2_u recovered to the right order",
   vc2$sigma2_u_bar > sigma2_u / 4 && vc2$sigma2_u_bar < sigma2_u * 4,
   sprintf("got %.4f, true %.4f", vc2$sigma2_u_bar, sigma2_u))
ok("sigma^2_job = sigma^2_u / 2", eq(vc2$sigma2_job, vc2$sigma2_u_bar / 2))

# >>> ITERATION 005 ONLY - TEMPORARY (see iteration-005-REVERT.md) <<<
cat("\n=== stratification adapts to the design (iter005) ===\n")

# Iteration 004 shape: no `relationship` column -> the standing six strata.
d4 <- data.frame(model = rep(c("sparse", "sparse_inf"), each = 3),
                 annotation_type = rep(c("none", "binary", "continuous"), 2),
                 stringsAsFactors = FALSE)
ok("iter004-shaped data yields the standing 6 strata",
   length(strata_for(d4)) == 6L, sprintf("got %d", length(strata_for(d4))))

# Iteration 005 shape: model constant, relationship varies -> 6 x 2 = 12.
d5 <- expand.grid(relationship = c("additive", "cooccur", "mixed",
                                   "nonmono", "threshold", "null"),
                  annotation_type = c("binary", "continuous"),
                  stringsAsFactors = FALSE)
d5$model <- "sparse"
st5 <- strata_for(d5)
ok("iter005-shaped data yields 12 strata, not 2",
   length(st5) == 12L, sprintf("got %d", length(st5)))

# The failure this guards against: stratifying Iteration 005 on (model,
# annotation_type) puts all twelve rows into two cells, averaging the five
# non-null relationships together. enrichment_fold does not rescue it - it is
# 5.4 on every non-null arm.
n_old <- length(unique(paste(d5$model, d5$annotation_type)))
ok("the OLD keys would have collapsed those 12 rows to 2 cells", n_old == 2L,
   sprintf("got %d", n_old))

# A stratum must select exactly its own rows.
d5$y <- seq_len(nrow(d5))
sub <- stratum_subset(d5, list(key = "nonmono_binary",
                               filter = list(relationship = "nonmono",
                                             annotation_type = "binary")))
ok("a relationship stratum selects exactly its own rows",
   !is.null(sub) && nrow(sub$data) == 1L &&
     sub$data$relationship == "nonmono" && sub$data$annotation_type == "binary")

# THE IMPORTANT ONE. If the constraint column is missing, returning the
# unfiltered frame would hand every sense script a SUPERSET of the stratum -
# silently pooling precisely what the stratification exists to separate. NULL
# (i.e. "absent, skipping") is the only safe answer.
noRel <- d5[, setdiff(names(d5), "relationship"), drop = FALSE]
sub2 <- stratum_subset(noRel, list(key = "nonmono_binary",
                                   filter = list(relationship = "nonmono",
                                                 annotation_type = "binary")))
ok("a filter on an ABSENT column returns NULL, never a superset",
   is.null(sub2), "silently pooling here would invert the iteration's result")

# Back-compat: the standing strata carry model/annot, not a filter list.
ok("legacy (model, annot) strata still subset correctly",
   {ss <- stratum_subset(cbind(d4, y = 1), STRATA[[2]])
    !is.null(ss) && nrow(ss$data) == 1L})
# >>> END ITERATION 005 TEMPORARY BLOCK <<<

cat(sprintf("\n---------------------------------------------\n%d passed, %d failed\n",
            PASS, FAIL))
quit(save = "no", status = if (FAIL > 0L) 1L else 0L)
