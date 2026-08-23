#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/make_test_fixture_iter005.R
#
# >>> ITERATION 005 ONLY - TEMPORARY. See docs/autoresearch/iteration-005-REVERT.md
#
# A synthetic benchmark tree in Iteration 005's shape: 6 relationships x 2
# annotation types = 12 rows, all model = "sparse", one region-size class.
#
# WHY A SEPARATE FIXTURE. Iteration 004's fixture varies (model,
# annotation_type) and so cannot detect the failure this exists to catch: an
# analysis that stratifies on those two columns puts all twelve Iteration 005
# rows into two cells and averages five relationships together. A fixture that
# cannot fail the way production fails is not a test.
#
# The skill model below is deliberately NOT flat across relationships: the
# linear-prior methods are strong on `additive`, lose most of that edge on the
# forms a linear model cannot express, and collapse to the annotation-blind
# baseline on `null`. If the pipeline pools relationships, that structure
# disappears into an average - which is exactly what the assertions check.
#
#   Rscript scripts/analysis/make_test_fixture_iter005.R <out_root>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x)) y else x
root <- args[1] %||% "/tmp/fmb_fixture_iter005"
set.seed(505)

bench <- file.path(root, "results", "benchmark")
dir.create(bench, recursive = TRUE, showWarnings = FALSE)

RELATIONSHIPS <- c("additive", "cooccur", "mixed", "nonmono", "threshold", "null")
ANNOT_TYPES   <- c("binary", "continuous")
S_LV   <- c(1L, 3L)
PHI_LV <- c(0.1, 0.4)
N_ITER <- 5L                       # 25 in production; 5 keeps the fixture fast
P_VEC  <- rep(80L, 4L)             # ONE size class, 4 draws - as Iteration 005
ENRICHMENT    <- 5.4
N_ANNOTATIONS <- 10L
N_INFORMATIVE <- 5L

METHODS <- c("susie", "beatrice",
             "polyfun_oracle", "polyfun_est", "polyfun_ldsc",
             "paintor", "sbayesrc", "funmap",
             "functional_beatrice", "fb_pooled", "fb_xregion")

# Baseline skill, on the additive arm where every prior is correctly specified.
SKILL_ADDITIVE <- c(susie = 0.00, beatrice = 0.05,
                    polyfun_oracle = 0.45, polyfun_est = 0.30, polyfun_ldsc = 0.26,
                    paintor = 0.24, sbayesrc = 0.22, funmap = 0.20,
                    functional_beatrice = 0.24, fb_pooled = 0.26, fb_xregion = 0.30)

# How much of that edge survives each relationship. Annotation-blind methods
# (susie, beatrice) are unaffected by construction - they never see annotations.
LINEAR   <- c("polyfun_oracle", "polyfun_est", "polyfun_ldsc",
              "paintor", "sbayesrc", "funmap", "fb_pooled")
LASSONET <- c("functional_beatrice", "fb_xregion")
RETENTION <- list(
  additive  = c(linear = 1.00, lassonet = 1.00),
  threshold = c(linear = 0.70, lassonet = 0.85),
  mixed     = c(linear = 0.55, lassonet = 0.90),
  cooccur   = c(linear = 0.35, lassonet = 0.80),
  nonmono   = c(linear = 0.20, lassonet = 0.75),
  null      = c(linear = 0.00, lassonet = 0.00)
)
# polyfun_oracle reads the STORED causal probabilities, so it is a true oracle
# under every relationship and keeps its edge throughout (except on null, where
# there is no signal in the annotations to have an edge from).
skill_for <- function(m, rel) {
  base <- SKILL_ADDITIVE[[m]]
  if (m %in% c("susie", "beatrice")) return(base)
  if (m == "polyfun_oracle") return(if (rel == "null") 0.05 else base)
  r <- RETENTION[[rel]]
  base * (if (m %in% LASSONET) r[["lassonet"]] else r[["linear"]])
}

rows <- list(); k <- 0L
for (rel in RELATIONSHIPS) for (at in ANNOT_TYPES) {
  k <- k + 1L
  rows[[k]] <- data.frame(
    job_id = k,
    label  = sprintf("rel%s_an%s", rel, if (at == "binary") "Binary" else "Cont"),
    model = "sparse", p_causal = NA_real_, annotation_type = at,
    enrichment_fold = if (rel == "null") NA_real_ else ENRICHMENT,
    enrichment_values = paste(
      if (rel == "null") rep(1, N_ANNOTATIONS)
      else c(rep(ENRICHMENT, N_INFORMATIVE),
             rep(1, N_ANNOTATIONS - N_INFORMATIVE)), collapse = "|"),
    relationship = rel, n_informative = N_INFORMATIVE,
    annotation_correlation = 0, n_ref = NA_integer_,
    n_annotations = N_ANNOTATIONS, n_regions = length(P_VEC),
    n = 5000L, n_iter = N_ITER,
    S_values   = paste(S_LV,   collapse = "|"),
    phi_values = paste(PHI_LV, collapse = "|"),
    p_values   = paste(P_VEC,  collapse = "|"),
    stringsAsFactors = FALSE)
}
grid <- do.call(rbind, rows)
grid_csv <- file.path(root, "params_grid.csv")
write.csv(grid, grid_csv, row.names = FALSE)

for (i in seq_len(nrow(grid))) {
  r   <- grid[i, ]
  rel <- r$relationship
  jd  <- file.path(bench, sprintf("job_%03d_%s", r$job_id, r$label))
  dir.create(jd, recursive = TRUE, showWarnings = FALSE)

  pg <- expand.grid(S = S_LV, phi = PHI_LV, iter = seq_len(N_ITER))
  scenarios <- vector("list", nrow(pg))
  region_effect <- rnorm(length(P_VEC), 0, 0.05)   # shared draw -> sigma^2_u

  for (sc in seq_len(nrow(pg))) {
    regs <- lapply(seq_along(P_VEC), function(rg) {
      p <- P_VEC[rg]
      list(truth = list(causal_indices = sort(sample.int(p, pg$S[sc])),
                        S = pg$S[sc], phi = pg$phi[sc], model = "sparse"))
    })
    scenarios[[sc]] <- list(scenario_id = sc, S = pg$S[sc], phi = pg$phi[sc],
                            p_causal = NA_real_, iter = pg$iter[sc],
                            model = "sparse", regions = regs)
  }
  saveRDS(list(scenarios = scenarios,
               genotypes = lapply(P_VEC, function(p) list(maf = runif(p, 0.05, 0.5))),
               params = list(p = P_VEC, n = 5000L)),
          file.path(jd, "sim.rds"))

  for (sc in seq_len(nrow(pg))) {
    sd <- file.path(jd, sprintf("scenario_%03d", sc))
    dir.create(sd, recursive = TRUE, showWarnings = FALSE)
    res <- list()
    for (m in METHODS) {
      sk <- skill_for(m, rel)
      fits <- lapply(seq_along(P_VEC), function(rg) {
        p  <- P_VEC[rg]
        ci <- scenarios[[sc]]$regions[[rg]]$truth$causal_indices
        sig <- (0.35 + sk) * pg$phi[sc] * 6 / (1 + 0.35 * pg$S[sc]) +
               region_effect[rg]
        z <- rnorm(p, 0, 1); z[ci] <- z[ci] + max(sig, 0.05) * 6
        pip <- exp(z) / sum(exp(z))
        list(pip = pip, scenario_id = 1L, region_id = rg, method = m)
      })
      res[[m]] <- list(results = fits, n_total = length(P_VEC), n_failed = 0L,
                       method_args = list(), total_runtime_seconds = 1)
    }
    res$methods_run       <- METHODS
    res$run_timestamp     <- "2026-08-23 00:00:00"
    res$simulation_params <- list(n = 5000L, p = P_VEC)
    saveRDS(res, file.path(sd, "results.rds"))
    saveRDS(list(), file.path(sd, "evaluation.rds"))
  }
}

cat(sprintf("iter005 fixture: %d rows x %d scenarios x %d regions x %d methods\n",
            nrow(grid), N_ITER * length(S_LV) * length(PHI_LV),
            length(P_VEC), length(METHODS)))
cat("bench root: ", bench, "\n", sep = "")
cat("grid csv  : ", grid_csv, "\n", sep = "")
