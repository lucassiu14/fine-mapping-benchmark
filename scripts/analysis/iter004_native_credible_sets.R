#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_native_credible_sets.R
#
# Score each method's OWN credible sets, for the comparison of set
# constructions. The uniform reconstruction used in the results section cuts
# the ranked PIP vector at 0.95 cumulative mass and yields exactly one set per
# fit; this scores what each method actually reported instead.
#
#   Rscript scripts/analysis/iter004_native_credible_sets.R \
#       [l3_rds] [credsets_dir] [piptail_dir] [out_rds]
#
# The four paths were hard-coded to results/iter004, which made the script
# unrunnable on any other iteration. They are arguments now, defaulting to
# Iteration 004 so existing invocations are unchanged.
#
# Truth comes from the PIP tails rather than sim.rds: extract_pip_tail.R retains
# every causal variant whatever its PIP, so the causal index set is complete
# there even though the PIP vector is truncated.
#
# Coverage is per REPORTED SET - a fit reporting three sets contributes three
# observations - because that is the quantity a method's own convention
# defines. Power is per fit: the causal variants captured by any of its sets,
# over the causal variants present.
# =============================================================================

# Coverage and power under each method's OWN credible sets, using causal indices
# recovered from the PIP tails (which retain every causal variant whatever its
# PIP, so the truth is complete there).
.a   <- commandArgs(TRUE)
L3F  <- if (length(.a) >= 1 && nzchar(.a[1])) .a[1] else "results/iter004/combined_scenario_metrics.rds"
CSD  <- if (length(.a) >= 2 && nzchar(.a[2])) .a[2] else "results/iter004/credsets"
PTD  <- if (length(.a) >= 3 && nzchar(.a[3])) .a[3] else "results/iter004/piptail"
OUT  <- if (length(.a) >= 4 && nzchar(.a[4])) .a[4] else "results/iter004/native_credible_sets.rds"
stopifnot(file.exists(L3F), dir.exists(CSD), dir.exists(PTD))

design <- unique(readRDS(L3F)[, c("job_dir","model","annotation_type")])
labs <- sub("^cs_", "", sub("\\.rds$", "", basename(
  list.files(CSD, pattern="^cs_.*\\.rds$"))))
stopifnot(length(labs) > 0L)
acc <- new.env(parent = emptyenv())
bump <- function(k, v) { o <- acc[[k]]; acc[[k]] <- if (is.null(o)) v else o + v }

for (lab in labs) {
  cs <- readRDS(file.path(CSD, paste0("cs_", lab, ".rds")))
  pt <- tryCatch(readRDS(file.path(PTD,
                                   paste0("piptail_", lab, ".rds"))), error=function(e) NULL)
  if (is.null(pt)) { message("no piptail for ", lab); next }
  truth <- new.env(parent = emptyenv())
  for (f in pt$tail)
    assign(paste(f$scenario_id, f$region_id, sep="\r"),
           f$idx[f$is_causal == 1L], envir = truth)
  d <- design[match(cs$job_dir, design$job_dir), ]
  key0 <- paste(d$model[1], d$annotation_type[1], sep="\r")
  for (r in cs$rows) {
    ci <- mget(paste(r$scenario_id, r$region_id, sep="\r"),
               envir = truth, ifnotfound = list(NULL))[[1]]
    if (is.null(ci) || !length(ci)) next
    k <- paste(key0, r$method, sep="\r")
    if (!isTRUE(r$reported) || !length(r$sets)) {
      bump(k, c(fits=1, nosets=1, sets=0, hits=0, size=0, cap=0, S=length(ci)))
      next
    }
    hits <- sum(vapply(r$sets, function(s) any(s %in% ci), TRUE))
    cap  <- length(unique(unlist(r$sets)[unlist(r$sets) %in% ci]))
    bump(k, c(fits=1, nosets=0, sets=length(r$sets), hits=hits,
              size=sum(vapply(r$sets, length, 1L)), cap=cap, S=length(ci)))
  }
}
keys <- ls(acc); p <- do.call(rbind, strsplit(keys, "\r", fixed=TRUE))
v <- do.call(rbind, lapply(keys, function(k) acc[[k]]))
out <- data.frame(model=p[,1], annotation_type=p[,2], method=p[,3], v,
                  stringsAsFactors=FALSE)
out$coverage_per_set <- out$hits / out$sets           # per reported set
out$power            <- out$cap  / out$S              # causal captured / causal present
out$mean_size        <- out$size / out$sets
out$sets_per_fit     <- out$sets / out$fits
out$pct_no_set       <- 100 * out$nosets / out$fits
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
saveRDS(out, OUT)
cat("wrote ", OUT, "  (", nrow(out), " rows )\n", sep="")
