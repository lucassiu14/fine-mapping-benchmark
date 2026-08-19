#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/qc_gate.R
#
# MACHINE-READABLE go/no-go on a finished benchmark run.
#
# qc_run.R is for a human to read. This exists so a script can DECIDE, because
# the Iteration 004 failure was not that the numbers were unavailable - it was
# that nobody looked at them until after the analysis had been built on top.
#
# Reports per-method failure rates split by model arm, then emits
#
#   GATE_WORST <method> <arm> <pct>
#   GATE_MAX_FAIL_PCT <pct>
#
# which the caller thresholds. Exits 0 regardless: the shell decides policy,
# this script only measures.
#
# STRUCTURAL ABSENCES ARE EXCLUDED, because they are correct behaviour and would
# otherwise pin the maximum at 100% forever:
#   - funmap on the `none` arm cannot run at all; it requires annotations.
# Silent all-NA output is NOT excluded and is reported separately - a method
# that returns NA without setting an error is the exact failure this project
# has been bitten by twice (CARMA for three iterations, finemap_inf at phi=0.6).
#
#   Rscript scripts/hpc/qc_gate.R <output_root> [--sample N]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
if (is.na(root) || !nzchar(root) || !dir.exists(root))
  stop("usage: Rscript scripts/hpc/qc_gate.R <output_root> [--sample N]", call. = FALSE)
sample_n <- 200L
if ("--sample" %in% args)
  sample_n <- suppressWarnings(as.integer(args[which(args == "--sample") + 1L]))
if (is.na(sample_n) || sample_n < 1L) sample_n <- 200L

jobs <- list.files(root, pattern = "^job_", full.names = TRUE)
if (!length(jobs)) stop("no job_* directories under ", root, call. = FALSE)

# Sample scenarios spread ACROSS rows rather than taking the first N, so a
# systematic failure confined to one arm cannot hide behind healthy neighbours.
set.seed(1L)
picks <- unlist(lapply(jobs, function(j) {
  sd <- list.files(j, pattern = "^scenario_", full.names = TRUE)
  if (!length(sd)) return(character(0))
  k <- max(1L, ceiling(sample_n / length(jobs)))
  sd[unique(round(seq(1, length(sd), length.out = min(k, length(sd)))))]
}))

tot <- new.env(parent = emptyenv())
add <- function(key, field, v) {
  cur <- if (exists(key, tot)) get(key, tot) else
           c(fits = 0, failed = 0, nafits = 0, secs = 0, nscen = 0)
  cur[[field]] <- cur[[field]] + v
  assign(key, cur, tot)
}

n_read <- 0L
for (sd in picks) {
  f <- file.path(sd, "results.rds")
  if (!file.exists(f)) next
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(r)) next
  n_read <- n_read + 1L
  jd  <- basename(dirname(sd))
  arm <- if (grepl("sparse_inf", jd)) "sparse_inf" else "sparse"
  ann <- if (grepl("anNone", jd)) "none" else if (grepl("anBinary", jd)) "binary" else "cont"
  for (m in names(r)) {
    v <- r[[m]]
    if (!is.list(v) || is.null(v$results)) next
    if (m == "funmap" && ann == "none") next          # structural, see header
    key <- paste(m, arm)
    nt <- if (is.null(v$n_total))  length(v$results) else v$n_total
    nf <- if (is.null(v$n_failed)) 0L                 else v$n_failed
    # Silent all-NA: returned a pip vector, declared no error, but nothing usable.
    na <- sum(vapply(v$results, function(x)
      isTRUE(!is.null(x$pip) && is.null(x$error) && all(!is.finite(x$pip))),
      logical(1)))
    add(key, "fits", nt); add(key, "failed", nf); add(key, "nafits", na)
    # Runtime is recorded per (method, scenario) by run_methods() as
    # total_runtime_seconds, covering all regions of that scenario. It exists
    # ONLY inside results.rds - no aggregate table carries it - which is why
    # the CARMA cost that justified dropping it survives merely as a code
    # comment. Accumulate it here so the number is recoverable from a normal QC
    # pass instead of requiring a bespoke dig through the raw output.
    rt <- if (is.null(v$total_runtime_seconds)) NA_real_ else
            suppressWarnings(as.numeric(v$total_runtime_seconds))
    if (length(rt) == 1L && is.finite(rt)) {
      add(key, "secs", rt); add(key, "nscen", 1)
    }
  }
}

keys <- sort(ls(tot))
cat(sprintf("scenarios read: %d\n\n", n_read))
cat(sprintf("%-22s %-11s %8s %8s %9s %9s\n",
            "method", "arm", "fits", "failed", "%failed", "silentNA"))
worst <- list(m = "none", arm = "-", pct = 0)
for (k in keys) {
  v <- get(k, tot); p <- strsplit(k, " ")[[1]]
  pct <- if (v[["fits"]] > 0) 100 * v[["failed"]] / v[["fits"]] else 0
  cat(sprintf("%-22s %-11s %8d %8d %8.1f%% %9d\n",
              p[1], p[2], v[["fits"]], v[["failed"]], pct, v[["nafits"]]))
  if (pct > worst$pct) worst <- list(m = p[1], arm = p[2], pct = pct)
}

# SILENT failure needs its own gate. %failed counts only fits that DECLARED an
# error; a method returning an all-NA posterior without setting one scores 0.0%
# and would sail through the gate above while contributing nothing. That is the
# exact shape of both failures this project has suffered - CARMA returned NA for
# three iterations unnoticed, and finemap_inf's phi=0.6 output passed as success.
#
# finemap_inf is exempt BY DEFAULT because its phi=0.6 all-NA output is an
# established, documented property of the method rather than a broken install
# (see docs/autoresearch/iteration-004.md). Exempting it keeps the gate usable;
# leaving it in would make every future run trip the alarm for a known cause,
# which is how alarms get ignored. Override with --silent-exempt "".
exempt <- c("finemap_inf")
if ("--silent-exempt" %in% args) {
  e <- args[which(args == "--silent-exempt") + 1L]
  exempt <- if (is.na(e) || !nzchar(e)) character(0) else strsplit(e, ",")[[1]]
}
ws <- list(m = "none", arm = "-", pct = 0)
for (k in keys) {
  v <- get(k, tot); p <- strsplit(k, " ")[[1]]
  if (p[1] %in% exempt) next
  pct <- if (v[["fits"]] > 0) 100 * v[["nafits"]] / v[["fits"]] else 0
  if (pct > ws$pct) ws <- list(m = p[1], arm = p[2], pct = pct)
}

# ---- runtime -------------------------------------------------------------
rt <- do.call(rbind, lapply(keys, function(k) {
  v <- get(k, tot); p <- strsplit(k, " ")[[1]]
  if (v[["nscen"]] == 0) return(NULL)
  data.frame(method = p[1], arm = p[2], secs = v[["secs"]],
             nscen = v[["nscen"]], stringsAsFactors = FALSE)
}))
if (!is.null(rt) && nrow(rt)) {
  agg <- tapply(rt$secs, rt$method, sum)
  nsc <- tapply(rt$nscen, rt$method, sum)
  per <- agg / pmax(nsc, 1)
  o <- order(-per)
  cat("\nRUNTIME per method (seconds per scenario = all regions of one scenario)\n\n")
  cat(sprintf("  %-22s %12s %10s %9s\n", "method", "s/scenario", "share", "scenarios"))
  tot_per <- sum(per)
  for (i in o)
    cat(sprintf("  %-22s %12.1f %9.1f%% %9d\n", names(per)[i], per[i],
                100 * per[i] / max(tot_per, 1e-9), nsc[i]))
  cat(sprintf("\n  %-22s %12.1f\n", "ALL METHODS", tot_per))
}

cat("\n")
if (length(exempt))
  cat(sprintf("silent-NA gate exempts: %s\n", paste(exempt, collapse = ", ")))
cat(sprintf("GATE_WORST %s %s %.1f\n", worst$m, worst$arm, worst$pct))
cat(sprintf("GATE_MAX_FAIL_PCT %.1f\n", worst$pct))
cat(sprintf("GATE_WORST_SILENT %s %s %.1f\n", ws$m, ws$arm, ws$pct))
cat(sprintf("GATE_MAX_SILENT_PCT %.1f\n", ws$pct))
