#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/inspect_deposit.R
#
# Look at what is in the staged deposit, without having to remember how any of
# it is shaped. Prints a readable report: what each file holds, what the key
# columns contain, and a few summaries that should look sensible if the data is
# right.
#
#   Rscript scripts/analysis/inspect_deposit.R [dir]
#
# dir defaults to zenodo/main-simulation.
#
# For Python, the two CSVs need no R at all:
#   import pandas as pd
#   fits = pd.read_csv("fit-metrics.csv.gz")        # reads the .gz directly
#   grid = pd.read_csv("simulation-grid.csv")
# =============================================================================

D <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "zenodo/main-simulation"
stopifnot(dir.exists(D))
rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")

rule("WHAT IS IN THE FOLDER")
f <- list.files(D, recursive = TRUE, full.names = TRUE)
sz <- file.info(f)$size
top <- sub(paste0("^", D, "/"), "", f)
grp <- ifelse(grepl("/", top), sub("/.*", "/", top), top)
agg <- tapply(sz, grp, sum); cnt <- table(grp)
for (g in names(sort(agg, decreasing = TRUE)))
  cat(sprintf("  %-36s %8.1f MB   %s\n", g, agg[[g]] / 1024^2,
              if (cnt[[g]] > 1) sprintf("%d files", cnt[[g]]) else ""))

# --- the design --------------------------------------------------------------
rule("simulation-grid.csv  -  the 45 parameter rows")
g <- read.csv(file.path(D, "simulation-grid.csv"), stringsAsFactors = FALSE)
cat("columns:", paste(names(g), collapse = ", "), "\n\n")
print(head(g[, c("file", "model", "p_causal", "annotation_type", "enrichment_fold")], 6),
      row.names = FALSE)
cat("  ... and", nrow(g) - 6, "more\n")

# --- the fit table -----------------------------------------------------------
rule("fit-metrics.rds  -  one row per fit")
L1 <- readRDS(file.path(D, "fit-metrics.rds"))
cat(sprintf("%s rows x %d columns\n\n", format(nrow(L1), big.mark = ","), ncol(L1)))

KEY <- c("model","p_causal","annotation_type","enrichment_fold","S","phi","iter",
         "region_id","method")
cat("the key, and what each column takes:\n")
for (k in KEY) {
  u <- unique(L1[[k]])
  cat(sprintf("  %-18s %3d values   %s\n", k, length(u),
              paste(head(sort(u, na.last = TRUE), 7), collapse = ", ")))
}
cat("\nis the key unique over every row: ",
    !anyDuplicated(do.call(paste, c(L1[KEY], sep = "\r"))), "\n", sep = "")

cat("\nall columns:\n")
print(matrix(c(names(L1), rep("", (3 - ncol(L1) %% 3) %% 3)), ncol = 3, byrow = TRUE),
      quote = FALSE)

cat("\none fit, in full:\n")
r <- L1[which(L1$method == "susie" & L1$S == 3 & !L1$failed)[1], ]
for (n in names(r)) cat(sprintf("  %-22s %s\n", n, format(r[[n]])))

# --- does it look right ------------------------------------------------------
rule("SUMMARIES THAT SHOULD LOOK SENSIBLE")
ok <- L1[!L1$failed & !L1$bad_pip, ]
cat("average precision by method, sparse model with binary annotations:\n")
s <- ok[ok$model == "sparse" & ok$annotation_type == "binary", ]
ap <- sort(tapply(s$ap, s$method, mean), decreasing = TRUE)
for (m in names(ap)) cat(sprintf("  %-22s %.3f\n", m, ap[[m]]))

cat("\nfailed or unusable fits, by method (only those with any):\n")
fr <- tapply(L1$failed | L1$bad_pip, L1$method, mean)
fr <- sort(fr[fr > 0], decreasing = TRUE)
for (m in names(fr)) cat(sprintf("  %-22s %5.2f%%\n", m, 100 * fr[[m]]))

cat("\nmethods reporting no credible set of their own:\n")
ns <- tapply(L1$native_n_sets, L1$method, function(z) all(z == 0, na.rm = TRUE))
cat(" ", paste(names(ns)[which(ns)], collapse = ", "), "\n")

# --- the long tables ---------------------------------------------------------
row1 <- g$file[1]
rule(sprintf("posterior-inclusion-probabilities/%s.rds", row1))
p <- readRDS(file.path(D, "posterior-inclusion-probabilities", paste0(row1, ".rds")))
cat(sprintf("%s rows x %d columns: %s\n\n", format(nrow(p), big.mark = ","), ncol(p),
            paste(names(p), collapse = ", ")))
print(head(p, 5), row.names = FALSE)
cat(sprintf("\npip: %.4f to %.4f, with %s below the 0.01 floor (all causal), %s NA\n",
            min(p$pip, na.rm = TRUE), max(p$pip, na.rm = TRUE),
            format(sum(p$pip < 0.01, na.rm = TRUE), big.mark = ","),
            format(sum(is.na(p$pip)), big.mark = ",")))

rule(sprintf("credible-sets/%s.rds", row1))
cs <- readRDS(file.path(D, "credible-sets", paste0(row1, ".rds")))
cat(sprintf("%s rows x %d columns: %s\n\n", format(nrow(cs), big.mark = ","), ncol(cs),
            paste(names(cs), collapse = ", ")))
print(head(cs, 5), row.names = FALSE)
cat("\nsets per fit, by method (mean):\n")
k <- do.call(paste, c(cs[c("S","phi","iter","region_id","method")], sep = "\r"))
nset <- tapply(cs$set, k, function(z) length(unique(z)))
mth  <- sub(".*\r", "", names(nset))
print(round(sort(tapply(nset, mth, mean), decreasing = TRUE), 2))

rule("scenario-runtimes.rds")
sr <- readRDS(file.path(D, "scenario-runtimes.rds"))
cat(sprintf("%s rows x %d columns: %s\n\n", format(nrow(sr), big.mark = ","), ncol(sr),
            paste(names(sr), collapse = ", ")))
cat("median seconds per scenario, the ten costliest methods:\n")
md <- sort(tapply(sr$runtime_seconds, sr$method, median, na.rm = TRUE), decreasing = TRUE)
for (m in head(names(md), 10)) cat(sprintf("  %-22s %8.1f\n", m, md[[m]]))

rule("DONE")
cat("Nothing above should surprise you. If it does, that is the point of it.\n")
