#!/usr/bin/env python
"""Build the annotation files PolyFun regresses on: base + quintile indicators.

WHY this exists. The Borzoi PCs are signed and continuous (chr22 range
[-841.8, 4471], all 102 columns carry negatives), and neither the Borzoi set
nor baselineLF ships a `base` column. That left the S-LDSC regressions with

  * no level term, so nothing absorbs the genome-wide mean per-SNP h2, and
  * negative column sums, i.e. a negative M -- which S-LDSC reads as a
    per-annotation SNP count (M_tot = sum(M) sets the regression weights via
    LD_Score_Regression._update_func, and the per-annotation M_j feed the h2
    and enrichment outputs).

So each PC is replaced by quintile indicators and an all-ones `base` is
added, which is the standard S-LDSC treatment of a continuous annotation and
exactly what the published baselineLD/baselineLF files contain (`base` there
is literally all ones -- it is a definition, not data to be sourced).

THE REFERENCE CATEGORY. All five quintile indicators of one PC are mutually
exclusive and sum to the all-ones vector, so including all five AND base is
exactly collinear and the regression is singular. The lowest quintile is
dropped as the reference category, leaving 4 per PC:

    sniff     102*4 + base                     = 409 annotations
    combined  102*4 + base + 187 baselineLF    = 596 annotations

baselineLF's own columns are kept exactly as published, including the 5 with
negative values -- those are its standardised continuous annotations, and
re-binning them would depart from the published baselineLF model. With `base`
present, M_tot is dominated by a genuine SNP count, so a few negative M_j only
affect those annotations' own reported h2, which this pipeline does not use.

BREAKPOINTS ARE GENOME-WIDE. S-LDSC is one genome-wide regression, so a
per-chromosome quintile would make the same annotation name mean a different
thing on each chromosome. Pass 1 pools each PC across all 22 chromosomes and
writes the 20/40/60/80th percentiles; pass 2 applies them. Pass 1 reads ONE PC
column at a time across the 22 files, so it holds ~70 MB rather than the
~7 GB the full matrix would take.

Quintiles are defined on the SNPs present in the annotation files. The
LD-scores are later computed on panel-and-annotation SNPs (98.5% of the
panel), so the two sets differ slightly; defining the breakpoints on the
annotation files keeps them independent of which panel is used.

  # pass 1, once
  python build_polyfun_annot.py breakpoints --annot-dir annot --out pc_quintiles.tsv
  # pass 2, per chromosome
  python build_polyfun_annot.py build --set sniff --chr 22 --breakpoints pc_quintiles.tsv
"""
import os, sys, argparse
import numpy as np
import pyarrow.parquet as pq
import pyarrow as pa

KEY = ["CHR", "SNP", "BP", "A1", "A2"]
QS = [0.2, 0.4, 0.6, 0.8]


def pc_columns(path):
    return [c for c in pq.ParquetFile(path).schema_arrow.names if "-PC" in c]


def cmd_breakpoints(a):
    files = [os.path.join(a.annot_dir, f"borzoi.{c}.annot.parquet") for c in range(1, 23)]
    for f in files:
        if not os.path.exists(f):
            sys.exit(f"FATAL missing {f}")
    pcs = pc_columns(files[0])
    for f in files[1:]:
        if pc_columns(f) != pcs:
            sys.exit(f"FATAL PC columns differ between {files[0]} and {f}")
    print(f"{len(pcs)} PCs across 22 chromosomes", flush=True)

    rows = []
    for i, pc in enumerate(pcs):
        # one column at a time: ~70 MB pooled, not the ~7 GB full matrix
        vals = np.concatenate([
            pq.read_table(f, columns=[pc]).column(0).to_numpy(zero_copy_only=False)
            for f in files])
        q = np.quantile(vals.astype(np.float64), QS)
        if len(np.unique(q)) < len(QS):
            print(f"  WARNING {pc}: tied breakpoints {q} - some bins will be empty")
        rows.append((pc, len(vals), *q))
        if (i + 1) % 20 == 0 or i + 1 == len(pcs):
            print(f"  {i+1}/{len(pcs)} PCs", flush=True)

    tmp = a.out + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("PC\tN\tq20\tq40\tq60\tq80\n")
        for pc, n, *q in rows:
            fh.write(pc + "\t" + str(n) + "\t" + "\t".join(repr(float(x)) for x in q) + "\n")
    os.replace(tmp, a.out)
    print(f"wrote {a.out}: {len(rows)} PCs, {rows[0][1]:,} SNPs pooled per PC")


def read_breakpoints(path):
    bp = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        assert hdr[:2] == ["PC", "N"], f"unexpected header: {hdr}"
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            bp[p[0]] = np.array([float(x) for x in p[2:6]])
    return bp


def cmd_build(a):
    bp = read_breakpoints(a.breakpoints)
    if a.set == "sniff":
        src = os.path.join(a.annot_dir, f"borzoi.{a.chr}.annot.parquet")
    else:
        src = os.path.join(a.combined_dir, f"combined.{a.chr}.annot.parquet")
    if not os.path.exists(src):
        sys.exit(f"FATAL missing {src}")
    out = os.path.join(a.out_dir, f"{a.prefix}.{a.chr}.annot.parquet")
    if os.path.exists(out) and not a.force and not os.path.islink(out):
        n = len(pq.ParquetFile(out).schema_arrow.names)
        print(f"chr{a.chr} ({a.set}): already built, {n} columns")
        return
    if os.path.islink(out):
        os.unlink(out)          # an older symlink to the raw annotations

    t = pq.read_table(src)
    names = list(t.schema.names)
    pcs = [c for c in names if "-PC" in c]
    extra = [c for c in names if c not in KEY and "-PC" not in c]
    missing = [p for p in pcs if p not in bp]
    if missing:
        sys.exit(f"FATAL no breakpoints for {missing[:5]}")

    arrays, out_names = [], []
    for k in KEY:
        arrays.append(t.column(k)); out_names.append(k)
    n = t.num_rows
    arrays.append(pa.array(np.ones(n, dtype=np.float32))); out_names.append("base")

    # float32 throughout: compute_ldscores takes df.values, and a homogeneous
    # dtype keeps that a float32 view rather than a float64 copy (which at 596
    # annotations x 1.65M rows is 7.9 GB instead of 3.9 GB).
    empty = []
    for pc in pcs:
        v = t.column(pc).to_numpy(zero_copy_only=False).astype(np.float64)
        q = bp[pc]
        # bin index 0..4 by the genome-wide quintile breakpoints; bin 0 is the
        # dropped reference category
        b = np.searchsorted(q, v, side="left")
        for k in range(1, 5):
            col = (b == k).astype(np.float32)
            if col.sum() == 0:
                empty.append(f"{pc}_q{k+1}")
            arrays.append(pa.array(col)); out_names.append(f"{pc}_q{k+1}")
    for c in extra:
        arrays.append(pa.array(t.column(c).to_numpy(zero_copy_only=False).astype(np.float32)))
        out_names.append(c)

    os.makedirs(a.out_dir, exist_ok=True)
    tmp = out + ".tmp"
    pq.write_table(pa.table(arrays, names=out_names), tmp)
    os.replace(tmp, out)
    n_ann = len(out_names) - len(KEY)
    print(f"chr{a.chr} ({a.set}): {len(pcs)} PCs -> {len(pcs)*4} indicators + base"
          + (f" + {len(extra)} baselineLF" if extra else "")
          + f" = {n_ann} annotations over {n:,} SNPs")
    if empty:
        print(f"  WARNING {len(empty)} all-zero indicator(s) on this chromosome, "
              f"e.g. {empty[:3]} - PolyFun's _check_variance will drop them")


ap = argparse.ArgumentParser()
sub = ap.add_subparsers(dest="cmd", required=True)
p1 = sub.add_parser("breakpoints")
p1.add_argument("--annot-dir", default="annot")
p1.add_argument("--out", default="pc_quintiles.tsv")
p2 = sub.add_parser("build")
p2.add_argument("--set", required=True, choices=["sniff", "combined"])
p2.add_argument("--chr", required=True, type=int)
p2.add_argument("--breakpoints", default="pc_quintiles.tsv")
p2.add_argument("--annot-dir", default="annot")
p2.add_argument("--combined-dir", default="annot_combined")
p2.add_argument("--out-dir", default=None)
p2.add_argument("--prefix", default="annot")
p2.add_argument("--force", action="store_true")
a = ap.parse_args()
if a.cmd == "breakpoints":
    cmd_breakpoints(a)
else:
    if a.out_dir is None:
        a.out_dir = f"ldscores_{a.set}"
    cmd_build(a)
