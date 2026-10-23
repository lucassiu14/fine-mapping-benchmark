#!/usr/bin/env python
"""Write the .l2.M files PolyFun requires beside each LD-score file.

compute_ldscores.py writes ONLY the LD-score parquet (`df_ldscores.to_parquet`),
but PolyFun also needs a .l2.M file per chromosome under the --ref-ld-chr
prefix. The authors' own example_data/ carries all three per chromosome:

    annotations.<chr>.annot.parquet
    annotations.<chr>.l2.ldscore.parquet
    annotations.<chr>.l2.M

and ldsc_polyfun/sumstats._read_M() reaches it via parse.M(), which does a
bare `open(prefix + '.l2.M').readline()` -- no parquet fallback, and no option
to derive it. polyfun.py's check_files() does NOT verify it, so a missing M
surfaces late, inside _read_M.

DEFINITION, taken from the authors' ldsc.py rather than invented:

    M      = np.sum(annot_matrix, axis=0)                  -> .l2.M
    M_5_50 = np.sum(annot_matrix[maf > 0.05, :], axis=0)   -> .l2.M_5_50

written as one tab-separated line, one value per annotation, in the same
column order as the LD-score file.

WHICH SNPs annot_matrix covers, traced through compute_ldscores.py: the
annotations are filtered to the panel, --allow-missing drops panel SNPs with
no annotation (keep_snps), PlinkBEDFile may drop more (geno_array.kept_snps),
then `df_annot = df_annot.loc[df_bim.index]`. Those surviving SNPs are exactly
the rows written to the LD-score parquet. So annot_matrix is recovered by
indexing the annotation file to the LD-score file's own rows -- no LD-score
recomputation needed, and no reimplementation of the filtering.

set_snpid_index is imported from PolyFun rather than reimplemented, so the
CHR.BP.A1s.A2s key (alleles alphabetised, indels forced A1-first) is theirs.

.l2.M_5_50 is deliberately NOT written: polyfun.py sets not_M_5_50 = True in
run_ldsc(), so _read_M always resolves '.l2.M' and never reads M_5_50.
Producing it would require per-SNP MAF from the plink panel, and writing a
guess at a file PolyFun never reads would be worse than omitting it.

  python make_ldscore_M.py --ldscores-dir ldscores_sniff --annot-prefix annot/borzoi
  python make_ldscore_M.py --ldscores-dir ldscores_combined \
                           --annot-prefix annot_combined/combined
"""
import os, sys, argparse
import numpy as np
import pyarrow.parquet as pq

ap = argparse.ArgumentParser()
ap.add_argument("--ldscores-dir", required=True)
ap.add_argument("--annot-prefix", required=True,
                help="path prefix; <prefix>.<chr>.annot.parquet")
ap.add_argument("--ldscore-prefix", default="annot",
                help="basename prefix inside --ldscores-dir (default: annot)")
ap.add_argument("--polyfun-dir", default=os.environ.get("POLYFUN_DIR", ""))
ap.add_argument("--force", action="store_true", help="rewrite existing .l2.M")
a = ap.parse_args()

if not a.polyfun_dir:
    sys.exit("set POLYFUN_DIR or pass --polyfun-dir")
sys.path.insert(0, a.polyfun_dir)
from polyfun_utils import set_snpid_index, SNP_COLUMNS   # noqa: E402

wrote = skipped = 0
base_seen = None
mins, maxs = [], []
for c in range(1, 23):
    ldf = os.path.join(a.ldscores_dir, f"{a.ldscore_prefix}.{c}.l2.ldscore.parquet")
    anf = f"{a.annot_prefix}.{c}.annot.parquet"
    out = os.path.join(a.ldscores_dir, f"{a.ldscore_prefix}.{c}.l2.M")
    for f in (ldf, anf):
        if not os.path.exists(f):
            sys.exit(f"FATAL missing input: {f}")

    ld = pq.read_table(ldf).to_pandas()
    acols = [x for x in ld.columns if x not in SNP_COLUMNS]   # order matters
    if os.path.exists(out) and not a.force:
        n = len(open(out).readline().split())
        if n == len(acols):
            skipped += 1
            continue
        print(f"  chr{c}: existing .l2.M has {n} values, expected {len(acols)} - rewriting")

    an = pq.read_table(anf).to_pandas()
    ld = set_snpid_index(ld)
    an = set_snpid_index(an)
    missing = ~ld.index.isin(an.index)
    if missing.any():
        sys.exit(f"FATAL chr{c}: {missing.sum()} LD-score SNPs absent from {anf}; "
                 "the LD-scores and annotations are not from the same build")
    absent = [x for x in acols if x not in an.columns]
    if absent:
        sys.exit(f"FATAL chr{c}: LD-score columns missing from the annotations: {absent[:5]}")

    M = an.loc[ld.index, acols].values.sum(axis=0, dtype=np.float64)
    assert M.shape == (len(acols),), f"chr{c}: M has shape {M.shape}"
    if not np.isfinite(M).all():
        sys.exit(f"FATAL chr{c}: M contains non-finite values")
    tmp = out + ".tmp"
    with open(tmp, "w") as f:
        f.write("\t".join(map(str, M)) + "\n")
    os.replace(tmp, out)
    wrote += 1
    mins.append(M.min()); maxs.append(M.max())
    if base_seen is None:
        base_seen = "base" in acols
    print(f"  chr{c}: {len(acols)} annotations over {len(ld):,} panel SNPs -> "
          f"{os.path.basename(out)}  (M in [{M.min():.4g}, {M.max():.4g}])", flush=True)

print(f"\n{wrote} written, {skipped} already correct")
if mins:
    print(f"M across chromosomes: min {min(mins):.4g}, max {max(maxs):.4g}")
    if min(mins) <= 0:
        print("\nNOTE  at least one annotation has a non-positive column sum. That is")
        print("      possible for CONTINUOUS annotations (the Borzoi PCs are signed)")
        print("      but M enters S-LDSC as the per-annotation SNP count, so a sum")
        print("      near zero makes that annotation's h2 share ill-defined.")
if base_seen is False:
    print("\nNOTE  there is NO 'base' annotation in this set. Standard S-LDSC always")
    print("      includes an all-ones 'base' column to absorb the genome-wide mean")
    print("      per-SNP heritability. Without it the regression has no such term.")
