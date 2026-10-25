#!/usr/bin/env python
"""Per block: annotation coverage AND the largest contiguous unannotated run.

Coverage alone is the wrong summary. ldEUR_11_0032 retains 23.8% of the GWAS
variants in its span, but they are not thinned at random -- the first 1.75 Mb
has ZERO retained variants and the last 1.2 Mb has 87-91%. A fine-mapper given
that block searches only the final 40% of the span, so a causal variant in the
missing part is not even a candidate and the credible set is confidently wrong
rather than uncertain.

A block at 90% coverage spread uniformly is fine. A block at 60% coverage where
the missing 40% is one contiguous hole is a truncated search space. This reports
both so the distinction is visible.

For each block: the GWAS variants in its span are marked retained or not (by
membership in snps.csv), and the longest run of consecutive unretained ones is
reported, with its base-pair extent.
"""
import os, csv, sys
import numpy as np, pandas as pd, pyarrow.parquet as pq

ROOT = os.environ.get("APP_ROOT", ".")
rows = list(csv.DictReader(open(os.path.join(ROOT, "blocks.csv"))))
ss = pq.read_table(os.path.join(ROOT, "chol_sumstats.parquet"), columns=["CHR", "BP"])
chrs, bps = ss.column("CHR").to_numpy(), ss.column("BP").to_numpy()
by_chr = {c: np.sort(bps[chrs == c]) for c in sorted({int(r["chrom"]) for r in rows})}

out = []
for r in rows:
    rid, c = r["region_id"], int(r["chrom"])
    lo, hi = int(r["block_start"]), int(r["block_end"])
    d = os.path.join(ROOT, "regions", rid)
    excluded = os.path.exists(os.path.join(d, "EXCLUDED"))
    f = os.path.join(d, "snps.csv")
    a = by_chr[c]
    span = a[(a >= lo) & (a <= hi)]
    if not os.path.exists(f) or len(span) == 0:
        out.append((rid, c, lo, hi, 0, len(span), 0.0, len(span),
                    (hi - lo) if len(span) else 0, excluded))
        continue
    kept = set(pd.read_csv(f)["BP"].astype(int).tolist())
    mark = np.fromiter((b in kept for b in span), dtype=bool, count=len(span))
    # longest run of consecutive GWAS variants with none retained
    best_n, best_bp, n, start = 0, 0, 0, None
    for i, ok in enumerate(mark):
        if not ok:
            if n == 0:
                start = i
            n += 1
            if n > best_n:
                best_n = n
                best_bp = int(span[i] - span[start])
        else:
            n = 0
    cov = mark.sum() / len(span)
    out.append((rid, c, lo, hi, int(mark.sum()), len(span), cov,
                best_n, best_bp, excluded))

kept = [o for o in out if not o[9]]
print(f"analysed blocks: {len(kept)}   excluded: {len(out) - len(kept)}")
gapbp = np.array([o[8] for o in kept], dtype=float)
cov = np.array([o[6] for o in kept])
print(f"coverage           : min {cov.min():.1%}  median {np.median(cov):.1%}")
print(f"largest gap (kb)   : median {np.median(gapbp)/1e3:,.0f}  "
      f"p90 {np.percentile(gapbp, 90)/1e3:,.0f}  max {gapbp.max()/1e3:,.0f}")
for thr in (250e3, 500e3, 1e6):
    print(f"  analysed blocks with a contiguous unannotated run > {thr/1e6:.2f} Mb : "
          f"{int((gapbp > thr).sum())}")
print("\nworst 15 analysed blocks by largest contiguous unannotated run:")
print(f"  {'block':<16}{'span Mb':>9}  {'cov':>6}  {'gap Mb':>7}  {'gap vars':>9}")
for rid, c, lo, hi, k, n, cv, gn, gb, _ in sorted(kept, key=lambda t: -t[8])[:15]:
    print(f"  {rid:<16}{(hi-lo)/1e6:>9.2f}  {cv:>5.1%}  {gb/1e6:>7.2f}  {gn:>9,}")
