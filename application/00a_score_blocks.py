"""Score every LDetect EUR block by its strongest association.

Input : gwfm_regions_ldetect_EUR.csv (ships with the benchmark package at
        inst/extdata/) and the munged sumstats parquet.
Output: blocks_scored.csv with max_chi2 and min_p per block, which
        00_make_blocks.py then thresholds.

The LDetect EUR blocks are Berisa & Pickrell (2016) approximately-independent
LD blocks for European ancestry, which is the right partition here because the
GWAS and the LD panel are both UK Biobank white-British.
"""
import os
import numpy as np, pandas as pd, pyarrow.parquet as pq

ROOT   = os.environ.get("APP_ROOT", os.path.dirname(os.path.abspath(__file__)))
BLOCKS = os.environ.get("LDETECT_CSV", os.path.join(ROOT, "gwfm_regions_ldetect_EUR.csv"))
SUMST  = os.path.join(ROOT, "chol_sumstats.parquet")
OUT    = os.path.join(ROOT, "blocks_scored.csv")

reg = pd.read_csv(BLOCKS)
ss  = pq.ParquetFile(SUMST).read(columns=["CHR", "BP", "Z"]).to_pandas()
ss["chi2"] = ss.Z.values ** 2

rows = []
for ch, g in reg.groupby("chrom"):
    s = ss[ss.CHR == ch].sort_values("BP")
    bp, chi2 = s.BP.values, s.chi2.values
    for r in g.itertuples():
        lo, hi = np.searchsorted(bp, [r.block_start, r.block_end], side="left")
        m = chi2[lo:hi].max() if hi > lo else np.nan
        rows.append(dict(region_id=r.region_id, chrom=int(ch),
                         start=int(r.start), end=int(r.end),
                         block_start=int(r.block_start), block_end=int(r.block_end),
                         notes=r.notes, max_chi2=m))

d = pd.DataFrame(rows)
# two-sided p from a 1-df chi-square, via the normal survival function in logs
from scipy.stats import norm
d["min_p"] = 2.0 * norm.sf(np.sqrt(d.max_chi2.values))
d.to_csv(OUT, index=False)
print(f"scored {len(d)} blocks -> {OUT}")
for thr in (1e-4, 1e-6, 5e-8):
    print(f"  min_p < {thr:g}: {(d.min_p < thr).sum()} blocks")
