"""Build the final block list: LDetect EUR blocks at p < 1e-6, each mapped to
one UKBB 3 Mb LD window.

A UKBB LD window is [k*1e6+1, (k+3)*1e6+1] for integer k, so windows step by
1 Mb and overlap by 2 Mb. Any block up to 2 Mb wide always fits inside one;
wider blocks may not, depending on phase. For those, we pick the window that
(a) contains the block's lead variant and (b) retains the most of the block,
and trim the block to that window. Measured cost over the 252 blocks: 6.83 Mb,
~34,700 variants, 1.58% of the span, with every lead variant retained.

The alternative - stitching two windows and zeroing pairs more than 2 Mb apart
- was rejected: it keeps ~35k more variants but makes the LD inexact for
exactly the long-range pairs that a 3 Mb locus is meant to model.
"""
import os, math
import numpy as np, pandas as pd, pyarrow.parquet as pq

ROOT   = os.environ.get("APP_ROOT", os.path.dirname(os.path.abspath(__file__)))
SCORED = os.path.join(ROOT, "blocks_scored.csv")
SUMST  = os.path.join(ROOT, "chol_sumstats.parquet")
OUT    = os.path.join(ROOT, "blocks.csv")
PTHR   = float(os.environ.get("P_THRESHOLD", 1e-6))
WIN_MB = 3

b = pd.read_csv(SCORED)
s = b[b.min_p < PTHR].copy().reset_index(drop=True)
print(f"{len(s)} blocks at p < {PTHR:g}")

ss = pq.ParquetFile(SUMST).read(columns=["CHR", "BP", "Z"]).to_pandas()

rec = []
for r in s.itertuples():
    x = ss[(ss.CHR == r.chrom) & (ss.BP >= r.block_start) & (ss.BP <= r.block_end)]
    i = x.Z.abs().idxmax()
    lead, leadz = int(x.BP[i]), float(abs(x.Z[i]))

    # Candidate windows are exactly those containing the LEAD variant: any
    # window [k, k+3] with k*1e6+1 <= lead <= (k+3)*1e6+1. Deriving the range
    # from the block edges instead fails whenever the block is wider than a
    # window, which is the case we most need to handle.
    best = None
    klo = max(0, math.ceil((lead - 1) / 1e6) - WIN_MB)
    khi = math.floor((lead - 1) / 1e6)
    for k in range(klo, khi + 1):
        ws, we = k * 10**6 + 1, (k + WIN_MB) * 10**6 + 1
        if not (ws <= lead <= we):
            continue
        lo, hi = max(ws, r.block_start), min(we, r.block_end)
        if hi <= lo:
            continue
        ov = hi - lo
        if best is None or ov > best[0]:
            best = (ov, k, lo, hi)
    assert best is not None, f"{r.region_id}: no window contains its lead variant"
    ov, k, lo, hi = best
    rec.append(dict(region_id=r.region_id, chrom=int(r.chrom),
                    block_start=int(lo), block_end=int(hi),
                    ld_region=f"chr{int(r.chrom)}_{k*10**6+1}_{(k+WIN_MB)*10**6+1}",
                    min_p=float(r.min_p), lead_bp=lead, lead_abs_z=round(leadz, 3),
                    orig_start=int(r.block_start), orig_end=int(r.block_end),
                    trimmed_bp=int((r.block_end - r.block_start) - ov)))

d = pd.DataFrame(rec)
d.to_csv(OUT, index=False)
tr = d[d.trimmed_bp > 0]
print(f"  blocks trimmed       : {len(tr)}")
print(f"  bp trimmed in total  : {tr.trimmed_bp.sum()/1e6:.2f} Mb "
      f"({tr.trimmed_bp.sum()/(d.orig_end-d.orig_start).sum()*100:.2f}% of span)")
print(f"  lead variants lost   : 0 (enforced by construction)")
print(f"  distinct LD windows  : {d.ld_region.nunique()}")
d[["ld_region"]].drop_duplicates().to_csv(os.path.join(ROOT, "ld_windows.txt"),
                                          index=False, header=False)
print(f"  wrote {OUT} and ld_windows.txt")
