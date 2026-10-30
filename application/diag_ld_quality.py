#!/usr/bin/env python3
"""Are the badly-conditioned blocks a prep bug or a genuine panel property?

Ten blocks came out of job 4487195 with min eig < -0.1, one at -0.6311, where
regularize_ld's own notes put real reference panels at -0.017 to -0.069. Two
candidate explanations, and they are distinguishable.

  A prep bug. prep_block.py does
        R = sparse.load_npz(pre + ".npz").toarray(); R = R + R.T
        np.fill_diagonal(R, 1.0)
    which is correct ONLY if the stored matrix is STRICTLY upper triangular.
    Any entry in the lower triangle is added to its mirror, doubling that
    off-diagonal - and a correlation above 1 makes a matrix violently
    non-PSD. The diagonal is fine either way, because fill_diagonal
    overwrites it, so this failure is invisible unless you look for it.

  A panel property. The Broad UKBB panels are stored sparsified: small
    correlations are zeroed and values are compressed. Zeroing entries
    destroys positive semi-definiteness on its own, and that is expected -
    it is why regularize_ld exists at all.

The test: max |off-diagonal| in the block matrix. Above 1.0 means the prep
bug; at or below 1.0 with a negative spectrum means sparsification, which is
the panel's own problem and not ours. Window triangularity is checked
directly, which settles it at the source.
"""
import argparse, os
import numpy as np
import pandas as pd
import scipy.sparse as sparse

BAD = ["ldEUR_17_0027", "ldEUR_03_0035", "ldEUR_15_0019", "ldEUR_04_0077",
       "ldEUR_20_0018", "ldEUR_07_0032", "ldEUR_17_0038", "ldEUR_01_0056",
       "ldEUR_04_0098", "ldEUR_14_0034"]
CONTROL = ["ldEUR_01_0034", "ldEUR_19_0031", "ldEUR_08_0054", "ldEUR_01_0067",
           "ldEUR_08_0045"]

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=os.environ.get("APP_ROOT", "."))
ap.add_argument("--ld-dir", default=None, help="default <root>/ld")
a = ap.parse_args()

blocks = pd.read_csv(os.path.join(a.root, "blocks.csv"))
win = dict(zip(blocks.region_id, blocks.ld_region))
ld_dir = a.ld_dir or os.path.join(a.root, "ld")

print(f"{'region':<16} {'p':>6} {'diag min':>9} {'diag max':>9} {'max|off|':>9} "
      f"{'n |r|>1':>8} {'n |r|>1.001':>12} {'NaN':>5} {'group':>8}")
rows = []
for grp, rids in (("BAD", BAD), ("control", CONTROL)):
    for rid in rids:
        f = os.path.join(a.root, "regions", rid, "locus_ld.npy")
        if not os.path.isfile(f):
            print(f"{rid:<16} (no locus_ld.npy)"); continue
        R = np.load(f).astype(np.float64)
        d = np.diag(R).copy()
        off = R.copy(); np.fill_diagonal(off, 0.0)
        mx = float(np.abs(off).max())
        n1 = int((np.abs(off) > 1.0).sum() // 2)
        n2 = int((np.abs(off) > 1.001).sum() // 2)
        print(f"{rid:<16} {R.shape[0]:>6} {d.min():>9.5f} {d.max():>9.5f} {mx:>9.5f} "
              f"{n1:>8} {n2:>12} {int(np.isnan(R).sum()):>5} {grp:>8}")
        rows.append((rid, grp, mx, n1))
        del R, off

print("\n=== window triangularity at source (settles prep bug vs panel) ===")
seen = {}
for rid in BAD + CONTROL:
    w = win.get(rid)
    if w is None or w in seen:
        continue
    p = os.path.join(ld_dir, f"{w}.npz")
    if not os.path.isfile(p):
        print(f"  {w:<28} (.npz not found under {ld_dir})"); seen[w] = 1; continue
    S = sparse.load_npz(p).tocoo()
    lower = int((S.row > S.col).sum())
    diag  = int((S.row == S.col).sum())
    upper = int((S.row < S.col).sum())
    mx = float(np.abs(S.data).max()) if S.nnz else 0.0
    print(f"  {w:<28} nnz={S.nnz:>10,}  upper={upper:>10,}  lower={lower:>8,}  "
          f"diag={diag:>7,}  max|v|={mx:.5f}  "
          f"{'STRICTLY UPPER (R+R.T is correct)' if lower == 0 else 'HAS LOWER ENTRIES -> R+R.T DOUBLES THEM'}")
    seen[w] = 1

bad_mx = [m for _, g, m, _ in rows if g == "BAD"]
if bad_mx and max(bad_mx) > 1.0:
    print("\nVERDICT: correlations exceed 1 -> prep bug in the R + R.T symmetrisation.")
else:
    print("\nVERDICT: no correlation exceeds 1. The negative spectrum is the panel's")
    print("own sparsification, not a prep bug - these blocks are genuinely far from")
    print("being correlation matrices and no choice of floor makes them trustworthy.")
