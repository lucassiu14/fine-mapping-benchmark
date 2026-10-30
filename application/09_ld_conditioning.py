#!/usr/bin/env python3
"""Per-block LD conditioning, recovered from an 03 log and written to results/.

Why this exists: the eigenvalue floor regularize_ld applies is the single
biggest determinant of how much a block's LD was altered, and it was only ever
visible in a job log. Pairing is positional and exact - joint_trainer's read_ld
calls regularize_ld once per region while loading, in manifest order, and
generate_credible_sets reuses the handed-over LD_t rather than re-regularising
- so the Nth "Projected LD" line is the Nth manifest row. The script refuses to
write if the counts disagree.

Flooring the spectrum at f caps every pairwise correlation near 1 - f: the
eigenvalues of [[1,r],[r,1]] are 1 +/- r, so lifting 1 - r to f bounds r at
1 - f. Measured on job 4487195: median floor 0.0305 (r ceiling 0.97), p90
0.0703 (0.93), max 0.6311 (0.37). So 200 of 254 blocks are capped below 0.98
and 10 below 0.90.

Two flags:
  distorted  floor > 0.05, so r is capped below 0.95. Credible sets widen and
             PIPs split across tagging variants; the co-inclusion penalty
             weakens (-11.8 to -2.2 nats at the extreme) but does NOT change
             sign, so this degrades resolution rather than inverting it.
  bad_ld     min eig < -0.1. regularize_ld's own notes put real reference
             panels at -0.017 to -0.069, so this is an order of magnitude
             outside the expected range and indicates the panel for that block
             is not close to a correlation matrix at all. No choice of floor
             makes such a block trustworthy - treat it as a data-quality
             problem, and note four LD windows arrived via the .npz2 404
             fallback.
"""
import argparse, os, re, sys
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=os.environ.get("APP_ROOT", "."))
ap.add_argument("--log", required=True, help="the 03_fb_xregion job log (logs/<id>.pbs-7.OU)")
ap.add_argument("--manifest", default=None, help="default <root>/manifest.tsv")
ap.add_argument("--distorted-floor", type=float, default=0.05)
ap.add_argument("--bad-ld-eig", type=float, default=-0.1)
a = ap.parse_args()

man = a.manifest or os.path.join(a.root, "manifest.tsv")
rows = pd.read_table(man)
rids = [t.rstrip("/").split("/")[-1] for t in rows["target"]]

pat = re.compile(r"min eigenvalue was ([-0-9.e+]+), floor ([-0-9.e+]+)")
hits = pat.findall(open(a.log, errors="replace").read())
if len(hits) != len(rids):
    sys.exit(f"REFUSING to write: {len(rids)} manifest rows but {len(hits)} floor "
             f"lines in {a.log}. The positional pairing is only valid when they "
             f"match - a block that took regularize_ld's Cholesky fast path "
             f"prints nothing, so a PD block would shift every later row.")

df = pd.DataFrame({"region_id": rids,
                   "min_eig": [float(m) for m, _ in hits],
                   "floor":   [float(f) for _, f in hits]})
df["r_ceiling"] = 1.0 - df["floor"]
df["distorted"] = df["floor"] > a.distorted_floor
df["bad_ld"]    = df["min_eig"] < a.bad_ld_eig

out = os.path.join(a.root, "results", "ld_conditioning.csv")
os.makedirs(os.path.dirname(out), exist_ok=True)
df.to_csv(out, index=False)

q = df["floor"].quantile([.5, .75, .9, .99]).to_dict()
print(f"{len(df)} blocks -> {out}")
print(f"  floor: median {q[0.5]:.4f} (r<={1-q[0.5]:.3f})  p75 {q[0.75]:.4f}  "
      f"p90 {q[0.9]:.4f} (r<={1-q[0.9]:.3f})  p99 {q[0.99]:.4f}  max {df.floor.max():.4f} "
      f"(r<={1-df.floor.max():.3f})")
print(f"  distorted (floor > {a.distorted_floor}): {int(df.distorted.sum())}")
print(f"  bad_ld (min eig < {a.bad_ld_eig}):       {int(df.bad_ld.sum())}")
if df.bad_ld.any():
    print("\n  blocks whose LD is too far from a correlation matrix to trust:")
    for _, r in df[df.bad_ld].sort_values("min_eig").iterrows():
        print(f"    {r.region_id:18s} min eig {r.min_eig:+.4f}  floor {r.floor:.4f}  "
              f"r capped at {r.r_ceiling:.3f}")
