#!/usr/bin/env python3
"""Count memo entries pinned at abf()'s 1e15 cap.

abf() in generate_credible_sets.py stores exp(logBF - mean_memo) clamped at
_CAP = 1e15, where mean_memo is the MAX log-BF over the *training* memo. The
credible-set sweep exists to evaluate configurations training never sampled,
so ex > 0 is the normal case and every configuration more than
log(1e15) = 34.54 nats better than the training best is stored as the
IDENTICAL value 1e15. Those configurations then carry equal weight in
calculate_pip, which turns the PIP into a count:

    pip_i = #{saturated k : i in k} / #{saturated k}

Signatures, all of which this script measures directly:
  n_cap   > 0        -> posterior mass is a count, not a posterior
  n_zero             -> ex < -745 underflowed; these read as pip == 0 even
                        though the variant WAS visited (do not call this
                        "unexplored")
  sum(PIP) integral  -> every saturated configuration has the same size
"""
import argparse, os, pickle, sys
import numpy as np

CAP = 1e15

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=os.environ.get("APP_ROOT", "."))
ap.add_argument("--limit", type=int, default=0, help="0 = every region")
a = ap.parse_args()

jd = os.path.join(a.root, "joint_out")
rids = sorted(d for d in os.listdir(jd) if os.path.isfile(os.path.join(jd, d, "res")))
if a.limit:
    rids = rids[: a.limit]

print(f"{'region':<16} {'n_memo':>8} {'n_cap':>7} {'%cap':>6} {'n_zero':>7} "
      f"{'sumPIP':>9} {'int?':>5} {'maxPIP':>7} {'capsizes':>12}")
hit = 0
for rid in rids:
    with open(os.path.join(jd, rid, "res"), "rb") as fh:
        memo = pickle.load(fh)["memo"]
    vals = {k: float(np.asarray(v).squeeze()) for k, v in memo.items()}
    n = len(vals)
    cap = [k for k, v in vals.items() if v >= CAP]
    zero = sum(1 for v in vals.values() if v == 0.0)

    tot = sum(vals.values())
    bp = max((max(k) for k in vals if k), default=-1) + 1
    pip = np.zeros(bp)
    for k, v in vals.items():
        for i in k:
            pip[i] += v
    pip = pip / tot if tot > 0 else pip
    s = pip.sum()
    integral = abs(s - round(s)) < 1e-3
    sizes = sorted({len(k) for k in cap})
    hit += bool(cap)
    print(f"{rid:<16} {n:>8} {len(cap):>7} {100*len(cap)/n:>5.1f}% {zero:>7} "
          f"{s:>9.3f} {str(integral):>5} {pip.max():>7.4f} "
          f"{str(sizes[:4]) if sizes else '-':>12}")

print(f"\n{hit}/{len(rids)} regions have >=1 configuration pinned at the 1e15 cap")
