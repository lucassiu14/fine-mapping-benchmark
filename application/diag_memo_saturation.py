#!/usr/bin/env python3
"""Where does fb_xregion's posterior mass actually sit?

The 1e15-cap hypothesis is REFUTED (n_cap == 0 across sampled regions), so
this measures the two things that are left.

1. Did the exhaustive sweep run?  add(..., range(bp)) in
   generate_credible_sets.main writes one memo entry per variant, so a sweep
   that ran leaves every index in [0, bp) present in some key. n_seen < bp
   means cond_stepwise_causal hit `return []` at its gamma_key gate (the only
   early exit before the first add) and the locus was never swept.

2. Whose configurations hold the mass?  The sweep can only write keys with
   |k| <= prior_causal + 1 = 11. Anything larger comes from the training
   sampler, where sparse_concrete sets the scale (~43-49). Mass sitting on
   |k| > 11 means the reported PIP is a variational selection frequency --
   how often the sampler nominated a variant -- not a posterior. A variant
   picked in every sample reads 1.0 whatever its Z score.
"""
import argparse, os, pickle, sys
import numpy as np
import pandas as pd

SWEEP_MAX = 11  # prior_causal (<=10) + the one index add() appends

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=os.environ.get("APP_ROOT", "."))
ap.add_argument("--limit", type=int, default=0, help="0 = every region")
a = ap.parse_args()

jd = os.path.join(a.root, "joint_out")
rids = sorted(d for d in os.listdir(jd)
               if os.path.isfile(os.path.join(jd, d, "res"))
               and os.path.isfile(os.path.join(jd, d, "pip.csv")))
if a.limit:
    rids = rids[: a.limit]

print(f"{'region':<16} {'bp':>7} {'n_seen':>7} {'swept':>6} {'n_memo':>7} "
      f"{'n_zero':>7} {'sizes(min/med/max)':>20} {'massBig':>8} {'sumPIP':>8} {'maxPIP':>7}")
n_unswept = n_big = 0
for rid in rids:
    with open(os.path.join(jd, rid, "res"), "rb") as fh:
        memo = pickle.load(fh)["memo"]
    bp = len(pd.read_csv(os.path.join(jd, rid, "pip.csv")))

    vals = {k: float(np.asarray(v).squeeze()) for k, v in memo.items()}
    nz = {k: v for k, v in vals.items() if v > 0.0}

    seen = set()
    for k in vals:
        seen.update(k)
    swept = len(seen) >= bp

    sizes = np.array([len(k) for k in nz]) if nz else np.array([0])
    tot = sum(nz.values())
    big = sum(v for k, v in nz.items() if len(k) > SWEEP_MAX)
    frac_big = big / tot if tot > 0 else 0.0

    pip = np.zeros(bp)
    for k, v in nz.items():
        for i in k:
            if i < bp:
                pip[i] += v
    pip = pip / tot if tot > 0 else pip

    n_unswept += (not swept)
    n_big += (frac_big > 0.5)
    print(f"{rid:<16} {bp:>7} {len(seen):>7} {str(swept):>6} {len(vals):>7} "
          f"{len(vals)-len(nz):>7} "
          f"{f'{sizes.min()}/{int(np.median(sizes))}/{sizes.max()}':>20} "
          f"{frac_big:>8.3f} {pip.sum():>8.3f} {pip.max():>7.4f}")

print(f"\n{n_unswept}/{len(rids)} regions were NEVER swept (gamma_key gate returned [])")
print(f"{n_big}/{len(rids)} regions carry >50% of mass on |k| > {SWEEP_MAX} "
      f"(training-sampler configurations, i.e. PIP = selection frequency)")
