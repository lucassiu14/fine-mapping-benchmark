#!/usr/bin/env python3
"""Did the credible-set sweep run, and did it change the answer?

generate_credible_sets.py LOADS <target>/res and never saves it back, so the
pickle on disk is the TRAINING memo, frozen before the sweep. pip.csv is
written after the sweep. Comparing the two is therefore an exact test:

  pip.csv == memo-derived PIP  ->  the sweep changed nothing. The only path
                                   to that is cond_stepwise_causal hitting
                                   `return []` at its gamma_key gate, which
                                   fires when max(training PIP) < gamma_key.
  pip.csv != memo-derived PIP  ->  the sweep ran.

`gate` reports max(training PIP) >= GAMMA_KEY, i.e. what cond_stepwise_causal
decides on. If gate is False the region is never swept and its PIPs are the
variational training posterior alone.

Refuted earlier hypotheses, recorded so they are not re-tried:
  - the 1e15 cap in abf(): n_cap == 0 across sampled regions
  - "unexplored" variants: pip == 0 also arises from exp() underflow
"""
import argparse, os, pickle
import numpy as np
import pandas as pd

GAMMA_KEY = 0.2   # beatrice_annot.py flags.DEFINE_float('gamma_key', 0.2, ...)

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

print(f"{'region':<16} {'bp':>7} {'preMax':>7} {'gate':>5} {'postMax':>8} "
      f"{'preSum':>8} {'postSum':>8} {'same':>5} {'post>0':>7} {'swept':>6}")
n_gate = n_same = 0
for rid in rids:
    with open(os.path.join(jd, rid, "res"), "rb") as fh:
        memo = pickle.load(fh)["memo"]
    post = pd.read_csv(os.path.join(jd, rid, "pip.csv"))["pip"].to_numpy(float)
    bp = len(post)

    vals = {k: float(np.asarray(v).squeeze()) for k, v in memo.items()}
    tot = sum(vals.values())
    pre = np.zeros(bp)
    for k, v in vals.items():
        for i in k:
            if i < bp:
                pre[i] += v
    pre = pre / tot if tot > 0 else pre

    gate = pre.max() >= GAMMA_KEY
    same = np.allclose(pre, post, atol=1e-9, rtol=0)
    n_gate += (not gate)
    n_same += same
    print(f"{rid:<16} {bp:>7} {pre.max():>7.4f} {str(gate):>5} {post.max():>8.4f} "
          f"{pre.sum():>8.3f} {post.sum():>8.3f} {str(same):>5} "
          f"{int((post > 0).sum()):>7} {str(not same):>6}")

print(f"\n{n_gate}/{len(rids)} regions fail the gamma_key={GAMMA_KEY} gate "
      f"(max training PIP < {GAMMA_KEY}) -> cond_stepwise_causal returns []")
print(f"{n_same}/{len(rids)} regions have pip.csv identical to the training "
      f"posterior -> the sweep changed nothing")
