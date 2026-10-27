#!/usr/bin/env python3
"""Classify every fb_xregion region by how its PIPs were actually produced.

generate_credible_sets.py loads <target>/res and never saves it back, so the
pickle is the TRAINING memo frozen before the sweep, while pip.csv is written
after it. That makes the mode exactly determinable:

  GATE      max(training PIP) < gamma_key, so cond_stepwise_causal returns []
            at its only early exit and no sweep ever runs. pip.csv IS the
            variational training posterior: PIP is a selection frequency over
            size-`sparse_concrete` samples, and |Z| never enters it.

  UNDERFLOW the gate passed and the sweep ran, but pip.csv is identical to the
            training posterior. Training latched onto large configurations
            with enormous log-BFs, so mean_memo is large and every smaller
            configuration the sweep proposes has ex < -745 and is stored as
            exactly 0.0. It cannot shift the normaliser. Same endpoint as GATE.

  SWEPT     the sweep changed the answer. Note cond_stepwise_causal calls
            add(..., start_set, range(bp)), which enumerates {s, i} for every
            i, so the seed s sits in every configuration it writes: postSum
            near an integer means all mass is on that one size and the seed's
            PIP approaches 1 by construction, not by evidence.

Refuted hypotheses, recorded so they are not re-tried:
  - the 1e15 cap in abf(): n_cap == 0 across sampled regions
  - pip == 0 as "never explored": exp() underflow produces it too
"""
import argparse, os, pickle
import numpy as np
import pandas as pd

GAMMA_KEY = 0.2    # beatrice_annot.py: flags.DEFINE_float('gamma_key', 0.2, ...)
CONF      = 0.95   # "confident" PIP
NULL_Z    = 2.0

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=os.environ.get("APP_ROOT", "."))
ap.add_argument("--joint-out", default=None,
                help="fb_xregion output dir; default <root>/joint_out. Point it at a\n                     preserved run (e.g. joint_out.pre_fix_d0fb022) to compare\n                     against the current one.")
ap.add_argument("--limit", type=int, default=0, help="0 = every region")
a = ap.parse_args()

jd = a.joint_out or os.path.join(a.root, "joint_out")
rids = sorted(d for d in os.listdir(jd)
               if os.path.isfile(os.path.join(jd, d, "res"))
               and os.path.isfile(os.path.join(jd, d, "pip.csv")))
if a.limit:
    rids = rids[: a.limit]

rows = []
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

    sp = os.path.join(a.root, "regions", rid, "snps.csv")
    z = pd.read_csv(sp)["Z"].to_numpy(float) if os.path.isfile(sp) else np.full(bp, np.nan)
    if len(z) != bp:
        z = np.full(bp, np.nan)

    same = np.allclose(pre, post, atol=1e-9, rtol=0)
    gate = pre.max() >= GAMMA_KEY
    mode = "GATE" if not gate else ("UNDERFLOW" if same else "SWEPT")

    j = int(np.argmax(post))
    s = post.sum()
    rows.append(dict(region=rid, bp=bp, mode=mode, postMax=post.max(),
                     postSum=s, integral=abs(s - round(s)) < 1e-3,
                     topZ=abs(z[j]), leadZ=np.nanmax(np.abs(z)),
                     leadPIP=post[int(np.nanargmax(np.abs(z)))] if np.isfinite(z).any() else np.nan))

df = pd.DataFrame(rows)
print(f"=== {len(df)} regions by mode ===")
for m, g in df.groupby("mode"):
    print(f"  {m:<10} {len(g):>4}  ({100*len(g)/len(df):.1f}%)   "
          f"median postSum {g.postSum.median():.3f}   "
          f"integral sum {int(g.integral.sum())}/{len(g)}   "
          f"postMax>={CONF}: {int((g.postMax >= CONF).sum())}")

conf = df[df.postMax >= CONF]
print(f"\n=== {len(conf)} regions with postMax >= {CONF}: |Z| of that variant ===")
for m, g in conf.groupby("mode"):
    bad = g[g.topZ < NULL_Z]
    print(f"  {m:<10} n={len(g):>3}  median |Z| {g.topZ.median():>6.2f}   "
          f"|Z|<{NULL_Z}: {len(bad)}/{len(g)}")

print(f"\n=== confident calls on null variants (postMax>={CONF}, |Z|<{NULL_Z}) ===")
bad = conf[conf.topZ < NULL_Z].sort_values("topZ")
if len(bad):
    print(bad[["region", "bp", "mode", "postMax", "postSum", "integral",
               "topZ", "leadZ", "leadPIP"]].to_string(index=False,
               float_format=lambda v: f"{v:.4g}"))
else:
    print("  none")

df.to_csv("/tmp/fb_modes.csv", index=False)
print(f"\nfull table -> /tmp/fb_modes.csv")
