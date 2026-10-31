#!/usr/bin/env python3
"""How did the four methods do on the cholesterol application?

There is no ground truth here, so nothing below is power or FDR. What can be
measured is resolution (how tightly each method concentrates its posterior),
concordance (whether methods agree on which variant), and whether confident
calls track the association signal at all.

Two stratifications are applied and NOT pooled, because pooling them would
average over cells that answer different questions:

  significance  A block whose strongest association is |Z| < 5.45 (p > 5e-8)
                has no genome-wide signal to fine-map. Every method is fitting
                noise there, so a method that declines to concentrate is
                behaving correctly, not timidly. Mixing those blocks with real
                signals makes a conservative method look weak and an
                aggressive one look powerful.

  LD quality    The adaptive eigenvalue floor caps pairwise correlations near
                1 - floor, so high-floor blocks have flattened LD: credible
                sets widen and PIPs split across tagging variants. `distorted`
                is floor > 0.05 (r capped below 0.95); `bad_ld` is min eig <
                -0.1, where the panel is too far from a correlation matrix to
                trust at all. Reported separately so a resolution difference
                driven by regularisation is not read as a difference between
                methods.

Credible sets are the smallest set of variants whose PIPs sum to >= 0.95. Where
a method's total PIP in a region is below 0.95 no such set exists, and that is
reported as a miss rather than silently dropped - a method can look
high-resolution purely by never forming a set.
"""
import argparse, os
import numpy as np
import pandas as pd

ZSIG = 5.45          # |Z| for p < 5e-8 two-sided
CONF = 0.95
NULLZ = 2.0

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=os.environ.get("APP_ROOT", "."))
a = ap.parse_args()
R = a.root

df = pd.read_parquet(os.path.join(R, "results", "pips.parquet"))
cond = pd.read_csv(os.path.join(R, "results", "ld_conditioning.csv"))
methods = sorted(df.method.unique())

# ---- per (region, method) summaries -------------------------------------
def cs_size(p):
    """Smallest count of variants whose PIP sums to >= CONF; NaN if total < CONF."""
    s = np.sort(p)[::-1]
    if s.sum() < CONF:
        return np.nan
    return int(np.searchsorted(np.cumsum(s), CONF) + 1)

recs = []
for (rid, meth), g in df.groupby(["region_id", "method"], sort=False):
    pip = g.pip.to_numpy(float)
    z = g.Z.to_numpy(float)
    j = int(np.argmax(pip))
    k = int(np.argmax(np.abs(z)))
    recs.append(dict(region_id=rid, method=meth, p=len(g),
                     sumPIP=pip.sum(), maxPIP=pip.max(),
                     cs=cs_size(pip),
                     top_snp=g.SNP.iloc[j], top_absZ=abs(z[j]),
                     lead_snp=g.SNP.iloc[k], lead_absZ=abs(z[k]),
                     lead_pip=pip[k], top_is_lead=(j == k)))
S = pd.DataFrame(recs).merge(cond[["region_id", "floor", "distorted", "bad_ld"]],
                             on="region_id", how="left")
S["sig"] = S.lead_absZ >= ZSIG

nsig = S[S.method == methods[0]].sig.sum()
print(f"{S.region_id.nunique()} regions, {len(methods)} methods, {len(df):,} rows")
print(f"regions with a genome-wide signal (|Z| >= {ZSIG}): {nsig}  "
      f"| without: {S[S.method==methods[0]].shape[0]-nsig}")
print(f"LD flags: distorted (floor>0.05) {int(cond.distorted.sum())}, "
      f"bad_ld (min eig<-0.1) {int(cond.bad_ld.sum())} of {len(cond)}\n")

def block(sub, title):
    if sub.empty:
        print(f"=== {title}: no regions ===\n"); return
    n = sub.region_id.nunique()
    print(f"=== {title}  (n={n} regions) ===")
    print(f"{'method':<18} {'med sumPIP':>11} {'med maxPIP':>11} {'maxPIP>=.95':>12} "
          f"{'med CS':>7} {'CS<=5':>7} {'no CS':>6} {'top=lead':>9} {'med lead PIP':>13}")
    for m in methods:
        g = sub[sub.method == m]
        cs = g.cs.dropna()
        print(f"{m:<18} {g.sumPIP.median():>11.3f} {g.maxPIP.median():>11.4f} "
              f"{int((g.maxPIP>=CONF).sum()):>12} "
              f"{(cs.median() if len(cs) else float('nan')):>7.0f} "
              f"{int((cs<=5).sum()):>7} {int(g.cs.isna().sum()):>6} "
              f"{int(g.top_is_lead.sum()):>9} {g.lead_pip.median():>13.4f}")
    print()

block(S, "ALL regions")
block(S[S.sig], f"regions WITH a signal (|Z| >= {ZSIG})")
block(S[~S.sig], f"regions WITHOUT a signal (|Z| < {ZSIG})")

# ---- do confident calls track the signal? ------------------------------
print(f"=== confident calls (PIP >= {CONF}): does the called variant have signal? ===")
print(f"{'method':<18} {'n calls':>8} {'med |Z|':>8} {'|Z|<2':>7} {'|Z|<5.45':>9} "
      f"{'in sig regions':>15}")
for m in methods:
    g = S[(S.method == m) & (S.maxPIP >= CONF)]
    if g.empty:
        print(f"{m:<18} {0:>8}"); continue
    print(f"{m:<18} {len(g):>8} {g.top_absZ.median():>8.2f} "
          f"{int((g.top_absZ<NULLZ).sum()):>7} {int((g.top_absZ<ZSIG).sum()):>9} "
          f"{int(g.sig.sum()):>15}")
print()

# ---- agreement ----------------------------------------------------------
print("=== pairwise agreement on the top-PIP variant (signal regions only) ===")
sig_ids = set(S[S.sig].region_id)
piv = (S[S.region_id.isin(sig_ids)]
       .pivot(index="region_id", columns="method", values="top_snp"))
print(f"{'pair':<40} {'same top variant':>17}")
for i, m1 in enumerate(methods):
    for m2 in methods[i+1:]:
        ok = (piv[m1] == piv[m2]).sum()
        print(f"{m1+' vs '+m2:<40} {ok:>8} / {len(piv):<6} ({100*ok/len(piv):.0f}%)")
print()

print("=== per-region Spearman of PIPs, median across signal regions ===")
rows = []
for rid, g in df[df.region_id.isin(sig_ids)].groupby("region_id", sort=False):
    w = g.pivot_table(index="SNP", columns="method", values="pip")
    for i, m1 in enumerate(methods):
        for m2 in methods[i+1:]:
            rows.append((f"{m1} vs {m2}", w[m1].corr(w[m2], method="spearman")))
sp = pd.DataFrame(rows, columns=["pair", "rho"]).groupby("pair").rho.median()
for k, v in sp.items():
    print(f"  {k:<40} {v:+.3f}")
print()

# ---- strong disagreement ------------------------------------------------
print(f"=== both confident (PIP >= {CONF}) but different variants, signal regions ===")
conf = S[(S.maxPIP >= CONF) & S.sig][["region_id", "method", "top_snp", "top_absZ"]]
cp = conf.pivot(index="region_id", columns="method", values="top_snp")
cz = conf.pivot(index="region_id", columns="method", values="top_absZ")
for i, m1 in enumerate(methods):
    for m2 in methods[i+1:]:
        both = cp[[m1, m2]].dropna()
        if both.empty:
            continue
        diff = both[both[m1] != both[m2]]
        print(f"  {m1} vs {m2}: both confident in {len(both)} regions, "
              f"disagree in {len(diff)}")
print()

# ---- LD-quality stratification -----------------------------------------
block(S[S.sig & ~S.distorted.astype(bool)], f"signal regions, LD NOT distorted (floor<=0.05)")
block(S[S.sig & S.distorted.astype(bool)],  f"signal regions, LD distorted (floor>0.05)")
if S.bad_ld.any():
    block(S[S.bad_ld.astype(bool)], "bad_ld regions (min eig < -0.1), signal or not")

S.to_csv(os.path.join(R, "results", "region_method_summary.csv"), index=False)
print(f"per-(region,method) table -> {os.path.join(R,'results','region_method_summary.csv')}")
