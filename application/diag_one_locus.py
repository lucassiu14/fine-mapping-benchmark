#!/usr/bin/env python3
"""Is the ~51-variant configuration genuinely better, or an artefact?

Reimplements abf() from generate_credible_sets.py exactly (minus the
mean_memo offset) so arbitrary configurations can be scored and compared on
one locus.

Self-check: reformat_memo normalises the training memo by its max, so the best
training configuration is stored as exp(0) = 1.0 and its log-BF IS mean_memo.
If logbf(top training config) reproduces mean_memo, the reimplementation is
faithful. That check is printed first -- read nothing below it if it fails.

The question this answers: abf stores exp(logBF - mean_memo), and float64
sends anything below -745 to exactly 0. If the 51-variant configuration beats
{lead} by far more than 745 nats, the preference is REAL and no amount of
log-space arithmetic will change the ranking -- the model genuinely wants 51
causal variants at this n, which is a calibration problem. If the gap is
smaller than that, or vanishes in float64, the arithmetic is the problem.

Decomposition is printed per term because they scale differently:
  quad    ~ O(n_sub), grows with every variant added
  logdet  ~ O(log n_sub) per variant
  prior   = k*log(p0) + (p-k)*log(1-p0), FIXED penalty ~log(p0) per variant
If quad/k exceeds |log(p0)| the posterior pays to add variants indefinitely,
until it hits sparse_concrete. That is the thing to look at.
"""
import argparse, os, pickle
import numpy as np
import pandas as pd
import torch

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=os.environ.get("APP_ROOT", "."))
ap.add_argument("--region", default="ldEUR_08_0045")
ap.add_argument("--n-sub", type=float, default=337000.0)   # 03_fb_xregion.pbs N_LIKELIHOOD
ap.add_argument("--sigma-sq", type=float, default=0.0899)  # beatrice_joint.py default
a = ap.parse_args()

rd = os.path.join(a.root, "regions", a.region)
zs = pd.read_table(os.path.join(rd, "locus.z"), sep=" ", header=None)
names = zs[0].to_numpy()
z_np = zs[1].to_numpy(float)
ld_np = np.load(os.path.join(rd, "locus_ld.npy"))
p = len(z_np)
print(f"region {a.region}: p={p}  n_sub={a.n_sub:g}  sigma_sq={a.sigma_sq}")

def regularize(A):
    """Verbatim regularize_ld: nearest PD CORRELATION matrix via eigenvalue clip."""
    A = (A + A.T) / 2
    try:
        np.linalg.cholesky(A - 1e-3 * np.eye(A.shape[0], dtype=A.dtype))
        return A, 0.0
    except np.linalg.LinAlgError:
        pass
    w, V = np.linalg.eigh(A)
    mn = w.min()
    w = np.clip(w, 1e-3, None)
    A2 = (V * w) @ V.T
    d = np.sqrt(np.diag(A2))
    A2 = A2 / np.outer(d, d)
    return A2, mn

def build(dtype):
    tz = torch.tensor(z_np, dtype=dtype)
    A, mn = regularize(ld_np.astype(np.float64 if dtype == torch.float64 else np.float32))
    tl = torch.tensor(A, dtype=dtype)
    S = torch.linalg.solve(tl, tz.unsqueeze(1))
    p0 = torch.full((p,), 1.0 / p, dtype=dtype)
    return tz, tl, S, p0, mn

def logbf(ind, tz, tl, S, p0, dtype):
    """abf() without the mean_memo offset. Returns (total, logdet, quad, prior)."""
    lp1 = torch.log(1 - p0)
    if len(ind) == 0:
        t = float(lp1.sum())
        return t, 0.0, 0.0, t
    ind = sorted(ind)
    w = a.n_sub * a.sigma_sq * torch.ones(len(ind), dtype=dtype)
    VU = tl[ind][:, ind] * w.unsqueeze(0)
    zU = (tz[ind] * w).unsqueeze(0)
    VS = torch.mm(tl[ind, :], S)
    eye = torch.eye(len(ind), dtype=dtype)
    quad = float(torch.mm(torch.mm(zU, torch.inverse(eye + VU)), VS) / 2)
    ldet = float(-torch.logdet(eye + VU) / 2)
    prior = lp1.clone(); prior[ind] = torch.log(p0[ind])
    pr = float(prior.sum())
    return ldet + quad + pr, ldet, quad, pr

with open(os.path.join(a.root, "joint_out", a.region, "res"), "rb") as fh:
    res = pickle.load(fh)
memo, mean_memo = res["memo"], float(np.asarray(res["mean_memo"]).squeeze())

vals = {k: float(np.asarray(v).squeeze()) for k, v in memo.items()}
top = max(vals, key=vals.get)
tot = sum(vals.values())
pip = np.zeros(p)
for k, v in vals.items():
    for i in k:
        if i < p:
            pip[i] += v
pip = pip / tot if tot > 0 else pip

lead = int(np.argmax(np.abs(z_np)))
seed = int(np.argmax(pip))
print(f"mean_memo (max training log-BF) = {mean_memo:.3f}")
print(f"top training config: |k|={len(top)}  stored value={vals[top]:.6g}")
print(f"lead  idx={lead} {names[lead]} |Z|={abs(z_np[lead]):.2f} trainPIP={pip[lead]:.4g}")
print(f"seed  idx={seed} {names[seed]} |Z|={abs(z_np[seed]):.2f} trainPIP={pip[seed]:.4g}")

by_z = list(np.argsort(-np.abs(z_np)))
cands = {
    "{} (empty)":            [],
    "{lead}":                [lead],
    "{seed}":                [seed],
    "{seed,lead}":           sorted({seed, lead}),
    "top10 by |Z|":          sorted(by_z[:10]),
    "top training config":   sorted(top),
    "topTrain cap to 10":    sorted(sorted(top, key=lambda i: -abs(z_np[i]))[:10]),
}

for dtype, nm in ((torch.float64, "float64"), (torch.float32, "float32")):
    tz, tl, S, p0, mn = build(dtype)
    print(f"\n=== {nm} (LD min eig before clip: "
          f"{'no shift needed' if mn == 0 else f'{mn:.4g}'}) ===")
    print(f"{'configuration':<22} {'|k|':>4} {'logBF':>12} {'-logdet/2':>11} "
          f"{'quad':>12} {'prior':>10} {'vs mean_memo':>13} {'exp()':>9} {'quad/|k|':>9}")
    for label, ind in cands.items():
        t, ldet, quad, pr = logbf(ind, tz, tl, S, p0, dtype)
        d = t - mean_memo
        ex = "0 (gone)" if d < -745 else f"{np.exp(min(d, 700.0)):.2e}"
        qk = f"{quad/len(ind):.2f}" if ind else "-"
        print(f"{label:<22} {len(ind):>4} {t:>12.2f} {ldet:>11.2f} {quad:>12.2f} "
              f"{pr:>10.2f} {d:>13.2f} {ex:>9} {qk:>9}")

tz, tl, S, p0, _ = build(torch.float64)
pen = float(np.log(1.0 / p) - np.log(1 - 1.0 / p))
print(f"\nper-variant prior penalty log(p0/(1-p0)) = {pen:.2f} nats")
print("marginal log-BF of growing the |Z|-ordered set (float64):")
prev = logbf([], tz, tl, S, p0, torch.float64)[0]
for k in (1, 2, 3, 5, 10, 20, 30, 40, 51):
    if k > p:
        break
    t = logbf(sorted(by_z[:k]), tz, tl, S, p0, torch.float64)[0]
    print(f"  |k|={k:>3}  logBF={t:>12.2f}   gain over previous row={t-prev:>10.2f}")
    prev = t
print("\nIf gain stays positive out to 51, the ABF pays to add variants at this n\n"
      "and the cap is what stops it: that is calibration, not arithmetic.")
