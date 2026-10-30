#!/usr/bin/env python3
"""Is co-inclusion of collinear variants still penalised under the new floor?

This is the soundness test for the adaptive eigenvalue floor, not an
effectiveness test. regularize_ld's own comment records why clipping replaced
a diagonal shift: the shift scaled every implied correlation by 1/(1+c), which
pushed strongly-correlated pairs across the ABF's crossover near r = 0.955 and
made the posterior FAVOUR including both members of a perfectly-tagged pair
(-11.1 nats at r = 0.996 became +8.2 nats at r = 0.931). Clipping at 1e-3 kept
correlations to within 0.0040 and preserved the penalty.

The adaptive floor distorts correlations more than that - 0.0132 to 0.0682
across blocks in job 4487195, and 0.0657 is what the rejected diagonal shift
cost. So the question is whether we have bought correct answers through a
mechanism that reintroduces the failure clipping was chosen to avoid.

For each strongly-correlated pair this prints

    delta = logBF({i,j}) - max(logBF({i}), logBF({j}))

NEGATIVE means co-inclusion is penalised, which is the correct behaviour: two
variants in near-perfect LD should not both be called causal. POSITIVE is the
diagonal-shift pathology. Compare the two floors side by side.
"""
import argparse, os
import numpy as np
import pandas as pd
import torch

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=os.environ.get("APP_ROOT", "."))
ap.add_argument("--region", default="ldEUR_08_0045")
ap.add_argument("--r-min", type=float, default=0.99, help="|r| threshold on the RAW LD")
ap.add_argument("--n-pairs", type=int, default=12)
ap.add_argument("--n-sub", type=float, default=float(os.environ.get("N_LIKELIHOOD", 373220)))
ap.add_argument("--sigma-sq", type=float, default=0.0899)
a = ap.parse_args()

rd = os.path.join(a.root, "regions", a.region)
zs = pd.read_table(os.path.join(rd, "locus.z"), sep=" ", header=None)
z_np = zs[1].to_numpy(float)
raw = np.load(os.path.join(rd, "locus_ld.npy")).astype(np.float64)
p = len(z_np)
print(f"region {a.region}: p={p}  n_sub={a.n_sub:g}  sigma_sq={a.sigma_sq}  "
      f"|r| > {a.r_min} on the RAW panel")

def regularize(A, floor):
    A = (A + A.T) / 2
    if floor is None or floor <= 1e-3:
        try:
            np.linalg.cholesky(A - 1e-3 * np.eye(A.shape[0], dtype=A.dtype))
            return A, 1e-3
        except np.linalg.LinAlgError:
            pass
    w, V = np.linalg.eigh(A)
    f = max(1e-3, abs(min(float(w.min()), 0.0))) if floor is None else float(floor)
    A2 = (V * np.clip(w, f, None)) @ V.T
    d = np.sqrt(np.diag(A2))
    return A2 / np.outer(d, d), f

def logbf(ind, tz, tl, S, p0):
    lp1 = torch.log(1 - p0)
    if len(ind) == 0:
        return float(lp1.sum())
    ind = sorted(ind)
    w = a.n_sub * a.sigma_sq * torch.ones(len(ind), dtype=torch.float64)
    VU = tl[ind][:, ind] * w.unsqueeze(0)
    zU = (tz[ind] * w).unsqueeze(0)
    VS = torch.mm(tl[ind, :], S)
    eye = torch.eye(len(ind), dtype=torch.float64)
    quad = float(torch.mm(torch.mm(zU, torch.inverse(eye + VU)), VS) / 2)
    ldet = float(-torch.logdet(eye + VU) / 2)
    prior = lp1.clone(); prior[ind] = torch.log(p0[ind])
    return ldet + quad + float(prior.sum())

# pairs from the RAW panel, so the selection does not depend on the floor
m = np.abs(raw) > a.r_min
np.fill_diagonal(m, False)
ii, jj = np.where(np.triu(m))
if len(ii) == 0:
    print(f"no pairs with |r| > {a.r_min} in this block; try --r-min 0.95")
    raise SystemExit(0)
# prefer pairs where at least one member carries signal: a penalty that only
# holds for noise pairs would not be reassuring
order = np.argsort(-np.maximum(np.abs(z_np[ii]), np.abs(z_np[jj])))
ii, jj = ii[order][: a.n_pairs], jj[order][: a.n_pairs]
print(f"{int(m.sum()//2)} such pairs in the block; showing {len(ii)} "
      f"with the largest max|Z|\n")

res = {}
for floor, nm in ((1e-3, "eig_floor=1e-3 (historical)"), (None, "adaptive floor")):
    A, f = regularize(raw.copy(), floor)
    tz = torch.tensor(z_np, dtype=torch.float64)
    tl = torch.tensor(A, dtype=torch.float64)
    S = torch.linalg.solve(tl, tz.unsqueeze(1))
    p0 = torch.full((p,), 1.0 / p, dtype=torch.float64)
    rows = []
    for i, j in zip(ii, jj):
        rows.append((abs(float(raw[i, j])), abs(float(A[i, j])),
                     abs(z_np[i]), abs(z_np[j]),
                     logbf([int(i)], tz, tl, S, p0),
                     logbf([int(j)], tz, tl, S, p0),
                     logbf([int(i), int(j)], tz, tl, S, p0)))
    res[nm] = (f, rows)

for nm, (f, rows) in res.items():
    print(f"=== {nm}  (floor applied {f:.4g}) ===")
    print(f"{'raw |r|':>8} {'reg |r|':>8} {'|Zi|':>7} {'|Zj|':>7} "
          f"{'logBF{i}':>10} {'logBF{j}':>10} {'logBF{i,j}':>11} {'delta':>9} {'penalised?':>11}")
    for rr, ar, zi, zj, bi, bj, bij in rows:
        d = bij - max(bi, bj)
        print(f"{rr:>8.4f} {ar:>8.4f} {zi:>7.2f} {zj:>7.2f} {bi:>10.2f} "
              f"{bj:>10.2f} {bij:>11.2f} {d:>9.2f} {str(d < 0):>11}")
    bad = sum(1 for r in rows if r[6] - max(r[4], r[5]) >= 0)
    print(f"  co-inclusion FAVOURED in {bad}/{len(rows)} pairs "
          f"({'PATHOLOGY' if bad else 'correct: all penalised'})\n")
