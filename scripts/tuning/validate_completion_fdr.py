#!/usr/bin/env python3
"""Does completing the sum fix FDR@0.95? Many loci, exact ground truth.

Separates two things the cluster cannot separate:

  ESTIMATOR error - PIPs normalised over the VISITED configurations instead of
                    all of them (trainer_annot.calculate_pip:38). Fixable by
                    completing the sum; that is joint_trainer.complete_posterior.
  MODEL error     - whatever FDR the EXACT posterior still has, once every
                    configuration is summed. No amount of completion touches it;
                    it is the prior and the likelihood.

For each simulated locus every configuration with |g| <= KMAX is enumerated, so
the exact posterior under FB's own model is known. We then build the partial
set a Concrete proposal would visit, and run the REAL complete_posterior on it.

Pooling across loci at PIP >= 0.95 gives FDR the same way the benchmark does
(evaluate.R:.fdr_power_curve pools variants across fits), so the three columns
are directly comparable:

    partial   -> what FB reports today
    completed -> what it reports after the fix
    exact     -> the ceiling: FB's model with a perfect sum
"""
import itertools
import os
import sys

import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..",
                                                "BEATRICE_annot_sparse")))
torch.set_default_dtype(torch.float64)

from scripts.joint_trainer import complete_posterior          # noqa: E402
from scripts.trainer_annot import calculate_pip, reformat_memo  # noqa: E402

KMAX = 3


def log_w(ind, z, R, n, sigma_sq, p0):
    ind = tuple(sorted(ind))
    base = float(np.sum(np.log(1.0 - p0)))
    if not ind:
        return base
    k = len(ind)
    w = n * sigma_sq
    M = np.eye(k) + w * R[np.ix_(ind, ind)]
    zk = z[list(ind)]
    _, logdet = np.linalg.slogdet(M)
    quad = 0.5 * w * zk @ np.linalg.solve(M, zk)
    lp = base + float(np.sum(np.log(p0[list(ind)]) - np.log(1.0 - p0[list(ind)])))
    return -0.5 * logdet + quad + lp


def one_locus(rng, m, n, sigma_sq, rho, n_causal):
    idx = np.arange(m)
    R = rho ** np.abs(idx[:, None] - idx[None, :]) + 1e-6 * np.eye(m)
    causal = np.sort(rng.choice(m, size=n_causal, replace=False))
    beta = np.zeros(m)
    beta[causal] = rng.normal(0, np.sqrt(sigma_sq), size=n_causal)
    z = np.sqrt(n) * R @ beta + np.linalg.cholesky(R) @ rng.normal(size=m)
    p0 = np.full(m, 1.0 / m)

    logs = {g: log_w(g, z, R, n, sigma_sq, p0)
            for k in range(KMAX + 1)
            for g in itertools.combinations(range(m), k)}
    mx = max(logs.values())
    full = {g: float(np.exp(v - mx)) for g, v in logs.items()}
    pip_exact = calculate_pip(full, m)

    # The Concrete proposal of trainer_annot.abf:420 - binary-concrete draws on
    # an importance vector, truncated to top-K_C and thresholded at eps.
    imp = 1.0 / (1.0 + np.exp(-(np.abs(z) - np.quantile(np.abs(z), 0.75))))
    imp = np.clip(imp, 1e-6, 1 - 1e-6)
    visited = {()}
    for _ in range(300):
        u = rng.uniform(1e-9, 1 - 1e-9, size=m)
        s = 1.0 / (1.0 + np.exp(-((np.log(imp) - np.log(1 - imp)
                                   + np.log(u) - np.log(1 - u)) / 0.5)))
        g = tuple(sorted(int(i) for i in np.argsort(-s)[:6] if s[i] > 0.5))
        if len(g) <= KMAX:
            visited.add(g)
    memo = {g: logs[g] for g in visited if g in logs}

    z_t = torch.tensor(z).reshape(m, 1)
    ld_t = torch.tensor(R)
    S_t = torch.linalg.solve(ld_t, z_t)
    p0_t = torch.tensor(p0)

    mean_memo = reformat_memo(memo, p0_t)
    pip_part = calculate_pip(memo, m)
    W_part = sum(memo.values())

    complete_posterior(memo, mean_memo, S_t, z_t, ld_t, n, sigma_sq, p0_t, m,
                       n_seed=8, width=32, max_rounds=4, tol=1e-3)
    pip_done = calculate_pip(memo, m)

    scale = float(np.exp(mx - mean_memo))
    W_full = sum(full.values()) * scale
    return dict(causal=set(int(c) for c in causal), n_causal=n_causal,
                exact=pip_exact, part=pip_part, done=pip_done,
                d_before=1.0 - W_part / W_full,
                d_after=1.0 - sum(memo.values()) / W_full)


def main():
    rng = np.random.default_rng(2026)
    rows = []
    for i in range(40):
        rows.append(one_locus(
            rng, m=16, n=5000, sigma_sq=0.05,
            rho=float(rng.choice([0.85, 0.9, 0.95])),
            n_causal=int(rng.choice([1, 2, 3]))))
        if (i + 1) % 10 == 0:
            print(f"  ... {i+1}/40 loci", flush=True)

    db = np.array([r["d_before"] for r in rows])
    da = np.array([r["d_after"] for r in rows])
    eb = np.array([np.abs(r["part"] - r["exact"]).max() for r in rows])
    ea = np.array([np.abs(r["done"] - r["exact"]).max() for r in rows])

    print(f"\nomitted mass delta   before {np.median(db):.4f}  "
          f"after {np.median(da):.4f}   (median over {len(rows)} loci)")
    print(f"max |PIP - exact|    before {np.median(eb):.4f}  "
          f"after {np.median(ea):.4f}")
    print(f"bound |err| <= delta holds on {int(np.sum(ea <= da + 1e-9))}"
          f"/{len(rows)} loci")

    print(f"\n{'threshold':>10} {'which':>10} {'nsel':>7} {'TP':>6} "
          f"{'FDR':>8} {'power':>8} {'mass ratio':>11}")
    for thr in (0.90, 0.95, 0.99):
        for key in ("part", "done", "exact"):
            nsel = tp = ncaus = 0
            mass = 0.0
            for r in rows:
                sel = np.where(r[key] >= thr)[0]
                nsel += len(sel)
                tp += sum(1 for j in sel if j in r["causal"])
                ncaus += r["n_causal"]
                mass += float(r[key].sum())
            fdr = (nsel - tp) / nsel if nsel else float("nan")
            lab = {"part": "partial", "done": "completed",
                   "exact": "exact"}[key]
            print(f"{thr:>10.2f} {lab:>10} {nsel:>7} {tp:>6} {fdr:>8.3f} "
                  f"{tp/ncaus:>8.3f} {mass/ncaus:>11.3f}")
        print()


if __name__ == "__main__":
    main()
