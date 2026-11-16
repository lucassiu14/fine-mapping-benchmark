#!/usr/bin/env python3
"""Verify joint_trainer.complete_posterior against exhaustive enumeration.

Runs the REAL torch code path (_abf_at / complete_posterior / reformat_memo /
calculate_pip), not a reimplementation, and checks that completing the visited
set drives the PIPs to the exact posterior computed by enumerating every
configuration. Also cross-validates the torch ABF against an independent numpy
evaluation of the same formula.

  cd BEATRICE_annot_sparse && python3 ../scripts/tuning/test_posterior_completion.py
"""
import itertools
import os
import sys

import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
BEAT = os.path.join(HERE, "..", "..", "BEATRICE_annot_sparse")
sys.path.insert(0, os.path.abspath(BEAT))

torch.set_default_dtype(torch.float64)

from scripts.joint_trainer import complete_posterior, _abf_at  # noqa: E402
from scripts.trainer_annot import calculate_pip, reformat_memo  # noqa: E402

RNG = np.random.default_rng(7)
FAIL = []


def chk(name, cond, extra=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}{(' - ' + extra) if extra else ''}")
    if not cond:
        FAIL.append(name)


def log_w_np(ind, z, R, n, sigma_sq, p0):
    """Independent numpy evaluation of the same ABF + Bernoulli prior."""
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


def main():
    m, n, sigma_sq, KMAX = 16, 5000, 0.05, 3
    idx = np.arange(m)
    R = 0.9 ** np.abs(idx[:, None] - idx[None, :]) + 1e-6 * np.eye(m)
    causal = np.array([4, 9])
    beta = np.zeros(m)
    beta[causal] = RNG.normal(0, np.sqrt(sigma_sq), size=len(causal))
    z_np = np.sqrt(n) * R @ beta + np.linalg.cholesky(R) @ RNG.normal(size=m)
    p0_np = np.full(m, 1.0 / m)

    # --- exact posterior over every configuration with |g| <= KMAX ----------
    allg = [g for k in range(KMAX + 1) for g in itertools.combinations(range(m), k)]
    logs = {g: log_w_np(g, z_np, R, n, sigma_sq, p0_np) for g in allg}
    mx = max(logs.values())
    full = {g: float(np.exp(v - mx)) for g, v in logs.items()}
    pip_exact = calculate_pip(full, m)

    # --- a deliberately incomplete visited set ------------------------------
    order = np.argsort(-np.abs(z_np))
    seedset = [tuple(sorted(int(i) for i in order[:j])) for j in (1, 2, 3)]
    seedset += [(int(order[0]),), (int(order[1]),),
                tuple(sorted((int(order[0]), int(order[2]))))]
    visited = {g for g in seedset if len(g) <= KMAX}
    memo_log = {g: logs[g] for g in visited}

    # torch tensors in the shapes joint_trainer uses: z is (m,1), S = LD^-1 z
    z_t = torch.tensor(z_np).reshape(m, 1)
    ld_t = torch.tensor(R)
    S_t = torch.linalg.solve(ld_t, z_t)
    p0_t = torch.tensor(p0_np)

    # --- cross-validate the torch ABF against numpy -------------------------
    probe = (int(order[0]), int(order[3]))
    tmp = {}
    _abf_at(S_t, z_t, ld_t, tmp, n, sigma_sq, p0_t, probe, 0.0)
    got = float(np.log(tmp[tuple(sorted(probe))]))
    want = log_w_np(probe, z_np, R, n, sigma_sq, p0_np)
    chk("torch ABF == numpy ABF", abs(got - want) < 1e-6,
        f"log w: torch {got:.6f} vs numpy {want:.6f}")

    # --- the PRODUCTION z shape ---------------------------------------------
    # A region dict stores the 1-D Z that `_load_region` returns, while
    # training passes Zc = Z.unsqueeze(1) (joint_trainer.py:381). This test
    # originally built z_t already 2-D, so it could not catch a caller that
    # passed the 1-D form - and that is precisely what shipped: every region
    # fit raised "self must be a matrix" inside _woodbury_pieces. Assert both
    # shapes are accepted and agree exactly.
    tmp1, tmp2 = {}, {}
    _abf_at(S_t, z_t.squeeze(-1), ld_t, tmp1, n, sigma_sq, p0_t, probe, 0.0)
    _abf_at(S_t, z_t, ld_t, tmp2, n, sigma_sq, p0_t, probe, 0.0)
    k = tuple(sorted(probe))
    chk("1-D and 2-D z give the same ABF",
        abs(tmp1[k] - tmp2[k]) <= 1e-12 * max(1.0, abs(tmp2[k])),
        f"1-D {tmp1[k]:.12g} vs 2-D {tmp2[k]:.12g}")

    # --- partial sum, then completion, through the real code ----------------
    mean_memo = reformat_memo(memo_log, p0_t)   # in place: log -> exp(.-m0)
    pip_part = calculate_pip(memo_log, m)
    W_part = sum(memo_log.values())

    hist = complete_posterior(memo_log, mean_memo, S_t,
                              z_t.squeeze(-1),      # the production shape
                              ld_t, n,
                              sigma_sq, p0_t, m,
                              n_seed=8, width=m, max_rounds=6, tol=1e-6)
    pip_done = calculate_pip(memo_log, m)

    # reformat_memo offsets by its OWN max, so put both on the enumeration's
    # scale before comparing total mass.
    scale = float(np.exp(mx - mean_memo))
    recovered = sum(memo_log.values()) / (sum(full.values()) * scale)

    err_before = float(np.abs(pip_part - pip_exact).max())
    err_after = float(np.abs(pip_done - pip_exact).max())

    print(f"\nvisited {len(visited)} configs, enumerated {len(full)}")
    print("completion rounds:")
    for h in hist:
        print(f"   round {h['round']}: +{h['added']:4d} evaluated, "
              f"{h['configs']:4d} total, mass +{h['mass_gained']:.6f}")
    print(f"\nrecovered mass fraction : {recovered:.4f}")
    print(f"max |PIP - exact| before: {err_before:.4f}")
    print(f"max |PIP - exact| after : {err_after:.4f}\n")

    chk("partial sum is biased", err_before > 0.05, f"max error {err_before:.4f}")
    # |PIP_M - PIP_exact| <= delta, the omitted mass fraction. Writing
    # A(j) = A_M(j) + a with 0 <= a <= W - W_M and delta = (W - W_M)/W,
    #   PIP_M - PIP_exact = [delta*A_M - (1-delta)*a] / (W*(1-delta)),
    # which is at most delta*PIP_M (at a = 0) and at least -delta (at a its
    # maximum). So the residual mass BOUNDS the PIP error - that is what makes
    # `mass_gained` a calibration bound rather than a progress message.
    delta = 1.0 - recovered
    chk("PIP error is bounded by the omitted mass", err_after <= delta + 1e-9,
        f"max error {err_after:.6f} <= delta {delta:.6f}")
    chk("completion reaches the exact posterior", err_after < 2e-3,
        f"max error {err_after:.6f}")
    chk("bias is upward at the top variant",
        pip_part[int(np.argmax(pip_part))] >= pip_exact[int(np.argmax(pip_part))],
        f"{pip_part[int(np.argmax(pip_part))]:.4f} vs exact "
        f"{pip_exact[int(np.argmax(pip_part))]:.4f}")
    chk("completion terminated on tolerance", hist[-1]["mass_gained"] < 1e-6
        or len(hist) == 6)
    chk("mass is non-decreasing", all(h["mass_gained"] >= -1e-12 for h in hist))

    print("\nOVERALL:", "all checks passed" if not FAIL else f"FAILURES: {FAIL}")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
