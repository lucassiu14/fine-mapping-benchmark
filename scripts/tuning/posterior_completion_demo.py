#!/usr/bin/env python3
"""Exact demonstration: FB's PIP estimator is a partial sum, and the omission
inflates PIPs upward. Runs locally in seconds, no cluster, no BEATRICE import.

THE ESTIMATOR, read from source. trainer_annot.calculate_pip:38 forms

    PIP_j = sum_{g in M, j in g} w(g) / sum_{g in M} w(g)

where M is `memo` - the set of configurations ever evaluated - and w(g) is the
exponentiated value stored by abf (trainer_annot.py:420, joint_trainer._abf:140,
generate_credible_sets.abf:183; all three are the same formula). With
S = LD^{-1} z (trainer_annot.py:957), VS = LD[ind,:] LD^{-1} z = z[ind], so

    log w(g) = -1/2 logdet(I + n s^2 R_gg)
               + 1/2 n s^2 z_g' (I + n s^2 R_gg)^{-1} z_g
               + sum_j log prior_j,     prior_j = p0_j if j in g else 1 - p0_j

i.e. the standard approximate Bayes factor times the independent-Bernoulli
configuration prior. The WEIGHTS ARE CORRECT. The defect is the index set: the
true posterior normalises over all 2^m configurations, the code normalises over
the visited ones.

CONSEQUENCE, proved here by exhaustive enumeration. Omitting a configuration
removes it from numerator AND denominator. For a variant j, configurations that
do NOT contain j are exactly the ones that would dilute PIP_j, and an
under-dispersed proposal misses those first. So the bias is UPWARD, and its
size is the omitted posterior mass. That is over-confidence, and it is why
FDR@0.95 has a floor no hyperparameter reached: `sparse_concrete`, `gamma` and
`temp` change WHICH configurations are visited, never that the sum is partial.

THE FIX is to enlarge M by neighbourhood moves (add/delete/swap) before calling
calculate_pip - the moves already exist in generate_credible_sets.py, where they
are used only AFTER the PIPs, to build credible sets. Enlarging M changes no
model, no prior and no likelihood: it adds terms that were always in the
definition of the posterior. This script measures what that recovers.
"""
import itertools
import numpy as np

RNG = np.random.default_rng(20261116)


def make_locus(m=18, n=5000, n_causal=2, sigma_sq=0.05, rho=0.92):
    """An LD block strong enough that competing configurations genuinely exist."""
    # AR(1) latent -> genotype correlation, the usual stand-in for an LD block.
    idx = np.arange(m)
    R = rho ** np.abs(idx[:, None] - idx[None, :])
    R = R + 1e-6 * np.eye(m)

    causal = np.sort(RNG.choice(m, size=n_causal, replace=False))
    beta = np.zeros(m)
    beta[causal] = RNG.normal(0.0, np.sqrt(sigma_sq), size=n_causal)

    # Marginal z under the standard summary-statistic model: z ~ N(sqrt(n) R b, R)
    L = np.linalg.cholesky(R)
    z = np.sqrt(n) * R @ beta + L @ RNG.normal(size=m)
    return z, R, causal, n, sigma_sq


def log_w(ind, z, R, n, sigma_sq, p0):
    """log[ ABF(ind) * p(ind) ] - mirrors generate_credible_sets.abf:183."""
    ind = tuple(sorted(ind))
    base = np.sum(np.log(1.0 - p0))
    if len(ind) == 0:
        return base
    k = len(ind)
    w = n * sigma_sq
    M = np.eye(k) + w * R[np.ix_(ind, ind)]
    zk = z[list(ind)]
    sign, logdet = np.linalg.slogdet(M)
    quad = 0.5 * w * zk @ np.linalg.solve(M, zk)
    logprior = base + np.sum(np.log(p0[list(ind)]) - np.log(1.0 - p0[list(ind)]))
    return -0.5 * logdet + quad + logprior


def pip_from(memo, m):
    """Exactly trainer_annot.calculate_pip:38, on whatever set it is given."""
    tot = sum(memo.values())
    pip = np.zeros(m)
    if tot <= 0:
        return pip
    for g, w in memo.items():
        for j in g:
            pip[j] += w
    return pip / tot


def enumerate_all(z, R, n, sigma_sq, p0, kmax):
    m = len(z)
    logs = {}
    for k in range(kmax + 1):
        for g in itertools.combinations(range(m), k):
            logs[g] = log_w(g, z, R, n, sigma_sq, p0)
    mx = max(logs.values())
    return {g: float(np.exp(v - mx)) for g, v in logs.items()}, mx


def concrete_visits(z, R, n, sigma_sq, p0, kmax, n_draws, K_C, eps, temp=0.5):
    """Mimic the Concrete/Gumbel proposal: binary-concrete draws on an
    importance vector, truncated to the top-K_C and thresholded at eps
    (trainer_annot.abf:420 does exactly this to form `ind`)."""
    m = len(z)
    # A trained head concentrates importance on the marginally strong variants.
    imp = 1.0 / (1.0 + np.exp(-(np.abs(z) - np.quantile(np.abs(z), 0.75))))
    imp = np.clip(imp, 1e-6, 1 - 1e-6)
    visited = set()
    for _ in range(n_draws):
        u = RNG.uniform(1e-9, 1 - 1e-9, size=m)
        logistic = np.log(u) - np.log(1 - u)
        s = 1.0 / (1.0 + np.exp(-((np.log(imp) - np.log(1 - imp) + logistic) / temp)))
        order = np.argsort(-s)[:K_C]
        ind = tuple(sorted(int(i) for i in order if s[i] > eps))
        if len(ind) <= kmax:
            visited.add(ind)
    visited.add(())
    return visited


def complete(memo, z, R, n, sigma_sq, p0, mx, kmax, n_seed, width, max_rounds, tol):
    """Neighbourhood expansion: add / delete / swap from the heaviest
    configurations, exactly the moves in generate_credible_sets.py:241-290."""
    m = len(z)
    history = []
    for rnd in range(max_rounds):
        before = sum(memo.values())
        pip = pip_from(memo, m)
        search = list(np.argsort(-pip)[:width])
        seeds = [g for g, _ in sorted(memo.items(), key=lambda kv: -kv[1])[:n_seed]]
        new = {}
        for g in seeds:
            cand = []
            for i in range(len(g)):                      # delete
                cand.append(tuple(g[:i] + g[i + 1:]))
            if len(g) < kmax:                            # add
                for i in search:
                    if i not in g:
                        cand.append(tuple(sorted(g + (int(i),))))
            for i in range(len(g)):                      # swap
                for j in search:
                    if j not in g:
                        cand.append(tuple(sorted(g[:i] + g[i + 1:] + (int(j),))))
            for c in cand:
                if c not in memo and c not in new and len(c) <= kmax:
                    new[c] = float(np.exp(log_w(c, z, R, n, sigma_sq, p0) - mx))
        memo.update(new)
        after = sum(memo.values())
        gained = (after - before) / after if after > 0 else 0.0
        history.append((rnd + 1, len(new), gained))
        if gained < tol:
            break
    return history


def main():
    KMAX = 3
    z, R, causal, n, sigma_sq = make_locus()
    m = len(z)
    p0 = np.full(m, 1.0 / m)          # published BEATRICE's prior: sum(p0) = 1

    full, mx = enumerate_all(z, R, n, sigma_sq, p0, KMAX)
    pip_exact = pip_from(full, m)

    visited = concrete_visits(z, R, n, sigma_sq, p0, KMAX,
                              n_draws=400, K_C=6, eps=0.5)
    memo = {g: full[g] for g in visited if g in full}
    pip_part = pip_from(memo, m)

    W_full = sum(full.values())
    W_part = sum(memo.values())
    delta = (W_full - W_part) / W_full

    print(f"\nlocus m={m}, causal={list(causal)}, |g|<={KMAX}")
    print(f"configurations: enumerated {len(full)}, visited {len(memo)}")
    print(f"omitted posterior mass  delta = {delta:.4f}")

    top = np.argsort(-pip_part)[:6]
    print(f"\n{'variant':>8} {'exact PIP':>10} {'partial':>9} {'inflation':>10} causal")
    for j in top:
        infl = pip_part[j] - pip_exact[j]
        print(f"{j:>8} {pip_exact[j]:>10.4f} {pip_part[j]:>9.4f} "
              f"{infl:>+10.4f} {'*' if j in causal else ''}")

    hist = complete(memo, z, R, n, sigma_sq, p0, mx, KMAX,
                    n_seed=8, width=12, max_rounds=6, tol=1e-3)
    pip_done = pip_from(memo, m)
    W_done = sum(memo.values())

    print("\ncompletion rounds (round, new configs, mass gained):")
    for r, k, g in hist:
        print(f"   round {r}: +{k:4d} configs, mass +{g:.4f}")
    print(f"\nconfigurations after completion: {len(memo)} / {len(full)}")
    print(f"recovered mass: {W_done / W_full:.4f} of the full posterior "
          f"(was {W_part / W_full:.4f})")

    err_before = np.abs(pip_part - pip_exact).max()
    err_after = np.abs(pip_done - pip_exact).max()
    print(f"\nmax |PIP - exact|  before {err_before:.4f}   after {err_after:.4f}")

    print(f"\n{'variant':>8} {'exact':>8} {'partial':>8} {'completed':>10} causal")
    for j in top:
        print(f"{j:>8} {pip_exact[j]:>8.4f} {pip_part[j]:>8.4f} "
              f"{pip_done[j]:>10.4f} {'*' if j in causal else ''}")

    # The calibration consequence, on this locus: what a PIP>=0.95 call means.
    for name, pp in (("partial", pip_part), ("completed", pip_done)):
        sel = np.where(pp >= 0.95)[0]
        if len(sel):
            tp = sum(1 for j in sel if j in causal)
            print(f"\n{name:>9}: {len(sel)} call(s) at PIP>=0.95, "
                  f"{tp} true -> realised precision {tp/len(sel):.3f}")
        else:
            print(f"\n{name:>9}: no call at PIP>=0.95")


if __name__ == "__main__":
    main()
