"""
joint_trainer.py  -  cross-region joint prior training for (Functional) BEATRICE.

Why. Functional BEATRICE fits its annotation prior (`LassoNetPrior`) from a
SINGLE locus. That is an underpowered estimation problem: one region supplies
too few causal variants to pin down what is a genome-wide quantity, so the prior
over-fits annotation noise and calibration suffers. This module trains ONE
shared prior head jointly across ALL regions given to it:

    L(phi, {psi_r}) = sum_r  ELBO_r( psi_r ; p0_r = f_phi(v_r) )

  - phi   = the ONE shared prior head (annotations -> per-SNP prior). It sees
            every region's evidence at every optimiser step.
  - psi_r = the region's own variational fine-mapping posterior (its own
            `network`). These stay SEPARATE: which SNPs are causal differs by
            region.
  - only p0 = f_phi(v) is shared; the likelihood (Z_r, LD_r) is per-region.

Each optimiser step accumulates the mean of the per-region losses, then does ONE
backward and ONE step, so phi's gradient is (1/R) sum_r d ELBO_r / d phi. This is
the JOINT scheme, not the sequential `--prior_weights` warm-start chain, which is
order-dependent and forgets earlier regions.

Regions are weighted EQUALLY, not by size: each locus is treated as one
exchangeable draw of the shared enrichment, so a large region does not get a
larger say in the annotation map.

Two swappable prior heads:
  - prior_head='lassonet' -> `LassoNetPrior` from trainer_annot.py
  - prior_head='linear'   -> `LinearPrior` below, a plain logistic map

The per-region ELBO maths (`_region_elbo`, `_abf`) is copied verbatim from
`finemapper_lassonet.train` / `.abf` in trainer_annot.py, so the joint and
single-region paths use identical likelihood, KL and regularisation terms and
differ only in the optimisation.
"""
import math
import os
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.autograd import Variable

from scripts.convert_to_gpu import gpu
from scripts.convert_to_gpu_and_tensor import gpu_t
from scripts.convert_to_gpu_scalar import gpu_ts
from scripts.convert_to_cpu import cpu
from scripts.device import rand_like_hostrng, device_report
from scripts.trainer_annot import (
    network, LassoNetPrior, regularize_ld, calculate_pip, reformat_memo,
    save_object,
)
import scripts.generate_credible_sets as gen_cred


# ---------------------------------------------------------------------------
# A simple SHARED logistic map on annotations.
# Same forward interface as LassoNetPrior (returns ..., imp_o = p0), so the joint
# loop is head-agnostic. No hidden layers, no L1 / hierarchy (proximal ops are
# no-ops), giving a funmap/PolyFun-style global-coefficient prior.
# ---------------------------------------------------------------------------
class LinearPrior(nn.Module):
    def __init__(self, m):
        super().__init__()
        self.m = m
        self.skip = nn.Linear(m, 2, bias=True)   # 2-logit softmax, matches LassoNet
        self.imp = []
        self.feature_importance = None

    def compute_feature_importance(self):
        theta = self.skip.weight.T               # (m, 2)
        contrast = theta[:, 1] - theta[:, 0]     # identifiable logit contrast
        self.feature_importance = torch.abs(contrast).detach().cpu().numpy()
        return torch.abs(contrast)

    # L1 / hierarchy are meaningless for a plain linear head -> no-ops so the joint
    # loop can call them unconditionally.
    def proximal_l1(self, lambda_l1, lr):
        return
    def apply_hierarchy_constraint(self):
        return
    def l1_penalty(self):
        return gpu_ts(0.0)

    def gumbel(self, alpha, t):
        u = (-torch.log(-torch.log(rand_like_hostrng(alpha) + 1e-10) + 1e-10) + alpha) / t
        return torch.nn.functional.softmax(u, dim=1)

    def forward(self, X, T, samples, p0_only=False):
        eps = 1e-7
        out = self.skip(X)                       # (bp, 2)
        # Same identified logit LassoNetPrior stashes, so anchor_prior_mass can
        # renormalise in log space for either head.
        self.last_logit = out[:, 1] - out[:, 0]
        imp = torch.exp(out)
        imp_o = imp[:, 1] / (torch.sum(imp, dim=1) + eps)
        if not p0_only:                          # both lines are device syncs
            self.imp = imp_o.detach().cpu().numpy()
            self.compute_feature_importance()
        if self.training and not p0_only:
            z_N = self.gumbel(torch.log(imp.repeat(samples, 1) + eps), T)
            bin_concrete = z_N[:, 1].reshape(samples, len(imp_o))
            return bin_concrete, None, None, imp_o
        return None, None, None, imp_o


# ---------------------------------------------------------------------------
# Per-region ELBO. Copied from finemapper_lassonet.train / .abf (trainer_annot.py
# lines ~554-668) so JOINT and single-region share identical math. `S` and `memo`
# are per-region (passed in); the shared prior head supplies p_0 = f_phi(v).
# Returns the differentiable region loss; populates `memo` in place (for gen_cred).
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Contracted Woodbury pieces (results.tex eq:contract).
#
# U = n * diag(sigma^2 c) restricted to the K_C columns in `ind`, so U has
# exactly K_C non-zeros: U[ind[k], k] = n * sigma^2 * c[ind[k]] =: w[k]. The
# original code built U as a DENSE (m_r x m_r) diagonal and sliced it, then
# formed sigma_inv = U inv V as a DENSE (m_r x m_r), only to contract it to a
# scalar. That is the 2 m_r^2 K_C term. Every product below is O(m_r K_C) or
# O(K_C^2) instead, and no m_r x m_r matrix is ever allocated.
#
#   V U        = LD[ind][:, ind] * w            (K, K)
#   z^T U      = z[ind]^T * w                   (1, K)
#   z^T Sig^-1 s = z^T s - (z^T U) inv (V s)
#
# Algebraically identical to the dense form; floating-point results differ in
# the last bits because the contraction order changes.
# ---------------------------------------------------------------------------
def _woodbury_pieces(S, z, ld, n_sub, sigma_sq, cc, ind):
    """Return (VU, zU, VS, zS) with no m_r x m_r intermediate."""
    w = n_sub * sigma_sq * cc[ind]            # (K,) the non-zero entries of U
    VU = ld[ind][:, ind] * w.unsqueeze(0)     # (K, K)  == V @ U
    zU = (z[ind].squeeze(-1) * w).unsqueeze(0)  # (1, K) == z^T @ U
    VS = torch.mm(ld[ind, :], S)              # (K, 1)  == V @ S
    zS = torch.mm(z.T, S)                     # (1, 1)  == z^T @ S
    return VU, zU, VS, zS


def _abf(S, z, ld, memo, n_sub, sigma_sq, cc, p0, K_C, eps):
    id_sort = np.argsort(cpu(cc).data.numpy())[::-1][:K_C]
    cc_t = cc[list(id_sort)]
    ind = sorted(id_sort[cpu(torch.where(cc_t > eps)[0]).data.numpy()])
    ind_m = tuple(ind)
    cc = gpu(torch.ones(len(z)))
    if len(ind) > 0:
        if ind_m in memo:
            return memo[ind_m]
        VU, zU, VS, _ = _woodbury_pieces(S, z, ld, n_sub, sigma_sq, cc, ind)
        eyeK = gpu(torch.eye(len(ind)))
        inv = torch.inverse(eyeK + VU)
        # _abf's sigma_inv is U inv V (no identity term), so the quadratic form
        # is (z^T U) inv (V S) with no z^T S.
        sigma2 = torch.mm(torch.mm(zU, inv), VS) / 2
        sigma = eyeK + VU
        prior = 1 - p0
        prior[ind] = p0[ind]
        res = -torch.logdet(sigma) / 2 + sigma2 + torch.sum(torch.log(prior))
        memo[ind_m] = float(res.detach().cpu().numpy().squeeze())
        return res
    return None


# ---------------------------------------------------------------------------
# Posterior completion
#
# `calculate_pip` (trainer_annot.py:38) forms
#
#     PIP_j = sum_{g in M, j in g} w(g) / sum_{g in M} w(g)
#
# over M = `memo`, the set of configurations the Concrete sampler happened to
# visit. The WEIGHTS are correct: with S = LD^{-1} z (trainer_annot.py:957),
# VS = LD[ind,:] LD^{-1} z = z[ind], so w(g) is the standard approximate Bayes
# factor times the independent-Bernoulli configuration prior. The INDEX SET is
# not: the posterior normalises over all 2^m configurations, this sum over the
# visited ones.
#
# Omitting a configuration drops it from the numerator AND the denominator.
#
# CORRECTED 2026-11-17. This comment used to argue that the omitted
# configurations are the ones that would DILUTE PIP_j, so the bias is upward /
# over-confident, and that the FDR@0.95 floor followed from it. The paired
# canary (jobs 4580/4581, 26 scenarios) refuted that: completing the sum raised
# realised mass in 26 of 26 scenarios. The mass M was missing therefore sits on
# LARGER configurations, so truncation biased the PIPs DOWNWARD, and completion
# cannot be the explanation for either the ece_hi=0.232 over-confidence or the
# FDR floor. Completion is still worth having on its own terms - mean dAP
# +0.0442 over those 26 scenarios, 19 wins / 6 losses - but as an accuracy fix,
# not an FDR fix.
#
# What the floor does follow from: sum_j PIP_j = E[|gamma| | gamma in M], and
# the trainer's draws have |gamma| = #{j in top-`sparse_concrete` : c_j >
# `gamma`} (:141-144, :395). The number of variants able to hold a high PIP is
# therefore set by `gamma`, which has been pinned at the Python default 0.1
# since it was first exposed (2026-11-10) and never searched. See the
# FB_SIZE_LOG diagnostic in generate_credible_sets.py, which measures the
# weighted size distribution the reported PIPs are formed from.
#
# Completing the sum changes NO model, NO prior and NO likelihood. It evaluates
# terms that were always in the definition of the posterior, using the same
# `w`. The moves are add / delete / swap - the neighbourhood of
# generate_credible_sets.py:241-290, which today runs only AFTER the PIPs, to
# build credible sets - i.e. FINEMAP's shotgun neighbourhood (Benner et al.
# 2016), applied where it belongs.
#
# `mass_gained` per round is reported so convergence is auditable: when a round
# discovers a negligible share of new mass, M covers the posterior and the PIPs
# are the posterior's, not the sampler's.
# ---------------------------------------------------------------------------
def _abf_at(S, z, ld, memo, n_sub, sigma_sq, p0, ind, mean_memo, _cap=10 ** 15):
    """ABF for an EXPLICIT configuration, on the post-`reformat_memo` scale.

    Identical formula to `_abf`; differs only in taking `ind` directly rather
    than deriving it from a Concrete draw, and in storing exp(res - mean_memo)
    instead of res - which is the scale the memo carries once `reformat_memo`
    has run, and the scale `generate_credible_sets.abf` writes, so the two
    interoperate on one memo.
    """
    key = tuple(sorted(int(i) for i in ind))
    if key in memo:
        return memo[key]
    if len(key) == 0:
        return None
    idx = list(key)
    # `_woodbury_pieces` computes zS = torch.mm(z.T, S) unconditionally, which
    # needs a MATRIX. The training path hands it Zc = Z.unsqueeze(1) (:381),
    # but a region dict stores the 1-D Z that `_load_region` returns, and
    # torch.mm on a 1-D z raises "self must be a matrix" on the first call -
    # failing every region fit. Normalise here so either shape is accepted.
    if z.dim() == 1:
        z = z.unsqueeze(-1)
    cc = gpu(torch.ones(len(z)))
    VU, zU, VS, _ = _woodbury_pieces(S, z, ld, n_sub, sigma_sq, cc, idx)
    eyeK = gpu(torch.eye(len(idx)))
    inv = torch.inverse(eyeK + VU)
    # Same contraction as _abf: sigma_inv is U inv V, so the quadratic form is
    # (z^T U) inv (V S) with no z^T S term.
    sigma2 = torch.mm(torch.mm(zU, inv), VS) / 2
    prior = 1 - p0
    prior[idx] = p0[idx]
    ex = (-torch.logdet(eyeK + VU) / 2 + sigma2
          + torch.sum(torch.log(prior)) - mean_memo)
    val = float(torch.clamp(torch.exp(torch.clamp(ex, max=math.log(_cap))),
                            max=_cap).detach().cpu().numpy().squeeze())
    memo[key] = val
    return val


def complete_posterior(memo, mean_memo, S, z, ld, n_sub, sigma_sq, p0, bp,
                       n_seed=8, width=32, max_rounds=4, tol=1e-3, kmax=None):
    """Enlarge `memo` by add/delete/swap until new mass falls below `tol`.

    Mutates `memo` in place and returns a per-round diagnostic. `kmax` defaults
    to the largest configuration the sampler itself reached, so completion
    finishes the sum over the configuration sizes the model already explores
    rather than widening the model's notion of sparsity.

    `width` restricts add/swap to the top-`width` variants by current PIP. A
    configuration that swaps in a variant with negligible marginal evidence
    carries negligible weight, so this costs little mass - and `mass_gained`
    reports what each round actually found, so the restriction is auditable
    rather than assumed.

    THE RESIDUAL MASS BOUNDS THE PIP ERROR. With W the total weight over all
    configurations, W_M the weight over the visited ones, delta = (W - W_M)/W,
    and A(j) = A_M(j) + a where 0 <= a <= W - W_M,

        PIP_M(j) - PIP_exact(j) = [delta*A_M(j) - (1 - delta)*a]
                                  / (W * (1 - delta))

    which is at most delta*PIP_M(j) (at a = 0, i.e. none of the omitted mass
    contains j) and at least -delta (at a maximal). So

        |PIP_M(j) - PIP_exact(j)| <= delta        for every j.

    The bias is therefore upward wherever the omitted configurations mostly
    exclude j - which is the common case, since an under-dispersed proposal
    misses the competing configurations first. `mass_gained` falling below the
    tolerance is evidence that delta is small, so it is a calibration bound and
    not merely a progress signal. Verified against exhaustive enumeration in
    scripts/tuning/test_posterior_completion.py.
    """
    sizes = [len(g) for g in memo if len(g) > 0]
    if not sizes:
        return []
    if z.dim() == 1:                    # see the note in `_abf_at`
        z = z.unsqueeze(-1)
    if kmax is None:
        kmax = max(sizes)
    history = []
    with torch.no_grad():
        for rnd in range(int(max_rounds)):
            before = float(sum(memo.values()))
            if before <= 0:
                break
            pip = calculate_pip(memo, bp)
            search = [int(i) for i in np.argsort(np.asarray(pip))[::-1][:width]]
            seeds = [g for g, _ in sorted(memo.items(),
                                          key=lambda kv: -kv[1])[:n_seed]]
            added = 0
            for g in seeds:
                cand = []
                for i in range(len(g)):                       # delete
                    cand.append(g[:i] + g[i + 1:])
                if len(g) < kmax:                             # add
                    for i in search:
                        if i not in g:
                            cand.append(tuple(sorted(g + (i,))))
                for i in range(len(g)):                       # swap
                    for j in search:
                        if j not in g:
                            cand.append(tuple(sorted(g[:i] + g[i + 1:] + (j,))))
                for c in cand:
                    if c and len(c) <= kmax and c not in memo:
                        _abf_at(S, z, ld, memo, n_sub, sigma_sq, p0, c,
                                mean_memo)
                        added += 1
            after = float(sum(memo.values()))
            gained = (after - before) / after if after > 0 else 0.0
            history.append({"round": rnd + 1, "added": added,
                            "configs": len(memo), "mass_gained": gained,
                            "mass_before": before, "mass_after": after})
            if gained < tol:
                break
    return history


def anchor_prior_mass(p0, mass, logit=None, eps=1e-7):
    """Rescale a region's prior so that sum(p0) == `mass`, differentiably.

    Why this exists. Published BEATRICE (sayangsep/Beatrice-Finemapping,
    scripts/trainer.py:386) does not learn the prior at all:

        p_0 = gpu_t(np.array([1/len(Z)]*len(Z)))

    i.e. p0_j = 1/m for every variant, so sum(p0) = 1 EXACTLY - a prior
    expectation of one causal variant per locus. `n_caus` does not change this;
    it only caps credible-set size in generate_credible_sets.py.

    Our LassoNetPrior instead emits p0_j = sigmoid(theta'x_j + f(x_j)) per
    variant, which fixes the SHAPE of the prior from annotations but leaves its
    SCALE unidentified: nothing in the ELBO ties sum(p0) to anything. The only
    term that touches the scale is lambda_reg * sum(p0**2), which pulls it
    toward 0 rather than toward a calibrated value. Measured on the v4 tuning
    grid, sum(PIP)/#causal ran 2.06 - 35.26 with a median of 5.73, and because
    `abf` uses a Bernoulli configuration prior,

        log p(config) = sum_{j in ind} log p0_j + sum_{j not in ind} log(1-p0_j)

    the posterior expected configuration size - hence sum_j PIP_j - is governed
    directly by sum_j p0_j. An unanchored scale therefore inflates every PIP,
    which is what drove worst-stratum FDR to 0.81.

    Renormalising restores the published scale while KEEPING the annotation
    head's contribution: the head now redistributes a fixed budget of prior
    mass across the locus instead of setting its total. With a flat logit this
    returns p0_j = mass/bp, which is published BEATRICE's 1/m when mass = 1, so
    mass = 1.0 is the faithful setting and larger values are a deliberate,
    documented relaxation for high-S strata.

    `mass <= 0` disables the anchor and returns p0 unchanged, which is the
    pre-existing (unanchored) behaviour.
    """
    if mass is None or mass <= 0:
        return p0
    if logit is not None:
        # The numerically sound form. softmax is invariant to a constant shift of
        # the logit - which is exactly the direction normalisation leaves
        # unidentified and therefore free to drift - and it subtracts the max
        # internally, so neither tail underflows.
        scaled = mass * torch.softmax(logit, dim=0)
    else:
        # Fallback for callers with no logit to hand. Renormalising the
        # probabilities is algebraically identical but NOT stable: see the note
        # above, and LassoNetPrior.forward's. Kept only so the function is total.
        tot = torch.sum(p0)
        tot_v = float(tot.detach())   # detach: `tot` carries requires_grad here
        if not (tot_v > 0.0) or tot_v != tot_v or tot_v == float("inf"):
            return p0
        scaled = p0 * (mass / tot)
    # log(1-p0) must stay finite, so cap just below 1. The cap only binds when a
    # single variant would take essentially the whole budget.
    return torch.clamp(scaled, eps, 1.0 - 1e-6)


def _region_elbo(model, p_0, S, Z, LD, temp, n_samples, sigma_sq,
                 n_sub, K_C, gamma, lambda_reg, memo):
    """One region's differentiable loss (lik + kl + reg), populating `memo`.
    Mirrors finemapper_lassonet.train's inner body; NO optimiser step here - the
    joint loop sums these across regions and steps ONCE.

    `p_0` is this region's slice of the SHARED prior head's output. The head is
    now evaluated ONCE per optimiser step on every region's annotations
    concatenated (see run_joint), rather than once per region: LassoNetPrior is
    purely per-variant - every reduction in it is over the 2-logit axis, none
    over the variant axis - so concatenating is exact and needs no padding."""
    eps = gpu_ts(1e-7)
    M = len(Z)
    c, c1, c2, imp = model(temp, n_samples)          # region posterior psi_r

    loss = gpu_ts(0.0)
    lik_loss = gpu_ts(0.0)
    Zc = Z.unsqueeze(1)
    valid_samples = 0
    best_u, best_ind = None, None

    for i in range(n_samples):
        u, ind = torch.topk(c[i], max(1, K_C))
        if ind.numel() == 0:
            continue
        valid_samples += 1
        K_C_eff = ind.numel()
        cc = c[i]
        if gamma > 0:
            # _abf's return value is DISCARDED - it runs purely to populate
            # `memo`, which is read after training by reformat_memo /
            # calculate_pip. It therefore contributes nothing to the gradient,
            # yet without no_grad it builds an autograd graph on every one of
            # its ~15,000 calls per scenario. The memo entry is already
            # detached, so this is numerically identical.
            with torch.no_grad():
                _abf(S, Zc, LD, memo, n_sub, sigma_sq, cc, p_0, K_C_eff, gamma)
        # Contracted Woodbury: no m_r x m_r matrix is formed. sigma_inv here
        # carries the identity term, so the quadratic form is
        #   z^T (I - U inv V) s = z^T s - (z^T U) inv (V s).
        VU, zU, VS, zS = _woodbury_pieces(S, Zc, LD, n_sub, sigma_sq, cc, ind)
        eyeKC = gpu(torch.eye(K_C_eff))
        inv = torch.inverse(eyeKC + VU)
        quad = zS - torch.mm(torch.mm(zU, inv), VS)
        sigma = eyeKC + VU + 1e-6 * eyeKC
        sigma2 = -quad / 2
        log_likelihood = -torch.logdet(sigma) / 2 + sigma2
        if not torch.isfinite(log_likelihood):
            continue
        lik_loss += -log_likelihood.squeeze()
        loss += -log_likelihood.squeeze()
        if (best_u is None) or (u.max() > best_u.max()):
            best_u, best_ind = u, ind

    if valid_samples == 0:
        best_ind = torch.topk(c[0], 1).indices
        valid_samples = 1

    x2 = imp[best_ind]
    x1 = p_0[best_ind]
    kl_loss = torch.sum(x2 * (torch.log(x2 + eps) - torch.log(x1 + eps))) + \
              torch.sum((1 - x2) * (torch.log(1 - x2 + eps) - torch.log(1 - x1 + eps)))
    reg_loss = lambda_reg * torch.sum(p_0 ** 2)

    region_loss = (loss / valid_samples) + kl_loss + reg_loss
    return region_loss, p_0


# ---------------------------------------------------------------------------
# The joint driver.
# ---------------------------------------------------------------------------
def _load_region(z_path, ld_path, annot_path):
    names = list(pd_read_col0(z_path))
    Z = gpu_t(pd_read(z_path)[:, 1].astype(float))
    if torch.max(Z) == torch.inf:
        Z[torch.where(Z == torch.inf)[0]] = 200
    LD = gpu(torch.from_numpy(read_ld(ld_path)))
    v = gpu_t(pd_read_annot(annot_path)[:, 1:].astype(float))
    if torch.max(v) == torch.inf:
        v[torch.where(v == torch.inf)[0]] = 200
    return names, Z, regularize_ld(LD), v


# small pandas readers kept local to avoid re-importing in hot code
import pandas as pd
def pd_read(path):
    return pd.read_table(path, sep=' ', header=None).to_numpy()
def pd_read_col0(path):
    return pd.read_table(path, sep=' ', header=None).to_numpy()[:, 0]
def pd_read_annot(path):
    return pd.read_table(path, sep=None, engine='python', header=None).to_numpy()


def read_ld(path):
    """Load an LD matrix, accepting .npy as well as whitespace text.

    The text path is the original format and still works. It is, however,
    unusable at scale: a block of p = 10,000 is an ~800 MB text file, and a
    genome-wide run over ~250 blocks would spend longer parsing ASCII than
    fitting the model. A float32 .npy of the same matrix is ~400 MB, memory-
    maps, and loads in milliseconds. Values are identical; only the encoding
    differs.

    Returned as float32 and C-contiguous so the caller can wrap it with
    torch.from_numpy directly. gpu_t() would be wrong here: it does
    `torch.tensor(x.astype(float)).float()`, round-tripping through float64,
    which transiently costs 4x the matrix for no benefit - about 1.6 GB for a
    single p = 10,000 block.
    """
    a = np.load(path) if str(path).endswith('.npy') else pd_read(path)
    return np.ascontiguousarray(a, dtype=np.float32)


def run_joint(options):
    """
    options keys:
      regions            : list of dicts {z, LD, annot, target, N}
      prior_head         : 'lassonet' | 'linear'
      NN                 : finemapping-net hidden dims  (default [10,200,10])
      prior_neural_network : LassoNet hidden dims       (default [10,200,10])
      max_iter, sigma_sq, sparse_concrete, MCMC_samples, temp_lower_bound,
      gamma, prior_regularisation, lambda_l1, hierarchy_M,
      n_causal, coverage_ths, selection_prob, key_thres, allow_duplicates,
      save_prior_weights (path or '')
    """
    torch.manual_seed(1)
    import random; random.seed(1)

    regions_cfg = options['regions']
    R = len(regions_cfg)
    if R == 0:
        raise ValueError("run_joint: no regions supplied")

    sigma_sq = options['sigma_sq']
    n_samples = options['MCMC_samples']
    n_epochs = options['max_iter']
    temp_lb = gpu_ts(options['temp_lower_bound'])
    gamma = options['gamma']
    lambda_reg = options['prior_regularisation']
    lambda_l1 = options.get('lambda_l1', 0.01)
    hierarchy_M = options.get('hierarchy_M', 10.0)
    NN = options['NN']
    prior_hidden = options['prior_neural_network']

    # --- load every region + build region-specific state --------------------
    reg = []
    m_annots = None
    for rc in regions_cfg:
        names, Z, LD, v = _load_region(rc['z'], rc['LD'], rc['annot'])
        if m_annots is None:
            m_annots = v.shape[1]
        elif v.shape[1] != m_annots:
            raise ValueError(f"annotation count differs across regions: "
                             f"{v.shape[1]} vs {m_annots}")
        model = gpu(network(len(Z), [1] + NN, 3, LD, Z))
        # solve(A, b) rather than inverse(A) @ b: same result, ~3x fewer flops
        # and better conditioned.
        S = torch.linalg.solve(LD, Z.unsqueeze(1))
        reg.append(dict(names=names, Z=Z, LD=LD, v=v, model=model, S=S,
                        memo={}, bp=len(Z), N=rc['N'], target=rc['target'],
                        z=rc['z'], LD_path=rc['LD']))

    # --- annotations of every region, concatenated ONCE ----------------------
    # The shared head is applied per variant, so evaluating it on the stacked
    # (sum_r m_r, d) matrix and slicing gives exactly the same p_0 as calling it
    # per region, with one dispatch per step instead of R. `offs` holds the
    # split points.
    v_all = torch.cat([r['v'] for r in reg], dim=0)
    offs, _acc = [], 0
    for r in reg:
        offs.append((_acc, _acc + r['bp'])); _acc += r['bp']

    # --- ONE shared prior head phi ------------------------------------------
    if options['prior_head'] == 'linear':
        prior_head = gpu(LinearPrior(m=m_annots))
        head_desc = "LinearPrior (shared logistic annotation prior)"
    else:
        prior_head = gpu(LassoNetPrior(m=m_annots, hidden_dims=prior_hidden,
                                       M=hierarchy_M))
        head_desc = (f"LassoNetPrior (shared, hidden={prior_hidden}, "
                     f"M={hierarchy_M}, one logit per variant)")
    print(f"[joint] {R} regions, {m_annots} annotations, head={head_desc}")
    print(f"[joint] {device_report()}")

    # --- ONE optimiser over all region posteriors + the shared prior --------
    params = list(prior_head.parameters())
    for r in reg:
        params += list(r['model'].parameters())
    opt = optim.Adam(params, lr=0.002, betas=(0.9, 0.999), weight_decay=0)
    scheduler = torch.optim.lr_scheduler.StepLR(opt, step_size=1000, gamma=0.5)

    prior_head.train()
    for r in reg:
        r['model'].train()

    # --- JOINT training loop: equal-weight mean over regions, ONE step ------
    from tqdm.auto import tqdm
    for n in tqdm(range(n_epochs + 1), desc="joint"):
        temp = torch.max(temp_lb, gpu_ts(np.exp(-0.0001 * n)))
        opt.zero_grad()
        total = gpu_ts(0.0)
        # ONE evaluation of the shared head for the whole scenario, then slice.
        _, _, _, p0_all = prior_head(v_all, temp, 1, p0_only=True)
        # The head sees every region's annotations concatenated, so the
        # logit must be sliced per region exactly as p0_all is: the
        # anchor is a PER-REGION budget, and a softmax over the whole
        # concatenation would instead split one budget across regions.
        logit_all = getattr(prior_head, 'last_logit', None)
        for i, r in enumerate(reg):
            K_C = min(r['bp'], options['sparse_concrete'])
            lo, hi = offs[i]
            region_loss, _ = _region_elbo(
                r['model'],
                anchor_prior_mass(p0_all[lo:hi], options.get('prior_mass', 0.0),
                                  logit=(logit_all[lo:hi]
                                         if logit_all is not None else None)),
                r['S'], r['Z'], r['LD'],
                temp, n_samples, gpu_ts(sigma_sq), r['N'], K_C, gamma,
                lambda_reg, r['memo'])
            total = total + region_loss / R          # EQUAL weight per region
        total.backward()
        opt.step()
        scheduler.step()
        # proximal L1 + hierarchy on the SHARED head, once per step (no-ops for linear)
        prior_head.proximal_l1(lambda_l1, opt.param_groups[0]['lr'])
        prior_head.apply_hierarchy_constraint()

    # --- per-region outputs: save res + generate credible sets --------------
    prior_head.eval()
    for r in reg:
        with torch.no_grad():
            _, _, _, p_0 = prior_head(r['v'], temp, 1)
            # This call is already per region, so last_logit needs no slicing.
            p_0 = anchor_prior_mass(p_0, options.get('prior_mass', 0.0),
                                    logit=getattr(prior_head, 'last_logit', None))
        mean_memo = reformat_memo(r['memo'], p_0)
        # Finish the partial sum BEFORE reading PIPs off it. `calculate_pip`
        # normalises over the visited configurations, which biases PIPs upward
        # by the omitted posterior mass; see complete_posterior above. Tolerance
        # 0 disables it, so the effect is measurable against the memo the
        # sampler alone produced.
        completion = []
        _tol = float(options.get('posterior_completion', 0) or 0)
        if _tol > 0:
            completion = complete_posterior(
                # Zc, not Z: `_woodbury_pieces` needs a matrix, and training
                # passes Z.unsqueeze(1) for exactly this reason (:381).
                r['memo'], mean_memo, r['S'], r['Z'].unsqueeze(1),
                r['LD'], r['N'],
                gpu_ts(sigma_sq), p_0, r['bp'],
                n_seed=int(options.get('completion_seeds', 8)),
                width=int(options.get('completion_width', 32)),
                max_rounds=int(options.get('completion_rounds', 4)),
                tol=_tol)
        if completion:
            # The ONLY stdout trace that completion ran. `res` is pickled into
            # a wrapper tempdir that R deletes, so without this line a silent
            # no-op is indistinguishable from "ran, but changed nothing" - the
            # ambiguity that wasted an earlier canary. delta_hat is the share
            # of the FINAL posterior mass that the sampler alone had MISSED,
            # i.e. the measured, production value of the delta that bounds the
            # PIP error above. last_gain is what the final round still found:
            # small means the enlarged sum has converged, so delta_hat is a
            # credible estimate of the omitted mass rather than a floor.
            _w0 = completion[0]['mass_before']
            _w1 = completion[-1]['mass_after']
            _d = (_w1 - _w0) / _w1 if _w1 > 0 else 0.0
            print(f"[completion] region={os.path.basename(r['target'])} "
                  f"rounds={len(completion)} "
                  f"configs={completion[-1]['configs']} "
                  f"delta_hat={_d:.4f} "
                  f"last_gain={completion[-1]['mass_gained']:.2e}",
                  flush=True)
        res = {'memo': r['memo'], 'mean_memo': mean_memo,
               'names': r['names'], 'pip': calculate_pip(r['memo'], r['bp']),
               'completion': completion,
               'feature_importance': prior_head.feature_importance}
        os.makedirs(r['target'], exist_ok=True)
        save_object(res, os.path.join(r['target'], 'res'))

        # Annotation importance for the SHARED head. Identical for every region
        # by construction (one head across all R), but written into each
        # region's target so the per-region wrapper parser picks it up without
        # needing to know about the joint layout. Same schema as the
        # single-locus path in trainer_annot.py so downstream analysis code
        # (wrapper -> extract_aux.R -> iter004_annotation_recovery.R) is
        # unchanged. Before this, the joint head's importances were computed and
        # then discarded with the wrapper's tempfile work_dir.
        if prior_head.feature_importance is not None:
            fi_df = pd.DataFrame({
                'annotation_index': range(len(prior_head.feature_importance)),
                'importance': prior_head.feature_importance
            })
            fi_df = fi_df.sort_values('importance', ascending=False)
            fi_df.to_csv(os.path.join(r['target'], 'feature_importance.csv'),
                         index=False)
        gen_opts = {
            'z': r['z'], 'LD': r['LD_path'], 'prior_location': '',
            'sigma_sq': sigma_sq, 'n_sub': r['N'], 'target': r['target'],
            'names': r['names'], 'key_thres': options['key_thres'],
            # hand over what we already hold: saves a second text parse of the
            # LD matrix, a second eigendecomposition and a second m_r x m_r
            # solve, per region.
            'Z_t': r['Z'], 'LD_t': r['LD'], 'S_t': r['S'],
            # The SAME p_0 the training memo was scored and normalised
            # with. Without it generate_credible_sets falls back to a
            # uniform prior and the annotations never reach the PIPs.
            'p0_t': p_0,
            # gen_cred recomputes the PIPs after its own sweep and writes
            # pip.csv, which is the only PIP anything downstream reads, so it
            # needs these to complete the sum at that point too.
            'posterior_completion': options.get('posterior_completion', 0),
            'completion_seeds': options.get('completion_seeds', 8),
            'completion_width': options.get('completion_width', 32),
            'completion_rounds': options.get('completion_rounds', 4),
            'n_causal': options['n_causal'], 'allow_duplicates': options['allow_duplicates'],
            'coverage_ths': options['coverage_ths'], 'selection_prob': options['selection_prob'],
            'purity': options.get('purity', 0.0),
        }
        gen_cred.main(gen_opts)

    # --- save the shared prior head (for inspection / reuse) ----------------
    if options.get('save_prior_weights'):
        ckpt = {'state_dict': prior_head.state_dict(),
                'n_annotations': m_annots, 'hidden_dims': prior_hidden,
                'hierarchy_M': hierarchy_M, 'prior_head': options['prior_head'],
                'feature_importance': prior_head.feature_importance}
        torch.save(ckpt, options['save_prior_weights'])
        print(f"[joint] saved shared prior head -> {options['save_prior_weights']}")
