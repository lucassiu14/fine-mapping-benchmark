#!/usr/bin/env python
"""Pinpoint why PolyFun's stage 3 (--compute-h2-bins) divides 0/0.

Stage 3 dies in ldsc_polyfun/regressions.py:348:

    prop = jk.RatioJackknife(cat / tot, ...)
    FloatingPointError: invalid value encountered in divide

with cat = M * coef and tot = sum(cat). "invalid value" (as opposed to
"divide by zero") means 0/0 or NaN/NaN, so either every coefficient is zero
-- the non-negative least squares collapsed -- or coef contains NaN, i.e. a
data problem. Those need opposite fixes, so measure which.

This loads EXACTLY what stage 3 loads, through PolyFun's own
sumstats._read_ld_sumstats, then reports the inputs and solves the regression
both unconstrained and non-negative so the two cases are distinguishable.

  python diag_polyfun_stage3.py --prefix pf_prior_sniff/chol \
      --w-ld-chr 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
"""
import os, sys, argparse
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--prefix", required=True, help="e.g. pf_prior_sniff/chol")
ap.add_argument("--w-ld-chr", required=True)
ap.add_argument("--polyfun-dir", default=os.environ.get("POLYFUN_DIR", ""))
ap.add_argument("--num-chr", type=int, default=22)
a = ap.parse_args()
if not a.polyfun_dir:
    sys.exit("set POLYFUN_DIR or pass --polyfun-dir")
sys.path.insert(0, a.polyfun_dir)

from ldsc_polyfun import sumstats                      # noqa: E402
from polyfun_utils import Logger                       # noqa: E402


class Args:
    pass


# Mirror what polyfun.py's run_ldsc(use_ridge=False) sets up.
args = Args()
args.h2 = "munged.parquet"
args.sumstats = args.h2
args.ref_ld = None
args.w_ld = None
args.ref_ld_chr = a.prefix + "."          # the BINS, as stage 3 uses
args.w_ld_chr = a.w_ld_chr
args.n_blocks = 2
args.M = None
args.not_M_5_50 = True
args.num_chr = a.num_chr
args.anno = None
args.chisq_max = None
args.no_intercept = False
args.two_step = None
args.keep_large = False

log = Logger()
print("loading, exactly as stage 3 does ...", flush=True)
M_annot, w_ld_cname, ref_ld_cnames, df_sumstats, novar_cols = \
    sumstats._read_ld_sumstats(args, log, args.h2)

print()
print("=== M (bin sizes, genome-wide) ===")
M = np.asarray(M_annot, dtype=np.float64).ravel()
print(f"  n_annot      : {len(M)}")
print(f"  min / max    : {M.min():,.0f} / {M.max():,.0f}")
print(f"  zeros        : {int((M == 0).sum())}")
print(f"  non-finite   : {int((~np.isfinite(M)).sum())}")
print(f"  M_tot        : {M.sum():,.0f}")

print()
print("=== design matrix (bin LD-scores on the regression SNPs) ===")
X = np.array(df_sumstats[list(ref_ld_cnames)], dtype=np.float64)
print(f"  shape        : {X.shape}")
print(f"  NaN / inf    : {int(np.isnan(X).sum())} / {int(np.isinf(X).sum())}")
print(f"  column means : min {X.mean(axis=0).min():.4g}  max {X.mean(axis=0).max():.4g}")
print(f"  all-zero cols: {int((np.abs(X).sum(axis=0) == 0).sum())}")
cn = np.linalg.cond(X.T @ X)
print(f"  cond(X'X)    : {cn:.4g}")

chisq = np.asarray(df_sumstats["Z"].values, dtype=np.float64) ** 2
N = np.asarray(df_sumstats["N"].values, dtype=np.float64)
w = np.asarray(df_sumstats[w_ld_cname].values, dtype=np.float64)
print()
print("=== response and weights ===")
print(f"  chi2  : mean {chisq.mean():.4g}  max {chisq.max():.4g}  NaN {int(np.isnan(chisq).sum())}")
print(f"  N     : min {N.min():,.0f}  mean {N.mean():,.0f}  max {N.max():,.0f}")
print(f"  w_ld  : min {w.min():.4g}  mean {w.mean():.4g}  NaN {int(np.isnan(w).sum())}"
      f"  non-positive {int((w <= 0).sum())}")

# The regression LDSC actually solves: (chi2 - 1) ~ N * X  (no weights here;
# the weighting only changes efficiency, not whether the solution collapses).
Nbar = N.mean()
Xs = (X * N[:, None]) / Nbar
y = chisq - 1.0

print()
print("=== unconstrained least squares ===")
coef_ols, *_ = np.linalg.lstsq(Xs, y, rcond=None)
coef_ols = coef_ols / Nbar
print(f"  coef   : min {coef_ols.min():.4g}  max {coef_ols.max():.4g}")
print(f"  negative: {int((coef_ols < 0).sum())} of {len(coef_ols)}")
print(f"  NaN    : {int(np.isnan(coef_ols).sum())}")
print(f"  tot = sum(M*coef) = {float(np.sum(M * coef_ols)):.6g}")

print()
print("=== non-negative least squares (what stage 3 uses, nn=True) ===")
from scipy.optimize import nnls                        # noqa: E402
coef_nn, _ = nnls(Xs, y)
coef_nn = coef_nn / Nbar
print(f"  coef   : min {coef_nn.min():.4g}  max {coef_nn.max():.4g}")
print(f"  exactly zero: {int((coef_nn == 0).sum())} of {len(coef_nn)}")
print(f"  NaN    : {int(np.isnan(coef_nn).sum())}")
tot_nn = float(np.sum(M * coef_nn))
print(f"  tot = sum(M*coef) = {tot_nn:.6g}")

print()
print("=== the two solvers PolyFun can use for the non-negative fit ===")
print("  (LstsqJackknifeSlow.__init__: --nnls-exact -> scipy nnls, otherwise")
print("   Lasso(alpha=1e-100, positive=True, precompute=xtx))")
import warnings
import sklearn
from sklearn.linear_model import Lasso
print(f"  sklearn {sklearn.__version__}")
with warnings.catch_warnings(record=True) as wl:
    warnings.simplefilter("always")
    try:
        xtx = Xs.T.dot(Xs)
        las = Lasso(alpha=1e-100, fit_intercept=False, precompute=xtx,
                    positive=True, max_iter=10000, random_state=0)
        c_las = las.fit(Xs, y).coef_ / Nbar
        print(f"  Lasso  : min {np.nanmin(c_las):.4g} max {np.nanmax(c_las):.4g} "
              f"NaN {int(np.isnan(c_las).sum())}  "
              f"tot {float(np.sum(M * c_las)):.6g}")
    except Exception as e:
        print(f"  Lasso  : RAISED {type(e).__name__}: {str(e)[:140]}")
    for w in wl[:3]:
        print(f"           warning: {w.category.__name__}: {str(w.message)[:120]}")
print(f"  scipy  : min {coef_nn.min():.4g} max {coef_nn.max():.4g} "
      f"NaN {int(np.isnan(coef_nn).sum())}  tot {tot_nn:.6g}")

print()
print("=== verdict ===")
if not np.isfinite(coef_ols).all() or np.isnan(X).any():
    print("  NaN in the inputs or the unconstrained fit -> a DATA problem,")
    print("  not an NNLS collapse.")
elif tot_nn == 0:
    print("  NNLS drove every coefficient to zero, so tot = 0 and cat/tot is 0/0.")
    print("  The non-negativity constraint is the binding problem.")
elif abs(tot_nn) < 1e-12:
    print(f"  tot is {tot_nn:.3g}: not exactly zero but numerically indistinguishable,")
    print("  which divides to NaN all the same.")
else:
    print(f"  tot = {tot_nn:.6g}, which is non-zero -- so the collapse happens inside")
    print("  the JACKKNIFE rather than in the full-sample fit. Stage 3 jackknifes over")
    print("  chromosome sets (evenodd_split=True); at least one delete-block must be")
    print("  producing an all-zero or NaN solution.")
