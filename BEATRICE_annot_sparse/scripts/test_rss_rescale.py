#!/usr/bin/env python3
"""Checks for `rss_rescale_z` -- the common-residual-scale correction.

Run:  python3 BEATRICE_annot_sparse/scripts/test_rss_rescale.py
"""
import os, sys, math
import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
torch.set_default_dtype(torch.float64)
from scripts.joint_trainer import rss_rescale_z

ok = 0


def check(name, cond, detail=""):
    global ok
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    assert cond, name
    ok += 1


print("rss_rescale_z")
N = 5000
rng = np.random.default_rng(0)
z = torch.tensor(rng.normal(0, 1, 400))
z[7] = 85.15                                    # a strong lead, as at phi=0.6

# 1. 'off' is a no-op, in every spelling
for m in ('off', '', None, 'none'):
    check(f"mode {m!r} is a no-op",
          torch.equal(rss_rescale_z(z, N, m), z))

# 2. 'null' reproduces susieR:::summary_stats_constructor EXACTLY:
#       adj <- (n-1)/(z^2 + n-2); z <- sqrt(adj)*z
adj = (N - 1.0) / (z.numpy() ** 2 + N - 2.0)
susie = np.sqrt(adj) * z.numpy()
got = rss_rescale_z(z, N, 'null').numpy()
check("mode 'null' == susieR's susie_rss adjustment",
      np.allclose(got, susie, rtol=0, atol=1e-13),
      f"max|diff| = {np.abs(got - susie).max():.2e}")

# 3. 'lead' leaves the largest |z| untouched -- no signal is discarded
lead = rss_rescale_z(z, N, 'lead')
i = int(torch.argmax(z.abs()))
check("mode 'lead' preserves max|z|",
      abs(float(lead[i]) - float(z[i])) < 1e-10,
      f"{float(z[i]):.4f} -> {float(lead[i]):.4f}")
# 'lead' references the SMALLEST se (the lead variant's), so every other
# variant's z is raised to that common scale while the lead is held fixed.
# 'null' references the null residual variance, the LARGEST, so it shrinks
# everything. Both are rescalings to one common S; the direction differs
# because the reference differs. What must hold is consistency, checked in (7).
check("mode 'lead' never decreases |z|, lead held fixed",
      bool(torch.all(lead.abs() >= z.abs() - 1e-12)))
# 'null' has its crossover at z^2 = 1, since (n-1)/(z^2+n-2) > 1 exactly when
# z^2 < 1: it shrinks the strong variants and very slightly lifts the weak ones.
nullz = rss_rescale_z(z, N, 'null')
big, sml = z.abs() > 1.0, z.abs() < 1.0
check("mode 'null' shrinks |z| where |z| > 1",
      bool(torch.all(nullz.abs()[big] <= z.abs()[big] + 1e-12)),
      f"n={int(big.sum())}")
check("mode 'null' lifts |z| where |z| < 1 (crossover at z^2=1)",
      bool(torch.all(nullz.abs()[sml] >= z.abs()[sml] - 1e-12)),
      f"n={int(sml.sum())}")

# 4. both shrink: monotone in |z|, and inert as z^2/n -> 0
small = torch.tensor([0.5, 1.0, 2.0])
for m in ('null', 'lead'):
    out = rss_rescale_z(small, 10 ** 8, m)
    check(f"mode {m!r} inert when z^2/n -> 0",
          torch.allclose(out, small, rtol=1e-7),
          f"max rel dev = {float((out / small - 1).abs().max()):.2e}")

# 5. the derivation: for a marginal regression on standardised variables,
#    se_j = 1/sqrt(n - 2 + z_j^2). The correction is z * se_j / se_ref.
se = 1.0 / np.sqrt(N - 2 + z.numpy() ** 2)
ref_lead = 1.0 / math.sqrt(N - 2 + float(z.numpy()[i]) ** 2)
check("'lead' == z * se_j / se[argmax|z|]",
      np.allclose(lead.numpy(), z.numpy() * se / ref_lead, rtol=0, atol=1e-11),
      f"max|diff| = {np.abs(lead.numpy() - z.numpy() * se / ref_lead).max():.2e}")

# 6. degenerate n is refused safely, bad mode is refused loudly
check("n <= 2 returns z unchanged", torch.equal(rss_rescale_z(z, 2, 'lead'), z))
try:
    rss_rescale_z(z, N, 'bogus'); bad = False
except ValueError:
    bad = True
check("an unknown mode raises ValueError", bad)

# 7. THE POINT: on the real scenario-7 locus the standardised conditional
#    residuals must come back to ~1 from 4.61.
SCR = "/private/tmp/claude-502/-Users-ls1020/81b16bdb-fc0f-4706-8fd7-d7709c33fd5d/scratchpad"
zp, ldp, cp = (os.path.join(SCR, f) for f in ("z.txt", "ld.txt", "causal.txt"))
if all(os.path.exists(p) for p in (zp, ldp, cp)):
    zr = np.loadtxt(zp); R = np.loadtxt(ldp)
    c = int(open(cp).read().strip()) - 1
    r = R[c]; m = np.abs(r) < 0.999
    sd = lambda v: float(((v[m] - r[m] * v[c]) / np.sqrt(1 - r[m] ** 2)).std())
    base = sd(zr)
    check("real locus: uncorrected residuals are badly over-dispersed",
          base > 4.0, f"sd = {base:.3f} (model assumes 1.0)")
    for m_ in ('lead', 'null'):
        s = sd(rss_rescale_z(torch.tensor(zr), 5000, m_).numpy())
        check(f"real locus: mode {m_!r} restores residual sd to ~1",
              s < 1.1, f"sd {base:.3f} -> {s:.3f}")
else:
    print("  [SKIP] real-locus check (scratch fixtures absent)")

print(f"\n{ok} checks passed")
