#!/usr/bin/env python
"""Stacked genome-wide figure: the GWAS on top, then each method's PIPs.

  python 08_plot_pips.py                      # genome-wide, all 254 regions
  python 08_plot_pips.py --region ldEUR_19_0031   # one locus, bp on the x-axis

Reads results/pips.parquet (one row per region x variant x method) and writes
results/fig_pips_genome.{pdf,png} or results/fig_pips_<region>.{pdf,png}.

Panels share an x-axis so a peak in the GWAS can be read straight down into
what each method did with it. The GWAS panel is -log10(p) recovered from Z via
norm.logsf, which stays finite where p underflows (the strongest locus here is
|Z| = 82.9, i.e. p ~ 1e-1500); it is capped for display and the cap is stated
on the panel.

Colours are Okabe-Ito, which is colourblind-safe: each method keeps its colour
across every figure, and the GWAS panel stays neutral so the methods' colours
mean "method" and nothing else.
"""
import os, sys, argparse
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.stats import norm

ROOT = os.environ.get("APP_ROOT", os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")

ap = argparse.ArgumentParser()
ap.add_argument("--region", default=None, help="plot one region instead of the genome")
ap.add_argument("--pip-floor", type=float, default=0.01,
                help="below this, PIPs are thinned for legibility (genome view only)")
ap.add_argument("--cap", type=float, default=60.0, help="-log10(p) display cap")
ap.add_argument("--provisional", default="polyfun_sniff,polyfun_combined",
                help="methods to mark provisional in the panel label")
a = ap.parse_args()

METHODS = ["fb_xregion", "funmap", "polyfun_sniff", "polyfun_combined"]
COLOUR = {"fb_xregion": "#0072B2", "funmap": "#D55E00",
          "polyfun_sniff": "#009E73", "polyfun_combined": "#CC79A7"}
LABEL = {"fb_xregion": "fb_xregion", "funmap": "funmap",
         "polyfun_sniff": "PolyFun (Borzoi only)",
         "polyfun_combined": "PolyFun (+ baselineLF)"}
PROV = set(x for x in a.provisional.split(",") if x)

df = pd.read_parquet(os.path.join(RES, "pips.parquet"))
have = [m for m in METHODS if m in set(df.method)]
if a.region:
    df = df[df.region_id == a.region]
    if df.empty:
        sys.exit(f"no rows for region {a.region}")

# -log10(p) from Z, finite where p underflows
one = df[df.method == have[0]][["region_id", "chrom", "BP", "SNP", "Z"]].copy()
one["logp"] = -(np.log10(2.0) + norm.logsf(np.abs(one.Z.values)) / np.log(10.0))

if a.region:
    one["x"] = one.BP.values
    xlabel = f"chr{int(one.chrom.iloc[0])} position (Mb)"
    bounds, centres, names = [], [], []
else:
    # concatenate chromosomes; offsets from the data actually plotted
    one = one.sort_values(["chrom", "BP"]).reset_index(drop=True)
    off, cum, bounds, centres, names = {}, 0.0, [], [], []
    for c in sorted(one.chrom.unique()):
        bp = one.loc[one.chrom == c, "BP"]
        off[c] = cum - bp.min()
        bounds.append((cum, cum + (bp.max() - bp.min())))
        centres.append(cum + (bp.max() - bp.min()) / 2.0)
        names.append(str(int(c)))
        cum += (bp.max() - bp.min()) + 2e7          # a visible gap between chromosomes
    one["x"] = one.BP.values + one.chrom.map(off).values
    xlabel = "chromosome"

xof = dict(zip(zip(one.region_id, one.SNP), one.x))

fig, axes = plt.subplots(len(have) + 1, 1, figsize=(13, 2.1 * (len(have) + 1)),
                         sharex=True, gridspec_kw=dict(hspace=0.16))
fig.patch.set_facecolor("white")

def shade(ax):
    for i, (lo, hi) in enumerate(bounds):
        if i % 2 == 0:
            ax.axvspan(lo, hi, color="#000000", alpha=0.035, lw=0, zorder=0)

# ---- panel 1: the GWAS -------------------------------------------------
ax = axes[0]
shade(ax)
capped = one.logp.clip(upper=a.cap)
n_cap = int((one.logp > a.cap).sum())
ax.scatter(one.x, capped, s=3, c="#444444", alpha=0.45, lw=0, rasterized=True)
ax.axhline(-np.log10(5e-8), color="#D55E00", lw=0.9, ls="--", zorder=3)
ax.text(0.004, 0.86, "GWAS  $-\\log_{10}p$" + (f"   (capped at {a.cap:.0f}; "
        f"{n_cap:,} variants above)" if n_cap else ""),
        transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")
ax.set_ylim(0, a.cap * 1.05)
ax.set_ylabel("$-\\log_{10}p$", fontsize=9)

# ---- one panel per method ----------------------------------------------
for ax, m in zip(axes[1:], have):
    shade(ax)
    d = df[df.method == m][["region_id", "SNP", "pip"]].copy()
    d["x"] = [xof.get(k, np.nan) for k in zip(d.region_id, d.SNP)]
    d = d.dropna(subset=["x"])
    if a.region:
        lo, hi = d, d.iloc[0:0]
    else:                      # thin the near-zero mass; keep every real signal
        hi = d[d.pip >= a.pip_floor]
        lo = d[d.pip < a.pip_floor].sample(frac=0.08, random_state=0)
    ax.scatter(lo.x, lo.pip, s=2, c=COLOUR[m], alpha=0.18, lw=0, rasterized=True)
    if len(hi):
        ax.scatter(hi.x, hi.pip, s=9, c=COLOUR[m], alpha=0.85, lw=0, rasterized=True)
    ax.axhline(0.95, color="#666666", lw=0.7, ls=":", zorder=3)
    ax.axhline(0.50, color="#999999", lw=0.7, ls=":", zorder=3)
    n95 = int((d.pip >= 0.95).sum()); n50 = int((d.pip >= 0.5).sum())
    tag = "   [provisional: 10-bin priors]" if m in PROV else ""
    ax.text(0.004, 0.88, f"{LABEL[m]}   PIP≥0.95: {n95:,}   ≥0.5: {n50:,}{tag}",
            transform=ax.transAxes, fontsize=9, fontweight="bold",
            color=COLOUR[m], va="top")
    ax.set_ylim(-0.03, 1.05)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_ylabel("PIP", fontsize=9)

for ax in axes:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.margins(x=0.004)
if a.region:
    axes[-1].xaxis.set_major_formatter(lambda v, p: f"{v/1e6:.2f}")
    axes[-1].xaxis.set_major_locator(MaxNLocator(8))
else:
    axes[-1].set_xticks(centres); axes[-1].set_xticklabels(names, fontsize=7)
axes[-1].set_xlabel(xlabel, fontsize=9)

title = (f"Cholesterol (UKBB) — {one.region_id.nunique()} LDetect blocks"
         if not a.region else f"{a.region}")
fig.suptitle(title, fontsize=11, fontweight="bold", x=0.5, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.985))

stem = f"fig_pips_{'genome' if not a.region else a.region}"
for ext in ("pdf", "png"):
    p = os.path.join(RES, f"{stem}.{ext}")
    fig.savefig(p, dpi=170, bbox_inches="tight", facecolor="white")
    print(f"wrote {p}")
