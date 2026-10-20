"""z-scores and PIPs for the three methods, two UKB cholesterol blocks."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

T = "/private/tmp/claude-502/-Users-ls1020/15a23a3d-176d-4bd7-84b8-4955967152e5/scratchpad/hpc_test"
OUT = "/private/tmp/claude-502/-Users-ls1020/15a23a3d-176d-4bd7-84b8-4955967152e5/scratchpad/app_pips.pdf"

# validated categorical palette, assigned in FIXED ORDER (never cycled)
C_FB, C_FM, C_PF = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, BG = "#0b0b0b", "#52514e", "#d8d7d2", "#fcfcfb"

REG = [("ldEUR_20_0038", "chr20", "chr20:62.19-62.97 Mb"),
       ("ldEUR_22_0024_part", "chr22", "chr22:50.00-51.24 Mb")]

def load(rid):
    s = pd.read_csv(f"{T}/regions/{rid}/snps.csv")
    fb = pd.read_csv(f"{T}/joint_out/{rid}/pip.csv")["pip"].to_numpy(float)
    fm = pd.read_csv(f"{T}/funmap_out/{rid}/pip.csv", header=None,
                     names=["SNP", "pip"])["pip"].to_numpy(float)
    pf_t = pd.read_csv(f"{T}/polyfun_out/{rid}.txt.gz", sep="\t")
    pf = s.SNP.map(dict(zip(pf_t.SNP, pf_t.PIP))).fillna(0).to_numpy(float)
    return s, fb, fm, pf

fig, axes = plt.subplots(4, 2, figsize=(13.5, 10.2), sharey="row",
                         gridspec_kw=dict(hspace=0.22, wspace=0.07, left=0.075,
                                          right=0.985, top=0.885, bottom=0.085))
fig.patch.set_facecolor(BG)

def style(ax, ylab=None):
    ax.set_facecolor(BG)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5, length=3)
    ax.grid(axis="y", color=GRID, lw=0.6, alpha=0.7); ax.set_axisbelow(True)
    if ylab: ax.set_ylabel(ylab, color=INK, fontsize=10)

for col, (rid, chrom, title) in enumerate(REG):
    s, fb, fm, pf = load(rid)
    x = s.BP.values / 1e6
    # --- z-scores -----------------------------------------------------------
    ax = axes[0, col]
    ax.scatter(x, s.Z.values, s=4, c=INK2, lw=0, alpha=0.55, rasterized=True)
    for t in (5.451, -5.451):
        ax.axhline(t, color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=1)
    style(ax, "GWAS z-score" if col == 0 else None)
    zm = np.abs(s.Z.values).max()
    ax.set_ylim(-zm * 1.18, zm * 1.18)
    if col == 0:
        ax.text(x.min(), 5.8, "p = 5e-8", fontsize=7.5, color=INK2, va="bottom")
    # --- PIP panels ---------------------------------------------------------
    for row, (v, c, lab) in enumerate(
            [(fb, C_FB, "fb_xregion"), (fm, C_FM, "funmap"), (pf, C_PF, "polyfun+SuSiE")], start=1):
        ax = axes[row, col]
        ax.scatter(x, v, s=4 + 110 * v ** 1.5, c=c, lw=0, alpha=0.85, rasterized=True)
        ax.axhline(0.5, color=INK2, lw=0.7, ls=(0, (4, 3)), zorder=1)
        style(ax, f"{lab}\nPIP" if col == 0 else None)
        ax.set_ylim(-0.05, 1.09)
        # direct-label the top variant: identity is never carried by colour alone
        i = int(np.argmax(v))
        if v[i] > 0.15:
            ax.annotate(f"{v[i]:.2f}", (x[i], v[i]), textcoords="offset points",
                        xytext=(0, 7), ha="center", fontsize=8, color=INK)
        ax.text(0.985, 0.88, f"Σ={v.sum():.2f}", transform=ax.transAxes,
                ha="right", va="top", fontsize=8, color=INK2)
    axes[3, col].set_xlabel(f"position on {chrom} (Mb)", color=INK, fontsize=9.5)
    axes[0, col].set_title(f"{title}    n = {len(s):,} variants",
                           color=INK, fontsize=10.5, pad=8)

fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=6, mfc=c, mec="none", label=l)
                    for c, l in [(C_FB, "fb_xregion"), (C_FM, "funmap"),
                                 (C_PF, "polyfun + SuSiE")]],
           loc="lower center", ncol=3, frameon=False, fontsize=9,
           bbox_to_anchor=(0.5, -0.002), labelcolor=INK)
fig.suptitle("UK Biobank cholesterol: GWAS z-scores and fine-mapping PIPs",
             x=0.075, ha="left", y=0.972, fontsize=13.5, color=INK)
fig.text(0.075, 0.935, "UKB in-sample LD (337K white-British) · Borzoi-102 annotations · "
         "LDetect EUR blocks · n = 373,220 · identical inputs to all three methods",
         ha="left", fontsize=9, color=INK2)
fig.savefig(OUT, dpi=190)
print("wrote", OUT)
for rid, _, title in REG:
    s, fb, fm, pf = load(rid)
    az = np.abs(s.Z.values); rk = (-az).argsort().argsort() + 1
    print(f"\n{title}  p={len(s)}  max|Z|={az.max():.2f}")
    for v, l in [(fb, "fb_xregion"), (fm, "funmap"), (pf, "polyfun")]:
        i = int(np.argmax(v))
        print(f"   {l:12s} sum {v.sum():6.2f}  max {v.max():.3f} at bp {int(s.BP[i]):,} "
              f"(|Z| {az[i]:.2f}, rank {rk[i]})")
