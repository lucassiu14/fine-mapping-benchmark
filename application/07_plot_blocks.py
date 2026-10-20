"""GWAS z-scores and every method's PIPs, for selected blocks.

    python 07_plot_blocks.py                 # the 4 strongest blocks
    python 07_plot_blocks.py ldEUR_01_0067 ldEUR_13_0062
    python 07_plot_blocks.py --n 6           # the 6 strongest

Writes results/block_pips.pdf. Blocks with no method output are skipped with
a note rather than failing, so this is usable while a run is still going.
"""
import os, sys, json
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.environ.get("APP_ROOT", os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results"); os.makedirs(RES, exist_ok=True)
OUT = os.path.join(RES, "block_pips.pdf")

# validated categorical palette, assigned in FIXED ORDER, never cycled
METHODS = [
    ("fb_xregion",       "#2a78d6", lambda r: (os.path.join(ROOT, "joint_out", r, "pip.csv"), "named")),
    ("funmap",           "#eb6834", lambda r: (os.path.join(ROOT, "funmap_out", r, "pip.csv"), "headerless")),
    ("polyfun (102)",    "#1baf7a", lambda r: (os.path.join(ROOT, "polyfun_out_sniff", f"{r}.txt.gz"), "finemapper")),
    ("polyfun (289)",    "#eda100", lambda r: (os.path.join(ROOT, "polyfun_out_combined", f"{r}.txt.gz"), "finemapper")),
]
INK, INK2, GRID, BG = "#0b0b0b", "#52514e", "#d8d7d2", "#fcfcfb"


def pick_blocks(argv):
    blocks = pd.read_csv(os.path.join(ROOT, "blocks.csv"))
    named = [a for a in argv if not a.startswith("-")]
    if named:
        miss = [a for a in named if a not in set(blocks.region_id)]
        if miss:
            sys.exit(f"not in blocks.csv: {miss}")
        return blocks[blocks.region_id.isin(named)]
    n = 4
    if "--n" in argv:
        n = int(argv[argv.index("--n") + 1])
    return blocks.nsmallest(n, "min_p")


def load(rid, snps):
    out = {}
    for name, colour, locate in METHODS:
        path, kind = locate(rid)
        if not os.path.exists(path):
            continue
        try:
            if kind == "headerless":            # funmap: SNP,PIP with no header
                v = pd.read_csv(path, header=None, names=["SNP", "pip"])["pip"].to_numpy(float)
            elif kind == "finemapper":          # PolyFun: its own TSV, keyed by SNP
                fm = pd.read_csv(path, sep="\t")
                col = "PIP" if "PIP" in fm.columns else "pip"
                v = snps.SNP.map(dict(zip(fm.SNP, fm[col]))).fillna(0.0).to_numpy(float)
            else:
                v = pd.read_csv(path)["pip"].to_numpy(float)
        except Exception as e:                  # a half-written file mid-run
            print(f"  {rid}: could not read {name} ({e})"); continue
        if len(v) != len(snps):
            print(f"  {rid}: {name} has {len(v)} PIPs for {len(snps)} variants - skipped")
            continue
        out[name] = (v, colour)
    return out


def main(argv):
    sel = pick_blocks(argv)
    panels = []
    for r in sel.itertuples():
        d = os.path.join(ROOT, "regions", r.region_id)
        if not os.path.exists(os.path.join(d, "snps.csv")):
            print(f"  {r.region_id}: not prepared yet - skipped"); continue
        snps = pd.read_csv(os.path.join(d, "snps.csv"))
        res = load(r.region_id, snps)
        if not res:
            print(f"  {r.region_id}: no method output yet - skipped"); continue
        panels.append((r, snps, res))
    if not panels:
        sys.exit("nothing to plot yet - no block has both prep and method output")

    names = [n for n, _, _ in METHODS if any(n in p[2] for p in panels)]
    nrow, ncol = 1 + len(names), len(panels)
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.6 * ncol, 2.4 * nrow + 0.8),
                             sharey="row", squeeze=False,
                             gridspec_kw=dict(hspace=0.25, wspace=0.07, left=0.085,
                                              right=0.985, top=0.88, bottom=0.09))
    fig.patch.set_facecolor(BG)

    def style(ax, ylab=None):
        ax.set_facecolor(BG)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
        for s in ("left", "bottom"): ax.spines[s].set_color(GRID)
        ax.tick_params(colors=INK2, labelsize=8.5, length=3)
        ax.grid(axis="y", color=GRID, lw=0.6, alpha=0.7); ax.set_axisbelow(True)
        if ylab: ax.set_ylabel(ylab, color=INK, fontsize=9.5)

    for col, (r, snps, res) in enumerate(panels):
        x = snps.BP.values / 1e6
        ax = axes[0][col]
        ax.scatter(x, snps.Z.values, s=4, c=INK2, lw=0, alpha=0.55, rasterized=True)
        for t in (5.451, -5.451):
            ax.axhline(t, color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=1)
        style(ax, "GWAS z-score" if col == 0 else None)
        zm = np.abs(snps.Z.values).max(); ax.set_ylim(-zm * 1.18, zm * 1.18)
        ax.set_title(f"{r.region_id}   chr{r.chrom}:{r.block_start/1e6:.2f}-{r.block_end/1e6:.2f} Mb"
                     f"   n = {len(snps):,}", color=INK, fontsize=10, pad=7)
        for row, nm in enumerate(names, start=1):
            ax = axes[row][col]
            if nm not in res:
                ax.text(0.5, 0.5, "not run", transform=ax.transAxes, ha="center",
                        va="center", fontsize=10, color=INK2, style="italic")
                style(ax, f"{nm}\nPIP" if col == 0 else None); ax.set_ylim(-0.05, 1.09)
                continue
            v, c = res[nm]
            ax.scatter(x, v, s=4 + 110 * v ** 1.5, c=c, lw=0, alpha=0.85, rasterized=True)
            ax.axhline(0.5, color=INK2, lw=0.7, ls=(0, (4, 3)), zorder=1)
            style(ax, f"{nm}\nPIP" if col == 0 else None); ax.set_ylim(-0.05, 1.09)
            i = int(np.argmax(v))
            if v[i] > 0.15:
                ax.annotate(f"{v[i]:.2f}", (x[i], v[i]), textcoords="offset points",
                            xytext=(0, 7), ha="center", fontsize=8, color=INK)
            ax.text(0.985, 0.88, f"Σ={v.sum():.2f}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=8, color=INK2)
        axes[-1][col].set_xlabel(f"position on chr{r.chrom} (Mb)", color=INK, fontsize=9)

    fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=6, mfc=c, mec="none", label=n)
                        for n, c, _ in METHODS if n in names],
               loc="lower center", ncol=len(names), frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, -0.002), labelcolor=INK)
    fig.suptitle("UK Biobank cholesterol: GWAS z-scores and fine-mapping PIPs",
                 x=0.085, ha="left", y=0.978, fontsize=13, color=INK)
    fig.text(0.085, 0.942, f"UKB in-sample LD · Borzoi-102 annotations · LDetect EUR blocks · "
             f"n = {os.environ.get('N_LIKELIHOOD', '373220')} · identical inputs to all methods",
             ha="left", fontsize=8.5, color=INK2)
    fig.savefig(OUT, dpi=190)
    print(f"\nwrote {OUT}")
    for r, snps, res in panels:
        az = np.abs(snps.Z.values); rk = (-az).argsort().argsort() + 1
        print(f"\n{r.region_id}  p={len(snps)}  max|Z|={az.max():.2f}")
        for nm, (v, _) in res.items():
            i = int(np.argmax(v))
            print(f"   {nm:16s} sum {v.sum():6.2f}  max {v.max():.3f} at bp {int(snps.BP[i]):,} "
                  f"(|Z| {az[i]:.2f}, rank {rk[i]})")


if __name__ == "__main__":
    main(sys.argv[1:])
