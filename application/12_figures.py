"""Figures for the UKB cholesterol application, from results/export/.

Runs LOCALLY on the bundle 11_export_results.py writes, not on the cluster:
  python 12_figures.py --export ~/Downloads/chol_results_v2/export --out figures/

Palette is Okabe-Ito, a published colourblind-safe qualitative set, rather
than a hand-picked one. Colour is reserved for the four METHODS; the three
ablations of our own code are greys with distinct line styles, so the figure
never implies that `beatrice` or `fb_memomean` is a method being recommended.
Panel y-axes are shared wherever panels are compared, including where that
leaves the ablation panels nearly empty -- that emptiness is the result.
"""
import argparse, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ap = argparse.ArgumentParser()
ap.add_argument("--export", required=True)
ap.add_argument("--out", default="figures")
a = ap.parse_args()
E, OUT = os.path.expanduser(a.export), os.path.expanduser(a.out)
os.makedirs(OUT, exist_ok=True)

# Okabe-Ito. Methods get hues in a fixed order; ablations get greys.
METHODS   = ["fb_xregion", "funmap", "polyfun_combined", "polyfun_sniff"]
ABLATIONS = ["fb_prefix", "fb_memomean", "beatrice"]
ORDER     = METHODS + ABLATIONS
COLOUR = {"fb_xregion": "#0072B2", "funmap": "#D55E00",
          "polyfun_combined": "#009E73", "polyfun_sniff": "#56B4E9",
          "fb_prefix": "#555555", "fb_memomean": "#888888", "beatrice": "#BBBBBB"}
DASH   = {m: (0, ()) for m in METHODS}
DASH.update({"fb_prefix": (0, (4, 1.5)), "fb_memomean": (0, (1, 1.5)),
             "beatrice": (0, (5, 1.5, 1, 1.5))})
LABEL = {"fb_xregion": "functional BEATRICE", "funmap": "Funmap",
         "polyfun_combined": "PolyFun (combined)", "polyfun_sniff": "PolyFun (Sniff)",
         "fb_prefix": "FB ablation: uniform prior + fixed clip",
         "fb_memomean": "FB ablation: mean memo offset",
         "beatrice": "BEATRICE (upstream, own defaults)"}
# Axis labels must stay distinct from one another: splitting LABEL on ":"
# rendered both FB ablations as "FB ablation".
SHORT = {"fb_xregion": "functional BEATRICE", "funmap": "Funmap",
         "polyfun_combined": "PolyFun (combined)", "polyfun_sniff": "PolyFun (Sniff)",
         "fb_prefix": "FB abl: uniform prior", "fb_memomean": "FB abl: mean offset",
         "beatrice": "BEATRICE (upstream)"}
INK, MUTED, GRID = "#1a1a1a", "#595959", "#d9d9d9"
CONF, ZSIG = 0.95, 5.45

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "pdf.fonttype": 42, "savefig.bbox": "tight", "figure.dpi": 150,
})

def recede(ax, grid="y"):
    ax.grid(axis=grid, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

wide = pd.read_csv(os.path.join(E, "pips_wide.csv.gz"))
S    = pd.read_csv(os.path.join(E, "region_method_summary.csv"))
val  = pd.read_csv(os.path.join(E, "validated_pips.csv"))
wide["absZ"] = wide.Z.abs()
sig = S[S.sig]
NSIG = sig.region_id.nunique()

# ---------------------------------------------------------------- figure 1
# Genome-wide PIP, one row per column. Only variants a method actually put
# mass on are drawn; everything else sits at 0 and would be an ink-filled
# baseline carrying no information.
chrom_len = wide.groupby("chrom").BP.max()
offset, run = {}, 0
for c in sorted(chrom_len.index):
    offset[c] = run; run += chrom_len[c] * 1.02
wide["gpos"] = wide.BP + wide.chrom.map(offset)

fig, axes = plt.subplots(len(ORDER), 1, figsize=(7.4, 10.2), sharex=True, sharey=True)
for ax, m in zip(axes, ORDER):
    for i, c in enumerate(sorted(offset)):          # alternating chromosome bands
        if i % 2 == 0:
            lo = offset[c]; ax.axvspan(lo, lo + chrom_len[c], color="#f4f4f4", lw=0, zorder=0)
    d = wide[wide[m] >= 0.01]
    hit = d[m] >= CONF
    ax.scatter(d.gpos[~hit], d[m][~hit], s=3.5, c=COLOUR[m], alpha=0.55,
               lw=0, zorder=2)
    # Confident calls get a ring in the surface colour so overlapping marks
    # stay countable.
    ax.scatter(d.gpos[hit], d[m][hit], s=18, facecolor=COLOUR[m],
               edgecolor="white", lw=0.5, zorder=3)
    ax.axhline(CONF, color=MUTED, lw=0.5, ls=(0, (3, 2)), zorder=1)
    ax.set_ylim(-0.04, 1.08); ax.set_yticks([0, 0.5, 0.95])
    ax.set_ylabel("PIP")
    # Labels go in the axes TITLE, not inside the panel: at 0.93 in axes
    # coordinates they sat on the 0.95 rule and on the densest band of marks.
    # "variants", not "calls": 10_analyse_pips.py counts REGIONS whose top
    # variant clears 0.95, which is a different and smaller number.
    n = int(hit.sum()); nb = int(d.region_id[hit].nunique())
    ax.set_title(LABEL[m], loc="left", pad=3, fontsize=8,
                 fontweight="bold" if m in METHODS else "normal")
    ax.set_title(f"{n} variant{'' if n == 1 else 's'} ≥ 0.95 "
                 f"in {nb} block{'' if nb == 1 else 's'}",
                 loc="right", pad=3, fontsize=7, color=MUTED)
    recede(ax, grid="y")
    if m == ABLATIONS[0]:
        ax.spines["top"].set_visible(True)
        ax.spines["top"].set_color(INK)
        ax.spines["top"].set_linewidth(0.9)
axes[-1].set_xticks([offset[c] + chrom_len[c] / 2 for c in sorted(offset)])
axes[-1].set_xticklabels([str(c) for c in sorted(offset)])
axes[-1].set_xlabel("chromosome")
fig.suptitle("Posterior inclusion probabilities across 254 LDetect blocks, "
             "UK Biobank total cholesterol", y=1.002, fontsize=10, ha="center")
fig.text(0.5, 0.986, "four methods (coloured) above the rule; three ablations of "
         "functional BEATRICE (grey) below it", ha="center", fontsize=7.5, color=MUTED)
fig.subplots_adjust(hspace=0.34)
fig.savefig(os.path.join(OUT, "fig1_pips_genome.pdf"))
fig.savefig(os.path.join(OUT, "fig1_pips_genome.png"), dpi=200)
plt.close(fig)
print("fig1_pips_genome")

# ---------------------------------------------------------------- figure 2
fig, axes = plt.subplots(2, 2, figsize=(7.6, 6.9))

# (a) localisation: what |Z| does the method's own top variant carry?
ax = axes[0, 0]
RAMP = ["#08519c", "#6baed6", "#deebf7"]      # ordered, one hue, dark = good
cats = ["top |Z| ≥ 0.99 × lead", "2 ≤ top |Z| < 0.99 × lead",
        "top |Z| < 2"]
y = np.arange(len(ORDER))[::-1]
for j, m in enumerate(ORDER):
    g = sig[sig.method == m]
    near = (g.top_absZ >= 0.99 * g.lead_absZ).sum()
    weak = (g.top_absZ < 2).sum()
    mid  = len(g) - near - weak
    left = 0
    for v, col in zip([near, mid, weak], RAMP):
        ax.barh(y[j], v, left=left, height=0.68, color=col, lw=0.6,
                edgecolor="white", zorder=2)
        left += v
    ax.text(NSIG + 2, y[j], f"{near}", va="center", ha="left",
            fontsize=7, color=INK)
ax.set_yticks(y); ax.set_yticklabels([SHORT[m] for m in ORDER])
ax.set_xlim(0, NSIG + 16); ax.set_xlabel(f"signal regions (n = {NSIG})")
ax.set_title("a  |Z| at each region's top-PIP variant", loc="left")
ax.legend(handles=[Patch(facecolor=c, label=l, edgecolor="white")
                   for c, l in zip(RAMP, cats)], loc="upper center",
          bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False,
          fontsize=6.2, handlelength=1.1, columnspacing=1.0)
recede(ax, grid="x")

# (b) every confident call, not just each region's best one. 10_analyse_pips
#     reports regions whose TOP variant clears 0.95, which cannot see a method
#     that calls six variants in one block.
ax = axes[0, 1]
conf = pd.read_csv(os.path.join(E, "confident_calls.csv"))
conf["absZ"] = conf.Z.abs()
bands = [("|Z| \u2265 5.45 (genome-wide)", lambda z: z >= ZSIG),
         ("2 \u2264 |Z| < 5.45", lambda z: (z >= 2) & (z < ZSIG)),
         ("|Z| < 2 (no evidence)", lambda z: z < 2)]
for j, m in enumerate(ORDER):
    g = conf[conf.method == m]
    if g.empty:
        ax.text(0.01, y[j], "no calls", va="center", fontsize=6.5, color=MUTED)
        continue
    left = 0.0
    for (lab, f), col in zip(bands, RAMP):
        v = float(f(g.absZ).sum()) / len(g)
        ax.barh(y[j], v, left=left, height=0.68, color=col, lw=0.6,
                edgecolor="white", zorder=2)
        left += v
    ax.text(1.02, y[j], f"n = {len(g)}", va="center", ha="left",
            fontsize=6.8, color=INK)
ax.set_yticks(y); ax.set_yticklabels([SHORT[m] for m in ORDER])
ax.set_xlim(0, 1.28); ax.set_xticks([0, 0.5, 1.0])
ax.set_xticklabels(["0", "50%", "100%"])
ax.set_xlabel("share of that method's confident calls")
ax.set_title("b  Every variant called at PIP \u2265 0.95", loc="left")
ax.legend(handles=[Patch(facecolor=c, label=l, edgecolor="white")
                   for c, (l, _) in zip(RAMP, bands)], loc="upper center",
          bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False,
          fontsize=6.2, handlelength=1.1, columnspacing=1.0)
recede(ax, grid="x")

# (c) sum(PIP): the posterior mean number of causal variants, against the
#     caps each model was configured with.
ax = axes[1, 0]
rng = np.random.default_rng(0)
for j, m in enumerate(ORDER):
    g = S[S.method == m].sumPIP.to_numpy()
    ax.scatter(g, y[j] + rng.uniform(-0.17, 0.17, len(g)), s=2.5,
               c=COLOUR[m], alpha=0.45, lw=0, zorder=2)
    ax.scatter(np.median(g), y[j], s=26, marker="|", c=INK, lw=1.1, zorder=3)
for xv, lab in [(10, "L / max-causal / n_caus = 10"),
                (51, "sparse_concrete = 51")]:
    ax.axvline(xv, color=MUTED, lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax.text(xv, len(ORDER) - 0.45, f"{lab} ", fontsize=5.8, color=MUTED,
            rotation=90, ha="center", va="top")
ax.set_xscale("log"); ax.set_yticks(y)
ax.set_yticklabels([SHORT[m] for m in ORDER])
ax.set_ylim(-0.7, len(ORDER) - 0.2)
ax.set_xlabel("sum(PIP) per block  —  posterior mean no. of causal variants")
ax.set_title("c  sum(PIP) per block, against each model's cap", loc="left")
recede(ax, grid="x")

# (d) credible-set size
ax = axes[1, 1]
for m in ORDER:
    g = sig[sig.method == m].cs.dropna().to_numpy()
    xs = np.sort(g)
    ax.step(np.clip(xs, 1, None), np.arange(1, len(xs) + 1) / len(xs),
            where="post", color=COLOUR[m], lw=1.6 if m in METHODS else 1.1,
            dashes=DASH[m][1] if DASH[m][1] else (None, None),
            label=SHORT[m], zorder=3 if m in METHODS else 2)
ax.set_xscale("log"); ax.set_xlabel("95% credible-set size (variants)")
ax.set_ylabel("fraction of signal regions")
ax.set_title("d  95% credible-set size", loc="left")
ax.legend(loc="lower right", frameon=False, fontsize=6.5)
recede(ax, grid="both")

fig.tight_layout(pad=1.1)
fig.savefig(os.path.join(OUT, "fig2_method_behaviour.pdf"))
fig.savefig(os.path.join(OUT, "fig2_method_behaviour.png"), dpi=200)
plt.close(fig)
print("fig2_method_behaviour")

# ---------------------------------------------------------------- figure 3
# Positive controls. The right-hand column is not decoration: a PIP is only
# interpretable next to how many variants in the block are indistinguishable
# from the control.
ok = val[val.status == "ok"].copy()
ok["key"] = ok.gene.str.replace(" / ", "/", regex=False) + "\n" + ok.rsid
Mcols = [f"pip_{m}" for m in ORDER]
Z = ok[Mcols].to_numpy(float)

fig, (ax, axn) = plt.subplots(1, 2, figsize=(7.8, 3.5),
                              gridspec_kw={"width_ratios": [1, 0.1], "wspace": 0.04})
# Sequential, one hue, light -> dark. PIP spans 20 orders of magnitude, so the
# colour is on log10 and the printed number carries the exact value.
im = ax.imshow(np.log10(np.clip(Z, 1e-21, 1)), cmap="Blues",
               vmin=-6, vmax=0, aspect="auto")
for i in range(Z.shape[0]):
    for j in range(Z.shape[1]):
        v = Z[i, j]
        s = f"{v:.3f}" if v >= 5e-4 else f"{v:.0e}".replace("e-0", "e-")
        ax.text(j, i, s, ha="center", va="center", fontsize=6.4,
                color="white" if np.log10(max(v, 1e-21)) > -1.4 else INK)
ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels([SHORT[m] for m in ORDER], rotation=32,
                   ha="right", fontsize=7)
for t, m in zip(ax.get_xticklabels(), ORDER):
    t.set_fontweight("bold" if m in METHODS else "normal")
ax.set_yticks(range(len(ok))); ax.set_yticklabels(ok.key, fontsize=7)
ax.axvline(len(METHODS) - 0.5, color=INK, lw=1.1)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_title("PIP on variants with variant-level functional evidence  "
             "(methods left of the rule, ablations right)", loc="left", fontsize=8.5)

# the identifiability column
axn.imshow(np.log10(ok.n_r2_ge_99.to_numpy(float)).reshape(-1, 1),
           cmap="Greys", vmin=-0.45, vmax=1.15, aspect="auto")
for i, (n, lead) in enumerate(zip(ok.n_r2_ge_99, ok.is_lead)):
    axn.text(0, i, f"{int(n)}", ha="center", va="center", fontsize=6.8,
             color="white" if n > 5 else INK)
axn.set_xticks([0]); axn.set_xticklabels(["variants at\nr² ≥ 0.99"],
                                         rotation=32, ha="right", fontsize=6.5)
axn.set_yticks([])
for s in axn.spines.values():
    s.set_visible(False)
fig.text(0.5, -0.20,
         "Where k variants sit at r\u00b2 \u2265 0.99, association data alone cannot "
         "exceed about 1/k on any one of them \u2014 SORT1 k = 4, LDLR k = 8.\n"
         "A PIP near 1 there is carried by the prior, not the likelihood. "
         "rs7412 and rs6511720 are not their block's lead variant; the other five are.",
         ha="center", va="top", fontsize=6.8, color=MUTED, linespacing=1.5)
fig.tight_layout(pad=0.9)
fig.savefig(os.path.join(OUT, "fig3_positive_controls.pdf"))
fig.savefig(os.path.join(OUT, "fig3_positive_controls.png"), dpi=200)
plt.close(fig)
print("fig3_positive_controls")
print("\nwrote to", OUT)
