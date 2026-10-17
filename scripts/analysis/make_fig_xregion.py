"""
Cross-region (joint) Functional BEATRICE architecture, for the paper.

Extends make_fig2_architecture.py from one region to R. The structure drawn
here is verified against BEATRICE_annot_sparse/scripts/joint_trainer.py:

  * ONE shared prior head G(.,psi).  joint_trainer.py:305 concatenates every
    region's annotation matrix into v_all, applies the head once, and slices
    the result back per region -- so the fan-in/fan-out is literal, not a
    schematic convenience.
  * A SEPARATE inference network per region (line 292), built inside the region
    loop with its own parameters added to the optimiser separately (line 325).
    Hence phi_r, but a bare psi.
  * The concrete temperature is computed once per epoch OUTSIDE the region loop
    (line 336), so lambda is global too -- it must NOT take an r.
  * The regions are coupled only through psi and the equal-weighted objective
    (1/R) sum_r L_r (line 348).  That collector is the whole content of the
    word "joint", so it is drawn explicitly.

Two regions are drawn in full with a vertical ellipsis between them; drawing
all R would make the figure unreadable without saying anything new.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, PathPatch
from matplotlib.path import Path

# ----------------------------------------------------------------- style ----
FONTDIR = "/Users/ls1020/Downloads/FB_poster_extract/Fonts"
for f in ["ImperialSansText-Regular.ttf", "ImperialSansText-Bold.ttf",
          "ImperialSansText-Medium.ttf", "ImperialSansText-Semibold.ttf"]:
    font_manager.fontManager.addfont(f"{FONTDIR}/{f}")
plt.rcParams["font.family"] = "Imperial Sans Text"
plt.rcParams["mathtext.fontset"] = "dejavusans"

BLACK   = "#000000"
VAR_ED  = "#3A3A3A"
NET_ED, NET_FC = "#0000CD", "#EDF0FC"
RND_ED, RND_FC = "#2E8B7A", "#E7F3F1"
GEN_ED, GEN_FC = "#B07D2B", "#F8F1E3"
AGG_ED, AGG_FC = "#8E44AD", "#F4ECF7"
PLATE_ED = "#9A9A9A"
NODE_IN, NODE_OUT, EDGE = "#5F5F5F", "#BFBFBF", "#B4B4B4"

import os
OUT = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- canvas ----
fig = plt.figure(figsize=(14.4, 6.9))
fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(-2, 218)
ax.set_ylim(-7, 98)
ax.set_aspect("equal")
ax.axis("off")

LW = 2.0

# ------------------------------------------------------------- geometry ----
Y1, YR = 72.0, 10.0          # the two drawn regions (gap left for the ellipsis)
YMID   = (Y1 + YR) / 2       # shared head sits between them
UP     = 17.0                # p_0r and Sigma_Xr ride this far above their row
Y_PRIOR = YMID + UP          # the head sits level with the midpoint of the two
                             # p_0r branches, not of the two region rows
GAP_Y   = (Y1 - 13.0 + YR + 22.0) / 2    # clear band between the plates

VAR_W, VAR_H = 7.5, 7.5
MOD_W, MOD_H = 23.0, 20.5
MOD_H_G      = 28.0          # no longer needs to span R separate inputs
PRC_W, PRC_H = 25.0, 12.5
GEN_W        = 28.0

X_V   = 7.0                  # v_1 ... v_R
X_VDOT = 19.0                # they merge here before entering G
X_G   = 36.0                 # shared prior network
X_BUS = 56.0                 # fan-out bus
X_P0  = 68.0                 # p_0r: outside the plate, right off the bus
X_Z   = 84.0
X_F   = 104.0
X_P   = 128.0
X_Q   = 152.0
X_CL  = 174.0
X_GEN = 196.0

PLATE_L, PLATE_R = 78.0, 214.0


def box(cx, cy, w, h, edge, face, r=1.6, lw=LW, z=2):
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle=f"round,pad=0,rounding_size={r}",
        linewidth=lw, edgecolor=edge, facecolor=face, zorder=z))


def arrow(p0, p1, lw=2.4, ms=22, color=BLACK, z=5):
    ax.add_patch(FancyArrowPatch(
        p0, p1, arrowstyle="-|>", mutation_scale=ms,
        linewidth=lw, color=color, zorder=z, shrinkA=0, shrinkB=0))


def line(p0, p1, lw=2.4, color=BLACK, z=5):
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, linewidth=lw,
            zorder=z, solid_capstyle="round")


def mathlabel(cx, cy, text, size=18, color=BLACK, z=6):
    ax.text(cx, cy, text, ha="center", va="center", fontsize=size,
            color=color, zorder=z)


def vdots(cx, cy, size=20, color="#555555"):
    """Three dots drawn as circles rather than mathtext.

    matplotlib's $\\vdots$ glyph carries about 0.22*fontsize of left side
    bearing, so ha="center" centres its bounding box and leaves the dots
    themselves visibly right of cx -- obvious as soon as they have to line up
    with a rule.  Circles sit exactly on cx."""
    r, gap = size * 0.018, size * 0.110
    for k in (-1, 0, 1):
        ax.add_patch(Circle((cx, cy + k * gap), r, facecolor=color,
                            edgecolor="none", zorder=6))


def _layer_coords(cx, cy, w, h, layers):
    xs = np.linspace(cx - w / 2, cx + w / 2, len(layers))
    coords = []
    for x, n in zip(xs, layers):
        ys = np.linspace(cy + h / 2, cy - h / 2, n) if n > 1 else [cy]
        coords.append([(x, y) for y in ys])
    return coords


def _nodes(coords, pruned=()):
    for li, layer in enumerate(coords):
        fc = NODE_IN if li == 0 else NODE_OUT
        for ni, (x, y) in enumerate(layer):
            faded = (li == 0 and ni in pruned)
            ax.add_patch(Circle(
                (x, y), 1.05,
                facecolor="white" if faded else fc,
                edgecolor=EDGE if faded else "white",
                linewidth=1.1 if faded else 0.8, zorder=4))


def network(cx, cy, w, h, layers=(3, 4, 3)):
    coords = _layer_coords(cx, cy, w, h, layers)
    for a, b in zip(coords[:-1], coords[1:]):
        for (x0, y0) in a:
            for (x1, y1) in b:
                ax.plot([x0, x1], [y0, y1], color=EDGE, linewidth=0.9,
                        zorder=3, solid_capstyle="round")
    _nodes(coords)


def lassonet(cx, cy, w, h, layers=(3, 4, 2), pruned=(1,)):
    coords = _layer_coords(cx, cy, w, h, layers)
    inp, out = coords[0], coords[-1]
    for li, (a, b) in enumerate(zip(coords[:-1], coords[1:])):
        for ni, (x0, y0) in enumerate(a):
            faded = (li == 0 and ni in pruned)
            for (x1, y1) in b:
                ax.plot([x0, x1], [y0, y1], color=EDGE, linewidth=0.9,
                        alpha=0.22 if faded else 1.0, zorder=3,
                        solid_capstyle="round")
    dip = cy - h / 2 - 2.6
    for ni, (x0, y0) in enumerate(inp):
        faded = ni in pruned
        for (x1, y1) in out:
            mx = (x0 + x1) / 2
            path = Path([(x0, y0), (mx, dip), (x1, y1)],
                        [Path.MOVETO, Path.CURVE3, Path.CURVE3])
            ax.add_patch(PathPatch(
                path, facecolor="none", edgecolor=NET_ED,
                linewidth=0.8 if faded else 1.7,
                alpha=0.28 if faded else 0.95, zorder=3.5))
    _nodes(coords, pruned=pruned)


# ======================================================== per-region plate ----
def region(Y, s):
    """One region's pipeline. `s` is the subscript string ('1' or 'R')."""
    # plate
    ax.add_patch(FancyBboxPatch(
        (PLATE_L, Y - 13.0), PLATE_R - PLATE_L, 35.0,
        boxstyle="round,pad=0,rounding_size=2.2",
        linewidth=1.6, edgecolor=PLATE_ED, facecolor="white",
        linestyle=(0, (5, 4)), zorder=1))
    ax.text((PLATE_L + PLATE_R) / 2, Y + 20.4, f"Region ${s}$", ha="center", va="top",
            fontsize=13, color=BLACK, zorder=6)

    box(X_Z, Y, VAR_W, VAR_H, VAR_ED, "white")
    mathlabel(X_Z, Y, rf"$\mathbf{{z}}_{{{s}}}$", 17)

    box(X_F, Y, MOD_W, MOD_H, NET_ED, NET_FC)
    mathlabel(X_F, Y + MOD_H / 2 - 3.6, rf"$\mathcal{{F}}(\cdot\,;\phi_{{{s}}})$", 16)
    network(X_F, Y - 2.6, 14.0, 7.2)

    box(X_P, Y, VAR_W, VAR_H, VAR_ED, "white")
    mathlabel(X_P, Y, rf"$\mathbf{{p}}_{{{s}}}$", 17)

    box(X_P0, Y + UP, VAR_W + 1.0, VAR_H, VAR_ED, "white")
    mathlabel(X_P0, Y + UP, rf"$\mathbf{{p}}_{{0}}^{{({s})}}$", 16)

    box(X_Q, Y, PRC_W, PRC_H, RND_ED, RND_FC)
    mathlabel(X_Q, Y, rf"$\mathbf{{c}}_{{{s}}}\sim q(\mathbf{{c}}_{{{s}}}\mid\mathbf{{p}}_{{{s}}},\lambda)$", 13)

    mathlabel(X_CL, Y, rf"$\mathbf{{c}}_{{{s}}}^{{(l)}}$", 17)

    box(X_GEN, Y + UP, VAR_W + 3.0, VAR_H, VAR_ED, "white")
    mathlabel(X_GEN, Y + UP, rf"$\widehat{{R}}_{{{s}}}$", 16)

    box(X_GEN, Y, GEN_W, PRC_H, GEN_ED, GEN_FC)
    mathlabel(X_GEN, Y,
              rf"$p(\mathbf{{z}}_{{{s}}}\mid\widehat{{R}}_{{{s}}},\mathbf{{c}}_{{{s}}}^{{(l)}})$", 14)

    arrow((X_Z + VAR_W / 2, Y), (X_F - MOD_W / 2, Y))
    arrow((X_F + MOD_W / 2, Y), (X_P - VAR_W / 2, Y))
    line((X_P0 + (VAR_W + 1.0) / 2, Y + UP), (X_P, Y + UP))
    arrow((X_P, Y + UP), (X_P, Y + VAR_H / 2))
    arrow((X_P + VAR_W / 2, Y), (X_Q - PRC_W / 2, Y))
    arrow((X_Q + PRC_W / 2, Y), (X_CL - 4.0, Y))
    arrow((X_CL + 4.0, Y), (X_GEN - GEN_W / 2, Y))
    arrow((X_GEN, Y + UP - VAR_H / 2), (X_GEN, Y + PRC_H / 2))


region(Y1, "1")
region(YR, "R")

# the regions between 1 and R, omitted
vdots((PLATE_L + PLATE_R) / 2, GAP_Y, 30)

# ============================================= shared prior network (psi) ----
HL = FancyBboxPatch((2.0, Y_PRIOR - 19.0), 48.0, 38.0,
                    boxstyle="round,pad=0,rounding_size=2.0",
                    facecolor=NET_FC, alpha=0.45, edgecolor=NET_ED,
                    linewidth=2.2, linestyle=(0, (6, 4)), zorder=0.5)
ax.add_patch(HL)

Y_V = [Y_PRIOR + 11.0, Y_PRIOR - 11.0]          # A_1 and A_R
for yv, s in zip(Y_V, ["1", "R"]):
    box(X_V, yv, VAR_W, VAR_H, VAR_ED, "white")
    mathlabel(X_V, yv, rf"$A_{{{s}}}$", 17)
    line((X_V + VAR_W / 2, yv), (X_VDOT, Y_PRIOR))      # merge into the junction
vdots(X_V, Y_PRIOR, 20)
ax.add_patch(Circle((X_VDOT, Y_PRIOR), 1.5, facecolor=BLACK, edgecolor="none",
                    zorder=6))
arrow((X_VDOT, Y_PRIOR), (X_G - MOD_W / 2, Y_PRIOR))

box(X_G, Y_PRIOR, MOD_W, MOD_H_G, NET_ED, NET_FC)
mathlabel(X_G, Y_PRIOR + MOD_H_G / 2 - 4.4, r"$\mathcal{G}(\cdot\,;\psi)$", 17)
lassonet(X_G, Y_PRIOR - 2.0, 14.0, 6.4)

# ---- fan-out: G -> p_0r, routed orthogonally so nothing crosses a pipeline --
line((X_G + MOD_W / 2, Y_PRIOR), (X_BUS, Y_PRIOR))
line((X_BUS, YR + UP), (X_BUS, GAP_Y - 6.5))   # break the bus level with the
line((X_BUS, GAP_Y + 6.5), (X_BUS, Y1 + UP))   # ellipsis between the plates
vdots(X_BUS, GAP_Y, 24)
for Y in (Y1, YR):
    arrow((X_BUS, Y + UP), (X_P0 - (VAR_W + 1.0) / 2, Y + UP))
ax.add_patch(Circle((X_BUS, Y_PRIOR), 1.5, facecolor=BLACK, edgecolor="none",
                    zorder=6))

fig.savefig(f"{OUT}/fig_xregion_architecture.pdf", format="pdf",
            bbox_inches="tight", pad_inches=0.06)
fig.savefig(f"{OUT}/fig_xregion_architecture.png", dpi=190,
            bbox_inches="tight", pad_inches=0.06)
print("written")
