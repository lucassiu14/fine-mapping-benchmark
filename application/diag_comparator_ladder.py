"""Audit the four-rung BEATRICE/fb_xregion ablation ladder.

06_collect.py asserts that each method's PIP vector has the SAME LENGTH as
snps.csv, but never that it is in the same ORDER. For our own fork that is
safe (prep writes locus.z and the model reads it row for row), but upstream
BEATRICE is a separate program reading a text LD file we materialise for it,
so its row order is an assumption, not a checked fact. If it were permuted,
every PIP would be attached to the wrong variant and the result would look
exactly like what we observe: rankings uncorrelated with |Z|.

This script closes that gap before the "upstream BEATRICE does not rank"
claim is written down anywhere. It reports, in order:

  [1] whether upstream's own variant_names column reproduces snps.csv order
      in all 254 analysed blocks (upstream writes the names it used, so this
      is its testimony, not ours). The names are NOT rsIDs: prep_block.py:172
      writes locus.z as "<region_id>_<BP>" and BEATRICE reads column 0 of that
      file back as its names, so the comparison has to be built from snps.csv
      BP. Comparing them to snps.csv SNP instead reports 254/254 mismatched
      no matter how well aligned the run is, which is exactly what the first
      version of this script did;
  [2] per-locus detail on the strongest loci: where each rung puts the lead
      variant, and the rank correlation between its PIPs and |Z|;
  [3] upstream's own credible_set.txt sizes, to confirm the gamma_key gate
      never fires rather than inferring it from the PIPs;
  [4] how much sum(PIP) varies across blocks, which is the tell that
      separates a saturated normaliser (near-constant) from a posterior that
      responds to the data.
"""
import os, glob
import numpy as np, pandas as pd

R = os.environ.get("APP_ROOT", os.getcwd())
S = pd.read_csv(os.path.join(R, "results", "region_method_summary.csv"))
rids = sorted(S.region_id.unique())

# The ladder, weakest rung first. Each is a directory of <rid>/pip.csv.
LADDER = [("beatrice",    "beatrice_out"),
          ("fb_memomean", "joint_out.memomean"),
          ("fb_prefix",   "joint_out.pre_fix_d0fb022"),
          ("fb_xregion",  "joint_out")]

def pip_of(sub, rid):
    d = pd.read_csv(os.path.join(R, sub, rid, "pip.csv"))
    col = next((c for c in d.columns if c.strip().lower() == "pip"), None)
    if col is None:                       # same fallback 06_collect.py uses
        col = [c for c in d.columns if pd.api.types.is_numeric_dtype(d[c])][-1]
    return d, d[col].to_numpy(float)

# --- [1] row-order agreement, all blocks ---------------------------------
# The label prep writes for variant i, which is what upstream echoes back.
bad_names, bad_index, unnamed, checked = [], [], 0, 0
for rid in rids:
    d, _ = pip_of("beatrice_out", rid)
    if "variant_names" not in d.columns:
        unnamed += 1; continue
    s = pd.read_csv(os.path.join(R, "regions", rid, "snps.csv"), usecols=["BP"])
    checked += 1
    expect = np.array([f"{rid}_{int(bp)}" for bp in s.BP.to_numpy()])
    if not np.array_equal(d.variant_names.astype(str).to_numpy(), expect):
        bad_names.append(rid)
    # Independently of the labels: upstream must emit one row per input
    # variant in input order, or the pip column cannot be zipped to snps.csv.
    if not np.array_equal(d.variant_index.to_numpy(), np.arange(len(s))):
        bad_index.append(rid)
print(f"[1] beatrice row order vs prep, {checked} blocks checked: "
      f"variant_names mismatched {len(bad_names)}, "
      f"variant_index not 0..p-1 in {len(bad_index)}, "
      f"no variant_names column {unnamed}")
for lbl, lst in (("names", bad_names), ("index", bad_index)):
    if lst:
        print(f"    first {lbl} mismatches:", ", ".join(lst[:5]))
        rid = lst[0]
        d, _ = pip_of("beatrice_out", rid)
        s = pd.read_csv(os.path.join(R, "regions", rid, "snps.csv"), usecols=["SNP", "BP"])
        print(f"    {rid} upstream says : {list(d.variant_names.astype(str)[:3])}")
        print(f"    {rid} prep expects  : {[f'{rid}_{int(b)}' for b in s.BP[:3]]}")
        print(f"    {rid} snps.csv SNP  : {list(s.SNP.astype(str)[:3])}")

# --- [2] the strongest loci, rung by rung --------------------------------
top = S[S.method == "fb_xregion"].nlargest(6, "lead_absZ")
for _, r in top.iterrows():
    rid = r.region_id
    s = pd.read_csv(os.path.join(R, "regions", rid, "snps.csv"), usecols=["SNP", "Z"])
    aZ = s.Z.abs().to_numpy()
    k = int(np.argmax(aZ))                        # the GWAS lead
    print(f"\n[2] {rid}  p={len(s)}  lead={s.SNP.iloc[k]}  |Z|={aZ[k]:.2f}")
    for m, sub in LADDER:
        _, p = pip_of(sub, rid)
        # Rank correlation against |Z|. A method may legitimately disagree
        # with |Z| on WHICH variant is causal -- that is what fine-mapping is
        # for -- but a whole-region rho near 0 means the PIPs carry no
        # information about the marginal evidence at all.
        rho = pd.Series(p).corr(pd.Series(aZ), method="spearman")
        print(f"    {m:12s} sum={p.sum():7.3f} max={p.max():.4f} "
              f"leadPIP={p[k]:.4g} leadRank={int((p > p[k]).sum()) + 1:>5d}/{len(p)} "
              f"rho(pip,|Z|)={rho:+.3f} |Z| at its top={aZ[int(np.argmax(p))]:.2f}")

# --- [3] upstream's own credible sets ------------------------------------
sz = [os.path.getsize(f) for f in glob.glob(os.path.join(R, "beatrice_out", "*", "credible_set.txt"))]
if sz:
    print(f"\n[3] beatrice credible_set.txt: {len(sz)} files, "
          f"{sum(1 for x in sz if x == 0)} zero-byte, largest {max(sz)} B")
else:
    print("\n[3] beatrice wrote no credible_set.txt at all")

# --- [4] does sum(PIP) respond to the data? ------------------------------
print("\n[4] sum(PIP) across the 254 blocks")
for m in sorted(S.method.unique()):
    g = S.loc[S.method == m, "sumPIP"]
    print(f"    {m:18s} min={g.min():8.3f} med={g.median():8.3f} "
          f"max={g.max():9.3f} sd={g.std():7.3f} cv={g.std()/g.mean():.3f}")


# --- [5] does the method localise, even when it cannot commit? -----------
# top_is_lead is a strict test: in a block where the lead has near-perfect LD
# partners, picking the partner is not an error worth distinguishing from
# picking the lead. So report, beside it, the |Z| the method's own top variant
# carries. A method that lands on a variant with ~the lead's |Z| has found the
# signal and merely cannot resolve within it; one whose top variant carries a
# small |Z| has not found it at all.
sig = S[S.sig]
print(f"\n[5] localisation across the {sig.region_id.nunique()} signal regions")
print(f"    {'method':18s} {'med top|Z|':>11s} {'med lead|Z|':>12s} "
      f"{'top|Z|>=.99 lead':>17s} {'top|Z|<2':>9s} {'top=lead':>9s}")
for m in sorted(sig.method.unique()):
    g = sig[sig.method == m]
    near = int((g.top_absZ >= 0.99 * g.lead_absZ).sum())
    print(f"    {m:18s} {g.top_absZ.median():>11.2f} {g.lead_absZ.median():>12.2f} "
          f"{near:>10d}/{len(g):<6d} {int((g.top_absZ < 2).sum()):>9d} "
          f"{int(g.top_is_lead.sum()):>9d}")
