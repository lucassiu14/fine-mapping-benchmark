"""Gather every method's PIPs into one long table, plus the annotation weights.

Output: results/pips.parquet with one row per (region, variant, method) and
        results/annotation_weights.csv comparing what each method learned.

There is no ground truth here - this is real GWAS data - so nothing in this
script computes power, FDR or calibration against a known causal set. What it
does produce is the per-region PIPs and credible sets, which is what the
comparison between selection (fb_xregion) and shrinkage (funmap) rests on.
"""
import os, glob, json
import numpy as np, pandas as pd

ROOT = os.environ.get("APP_ROOT", os.path.dirname(os.path.abspath(__file__)))
RES  = os.path.join(ROOT, "results"); os.makedirs(RES, exist_ok=True)
blocks = pd.read_csv(os.path.join(ROOT, "blocks.csv"))

rows, missing = [], []
for rid in blocks.region_id:
    d = os.path.join(ROOT, "regions", rid)
    if not os.path.exists(os.path.join(d, "DONE")):
        missing.append((rid, "prep")); continue
    snps = pd.read_csv(os.path.join(d, "snps.csv"))
    base = dict(region_id=rid, chrom=int(snps.BP.notna() and blocks.loc[blocks.region_id == rid, "chrom"].iloc[0]))

    src = {
        "fb_xregion":    (os.path.join(ROOT, "joint_out", rid, "pip.csv"), "csv_named"),
        "funmap":        (os.path.join(ROOT, "funmap_out", rid, "pip.csv"), "csv_headerless"),
        "polyfun_susie": (os.path.join(ROOT, "polyfun_out", f"{rid}.txt.gz"), "finemapper"),
    }
    for meth, (path, kind) in src.items():
        if not os.path.exists(path):
            missing.append((rid, meth)); continue
        if kind == "csv_headerless":          # funmap writes SNP,PIP with no header
            p = pd.read_csv(path, header=None, names=["SNP", "pip"])["pip"].to_numpy(float)
        elif kind == "finemapper":
            # PolyFun's finemapper writes its own TSV, keyed by SNP and NOT
            # necessarily in our row order or covering every variant, so align
            # on SNP rather than assuming positions line up.
            fm = pd.read_csv(path, sep="\t")
            col = "PIP" if "PIP" in fm.columns else "pip"
            p = snps.SNP.map(dict(zip(fm.SNP, fm[col]))).fillna(0.0).to_numpy(float)
        else:
            p = pd.read_csv(path)["pip"].to_numpy(float)
        if len(p) != len(snps):
            missing.append((rid, f"{meth}:length {len(p)} != {len(snps)}")); continue
        rows.append(pd.DataFrame(dict(
            region_id=rid, chrom=base["chrom"], BP=snps.BP.values, SNP=snps.SNP.values,
            Z=snps.Z.values, method=meth, pip=p)))

if rows:
    df = pd.concat(rows, ignore_index=True)
    df.to_parquet(os.path.join(RES, "pips.parquet"), index=False)
    print(f"wrote {len(df):,} rows for {df.region_id.nunique()} regions, "
          f"methods: {sorted(df.method.unique())}")
    s = (df.groupby(["method", "region_id"]).pip.agg(["sum", "max"])
           .groupby("method").agg(mean_sumPIP=("sum", "mean"), mean_maxPIP=("max", "mean"),
                                  regions=("sum", "size")))
    print(s.round(3).to_string())
else:
    print("no results found")

if missing:
    print(f"\n{len(missing)} missing outputs (first 20):")
    for rid, why in missing[:20]:
        print(f"  {rid:20s} {why}")

# annotation weights: fb_xregion writes feature_importance.csv per region
# (identical across regions - one shared head). Funmap writes its weights
# into sets.txt / its own output; collected here if present.
fi = glob.glob(os.path.join(ROOT, "joint_out", "*", "feature_importance.csv"))
if fi:
    w = pd.read_csv(fi[0])
    names = [l.strip() for l in open(os.path.join(ROOT, "annot_names.txt"))] \
            if os.path.exists(os.path.join(ROOT, "annot_names.txt")) else None
    if names and len(names) == len(w):
        w["annotation"] = [names[i] for i in w.annotation_index]
    w.to_csv(os.path.join(RES, "annotation_weights.csv"), index=False)
    print(f"\nfb_xregion annotation importances -> results/annotation_weights.csv")
    print(w.head(10).to_string(index=False))
