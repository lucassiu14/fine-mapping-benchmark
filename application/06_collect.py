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

chrom_of = dict(zip(blocks.region_id, blocks.chrom))

rows, missing, excluded = [], [], []
for rid in blocks.region_id:
    d = os.path.join(ROOT, "regions", rid)
    # A block excluded by the unannotated-gap rule is not a missing output:
    # it is absent from every manifest by construction and no method wrote it.
    if os.path.exists(os.path.join(d, "EXCLUDED")):
        excluded.append(rid); continue
    if not os.path.exists(os.path.join(d, "DONE")):
        missing.append((rid, "prep")); continue
    snps = pd.read_csv(os.path.join(d, "snps.csv"))
    # `int(snps.BP.notna() and ...)` was here and raises
    # "truth value of a Series is ambiguous" on every block -- this script had
    # never been run. The chromosome comes from blocks.csv, not from the SNPs.
    base = dict(region_id=rid, chrom=int(chrom_of[rid]))

    src = {
        "fb_xregion":    (os.path.join(ROOT, "joint_out", rid, "pip.csv"), "csv_named"),
        "funmap":        (os.path.join(ROOT, "funmap_out", rid, "pip.csv"), "csv_headerless"),
        "polyfun_sniff": (os.path.join(ROOT, "polyfun_out_sniff", f"{rid}.txt.gz"), "finemapper"),
        "polyfun_combined": (os.path.join(ROOT, "polyfun_out_combined", f"{rid}.txt.gz"), "finemapper"),
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

    # fb_xregion QC. sum(PIP) is the posterior mean number of causal variants,
    # so a value above the model's own prior cap (n_caus, default 10) is
    # self-contradictory and the region's PIPs are not a posterior.
    #
    # It happens when the variational posterior saturates its sparsity cap
    # (sparse_concrete, default 51). mean_memo is then the log-BF of a
    # ~51-variant configuration, so every configuration the credible-set sweep
    # proposes lands more than 745 nats below it, is stored as exactly 0.0 by
    # abf(), and cannot shift the normaliser. The sweep runs and is inert, and
    # the surviving PIP is a variational selection frequency: 1.0 for whatever
    # the sampler always picked, 0 for the true lead, with |Z| never entering.
    # Measured on the UKBB cholesterol run: 14 of 254 regions, all 14 placing
    # PIP 1.0 on a variant with |Z| < 2 while the lead -- |Z| up to 13.05 --
    # got PIP 0. The likelihood gain per added variant scales with n while the
    # per-variant prior penalty is fixed near log(p0), so this binds at
    # n ~ 3e5 and not at the n ~ 1e3-1e4 the default was tuned on.
    #
    # Flagged, never dropped: which regions to exclude is the analyst's call.
    N_CAUS = 10
    sumpip = df.groupby(["method", "region_id"]).pip.transform("sum")
    df["valid"] = ~((df.method == "fb_xregion") & (sumpip > N_CAUS))

    bad = (df[~df.valid].groupby("region_id").pip.agg(["sum", "max"])
             .sort_values("sum", ascending=False))
    if len(bad):
        print(f"\nQC: {len(bad)} fb_xregion regions have sum(PIP) > n_caus={N_CAUS} "
              f"and are marked valid=False (not dropped):")
        for rid, r in bad.iterrows():
            print(f"  {rid:20s} sum(PIP)={r['sum']:8.2f}  max(PIP)={r['max']:.4f}")

    df.to_parquet(os.path.join(RES, "pips.parquet"), index=False)
    print(f"wrote {len(df):,} rows for {df.region_id.nunique()} regions, "
          f"methods: {sorted(df.method.unique())}")
    s = (df.groupby(["method", "region_id"]).pip.agg(["sum", "max"])
           .groupby("method").agg(mean_sumPIP=("sum", "mean"), mean_maxPIP=("max", "mean"),
                                  regions=("sum", "size")))
    print(s.round(3).to_string())
else:
    print("no results found")

if excluded:
    print(f"\n{len(excluded)} blocks excluded by the unannotated-gap rule, not collected:")
    print("  " + ", ".join(excluded))
if missing:
    print(f"\n{len(missing)} MISSING outputs among the analysed blocks (first 20):")
    for rid, why in missing[:20]:
        print(f"  {rid:20s} {why}")

# annotation weights: fb_xregion writes feature_importance.csv per region
# (identical across regions - one shared head). Funmap writes its weights
# into sets.txt / its own output; collected here if present.
fi = glob.glob(os.path.join(ROOT, "joint_out", "*", "feature_importance.csv"))
if fi:
    w = pd.read_csv(fi[0])
    # Annotation names, in the column order the methods saw. annot_names.txt
    # was only ever READ by this script and never written by anything, so the
    # importances would have come out as bare indices and the pre-registered
    # RNAGTEX-PC5 check would have been impossible. prep_block.py writes the
    # names as the header of funmap.anno.txt (SNP + the 102 PCs, in acols
    # order), so take them from there, and check every block agrees.
    names = None
    npath = os.path.join(ROOT, "annot_names.txt")
    if os.path.exists(npath):
        names = [l.strip() for l in open(npath) if l.strip()]
    else:
        heads = glob.glob(os.path.join(ROOT, "regions", "*", "funmap.anno.txt"))
        seen = {}
        for h in heads[:40]:                     # a sample is enough to detect disagreement
            with open(h) as fh:
                cols = tuple(fh.readline().rstrip("\n").split("\t")[1:])
            seen.setdefault(cols, []).append(os.path.basename(os.path.dirname(h)))
        if len(seen) == 1:
            names = list(next(iter(seen)))
            with open(npath, "w") as fh:
                fh.write("\n".join(names) + "\n")
            print(f"\nderived {len(names)} annotation names from funmap.anno.txt "
                  f"-> annot_names.txt")
        elif len(seen) > 1:
            print(f"\nWARNING annotation column order differs between blocks "
                  f"({len(seen)} distinct orders); importances left unnamed")
    if names and len(names) == len(w):
        w["annotation"] = [names[i] for i in w.annotation_index]
    w.to_csv(os.path.join(RES, "annotation_weights.csv"), index=False)
    print(f"\nfb_xregion annotation importances -> results/annotation_weights.csv")
    print(w.head(10).to_string(index=False))
