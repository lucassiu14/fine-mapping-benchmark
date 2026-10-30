#!/usr/bin/env python3
"""What do the four methods say at known cholesterol variants?

The rsIDs are looked up in results/pips.parquet, whose SNP column is the UKBB
LD panel's own `rsid` (prep_block.py renames it, 06_collect.py passes it
through). So chromosome, position and Z come from the data; nothing here
depends on coordinates recalled from outside it. The GENE LABELS are mine and
are not verified by this script - they are conventional assignments for these
variants, printed only to make the table readable.

Reported per variant: each method's PIP, that method's maximum PIP anywhere in
the same region, and the variant's PIP rank within the region. A method can
agree on a lead without being confident, and can be confident about the wrong
variant, so rank and max are both needed to read the result.
"""
import argparse, os
import numpy as np
import pandas as pd

DEFAULT = {
    "rs11591147": "PCSK9 (R46L)",
    "rs12740374": "SORT1/CELSR2 1p13",
    "rs429358":   "APOE e4",
    "rs7412":     "APOE e2",
    "rs6511720":  "LDLR",
}

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=os.environ.get("APP_ROOT", "."))
ap.add_argument("--pips", default=None, help="default <root>/results/pips.parquet")
ap.add_argument("--snps", default=None,
                help="comma-separated rsIDs; default is a small canonical set")
a = ap.parse_args()

pq = a.pips or os.path.join(a.root, "results", "pips.parquet")
df = pd.read_parquet(pq)
want = ([s.strip() for s in a.snps.split(",")] if a.snps else list(DEFAULT))
labels = DEFAULT if not a.snps else {s: "" for s in want}

print(f"{len(df):,} rows, {df.region_id.nunique()} regions, "
      f"methods: {sorted(df.method.unique())}\n")

for rs in want:
    hit = df[df.SNP == rs]
    if hit.empty:
        print(f"{rs:<12} {labels.get(rs,''):<22} NOT PRESENT in the analysed blocks")
        print("   (excluded by the unannotated-gap rule, absent from the LD panel,")
        print("    or not matched to the sumstats by position+alleles)\n")
        continue
    rid = hit.region_id.iloc[0]
    reg = df[df.region_id == rid]
    r0 = hit.iloc[0]
    print(f"{rs:<12} {labels.get(rs,''):<22} chr{int(r0.chrom)}:{int(r0.BP)}  "
          f"Z={r0.Z:+.2f}  region {rid} ({reg.SNP.nunique():,} variants)")
    # the region's strongest association, for context on what "the lead" is
    lead = reg.loc[reg.Z.abs().idxmax()]
    print(f"{'':14}region lead by |Z|: {lead.SNP} at chr{int(lead.chrom)}:"
          f"{int(lead.BP)} Z={lead.Z:+.2f}")
    print(f"{'':14}{'method':<18} {'PIP':>8} {'region max':>11} {'rank':>6} {'of':>7}")
    for meth in sorted(df.method.unique()):
        mr = reg[reg.method == meth]
        v = mr[mr.SNP == rs]
        if v.empty:
            print(f"{'':14}{meth:<18} {'absent':>8}")
            continue
        pip = float(v.pip.iloc[0])
        rank = int((mr.pip > pip).sum()) + 1
        print(f"{'':14}{meth:<18} {pip:>8.4f} {mr.pip.max():>11.4f} "
              f"{rank:>6} {len(mr):>7}")
    print()
