"""Build the small local bundle: everything needed to plot and to run the
positive-control check, without moving the 15.4M-row pips.parquet.

Writes results/export/. Nothing here recomputes a PIP -- it only selects,
reshapes and annotates what 06_collect.py already produced.

  region_method_summary.csv  per (region, method), copied
  annotation_weights.csv     the 102 LassoNet skip coefficients, copied
  ld_conditioning.csv        per-block min eigenvalue / floor / r-ceiling, copied
  pips_wide.csv.gz           one row per variant, one column per method, for
                             variants that any method or the data makes
                             interesting (see PIP_KEEP / ZKEEP)
  validated_pips.csv         every method's PIP and within-block RANK on each
                             curated positive control, with the LD context
                             needed to say whether it is identifiable at all
  confident_calls.csv        every (method, variant) with PIP >= 0.95, for
                             the coding-consequence enrichment check

The two-pass structure is deliberate. A one-pass row filter keeps the row that
triggered it and drops the SAME variant's rows for the other methods, which
silently turns "funmap gave it 0.0003" into a missing value. Pass 1 collects
the qualifying variant ids and the blocks holding a positive control; pass 2
keeps every method's row for those variants, and every row of those blocks so
the ranks below are exact rather than ranks within a thinned table.
"""
import os
import numpy as np, pandas as pd
import pyarrow as pa, pyarrow.parquet as pq, pyarrow.compute as pc

R   = os.environ.get("APP_ROOT", os.getcwd())
RES = os.path.join(R, "results")
OUT = os.path.join(RES, "export"); os.makedirs(OUT, exist_ok=True)
PQ  = os.path.join(RES, "pips.parquet")

# PIP_KEEP is well below any threshold we report on, so nothing that matters
# is thinned away; ZKEEP keeps the marginal evidence visible even where no
# method put mass on the variant.
PIP_KEEP, ZKEEP, CONF = 0.01, 4.0, 0.95

val  = pd.read_csv(os.path.join(R, "validated_variants.csv"))
want = pa.array(sorted(set(val.rsid.astype(str))))
print(f"positive controls: {len(val)}")

# ---- pass 1: which variants qualify, and which blocks hold a control ----
pf = pq.ParquetFile(PQ)
keep_snps, val_regions = set(), set()
for batch in pf.iter_batches(batch_size=2_000_000, columns=["region_id", "SNP", "Z", "pip"]):
    t = pa.Table.from_batches([batch])
    hit = pc.is_in(t["SNP"], value_set=want)
    m = pc.or_(pc.greater_equal(t["pip"], PIP_KEEP),
               pc.greater_equal(pc.abs(t["Z"]), ZKEEP))
    keep_snps.update(t.filter(pc.or_(m, hit)).column("SNP").to_pylist())
    val_regions.update(t.filter(hit).column("region_id").to_pylist())
print(f"pass 1: {len(keep_snps):,} variants qualify; "
      f"{len(val_regions)} blocks hold a positive control")

# ---- pass 2 -------------------------------------------------------------
vs, vr = pa.array(sorted(keep_snps)), pa.array(sorted(val_regions))
parts = []
for batch in pf.iter_batches(batch_size=2_000_000):
    t = pa.Table.from_batches([batch])
    parts.append(t.filter(pc.or_(pc.is_in(t["SNP"], value_set=vs),
                                 pc.is_in(t["region_id"], value_set=vr))))
full = pa.concat_tables(parts).to_pandas()
methods = sorted(full.method.unique())
print(f"pass 2: {len(full):,} rows, {len(methods)} methods")

# The thinned table is the qualifying variants only; the control blocks were
# pulled in whole purely so the ranks below are exact.
df = full[full.SNP.isin(keep_snps)]

wide = (df.pivot_table(index=["region_id", "chrom", "BP", "SNP", "Z"],
                       columns="method", values="pip").reset_index())
wide.columns.name = None
wide.to_csv(os.path.join(OUT, "pips_wide.csv.gz"), index=False,
            compression="gzip", float_format="%.6g")
print(f"pips_wide.csv.gz: {len(wide):,} variants x {len(methods)} methods")

conf = df[df.pip >= CONF][["region_id", "chrom", "BP", "SNP", "Z", "method", "pip"]]
conf.to_csv(os.path.join(OUT, "confident_calls.csv"), index=False)
print(f"confident_calls.csv: {len(conf):,} rows; " +
      ", ".join(f"{m} {int((conf.method == m).sum())}" for m in methods))

# ---- the positive controls, with LD context -----------------------------
# A PIP on a validated variant is only interpretable next to how many other
# variants in the block are indistinguishable from it. rs6511720 sits in
# complete LD with three other intron-1 SNPs, so no method can single it out
# and a PIP near 1/k there is the CORRECT answer, not a miss.
rows = []
for _, v in val.iterrows():
    g = full[full.SNP == v.rsid]
    base = dict(rsid=v.rsid, gene=v.gene, tier=v.tier, identifiable=v.identifiable)
    if g.empty:
        rows.append(dict(base, status="not in the analysed blocks")); continue
    rid  = g.region_id.iloc[0]
    snps = pd.read_csv(os.path.join(R, "regions", rid, "snps.csv"))
    i  = int(np.flatnonzero(snps.SNP.values == v.rsid)[0])
    aZ = snps.Z.abs().to_numpy(); k = int(np.argmax(aZ))
    Rld = np.load(os.path.join(R, "regions", rid, "locus_ld.npy"), mmap_mode="r")
    r2 = np.asarray(Rld[i], dtype=np.float64) ** 2      # one row, not 578 MB
    rec = dict(base, status="ok", region_id=rid, p=len(snps), absZ=float(aZ[i]),
               is_lead=bool(i == k), lead_snp=snps.SNP.iloc[k],
               lead_absZ=float(aZ[k]), r2_with_lead=float(r2[k]),
               n_r2_ge_99=int((r2 >= 0.99).sum()),
               n_r2_ge_95=int((r2 >= 0.95).sum()))
    blk = full[full.region_id == rid]
    for m in methods:
        h = blk[blk.method == m]
        p = float(h.loc[h.SNP == v.rsid, "pip"].iloc[0])
        rec[f"pip_{m}"]  = p
        rec[f"rank_{m}"] = int((h.pip > p).sum()) + 1
    rows.append(rec)
vp = pd.DataFrame(rows)
vp.to_csv(os.path.join(OUT, "validated_pips.csv"), index=False)
print(f"validated_pips.csv: {len(vp)} controls, "
      f"{int((vp.status == 'ok').sum())} present in the analysed blocks")

for f in ["region_method_summary.csv", "annotation_weights.csv", "ld_conditioning.csv"]:
    src = os.path.join(RES, f)
    if os.path.exists(src):
        pd.read_csv(src).to_csv(os.path.join(OUT, f), index=False)

print(f"\nexport written to {OUT}")
for f in sorted(os.listdir(OUT)):
    print(f"  {os.path.getsize(os.path.join(OUT, f))/1e6:8.2f} MB  {f}")
