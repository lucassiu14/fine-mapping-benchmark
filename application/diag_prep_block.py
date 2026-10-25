#!/usr/bin/env python
"""Report where a block's variants vanish in prep_block.py.

Six of the 265 blocks prepped to p = 0 and still got a DONE marker, because
every assertion in prep_block.py passes trivially at zero rows:

  ldEUR_06_0020  chr6_24000001_27000001   lead|Z| 10.9    (.npz)
  ldEUR_08_0012  chr8_7000001_10000001    lead|Z| 18.053  (.npz)
  ldEUR_08_0013  chr8_8000001_11000001    lead|Z| 5.522   (.npz2)
  ldEUR_11_0029  chr11_47000001_50000001  lead|Z| 7.197   (.npz2)
  ldEUR_11_0030  chr11_50000001_53000001  lead|Z| 4.959   (.npz2)
  ldEUR_11_0031  chr11_53000001_56000001  lead|Z| 5.249   (.npz2)

Four are the windows published as .npz2; the other two are chr8p23.1 (the
polymorphic inversion / defensin region) and the approach to the MHC.

This replays prep_block.py's filters in order, reusing ITS helper functions
rather than reimplementing them, and prints how many variants survive each
step. It reads only the .npz's stored `shape` array, never densifying the
matrix, so it is safe to run anywhere.

  python diag_prep_block.py ldEUR_08_0012 ldEUR_08_0013
  python diag_prep_block.py --zero          # every block currently at p == 0
"""
import os, sys
import numpy as np, pandas as pd, pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import prep_block as pb                                      # noqa: E402

blocks = pd.read_csv(pb.BLOCKS)

args = [a for a in sys.argv[1:] if not a.startswith("-")]
if "--zero" in sys.argv:
    args = []
    for _, r in blocks.iterrows():
        f = os.path.join(pb.OUT, str(r.region_id), "locus_ld.npy")
        if os.path.exists(f) and int(np.load(f, mmap_mode="r").shape[0]) == 0:
            args.append(str(r.region_id))
    print(f"blocks at p == 0: {args}\n")
if not args:
    sys.exit("name at least one region_id, or pass --zero")

# sumstats once
ss_all = pq.ParquetFile(pb.SUMST).read(
    columns=["SNP", "CHR", "BP", "A1", "A2", "MAF", "N", "Z"]).to_pandas()

for rid in args:
    sel = blocks[blocks.region_id == rid]
    if sel.empty:
        print(f"!! {rid} is not in blocks.csv"); continue
    r = sel.iloc[0]
    ch = int(r.chrom)
    pre = os.path.join(pb.LD_DIR, str(r.ld_region))
    print("=" * 72)
    print(f"{rid}  chr{ch}:{int(r.block_start):,}-{int(r.block_end):,}"
          f"  window {r.ld_region}  lead|Z| {r.lead_abs_z}")

    # --- the LD panel's variant list and the matrix it indexes ------------
    ld0 = pd.read_table(pre + ".gz", sep=r"\s+")
    print(f"  .gz columns            : {list(ld0.columns)[:10]}")
    print(f"  .gz rows               : {len(ld0):,}")
    with np.load(pre + ".npz") as z:
        keys = list(z.keys())
        shape = tuple(int(x) for x in z["shape"]) if "shape" in z else None
    print(f"  .npz keys              : {keys}")
    print(f"  .npz shape             : {shape}")
    if shape and shape[0] != len(ld0):
        print(f"  ** MISMATCH: the matrix is {shape[0]} x {shape[1]} but the .gz lists "
              f"{len(ld0)} variants. Row indices derived from the .gz do not address "
              f"this matrix.")

    ld = ld0.rename(columns={"rsid": "SNP", "chromosome": "CHR", "position": "BP",
                             "allele1": "A1", "allele2": "A2"})
    missing = [c for c in ("SNP", "CHR", "BP", "A1", "A2") if c not in ld.columns]
    if missing:
        print(f"  ** the .gz lacks the expected columns after renaming: {missing}")
        continue

    in_range = (ld.BP >= r.block_start) & (ld.BP <= r.block_end)
    bi = pb.biallelic_snv(ld)
    print(f"  in the block's bp range: {int(in_range.sum()):,}")
    print(f"  ... and biallelic SNV  : {int((in_range & bi).sum()):,}")
    keep = in_range & bi
    ld = ld.loc[keep].reset_index(drop=True)
    ld["row"] = np.where(keep.values)[0]
    n_before = len(ld)
    ld = ld.drop_duplicates("BP", keep=False)
    print(f"  ... after dropping duplicate positions (keep=False): {len(ld):,}"
          f"   (removed {n_before - len(ld):,})")
    if len(ld) == 0:
        print("  -> VANISHED at the LD-panel stage"); continue

    ss = ss_all[(ss_all.CHR == ch) & pb.biallelic_snv(ss_all)]
    ss = ss.drop_duplicates("BP", keep=False)
    in_blk = ss[(ss.BP >= r.block_start) & (ss.BP <= r.block_end)]
    print(f"  sumstats on chr{ch}      : {len(ss):,}  (in this block: {len(in_blk):,})")

    ld = ld.assign(k=pb.unordkey(ld))
    ss = ss.assign(k=pb.unordkey(ss))
    m = ld.merge(ss[["k", "A1", "A2", "Z", "N", "MAF"]], on="k", how="inner",
                 suffixes=("_ld", "_ss")).drop_duplicates("k")
    print(f"  ... merged with sumstats (unordered allele key): {len(m):,}")
    if len(m) == 0:
        ldk, ssk = set(ld.k), set(in_blk.assign(k=pb.unordkey(in_blk)).k)
        print(f"  -> VANISHED at the sumstats merge. Key overlap in-block: "
              f"{len(ldk & ssk)}")
        print(f"     example LD keys      : {sorted(ldk)[:3]}")
        print(f"     example sumstat keys : {sorted(ssk)[:3]}")
        continue

    an = pq.ParquetFile(os.path.join(pb.AN_DIR, f"borzoi.{ch}.annot.parquet")).read().to_pandas()
    an = an[an.CHR == ch]
    an["ko"] = pb.ordkey(an, "A1", "A2")
    m["ko"] = pb.ordkey(m, "A1_ld", "A2_ld")
    acols = [c for c in an.columns if "-PC" in c]
    m2 = m.merge(an[["ko"] + acols].drop_duplicates("ko"), on="ko", how="inner")
    print(f"  ... merged with annotations (ordered allele key): {len(m2):,}")
    if len(m2) == 0:
        print("  -> VANISHED at the annotation merge")
        print(f"     example merged keys    : {list(m.ko[:3])}")
        print(f"     example annotation keys: {list(an.ko[:3])}")
        continue
    print(f"  -> final p = {len(m2):,}")
