"""Build fb_xregion / funmap / polyfun+susie inputs for ONE LDetect block.

Usage:  python prep_block.py <row_index>        # 0-based row of blocks.csv

Per block: match the UKBB LD panel variants to the cholesterol sumstats on
position + unordered allele pair, harmonise Z to the LD panel's allele coding,
join the Borzoi-102 annotations on the LD panel's allele ORDER, take the LD
submatrix, standardise annotations, and write the per-method inputs.

Three correctness rules are enforced here, each from a bug that actually bit:

 1. The LD row index is carried through every merge. It is NEVER re-derived
    from position: position is not unique, and a position-keyed lookup hands
    all alleles at a site one shared LD row, producing r = 1 pairs with
    discordant z - the worst possible RSS input.
 2. Annotations are matched on the ORDERED allele pair. Borzoi lists some
    indels in both orders with different values, so an unordered key
    silently duplicates rows.
 3. Only biallelic SNVs are kept and duplicate positions are dropped, which
    is what Funmap's own pipeline does, and removes the problem at source.
"""
import os, sys, json
import numpy as np, pandas as pd, scipy.sparse as sparse, pyarrow.parquet as pq

ROOT   = os.environ.get("APP_ROOT", os.path.dirname(os.path.abspath(__file__)))
BLOCKS = os.path.join(ROOT, "blocks.csv")
SUMST  = os.path.join(ROOT, "chol_sumstats.parquet")
LD_DIR = os.path.join(ROOT, "ld")
AN_DIR = os.path.join(ROOT, "annot")
OUT    = os.path.join(ROOT, "regions")
N_LIK  = int(os.environ.get("N_LIKELIHOOD", 373220))   # GWAS n; the LD panel (~337k) is a subset


def ordkey(df, a1, a2, bp="BP"):
    return df[bp].astype(str) + "_" + df[a1].str.upper() + "_" + df[a2].str.upper()


def unordkey(df, a1="A1", a2="A2", bp="BP"):
    a, b = df[a1].str.upper(), df[a2].str.upper()
    lo = np.where(a < b, a, b); hi = np.where(a < b, b, a)
    return df[bp].astype(str) + "_" + lo + "_" + hi


def biallelic_snv(df, a1="A1", a2="A2"):
    return (df[a1].str.len() == 1) & (df[a2].str.len() == 1)


def main(i):
    blocks = pd.read_csv(BLOCKS)
    r   = blocks.iloc[i]
    ch  = int(r.chrom)
    tag = str(r.region_id)
    d   = os.path.join(OUT, tag)
    if os.path.exists(os.path.join(d, "DONE")):
        print(f"{tag}: already done"); return
    os.makedirs(d, exist_ok=True)

    pre = os.path.join(LD_DIR, str(r.ld_region))
    ld  = pd.read_table(pre + ".gz", sep=r"\s+").rename(columns={
        "rsid": "SNP", "chromosome": "CHR", "position": "BP",
        "allele1": "A1", "allele2": "A2"})

    keep   = (ld.BP >= r.block_start) & (ld.BP <= r.block_end) & biallelic_snv(ld)
    ld_idx = np.where(keep.values)[0]          # row/col index into the LD matrix
    ld     = ld.loc[keep].reset_index(drop=True)
    ld["row"] = ld_idx
    ld = ld.drop_duplicates("BP", keep=False)  # ambiguous sites removed entirely

    ss = pq.ParquetFile(SUMST).read(
        columns=["SNP", "CHR", "BP", "A1", "A2", "MAF", "N", "Z"]).to_pandas()
    ss = ss[(ss.CHR == ch) & biallelic_snv(ss)]
    ss = ss.drop_duplicates("BP", keep=False)

    ld["k"] = unordkey(ld); ss = ss.assign(k=unordkey(ss))
    m = ld.merge(ss[["k", "A1", "A2", "Z", "N", "MAF"]], on="k", how="inner",
                 suffixes=("_ld", "_ss")).drop_duplicates("k")
    flip = m.A1_ld.str.upper() != m.A1_ss.str.upper()
    m.loc[flip, "Z"] = -m.loc[flip, "Z"]       # harmonise Z to the LD allele coding

    an = pq.ParquetFile(os.path.join(AN_DIR, f"borzoi.{ch}.annot.parquet")).read().to_pandas()
    an = an[an.CHR == ch]
    an["ko"] = ordkey(an, "A1", "A2")
    acols = [c for c in an.columns if "-PC" in c]
    m["ko"] = ordkey(m, "A1_ld", "A2_ld")
    m = m.merge(an[["ko"] + acols].drop_duplicates("ko"), on="ko", how="inner")

    rows = m.row.values
    assert len(set(rows)) == len(rows), f"{tag}: two retained variants share one LD row"
    R = sparse.load_npz(pre + ".npz").toarray(); R = R + R.T
    np.fill_diagonal(R, 1.0)
    Rs = R[np.ix_(rows, rows)].astype(np.float32)
    del R

    A = m[acols].values.astype(np.float64)
    sd = A.std(0); A = (A - A.mean(0)) / np.where(sd > 0, sd, 1.0)
    A = np.clip(A, -10, 10)

    names = [f"{tag}_{bp}" for bp in m.BP.values]
    # fb_xregion / BEATRICE inputs
    with open(os.path.join(d, "locus.z"), "w") as f:
        for nm, z in zip(names, m.Z.values): f.write(f"{nm} {z:.6f}\n")
    with open(os.path.join(d, "locus.annot"), "w") as f:
        for nm, row in zip(names, A):
            f.write(nm + "\t" + "\t".join(f"{v:.5f}" for v in row) + "\n")
    np.save(os.path.join(d, "locus_ld.npy"), Rs)      # .npy, NOT 800 MB of text
    # funmap inputs, in the format UPSTREAM funmap_cmd.py parses (we run the
    # authors' code unmodified, so we match it rather than patching its readers):
    #   zdir : pd.read_csv(sep="\t", header=None)  -> [name, z]
    #   Adir : pd.read_csv(sep="\t")               -> HEADER, col 0 = SNP id
    #   Rdir : pd.read_csv(sep="\s+", header=None) -> written per-task to
    #          $TMPDIR from locus_ld.npy, since the text form is ~800 MB.
    pd.DataFrame({"SNP": names, "Z": m.Z.values}).to_csv(
        os.path.join(d, "funmap.z.txt"), sep="\t", index=False, header=False)
    ann = pd.DataFrame(A, columns=acols)
    ann.insert(0, "SNP", names)
    ann.to_csv(os.path.join(d, "funmap.anno.txt"), sep="\t", index=False,
               float_format="%.5f")
    m[["SNP", "BP", "A1_ld", "A2_ld", "Z", "MAF", "row"]].to_csv(
        os.path.join(d, "snps.csv"), index=False)

    meta = dict(tag=tag, chrom=ch, p=int(len(m)), n_likelihood=N_LIK,
                block=f"{int(r.block_start)}-{int(r.block_end)}",
                max_abs_z=float(np.abs(m.Z).max()), flipped=int(flip.sum()),
                min_p=float(r.min_p), n_annot=len(acols))
    json.dump(meta, open(os.path.join(d, "meta.json"), "w"), indent=1)
    open(os.path.join(d, "DONE"), "w").close()
    print(f"{tag}: p={len(m)} max|Z|={meta['max_abs_z']:.2f} flips={meta['flipped']}")


if __name__ == "__main__":
    main(int(sys.argv[1]))
