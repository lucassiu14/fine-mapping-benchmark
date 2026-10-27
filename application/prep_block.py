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
# THE BLOCK-SET RULE. A block is analysed only if it contains no contiguous
# stretch longer than this in which GWAS variants exist but none are
# annotated. Default 100 kb.
#
# Why a gap rather than a retained fraction: the variants Borzoi does not
# annotate are present in both the GWAS and the LD panel, so dropping them
# hands the methods a locus with a contiguous chunk missing, in which the
# retained variants flanking the hole look mutually LD-linked without the
# intervening ones. If the causal variant lies in the hole, its proxies inside
# the retained set carry the signal and the credible set excludes the truth
# CONFIDENTLY. A locus thinned uniformly to 60% is fine; a locus missing a
# contiguous 1 Mb is a truncated search space, and the retained fraction does
# not distinguish the two.
#
# Measured over the 265 blocks: 254 have a largest gap of ~0 (at most a
# handful of variants), and 11 have holes from 0.18 Mb to the entire span.
# Nothing lies between 0.00 and 0.18 Mb, so 100 kb sits in an empty interval.
# The holes are all at three regions Borzoi does not cover at all: the start
# of the extended MHC, the chr8p23.1 inversion, and the chr11 centromere.
#
# This also replaced an earlier retained-fraction rule, which kept three
# blocks whose holes were 1.00, 0.72 and 0.18 Mb because their overall
# coverage was 58-72%.
MAX_GAP_BP = int(float(os.environ.get("PREP_MAX_ANNOT_GAP_BP", "100000")))


def longest_unannotated_gap(span_bp, retained_bp):
    """Longest run of consecutive GWAS variants in a span with none retained.

    Returns (n_variants, bp_extent). Shared with apply_block_filter.py so the
    rule has one definition; a second copy would be free to drift.
    """
    span_bp = np.asarray(span_bp)
    if span_bp.size == 0:
        return 0, 0
    kept = np.isin(span_bp, np.asarray(retained_bp))
    best_n = best_bp = n = start = 0
    for i, ok in enumerate(kept):
        if ok:
            n = 0
            continue
        if n == 0:
            start = i
        n += 1
        if n > best_n:
            best_n, best_bp = n, int(span_bp[i] - span_bp[start])
    return best_n, best_bp


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
    # GWAS variants in the block's span, before any of prep's own filtering,
    # so the gap rule means the same thing for every block
    span_mask = ((ss.CHR == ch) & (ss.BP >= r.block_start) & (ss.BP <= r.block_end))
    ss_span_bp = np.sort(ss.loc[span_mask, "BP"].values)
    n_ss_span = int(span_mask.sum())
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

    # Apply the rule HERE, before the O(p^2) LD extraction below, so an
    # excluded block costs nothing and no degenerate inputs are written.
    cov = (len(m) / n_ss_span) if n_ss_span else 0.0
    gap_n, gap_bp = longest_unannotated_gap(ss_span_bp, m.BP.values)
    if gap_bp > MAX_GAP_BP:
        json.dump(dict(region_id=tag, chrom=ch,
                       block_start=int(r.block_start), block_end=int(r.block_end),
                       p_retained=int(len(m)), gwas_variants_in_span=n_ss_span,
                       annotation_coverage=round(cov, 6),
                       largest_unannotated_gap_bp=gap_bp,
                       largest_unannotated_gap_variants=gap_n,
                       max_unannotated_gap_bp=MAX_GAP_BP,
                       reason="a contiguous stretch of this block has GWAS and "
                              "LD-panel variants but no Borzoi annotations, so "
                              "the analysed locus would be a truncated search "
                              "space rather than a thinned one"),
                  open(os.path.join(d, "EXCLUDED"), "w"), indent=2)
        print(f"{tag}: EXCLUDED - largest unannotated gap {gap_bp/1e6:.2f} Mb "
              f"({gap_n:,} variants) > {MAX_GAP_BP/1e6:.2f} Mb; "
              f"retained {len(m)} of {n_ss_span} ({cov:.1%})")
        return

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
                min_p=float(r.min_p), n_annot=len(acols),
                gwas_variants_in_span=n_ss_span, annotation_coverage=round(cov, 6),
                largest_unannotated_gap_bp=gap_bp)
    json.dump(meta, open(os.path.join(d, "meta.json"), "w"), indent=1)
    # The annotation names in the exact column order every method is handed.
    # 06_collect.py needs these to label fb_xregion's feature importances --
    # without them they are bare indices and a named prediction such as
    # RNAGTEX-PC5 cannot be checked. Written once, at APP_ROOT.
    npath = os.path.join(ROOT, "annot_names.txt")
    if not os.path.exists(npath):
        with open(npath, "w") as fh:
            fh.write("\n".join(acols) + "\n")
    open(os.path.join(d, "DONE"), "w").close()
    print(f"{tag}: p={len(m)} max|Z|={meta['max_abs_z']:.2f} flips={meta['flipped']}")


if __name__ == "__main__":
    main(int(sys.argv[1]))
