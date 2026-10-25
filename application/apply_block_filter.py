#!/usr/bin/env python
"""Audit annotation coverage per block, and exclude blocks that fall short.

THE RULE (the user's decision, 2026-10-25): a block enters the analysis only
if it retains at least MIN_COVERAGE of the GWAS variants lying in its span
after prep's annotation merge. Default 0.25.

Why coverage rather than a size floor: the cause is annotation coverage, so
the criterion names the cause. It also generalises -- on another GWAS or
annotation set a small but fully annotated block should not be dropped for
being small.

What it excludes here, measured: Borzoi has NO annotations in three regions,
so six blocks retained zero variants and a seventh retained 176 of 9,011.

  ldEUR_06_0020  chr6:25.68-26.79 Mb   0 / 6,632    start of the extended MHC
  ldEUR_08_0012  chr8:9.15-9.64 Mb     0 / 3,903    chr8p23.1 inversion
  ldEUR_08_0013  chr8:9.64-10.46 Mb    0 / 6,617    chr8p23.1 inversion
  ldEUR_11_0029  chr11:47.0-49.87 Mb   0 / 13,842   centromere
  ldEUR_11_0030  chr11:50.0-53.0 Mb    0 / 9,182    centromere
  ldEUR_11_0031  chr11:54.70-55.08 Mb  0 / 2,338    centromere
  ldEUR_08_0011  chr8:7.15-9.15 Mb   176 / 9,011    abuts the chr8p23.1 gap

Verified not to be an artifact of this pipeline: the .orig backups of the
rewritten chr6/8/11 annotation files have identical coverage (0 in range,
same row counts), so fix_annot_index.py did not cause it, and the .npz
shapes match their .gz row counts exactly, so the LD indexing is sound.

Healthy blocks retain 70-90%, and the gap below them is stark -- the next
smallest p after 176 is 2,759 -- so 0.25 sits far from any real block.

An excluded block gets an EXCLUDED marker carrying the numbers, in place of
DONE. Every stage keys on DONE, so this removes it from the manifests; the
marker is what lets preflight and the stages tell a deliberate exclusion
apart from a prep that failed.

  python apply_block_filter.py            # audit only
  python apply_block_filter.py --apply    # write the EXCLUDED markers
"""
import os, sys, csv, json
import numpy as np
import pyarrow.parquet as pq

ROOT = os.environ.get("APP_ROOT", os.path.dirname(os.path.abspath(__file__)))
MIN_COV = float(os.environ.get("PREP_MIN_ANNOT_COVERAGE", "0.25"))
APPLY = "--apply" in sys.argv

rows = list(csv.DictReader(open(os.path.join(ROOT, "blocks.csv"))))
ss = pq.read_table(os.path.join(ROOT, "chol_sumstats.parquet"), columns=["CHR", "BP"])
chrs = ss.column("CHR").to_numpy()
bps = ss.column("BP").to_numpy()
by_chr = {c: np.sort(bps[chrs == c]) for c in sorted({int(r["chrom"]) for r in rows})}

print(f"minimum annotation coverage: {MIN_COV:.0%}"
      f"   (PREP_MIN_ANNOT_COVERAGE)")
audit = []
for r in rows:
    rid, c = r["region_id"], int(r["chrom"])
    lo, hi = int(r["block_start"]), int(r["block_end"])
    d = os.path.join(ROOT, "regions", rid)
    f = os.path.join(d, "locus_ld.npy")
    if not os.path.exists(f):
        audit.append((rid, c, lo, hi, None, 0, None)); continue
    p = int(np.load(f, mmap_mode="r").shape[0])
    a = by_chr[c]
    n_ss = int(np.searchsorted(a, hi, "right") - np.searchsorted(a, lo, "left"))
    cov = (p / n_ss) if n_ss else 0.0
    audit.append((rid, c, lo, hi, p, n_ss, cov))

ok = [t for t in audit if t[6] is not None and t[6] >= MIN_COV]
low = [t for t in audit if t[6] is not None and t[6] < MIN_COV]
absent = [t for t in audit if t[6] is None]

print(f"  blocks in blocks.csv : {len(audit)}")
print(f"  meeting the rule     : {len(ok)}")
print(f"  below the rule       : {len(low)}")
if absent:
    print(f"  not prepped at all   : {len(absent)}  {[t[0] for t in absent][:5]}")
covs = np.array([t[6] for t in ok])
if len(covs):
    print(f"  coverage among those kept: min {covs.min():.1%}  "
          f"median {np.median(covs):.1%}  max {covs.max():.1%}")

print("\nbelow the rule:")
for rid, c, lo, hi, p, n_ss, cov in sorted(low, key=lambda t: t[6]):
    print(f"  {rid}  chr{c}:{lo:,}-{hi:,}  p={p:<6} of {n_ss:<7} GWAS variants"
          f"  = {cov:.2%}")

# closest kept block, so it is visible whether the threshold is near anything
if len(covs):
    nearest = min(ok, key=lambda t: t[6])
    print(f"\nclosest KEPT block to the threshold:")
    print(f"  {nearest[0]}  p={nearest[4]:,} of {nearest[5]:,} = {nearest[6]:.1%}")

if not APPLY:
    print("\naudit only; re-run with --apply to write the EXCLUDED markers")
    sys.exit(0)

print()
for rid, c, lo, hi, p, n_ss, cov in low:
    d = os.path.join(ROOT, "regions", rid)
    done, excl = os.path.join(d, "DONE"), os.path.join(d, "EXCLUDED")
    payload = {"region_id": rid, "chrom": c, "block_start": lo, "block_end": hi,
               "p_retained": p, "gwas_variants_in_span": n_ss,
               "annotation_coverage": round(cov, 6),
               "min_annotation_coverage": MIN_COV,
               "reason": "annotation coverage below the minimum; Borzoi has no "
                         "annotations covering this span"}
    with open(excl, "w") as fh:
        json.dump(payload, fh, indent=2)
        fh.write("\n")
    if os.path.exists(done):
        os.remove(done)
    print(f"  excluded {rid}  ({cov:.2%})")

# and clear any stale EXCLUDED on blocks that now pass
restored = 0
for rid, c, lo, hi, p, n_ss, cov in ok:
    excl = os.path.join(ROOT, "regions", rid, "EXCLUDED")
    if os.path.exists(excl):
        os.remove(excl); restored += 1
        print(f"  cleared a stale EXCLUDED on {rid} ({cov:.1%})")

print(f"\n{len(low)} excluded, {len(ok)} analysed"
      + (f", {restored} stale marker(s) cleared" if restored else ""))
