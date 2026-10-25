#!/usr/bin/env python
"""Audit, and apply, the block-set rule.

THE RULE (the user's decision, 2026-10-25): a block is analysed only if it
contains no contiguous stretch longer than PREP_MAX_ANNOT_GAP_BP (default
100 kb) in which GWAS variants exist but none are annotated.

Why a gap and not a retained fraction. The variants Borzoi does not annotate
are present in both the GWAS and the LD panel; prep drops them only for want
of annotations. So the methods receive a locus with a contiguous chunk
missing, in which the retained variants flanking the hole appear mutually
LD-linked without the intervening ones. If the causal variant lies in the
hole, its proxies inside the retained set carry the signal and the credible
set excludes the truth CONFIDENTLY -- the same class of failure as feeding an
RSS model r = 1 pairs with discordant z. A locus thinned uniformly to 60% is
fine; a locus missing a contiguous 1 Mb is a truncated search space, and the
retained fraction does not tell the two apart. It kept three blocks whose
holes were 1.00, 0.72 and 0.18 Mb because their overall coverage read
58-72%.

Measured over the 265 blocks: 254 have a largest gap of ~0 (a handful of
variants at most) and 11 have holes from 0.18 Mb up to the entire span.
Nothing falls between 0.00 and 0.18 Mb, so 100 kb sits in an empty interval
and the threshold is not doing fine-grained work. All the holes lie at three
regions Borzoi does not cover: the start of the extended MHC, the chr8p23.1
inversion, and the chr11 centromere.

An excluded block gets an EXCLUDED marker holding its numbers and reason, in
place of DONE. Every stage keys on DONE, so the manifests shed it; the marker
is what lets preflight and the stages tell a deliberate exclusion from a prep
that failed.

  python apply_block_filter.py            # audit only
  python apply_block_filter.py --apply    # write the EXCLUDED markers
"""
import os, sys, csv, json
import numpy as np, pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prep_block import longest_unannotated_gap, MAX_GAP_BP     # one definition

ROOT = os.environ.get("APP_ROOT", os.path.dirname(os.path.abspath(__file__)))
APPLY = "--apply" in sys.argv

rows = list(csv.DictReader(open(os.path.join(ROOT, "blocks.csv"))))
ss = pq.read_table(os.path.join(ROOT, "chol_sumstats.parquet"), columns=["CHR", "BP"])
chrs, bps = ss.column("CHR").to_numpy(), ss.column("BP").to_numpy()
by_chr = {c: np.sort(bps[chrs == c]) for c in sorted({int(r["chrom"]) for r in rows})}

print(f"maximum contiguous unannotated gap: {MAX_GAP_BP/1e3:,.0f} kb"
      f"   (PREP_MAX_ANNOT_GAP_BP)")

audit = []
for r in rows:
    rid, c = r["region_id"], int(r["chrom"])
    lo, hi = int(r["block_start"]), int(r["block_end"])
    d = os.path.join(ROOT, "regions", rid)
    a = by_chr[c]
    span = a[(a >= lo) & (a <= hi)]
    snps = os.path.join(d, "snps.csv")
    if os.path.exists(snps):
        kept_bp = pd.read_csv(snps)["BP"].astype(np.int64).values
    else:
        kept_bp = np.array([], dtype=np.int64)     # never prepped, or p == 0
    gap_n, gap_bp = longest_unannotated_gap(span, kept_bp)
    cov = (len(kept_bp) / len(span)) if len(span) else 0.0
    audit.append((rid, c, lo, hi, len(kept_bp), len(span), cov, gap_n, gap_bp))

ok = [t for t in audit if t[8] <= MAX_GAP_BP]
bad = [t for t in audit if t[8] > MAX_GAP_BP]
gaps = np.array([t[8] for t in ok], dtype=float)
print(f"  blocks in blocks.csv : {len(audit)}")
print(f"  analysed             : {len(ok)}")
print(f"  excluded             : {len(bad)}")
if len(gaps):
    print(f"  largest gap among those analysed: median {np.median(gaps)/1e3:,.1f} kb, "
          f"max {gaps.max()/1e3:,.1f} kb")

print("\nexcluded:")
for rid, c, lo, hi, k, n, cov, gn, gb in sorted(bad, key=lambda t: -t[8]):
    print(f"  {rid:<16} chr{c}:{lo:,}-{hi:,}  gap {gb/1e6:5.2f} Mb "
          f"({gn:>6,} variants)  retained {k:>6,}/{n:<6,} = {cov:5.1%}")
if len(gaps):
    worst = max(ok, key=lambda t: t[8])
    print(f"\nclosest ANALYSED block to the threshold:")
    print(f"  {worst[0]}  gap {worst[8]/1e3:,.1f} kb ({worst[7]:,} variants), "
          f"retained {worst[6]:.1%}")

if not APPLY:
    print("\naudit only; re-run with --apply to write the EXCLUDED markers")
    sys.exit(0)

print()
for rid, c, lo, hi, k, n, cov, gn, gb in bad:
    d = os.path.join(ROOT, "regions", rid)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "EXCLUDED"), "w") as fh:
        json.dump({"region_id": rid, "chrom": c, "block_start": lo, "block_end": hi,
                   "p_retained": k, "gwas_variants_in_span": n,
                   "annotation_coverage": round(cov, 6),
                   "largest_unannotated_gap_bp": gb,
                   "largest_unannotated_gap_variants": gn,
                   "max_unannotated_gap_bp": MAX_GAP_BP,
                   "reason": "a contiguous stretch of this block has GWAS and "
                             "LD-panel variants but no Borzoi annotations, so the "
                             "analysed locus would be a truncated search space "
                             "rather than a thinned one"}, fh, indent=2)
        fh.write("\n")
    done = os.path.join(d, "DONE")
    if os.path.exists(done):
        os.remove(done)
    print(f"  excluded {rid}  gap {gb/1e6:.2f} Mb")

restored = 0
for rid, *_ in ok:
    excl = os.path.join(ROOT, "regions", rid, "EXCLUDED")
    if os.path.exists(excl):
        os.remove(excl); restored += 1
        print(f"  cleared a stale EXCLUDED on {rid}")
if restored:
    print(f"  NOTE {restored} block(s) now pass the rule but have no DONE marker; "
          f"re-run 02_prep to rebuild them")
print(f"\n{len(bad)} excluded, {len(ok)} analysed")
