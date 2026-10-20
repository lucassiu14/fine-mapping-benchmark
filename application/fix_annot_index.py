#!/usr/bin/env python
"""Drop stray non-annotation columns from the Borzoi annotation parquet files.

chr6, chr8 and chr11 of the Sniff Borzoi release carry an extra `index`
column: a pandas artifact from a to_parquet() call without index=False. The
other 19 chromosomes do not have it. All 22 carry exactly 102 "-PC" columns,
so the annotation set itself is consistent -- only the serialisation differs.

Why it matters, and to whom:

  prep_block.py  selects annotations with [c for c in an.columns if "-PC" in c],
                 so fb_xregion and funmap never see it. Unaffected.
  00f            passes the whole parquet to compute_ldscores.py, which would
                 treat `index` as a 103rd annotation on those three
                 chromosomes.
  00b2           excludes only CHR/SNP/BP/A1/A2, so `index` would be carried
                 into the combined set as a Borzoi annotation.
  00g            is a GENOME-WIDE S-LDSC regression, so an annotation set that
                 differs by chromosome contaminates the whole prior, not just
                 the three chromosomes.

This rewrites only the affected files, keeping the original alongside as
<name>.orig, and deletes the derived files computed from them so they are
rebuilt. It is idempotent: a second run reports nothing to do.

  python fix_annot_index.py            # report only
  python fix_annot_index.py --apply    # rewrite and invalidate
"""
import os, sys
import pyarrow.parquet as pq
import pandas as pd

ROOT = os.environ.get("APP_ROOT", ".")
KEY = {"CHR", "SNP", "BP", "A1", "A2", "CM", "MAF", "base"}
APPLY = "--apply" in sys.argv

def extras(cols):
    return [c for c in cols if c not in KEY and "-PC" not in c]

bad = {}
for c in range(1, 23):
    f = os.path.join(ROOT, "annot", f"borzoi.{c}.annot.parquet")
    if not os.path.exists(f):
        print(f"  chr{c}: MISSING {f}"); continue
    names = list(pq.ParquetFile(f).schema_arrow.names)
    npc = len([x for x in names if "-PC" in x])
    if npc != 102:
        print(f"  chr{c}: FATAL {npc} '-PC' columns, expected 102. "
              f"Not a serialisation artifact - stopping.")
        sys.exit(2)
    e = extras(names)
    if e:
        bad[c] = e

if not bad:
    print("all 22 Borzoi files carry exactly 102 '-PC' columns and no extras")
    sys.exit(0)

print(f"{len(bad)} chromosomes carry extra non-annotation columns:")
for c, e in sorted(bad.items()):
    print(f"  chr{c:<3} {e}")

if not APPLY:
    print("\nre-run with --apply to rewrite these and invalidate what came from them")
    sys.exit(0)

print()
for c, e in sorted(bad.items()):
    f = os.path.join(ROOT, "annot", f"borzoi.{c}.annot.parquet")
    orig = f + ".orig"
    d = pq.ParquetFile(f).read().to_pandas()
    before = list(d.columns)
    d = d.drop(columns=e)
    assert len([x for x in d.columns if "-PC" in x]) == 102, f"chr{c}: lost a PC"
    assert len(d.columns) == len(before) - len(e)
    if not os.path.exists(orig):
        os.rename(f, orig)
    d.to_parquet(f, index=False)
    print(f"  chr{c}: dropped {e}, {len(before)} -> {len(d.columns)} columns, "
          f"{len(d):,} rows (original kept as {os.path.basename(orig)})")

# Derived files built from the un-normalised input. The prior is genome-wide,
# so it is invalidated by ANY affected chromosome, not only those three.
print("\ninvalidating derived files:")
removed = 0
for c in sorted(bad):
    for q in (os.path.join(ROOT, "annot_combined", f"combined.{c}.annot.parquet"),
              os.path.join(ROOT, "ldscores_sniff", f"annot.{c}.l2.ldscore.parquet"),
              os.path.join(ROOT, "ldscores_combined", f"annot.{c}.l2.ldscore.parquet")):
        if os.path.exists(q):
            os.remove(q); removed += 1
            print(f"  removed {os.path.relpath(q, ROOT)}")
for s in ("sniff", "combined"):
    d = os.path.join(ROOT, f"pf_prior_{s}")
    if os.path.isdir(d):
        for fn in sorted(os.listdir(d)):
            os.remove(os.path.join(d, fn)); removed += 1
        if os.listdir(d) == []:
            print(f"  cleared pf_prior_{s}/ (S-LDSC is genome-wide)")
if removed == 0:
    print("  nothing derived had been built yet")
print("\ndone. Resubmit 00f and 00g; finished clean chromosomes are skipped.")
