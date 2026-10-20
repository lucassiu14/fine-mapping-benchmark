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
    orig, tmp = f + ".orig", f + ".tmp"

    # Work on Arrow row groups, not a pandas round-trip. Two reasons:
    #
    #  1. `index` is recorded in the parquet footer's pandas metadata, so
    #     to_pandas() restores it as the DataFrame INDEX rather than a
    #     column. d.drop(columns=["index"]) therefore raises KeyError even
    #     though the column is plainly there in schema_arrow.names.
    #  2. these files reach 1.37 GB (chr6); reading one whole into pandas on
    #     a login node is inconsiderate at best. One row group at a time
    #     keeps the footprint small.
    #
    # replace_schema_metadata(None) drops the stale pandas metadata, which
    # still names `index` as the index column and would otherwise be carried
    # into the rewritten file.
    pfh = pq.ParquetFile(f)
    keep = [x for x in pfh.schema_arrow.names if x not in e]
    n_pc = len([x for x in keep if "-PC" in x])
    assert n_pc == 102, f"chr{c}: would keep {n_pc} PCs, not 102 - refusing"

    rows, writer = 0, None
    try:
        for i in range(pfh.num_row_groups):
            t = pfh.read_row_group(i, columns=keep).replace_schema_metadata(None)
            if writer is None:
                writer = pq.ParquetWriter(tmp, t.schema)
            writer.write_table(t)
            rows += t.num_rows
    finally:
        if writer is not None:
            writer.close()

    # Verify the rewrite BEFORE displacing the original.
    chk = pq.ParquetFile(tmp)
    got = list(chk.schema_arrow.names)
    assert chk.metadata.num_rows == pfh.metadata.num_rows, (
        f"chr{c}: row count changed {pfh.metadata.num_rows} -> {chk.metadata.num_rows}")
    assert len([x for x in got if "-PC" in x]) == 102, f"chr{c}: PC count wrong after rewrite"
    survived = [x for x in got if x not in KEY and "-PC" not in x]
    assert not survived, f"chr{c}: extras survived the rewrite: {survived}"
    del chk

    if not os.path.exists(orig):
        os.rename(f, orig)
    os.replace(tmp, f)
    print(f"  chr{c}: dropped {e}, {len(pfh.schema_arrow.names)} -> {len(got)} columns, "
          f"{rows:,} rows (original kept as {os.path.basename(orig)})")

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
