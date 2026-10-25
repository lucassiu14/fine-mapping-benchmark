"""Check a stage's inputs are complete before committing CPU to it.

    python preflight.py <stage>      # data | prep | methods | collect

Exits non-zero and names what is missing. The point is to fail in seconds
rather than discover a third of the way through a 265-task array that a
download silently truncated.
"""
import os, sys, csv, gzip

ROOT = os.environ.get("APP_ROOT", os.path.dirname(os.path.abspath(__file__)))
B = os.path.join(ROOT, "blocks.csv")
ok = True


fails = 0


def fail(msg):
    global ok, fails
    ok = False
    fails += 1
    print(f"  MISSING  {msg}")


def blocks():
    with open(B) as f:
        return list(csv.DictReader(f))


def check_data():
    print("inputs for 02_prep:")
    rows = blocks()
    # Report each group's result independently. A single failing group used to
    # suppress the one summary line at the end, so a run that failed on the LD
    # files said nothing at all about the annotations - which is exactly when
    # you want to know they are fine.
    f0 = fails
    s = os.path.join(ROOT, "chol_sumstats.parquet")
    if not os.path.exists(s) or os.path.getsize(s) < 400_000_000:
        fail(f"chol_sumstats.parquet (expect ~526 MB, got "
             f"{os.path.getsize(s)/1e6:.0f} MB)" if os.path.exists(s)
             else "chol_sumstats.parquet")
    elif fails == f0:
        print(f"  OK       sumstats, {os.path.getsize(s)/1e6:.0f} MB")
    f0 = fails
    # A truncated download is the failure mode that actually happens, and a
    # size floor does not catch it: an interrupted transfer leaves a file that
    # exists and is far over any floor worth setting. Both formats here keep
    # their index at the END of the file (zip central directory, parquet
    # footer), which is precisely what truncation destroys - so open that
    # index. Neither read touches the bulk data.
    missing_ld, broken_ld = [], []
    for w in sorted({r["ld_region"] for r in rows}):
        for ext in ("npz", "gz"):
            q = os.path.join(ROOT, "ld", f"{w}.{ext}")
            if not os.path.exists(q) or os.path.getsize(q) < 1000:
                missing_ld.append(f"{w}.{ext}")
            elif ext == "npz":
                try:
                    import zipfile
                    zipfile.ZipFile(q).namelist()
                except Exception as e:
                    broken_ld.append(f"{w}.npz ({type(e).__name__})")
    if missing_ld:
        fail(f"{len(missing_ld)} LD files, e.g. {missing_ld[:3]}")
    if broken_ld:
        fail(f"{len(broken_ld)} truncated/corrupt LD files, e.g. {broken_ld[:3]}")
    if fails == f0:
        print(f"  OK       {len({r['ld_region'] for r in rows})} LD windows")
    f0 = fails

    # Annotations. Nothing else opens these before 00f hands them to
    # compute_ldscores.py (sniff) or 00b2 merges them (combined), where a bad
    # file surfaces as an obscure error deep inside someone else's code. Also
    # report the annotation column count: a count that differs from 102 means
    # the set is not what fb_xregion and funmap are given, which matters more
    # than truncation because it would not crash anything.
    chroms = sorted({int(r["chrom"]) for r in rows})
    missing_an, broken_an, ncols, stray_an = [], [], {}, {}
    META = {"CHR", "SNP", "BP", "A1", "A2", "CM", "MAF", "base"}
    try:
        import pyarrow.parquet as pq
    except ImportError:
        pq = None
    for c in chroms:
        q = os.path.join(ROOT, "annot", f"borzoi.{c}.annot.parquet")
        if not os.path.exists(q):
            missing_an.append(c)
        elif pq is not None:
            try:
                names = list(pq.ParquetFile(q).schema_arrow.names)
                ncols.setdefault(len([x for x in names if "-PC" in x]), []).append(c)
                st = [x for x in names if x not in META and "-PC" not in x]
                if st:
                    stray_an[c] = st
            except Exception as e:
                broken_an.append(f"chr{c} ({type(e).__name__})")
    if missing_an:
        fail(f"annotations for chromosomes {missing_an}")
    if broken_an:
        fail(f"{len(broken_an)} unreadable annotation files: {broken_an[:5]}")
    # A stray column is the dangerous case: it does not crash anything, it
    # just means PolyFun regresses on one more annotation on some
    # chromosomes than others. chr6/8/11 of the Sniff release ship a pandas
    # `index` column; fb_xregion and funmap are insulated because
    # prep_block.py selects on "-PC", but 00f passes the whole file through.
    if stray_an:
        fail("non-annotation columns in the Borzoi files - run "
             "'python fix_annot_index.py --apply': "
             + "; ".join(f"chr{c} {v}" for c, v in sorted(stray_an.items())))
    if len(ncols) > 1:
        fail("the number of '-PC' annotations differs by chromosome: "
             + ", ".join(f"{n} PCs on chr {cs}" for n, cs in sorted(ncols.items())))
    elif ncols and list(ncols)[0] != 102:
        fail(f"expected 102 Borzoi PCs, found {list(ncols)[0]}")
    if fails == f0:
        print(f"  OK       {len(chroms)} annotation files"
              + (f", {list(ncols)[0]} PCs each" if len(ncols) == 1 else ""))


def check_prep():
    print("prepared blocks (inputs for 03/04/05):")
    rows = blocks()
    bad, excluded = [], []
    for r in rows:
        d = os.path.join(ROOT, "regions", r["region_id"])
        # A deliberate exclusion is not a missing prep. apply_block_filter.py
        # and prep_block.py both write EXCLUDED, carrying the numbers, when a
        # block falls below the annotation-coverage rule.
        if os.path.exists(os.path.join(d, "EXCLUDED")):
            excluded.append(r["region_id"]); continue
        if not os.path.exists(os.path.join(d, "DONE")):
            bad.append(r["region_id"]); continue
        for f in ("locus.z", "locus_ld.npy", "locus.annot", "funmap.z.txt",
                  "funmap.anno.txt", "snps.csv"):
            if not os.path.exists(os.path.join(d, f)):
                bad.append(f"{r['region_id']}/{f}")
    if bad:
        fail(f"{len(bad)} of {len(rows)} blocks incomplete, e.g. {bad[:3]}")
    else:
        print(f"  OK  all {len(rows) - len(excluded)} analysed blocks prepared")
    if excluded:
        print(f"  --  {len(excluded)} excluded by the annotation-coverage rule: "
              + ", ".join(excluded[:8]) + ("" if len(excluded) <= 8 else ", ..."))
    # The PolyFun prior, computed by 00g from our own annotations. Needed by
    # 05 only. Note S-LDSC is genome-wide, so all 22 must be present even
    # though only 22 chromosomes carry selected blocks.
    for st in ("sniff", "combined"):
        miss = [c for c in range(1, 23)
                if not os.path.exists(os.path.join(ROOT, f"pf_prior_{st}", f"chol.{c}.snpvar_constrained.gz"))]
        if miss:
            fail(f"PolyFun prior ({st}) missing for chromosomes {miss} - run "
                 f"00d/00e/00f/00g with ANNOT_SET={st} (needed by 05 only)")
        else:
            print(f"  OK  PolyFun prior ({st}) computed for all 22 chromosomes")


def check_methods():
    print("method outputs (inputs for 06_collect):")
    rows = blocks()
    for name, path, label in [
        ("fb_xregion", lambda r: os.path.join(ROOT, "joint_out", r["region_id"], "pip.csv"), "03"),
        ("funmap", lambda r: os.path.join(ROOT, "funmap_out", r["region_id"], "pip.csv"), "04"),
        ("polyfun_sniff", lambda r: os.path.join(ROOT, "polyfun_out_sniff", f"{r['region_id']}.txt.gz"), "05/sniff"),
        ("polyfun_combined", lambda r: os.path.join(ROOT, "polyfun_out_combined", f"{r['region_id']}.txt.gz"), "05/combined"),
    ]:
        miss = [r["region_id"] for r in rows if not os.path.exists(path(r))]
        if miss:
            # this stage gates 06_collect, so an incomplete method is a
            # failure, not just a progress note - collecting now would
            # silently produce a table missing a third of its rows
            fail(f"{name}: {len(miss)} of {len(rows)} blocks missing from "
                 f"stage {label}, e.g. {miss[:2]}")
        else:
            print(f"  {name:14s} {len(rows):4d}/{len(rows)} done")


def check_collect():
    p = os.path.join(ROOT, "results", "pips.parquet")
    print("collected results:")
    print(f"  {'OK  ' if os.path.exists(p) else 'MISSING  '}{p}")
    if not os.path.exists(p):
        fail("results/pips.parquet")


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "data"
    print(f"[preflight: {stage}]  APP_ROOT={ROOT}")
    {"data": check_data, "prep": check_prep,
     "methods": check_methods, "collect": check_collect}[stage]()
    if not ok:
        print("\npreflight FAILED - fix the above before submitting")
        sys.exit(1)
    print("\npreflight passed")
