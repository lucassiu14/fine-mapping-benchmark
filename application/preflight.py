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


def fail(msg):
    global ok
    ok = False
    print(f"  MISSING  {msg}")


def blocks():
    with open(B) as f:
        return list(csv.DictReader(f))


def check_data():
    print("inputs for 02_prep:")
    rows = blocks()
    s = os.path.join(ROOT, "chol_sumstats.parquet")
    if not os.path.exists(s) or os.path.getsize(s) < 400_000_000:
        fail(f"chol_sumstats.parquet (expect ~526 MB, got "
             f"{os.path.getsize(s)/1e6:.0f} MB)" if os.path.exists(s)
             else "chol_sumstats.parquet")
    missing_ld = []
    for w in sorted({r["ld_region"] for r in rows}):
        for ext in ("npz", "gz"):
            p = os.path.join(ROOT, "ld", f"{w}.{ext}")
            # a truncated .npz is the failure mode that actually happens
            if not os.path.exists(p) or os.path.getsize(p) < 1000:
                missing_ld.append(f"{w}.{ext}")
    if missing_ld:
        fail(f"{len(missing_ld)} LD files, e.g. {missing_ld[:3]}")
    missing_an = [c for c in sorted({int(r["chrom"]) for r in rows})
                  if not os.path.exists(os.path.join(ROOT, "annot", f"borzoi.{c}.annot.parquet"))]
    if missing_an:
        fail(f"annotations for chromosomes {missing_an}")
    if ok:
        print(f"  OK  sumstats, {len({r['ld_region'] for r in rows})} LD windows, "
              f"{len({r['chrom'] for r in rows})} annotation files")


def check_prep():
    print("prepared blocks (inputs for 03/04/05):")
    rows = blocks()
    bad = []
    for r in rows:
        d = os.path.join(ROOT, "regions", r["region_id"])
        if not os.path.exists(os.path.join(d, "DONE")):
            bad.append(r["region_id"]); continue
        for f in ("locus.z", "locus_ld.npy", "locus.annot", "funmap.z.txt",
                  "funmap.anno.txt", "snps.csv"):
            if not os.path.exists(os.path.join(d, f)):
                bad.append(f"{r['region_id']}/{f}")
    if bad:
        fail(f"{len(bad)} of {len(rows)} blocks incomplete, e.g. {bad[:3]}")
    else:
        print(f"  OK  all {len(rows)} blocks prepared")
    # The PolyFun prior, computed by 00g from our own annotations. Needed by
    # 05 only. Note S-LDSC is genome-wide, so all 22 must be present even
    # though only 22 chromosomes carry selected blocks.
    pd_ = os.environ.get("PF_PRIOR_DIR", "pf_prior")
    miss = [c for c in range(1, 23)
            if not os.path.exists(os.path.join(ROOT, pd_, f"chol.{c}.snpvar_constrained.gz"))]
    if miss:
        fail(f"PolyFun prior missing for chromosomes {miss} - run 00d/00e/00f/00g "
             f"(needed by 05 only)")
    else:
        print("  OK  PolyFun prior computed for all 22 chromosomes")


def check_methods():
    print("method outputs (inputs for 06_collect):")
    rows = blocks()
    for name, path, label in [
        ("fb_xregion", lambda r: os.path.join(ROOT, "joint_out", r["region_id"], "pip.csv"), "03"),
        ("funmap", lambda r: os.path.join(ROOT, "funmap_out", r["region_id"], "pip.csv"), "04"),
        ("polyfun_susie", lambda r: os.path.join(ROOT, "polyfun_out", f"{r['region_id']}.txt.gz"), "05"),
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
