#!/usr/bin/env python3
"""Apply the selection rule to the jointly re-scored candidates and name ONE config.

The rule (user's decision, 2026-11-19), over candidates scored on all 27 cells
by submit_fb_zfix_candidates.sh:

  maximise   min over the 9 SPARSE cells of AP
  subject to (a) mass feasible on BOTH strata, sum(PIP)/#causal in [0.5, 2.0]
             (b) no regression vs the DEPLOYED config on sparse_inf, in
                 stratum AP and in worst-sparse_inf-cell FDR@0.95

Why this rule rather than "best worst cell over all 27": sparse_inf's cells are
the minimum in every trial, so worst-cell-over-everything is a sparse_inf-only
selection, and no setting of the hyperparameters makes sparse_inf calibrated -
measured on its own front, cutting worst-cell FDR 0.62 -> 0.41 costs 58% of mean
AP. Selecting on sparse keeps the regime the model is actually specified for,
and (b) forbids buying sparse gains out of sparse_inf.

(a) comes free: optuna_fb records `mass_feasible` as the WORST stratum's window
violation (optuna_fb.py:500-503), so on a two-stratum canary it is already
"both strata inside the window".

TOLERANCE. Each cell is 10 regions at one seed, so differences of ~0.01-0.02 in
AP are Monte Carlo noise, not signal. Constraint (b) is therefore applied with a
tolerance, and the table prints each candidate's MARGIN so a near-miss is
visible rather than silently disqualifying. The verdict is reported at both
tol=0 and the chosen tolerance; if they disagree, that is said out loud and the
choice is the user's, not this script's.

Usage:
  python3 scripts/tuning/select_fb_config.py                    # default globs
  python3 scripts/tuning/select_fb_config.py --tol 0.0
  python3 scripts/tuning/select_fb_config.py --candidates 'tuning/fb_sel_*_canary.journal' \
      --baseline tuning/fb_zfix01_canary.journal
"""
import argparse
import glob
import os
import sys

import optuna

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from optuna_fb import make_storage  # noqa: E402

optuna.logging.set_verbosity(optuna.logging.WARNING)

PARAM_ORDER = ["sigma_sq", "n_caus", "sparse_concrete", "prior_regularisation",
               "lambda_l1", "hierarchy_M", "prior_mass", "gamma", "gamma_key"]


def load(path):
    """The one COMPLETE trial in a pinned-config journal, as a flat dict.

    A canary journal holds a single trial, but a resubmitted or crashed job can
    leave extra RUNNING/FAIL entries, so COMPLETE is filtered for explicitly
    and the LAST one wins (a re-run supersedes).
    """
    storage = make_storage(path)
    out = None
    for s in optuna.get_all_study_summaries(storage=storage):
        st = optuna.load_study(study_name=s.study_name, storage=storage)
        done = [t for t in st.trials if t.state == optuna.trial.TrialState.COMPLETE]
        if not done:
            continue
        t = done[-1]
        a = dict(t.user_attrs)
        # --fix holds parameters OUT of the search space, so t.params is empty
        # for a pinned trial and the values live under 'fixed_params'.
        a["_params"] = dict(a.get("fixed_params") or {}) or dict(t.params)
        a["_n_complete"] = len(done)
        out = a
    return out


def cells_of(a, stratum):
    """{S: (ap, fdr)} for one stratum's cells. Keys look like 'sparse|S=3'."""
    got = {}
    for k in a:
        if k.startswith("ap_") and "|S=" in k:
            lab, _, s = k[3:].partition("|S=")
            if lab != stratum:
                continue
            got[int(s)] = (a[k], a.get(f"fdr_at_95_{lab}|S={s}", float("nan")))
    return dict(sorted(got.items()))


def summarise(a):
    sp, si = cells_of(a, "sparse"), cells_of(a, "sparse_inf")
    f = lambda xs: [v for v in xs if v == v]                       # noqa: E731
    sp_ap = f([v[0] for v in sp.values()])
    si_ap = f([v[0] for v in si.values()])
    si_fd = f([v[1] for v in si.values()])
    sp_fd = f([v[1] for v in sp.values()])
    return {
        "worst_sparse_cell_ap": min(sp_ap) if sp_ap else float("nan"),
        "worst_sparse_cell_fdr": max(sp_fd) if sp_fd else float("nan"),
        "worst_sinf_cell_ap": min(si_ap) if si_ap else float("nan"),
        "worst_sinf_cell_fdr": max(si_fd) if si_fd else float("nan"),
        "ap_sparse": a.get("ap_sparse", float("nan")),
        "ap_sparse_inf": a.get("ap_sparse_inf", float("nan")),
        "fdr_sparse": a.get("fdr_at_95_sparse", float("nan")),
        "fdr_sparse_inf": a.get("fdr_at_95_sparse_inf", float("nan")),
        "mass_sparse": a.get("total_mass_ratio_sparse", float("nan")),
        "mass_sparse_inf": a.get("total_mass_ratio_sparse_inf", float("nan")),
        "ece_sparse": a.get("ece_hi_sparse", float("nan")),
        "ece_sparse_inf": a.get("ece_hi_sparse_inf", float("nan")),
        "feasible": a.get("mass_feasible"),
        "cells_sparse": sp,
        "cells_sinf": si,
    }


def verdict(c, b, tol):
    """Disqualification reasons under constraints (a) and (b). Empty == passes."""
    why = []
    if c["feasible"] is not True:
        why.append(f"mass infeasible (sparse {c['mass_sparse']:.2f}, "
                   f"sparse_inf {c['mass_sparse_inf']:.2f}; window [0.5, 2.0])")
    if c["ap_sparse_inf"] < b["ap_sparse_inf"] - tol:
        why.append(f"sparse_inf AP regressed {b['ap_sparse_inf']:.3f} -> "
                   f"{c['ap_sparse_inf']:.3f}")
    if c["worst_sinf_cell_fdr"] > b["worst_sinf_cell_fdr"] + tol:
        why.append(f"sparse_inf worst-cell FDR regressed "
                   f"{b['worst_sinf_cell_fdr']:.3f} -> {c['worst_sinf_cell_fdr']:.3f}")
    return why


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", default="tuning/fb_sel_*_canary.journal")
    ap.add_argument("--baseline", default="tuning/fb_zfix01_canary.journal",
                    help="the DEPLOYED config scored on the same 27 cells")
    ap.add_argument("--tol", type=float, default=0.02,
                    help="Monte Carlo slack on constraint (b); 10 regions/cell, "
                         "one seed. Default 0.02.")
    args = ap.parse_args()

    bl = load(args.baseline)
    if bl is None:
        sys.exit(f"no COMPLETE trial in the baseline {args.baseline}")
    b = summarise(bl)
    print(f"BASELINE (deployed config) {args.baseline}")
    print(f"  sparse      AP {b['ap_sparse']:.3f}  FDR@.95 {b['fdr_sparse']:.3f}  "
          f"ece {b['ece_sparse']:.3f}  mass {b['mass_sparse']:.3f}  "
          f"worst cell AP {b['worst_sparse_cell_ap']:.3f}")
    print(f"  sparse_inf  AP {b['ap_sparse_inf']:.3f}  FDR@.95 {b['fdr_sparse_inf']:.3f}  "
          f"ece {b['ece_sparse_inf']:.3f}  mass {b['mass_sparse_inf']:.3f}  "
          f"worst cell FDR {b['worst_sinf_cell_fdr']:.3f}")
    if b["feasible"] is not True:
        print("  NOTE: the baseline itself is mass-INFEASIBLE, so constraint (a)")
        print("        alone guarantees the selected config improves on it.")

    rows = []
    for path in sorted(glob.glob(args.candidates)):
        tag = os.path.basename(path).replace("fb_sel_", "").replace("_canary.journal", "")
        a = load(path)
        if a is None:
            print(f"\n[{tag}] no COMPLETE trial yet - still running or failed")
            continue
        rows.append((tag, summarise(a), a["_params"]))

    if not rows:
        sys.exit("\nno candidate has completed yet.")

    print(f"\n{'cand':<6} {'wSPcellAP':>9} {'spAP':>6} {'spFDR':>6} {'spMass':>7} "
          f"{'siAP':>6} {'wSIcFDR':>8} {'siMass':>7} {'sigma_sq':>9}  verdict")
    survivors = {}
    for tol in (0.0, args.tol):
        survivors[tol] = []
    for tag, c, p in sorted(rows, key=lambda r: -r[1]["worst_sparse_cell_ap"]):
        why = verdict(c, b, args.tol)
        for tol in survivors:
            if not verdict(c, b, tol):
                survivors[tol].append((tag, c, p))
        print(f"{tag:<6} {c['worst_sparse_cell_ap']:>9.3f} {c['ap_sparse']:>6.3f} "
              f"{c['fdr_sparse']:>6.3f} {c['mass_sparse']:>7.3f} "
              f"{c['ap_sparse_inf']:>6.3f} {c['worst_sinf_cell_fdr']:>8.3f} "
              f"{c['mass_sparse_inf']:>7.3f} {p.get('sigma_sq', float('nan')):>9.4f}  "
              + ("PASS" if not why else "; ".join(why)))

    best = {tol: (max(v, key=lambda r: r[1]["worst_sparse_cell_ap"]) if v else None)
            for tol, v in survivors.items()}
    print(f"\n{len(survivors[args.tol])}/{len(rows)} candidates pass at "
          f"tol={args.tol}; {len(survivors[0.0])}/{len(rows)} at tol=0")
    if best[0.0] and best[args.tol] and best[0.0][0] != best[args.tol][0]:
        print(f"  !! THE TOLERANCE CHANGES THE WINNER: {best[0.0][0]} at tol=0 vs "
              f"{best[args.tol][0]} at tol={args.tol}. That is a judgement call,")
        print("     not a result - take it to the user before deploying either.")
    if not best[args.tol]:
        print("\nNO CANDIDATE PASSES. That is a finding: the sparse front cannot be")
        print("reached without giving up sparse_inf. Either deploy the best")
        print("feasible point and state the sparse_inf regression, or relax (b).")
        return

    tag, c, p = best[args.tol]
    print(f"\n=== SELECTED: {tag} ===")
    print(f"  worst sparse cell AP {c['worst_sparse_cell_ap']:.3f} "
          f"(baseline {b['worst_sparse_cell_ap']:.3f})")
    print(f"  sparse      AP {c['ap_sparse']:.3f}  FDR@.95 {c['fdr_sparse']:.3f}  "
          f"ece {c['ece_sparse']:.3f}  mass {c['mass_sparse']:.3f}")
    print(f"  sparse_inf  AP {c['ap_sparse_inf']:.3f}  FDR@.95 {c['fdr_sparse_inf']:.3f}  "
          f"ece {c['ece_sparse_inf']:.3f}  mass {c['mass_sparse_inf']:.3f}")
    print("  per-cell sparse     " + "  ".join(
        f"S={s}:AP={v[0]:.3f}/FDR={v[1]:.3f}" for s, v in c["cells_sparse"].items()))
    print("  per-cell sparse_inf " + "  ".join(
        f"S={s}:AP={v[0]:.3f}/FDR={v[1]:.3f}" for s, v in c["cells_sinf"].items()))
    print("\n  deploy into scripts/hpc/run_benchmark_job.R, METHOD_ARGS$functional_beatrice")
    print("  (prior_mass goes on .FBj, which every joint arm inherits):")
    for k in PARAM_ORDER:
        if k in p:
            v = p[k]
            print(f"    {k:<22} = {int(v) if k in ('n_caus', 'sparse_concrete') else v}")
    print("\n  gamma and gamma_key are NOT currently in run_benchmark_job.R - every")
    print("  run so far used the python defaults (0.1, 0.2). If the selected values")
    print("  differ, they must be ADDED or the deployed model is not the tuned one.")


if __name__ == "__main__":
    main()
