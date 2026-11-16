#!/usr/bin/env python3
"""Is FB's FDR a function of its posterior mass? Scan EVERY trial, not the front.

Why this exists. Every read we have done of the v5 studies went through
`study.best_trials`, which under constrained NSGA-II holds only FEASIBLE trials
and orders them by the objectives - AP first. A trial that drove posterior mass
near 1 and bought a good FDR at the cost of some AP is dominated on AP, so it
never appears on the front and has been invisible in every report so far.

Each trial already recorded, per stratum (optuna_fb.py:382-404):
    ap_<lab>, fdr_at_95_<lab>, ece_hi_<lab>, nsel_at_95_<lab>,
    total_mass_ratio_<lab>
so the question "does mass -> 1 buy an acceptable FDR" is answerable from the
journals with no new compute.

Under the independent-Bernoulli configuration prior, sum_j PIP_j = E[|gamma| | z],
so total_mass_ratio = sum(PIP)/#causal is the posterior's own claim about how
many causal variants it found, divided by how many there were. A calibrated
posterior sits at 1. If FDR tracks that ratio, the defect is the prior's scale
and the fix is to make it adaptive. If FDR stays bad even at ratio 1, the
defect is elsewhere (over-confident variational posterior, or the missing
infinitesimal component) and no amount of mass control will reach it.

Usage:
  python3 scripts/tuning/fdr_vs_mass.py tuning/*.journal
  python3 scripts/tuning/fdr_vs_mass.py --min-nsel 50 tuning/fb_retune_sparse_v5.journal
"""
import argparse
import os
import sys

import optuna

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from optuna_fb import make_storage  # noqa: E402

optuna.logging.set_verbosity(optuna.logging.WARNING)

# Params worth printing beside a good-FDR configuration.
SHOW = ("prior_mass", "sigma_sq", "n_caus", "sparse_concrete",
        "prior_regularisation", "gamma", "lambda_l1", "hierarchy_M")


def strata_of(attrs):
    """Stratum labels, excluding cell keys ('sparse|S=3') and running means."""
    out = []
    for k in attrs:
        if k.startswith("total_mass_ratio_"):
            lab = k[len("total_mass_ratio_"):]
            if "|" in lab or "=" in lab or lab.startswith("after_"):
                continue
            out.append(lab)
    return sorted(out)


def spearman(xs, ys):
    """Rank correlation, implemented here to avoid a scipy dependency."""
    n = len(xs)
    if n < 3:
        return float("nan")

    def ranks(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:                      # average ranks within ties
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return num / (dx * dy) if dx and dy else float("nan")


def collect(paths):
    rows = []
    for path in paths:
        if not os.path.exists(path):
            print(f"  (missing: {path})")
            continue
        storage = make_storage(path)
        for s in optuna.get_all_study_summaries(storage=storage):
            st = optuna.load_study(study_name=s.study_name, storage=storage)
            for t in st.trials:
                if t.state != optuna.trial.TrialState.COMPLETE:
                    continue
                a = t.user_attrs
                params = dict(a.get("fixed_params") or {})
                params.update(t.params)
                for lab in strata_of(a):
                    mass = a.get(f"total_mass_ratio_{lab}")
                    fdr = a.get(f"fdr_at_95_{lab}")
                    if mass is None or fdr is None:
                        continue
                    if mass != mass or fdr != fdr:        # NaN
                        continue
                    rows.append({
                        "study": s.study_name, "trial": t.number, "stratum": lab,
                        "mass": float(mass), "fdr": float(fdr),
                        "ap": a.get(f"ap_{lab}", float("nan")),
                        "ece": a.get(f"ece_hi_{lab}", float("nan")),
                        "nsel": a.get(f"nsel_at_95_{lab}", 0) or 0,
                        "params": params,
                    })
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("journals", nargs="+")
    p.add_argument("--min-nsel", type=int, default=20,
                   help="drop trial-strata selecting fewer than this many "
                        "variants at PIP>=0.95; an FDR over 3 selections is noise")
    args = p.parse_args()

    rows = collect(args.journals)
    if not rows:
        print("no completed trials with both mass and FDR recorded")
        return

    kept = [r for r in rows if r["nsel"] >= args.min_nsel]
    print(f"\n{len(rows)} trial-strata recorded; {len(kept)} with "
          f"nsel_at_95 >= {args.min_nsel}\n")

    for lab in sorted({r["stratum"] for r in kept}):
        d = [r for r in kept if r["stratum"] == lab]
        if not d:
            continue
        print("=" * 72)
        print(f"STRATUM {lab}   ({len(d)} trials)")
        print("=" * 72)

        rho = spearman([r["mass"] for r in d], [r["fdr"] for r in d])
        print(f"  Spearman(mass ratio, FDR@0.95) = {rho:+.3f}")

        # FDR against the mass ratio, binned. A calibrated posterior sits at 1.
        edges = [(0, 0.8), (0.8, 1.25), (1.25, 2), (2, 3), (3, 5), (5, 1e9)]
        print(f"\n  {'mass ratio':>14} {'n':>5} {'median FDR':>11} "
              f"{'min FDR':>9} {'median AP':>10}")
        for lo, hi in edges:
            b = [r for r in d if lo <= r["mass"] < hi]
            if not b:
                continue
            f = sorted(r["fdr"] for r in b)
            a = sorted(r["ap"] for r in b)
            med = f[len(f) // 2]
            meda = a[len(a) // 2]
            name = f"{lo:g}-{hi:g}" if hi < 1e9 else f">{lo:g}"
            print(f"  {name:>14} {len(b):>5} {med:>11.3f} {f[0]:>9.3f} "
                  f"{meda:>10.3f}")

        # The direct question: did ANY configuration reach an acceptable FDR?
        print()
        for thr in (0.05, 0.10, 0.20, 0.30):
            ok = [r for r in d if r["fdr"] <= thr]
            if ok:
                best = max(ok, key=lambda r: r["ap"])
                print(f"  FDR <= {thr:.2f}: {len(ok):>3} trials; "
                      f"best AP among them {best['ap']:.3f} "
                      f"(mass {best['mass']:.2f}, nsel {best['nsel']})")
            else:
                print(f"  FDR <= {thr:.2f}: none")

        # What the lowest-FDR configurations actually look like.
        print("\n  Ten lowest-FDR trials:")
        print(f"  {'trial':>6} {'FDR':>7} {'AP':>6} {'mass':>6} {'ece':>6} "
              f"{'nsel':>7}  params")
        for r in sorted(d, key=lambda r: r["fdr"])[:10]:
            ps = " ".join(f"{k}={r['params'][k]:.4g}"
                          if isinstance(r["params"].get(k), float)
                          else f"{k}={r['params'][k]}"
                          for k in SHOW if k in r["params"])
            print(f"  {r['trial']:>6} {r['fdr']:>7.3f} {r['ap']:>6.3f} "
                  f"{r['mass']:>6.2f} {r['ece']:>6.3f} {r['nsel']:>7} {ps}")
        print()


if __name__ == "__main__":
    main()
