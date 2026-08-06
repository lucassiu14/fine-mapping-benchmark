#!/usr/bin/env python
"""
scripts/tuning/sweep_from_study.py

Slice a finished Optuna study by ONE parameter and report what each slice can
achieve on every objective.

Motivation. `sparse_concrete` (K_C) is not a modelling prior: it is the size of
the working set retained when evaluating the likelihood, so shrinking it
TRUNCATES the likelihood rather than merely making it cheaper. Study
fb_insample_v2 put every Pareto point at K_C = 263-269 (its ceiling is 300) and
the production benchmark runs at 50, so we need to know what K_C actually costs
and buys before choosing a production value.

We do not need a fresh sweep to answer that. The search sampled K_C
log-uniformly over [1, 300], so it has already evaluated the whole range many
times over with every other parameter varying around it. Binning the completed
trials recovers the sweep for free.

Read the output as an ACHIEVABILITY frontier, not an average: within a bin the
other hyperparameters differ, so "best ece_hi in this bin" answers "how well
calibrated can the model be at this K_C", which is the question that matters
when deciding where to fix it.

Usage:
  python scripts/tuning/sweep_from_study.py \
      --study fb_insample_v2 --storage tuning/fb_insample_v2.journal \
      --objective triple --param sparse_concrete
"""
import argparse
import sys

import optuna

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from optuna_fb import make_storage                      # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--study", required=True)
    ap.add_argument("--storage", required=True)
    ap.add_argument("--objective", choices=["multi", "scalar", "triple"],
                    default="triple")
    ap.add_argument("--param", default="sparse_concrete")
    ap.add_argument("--edges", default="1,5,15,40,60,100,160,240,300",
                    help="bin edges, comma separated")
    ap.add_argument("--min-ap", type=float, default=None,
                    help="also list the best-CALIBRATED trials (mass ratio closest "
                         "to 1.0) subject to mean_ap >= this floor. Mass ratio is "
                         "not one of the objectives, so well-calibrated trials can "
                         "be dominated on the recorded front while still being the "
                         "configuration you want.")
    ap.add_argument("--top", type=int, default=12,
                    help="how many to list for --min-ap")
    args = ap.parse_args()

    storage = make_storage(args.storage)
    kwargs = dict(study_name=args.study, storage=storage, load_if_exists=True)
    if args.objective == "triple":
        kwargs["directions"] = ["maximize", "minimize", "minimize"]
    elif args.objective == "multi":
        kwargs["directions"] = ["maximize", "minimize"]
    else:
        kwargs["direction"] = "maximize"
    study = optuna.create_study(**kwargs)

    complete = [t for t in study.get_trials(deepcopy=False)
                if t.state == optuna.trial.TrialState.COMPLETE
                and args.param in t.params]
    if not complete:
        print(f"no completed trials carry parameter '{args.param}'")
        return 1

    edges = [float(x) for x in args.edges.split(",")]

    def attr(t, k):
        v = t.user_attrs.get(k)
        return v if isinstance(v, float) and v == v else None

    print(f"\n=== study '{args.study}': {len(complete)} completed trials, "
          f"sliced by {args.param} ===")
    print("  'best' = the most favourable value ACHIEVED in the bin, i.e. what is")
    print("  attainable at that setting with the other parameters free.\n")
    hdr = (f"{args.param:>16} {'n':>5} {'bestAP':>8} {'medAP':>8} "
           f"{'best_ece':>9} {'med_ece':>9} {'mass@bestAP':>12} {'closest_mass':>13}")
    print(hdr)
    print("-" * len(hdr))

    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = [t for t in complete if lo <= float(t.params[args.param]) < hi]
        if not sel:
            continue
        aps = [a for a in (attr(t, "mean_ap") for t in sel) if a is not None]
        eces = [e for e in (attr(t, "mean_ece_hi") for t in sel) if e is not None]
        if not aps:
            continue
        best = max(sel, key=lambda t: attr(t, "mean_ap") or -1)
        mass_at_best = attr(best, "mean_total_mass_ratio")
        masses = [m for m in (attr(t, "mean_total_mass_ratio") for t in sel)
                  if m is not None]
        # Calibration target is a mass ratio of exactly 1.0; report the closest
        # value any configuration in this bin reached.
        closest = min(masses, key=lambda m: abs(m - 1.0)) if masses else float("nan")
        aps_s = sorted(aps); eces_s = sorted(eces) if eces else [float("nan")]
        print(f"{f'[{lo:g},{hi:g})':>16} {len(sel):>5} "
              f"{max(aps):>8.4f} {aps_s[len(aps_s)//2]:>8.4f} "
              f"{(min(eces) if eces else float('nan')):>9.4f} "
              f"{eces_s[len(eces_s)//2]:>9.4f} "
              f"{(mass_at_best if mass_at_best is not None else float('nan')):>12.2f} "
              f"{closest:>13.2f}")

    if args.min_ap is not None:
        elig = [t for t in complete
                if (attr(t, "mean_ap") or -1) >= args.min_ap
                and attr(t, "mean_total_mass_ratio") is not None]
        print(f"\n=== best-calibrated trials with AP >= {args.min_ap} "
              f"({len(elig)} eligible of {len(complete)}) ===")
        if not elig:
            print("  none - lower --min-ap")
        else:
            print(f"{'AP':>8} {'mass':>7} {'ece_hi':>8} {'reliab':>7} "
                  f"{'fdr@95':>8} {'trial':>7}  params")
            elig.sort(key=lambda t: abs(attr(t, "mean_total_mass_ratio") - 1.0))
            for t in elig[:args.top]:
                f95 = attr(t, "mean_fdr_at_95")
                print(f"{attr(t,'mean_ap'):>8.4f} "
                      f"{attr(t,'mean_total_mass_ratio'):>7.2f} "
                      f"{(attr(t,'mean_ece_hi') or float('nan')):>8.4f} "
                      f"{(attr(t,'mean_hi_pip_reliab') or float('nan')):>7.2f} "
                      f"{(f95 if f95 is not None else float('nan')):>8.4f} "
                      f"{t.number:>7}  {t.params}")

    print("\nRead the last two columns together with bestAP: if a low-K_C bin")
    print("reaches a mass ratio near 1.0 while a high-K_C bin cannot, then K_C is")
    print("controlling calibration and not only runtime, and the production value")
    print("must be chosen on that basis rather than on cost alone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
