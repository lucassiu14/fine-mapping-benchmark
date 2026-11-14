#!/usr/bin/env python
"""Why is no trial mass-feasible?

  python scripts/tuning/inspect_mass.py --storage tuning/<study>.journal \
      --study <name> [--lo 0.5] [--hi 2.0]

Reads a finished or running Optuna study and answers one question: can the
post-fix model reach sum(PIP)/#causal ~ 1 anywhere inside the published
parameter space, and if so which knob gets it there?

Prints
  * the mass distribution over completed trials, against the feasible window;
  * rank correlation of mass with every parameter, so the lever is identified
    rather than assumed (v4 was launched on the argument that `gamma` is the
    remaining lever - this is the test of that claim);
  * the lowest-mass trials in full, which is the configuration to extend
    toward if the window is merely out of reach rather than unreachable;
  * mass in the top vs bottom third of `gamma`, holding nothing else fixed -
    crude, but enough to see whether raising gamma moves mass at all.

Read-only. Touches no study state.
"""
import argparse, importlib.util as iu, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load_storage(path):
    spec = iu.spec_from_file_location("ofb", os.path.join(HERE, "optuna_fb.py"))
    mod = iu.module_from_spec(spec)
    argv, sys.argv = sys.argv, ["optuna_fb"]
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.argv = argv
    return mod.make_storage(path)


def spearman(xs, ys):
    """Rank correlation, no scipy. Ties get averaged ranks."""
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    if len(xs) < 3:
        return float("nan")
    rx, ry = ranks(xs), ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return num / (dx * dy) if dx > 0 and dy > 0 else float("nan")


def quantile(sorted_vals, q):
    if not sorted_vals:
        return float("nan")
    i = q * (len(sorted_vals) - 1)
    lo, hi = int(i), min(int(i) + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (i - lo) * (sorted_vals[hi] - sorted_vals[lo])


def main():
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    ap = argparse.ArgumentParser()
    ap.add_argument("--storage", required=True)
    ap.add_argument("--study", required=True)
    ap.add_argument("--lo", type=float, default=0.5)
    ap.add_argument("--hi", type=float, default=2.0)
    ap.add_argument("--show", type=int, default=8)
    a = ap.parse_args()

    study = optuna.load_study(study_name=a.study, storage=load_storage(a.storage))
    done = [t for t in study.get_trials(deepcopy=False)
            if t.state == optuna.trial.TrialState.COMPLETE]
    rows = [(t, t.user_attrs.get("mean_total_mass_ratio")) for t in done]
    rows = [(t, m) for t, m in rows if isinstance(m, float) and m == m]
    if not rows:
        print("no completed trial carries mean_total_mass_ratio")
        return
    print(f"=== {a.study}: {len(rows)} completed trials with a mass ratio ===")

    ms = sorted(m for _, m in rows)
    print(f"  feasible window [{a.lo}, {a.hi}]")
    print("  mass ratio:  min %.2f  q10 %.2f  median %.2f  q90 %.2f  max %.2f"
          % (ms[0], quantile(ms, .10), quantile(ms, .50), quantile(ms, .90), ms[-1]))
    n_in = sum(1 for m in ms if a.lo <= m <= a.hi)
    n_lo = sum(1 for m in ms if m < a.lo)
    print(f"  below window: {n_lo}   inside: {n_in}   above: {len(ms) - n_in - n_lo}")
    if n_in == 0:
        gap = ms[0] - a.hi
        print(f"  CLOSEST trial is {ms[0]:.2f}, i.e. {gap:+.2f} past the ceiling"
              if ms[0] > a.hi else
              f"  CLOSEST trial is {ms[0]:.2f}, below the floor")

    names = sorted({k for t, _ in rows for k in t.params})
    print("\n  rank correlation of mass with each parameter"
          "  (negative = raising it REDUCES mass):")
    for nm in names:
        sub = [(t.params[nm], m) for t, m in rows if nm in t.params]
        if len(sub) >= 3:
            r = spearman([x for x, _ in sub], [y for _, y in sub])
            mark = "   <- lowers mass" if r == r and r <= -0.3 else ""
            print(f"    {nm:<22} {r:+.2f}  (n={len(sub)}){mark}")

    print(f"\n  the {a.show} LOWEST-mass trials:")
    for t, m in sorted(rows, key=lambda r: r[1])[:a.show]:
        ua = t.user_attrs
        print(f"    mass={m:.2f}  minAP={ua.get('min_cell_ap', float('nan')):.4f}"
              f"  meanAP={ua.get('mean_ap', float('nan')):.4f}"
              f"  worstFDR={ua.get('worst_cell_fdr_at_95', float('nan')):.3f}"
              f"  #{t.number}")
        print(f"      {t.params}")

    gs = [(t.params["gamma"], m) for t, m in rows if "gamma" in t.params]
    if len(gs) >= 9:
        gs.sort()
        k = len(gs) // 3
        lo_third = [m for _, m in gs[:k]]
        hi_third = [m for _, m in gs[-k:]]
        print(f"\n  gamma's effect on mass (crude, nothing held fixed):")
        print(f"    lowest  third of gamma ({gs[0][0]:.3f}-{gs[k-1][0]:.3f}): "
              f"mean mass {sum(lo_third)/len(lo_third):.2f}")
        print(f"    highest third of gamma ({gs[-k][0]:.3f}-{gs[-1][0]:.3f}): "
              f"mean mass {sum(hi_third)/len(hi_third):.2f}")
        print(f"    gamma range explored: {gs[0][0]:.4f} to {gs[-1][0]:.4f}")


if __name__ == "__main__":
    main()
