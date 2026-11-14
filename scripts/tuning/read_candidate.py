#!/usr/bin/env python3
"""Print a pinned-configuration (CANARY) trial's recorded metrics.

Why this exists: `optuna_fb.py --report` prints its headline block from
`study.best_trials`, which under constrained NSGA-II contains only FEASIBLE
trials. A candidate scored across both strata will normally violate the
posterior-mass constraint, because `sparse_inf` mass runs well above any
ceiling that suits `sparse` - the denominator counts only the sparse causals.
The trial still COMPLETEs and still records everything as user attributes; the
Pareto front is simply empty, which is the same symptom v4 showed
(optuna_fb.py:392, "0/76 trials feasible").

One configuration has to serve every setup, so the number that matters is the
WORST stratum, not the mean. That is what this prints first.

Usage:
  python3 scripts/tuning/read_candidate.py tuning/fb_cand01.journal [...]
"""
import sys

import optuna
from optuna.storages import JournalStorage
from optuna.storages.journal import JournalFileBackend

optuna.logging.set_verbosity(optuna.logging.WARNING)


def strata_of(attrs, prefix):
    """Stratum labels for which `<prefix>_<label>` was recorded."""
    out = []
    for k in attrs:
        if k.startswith(prefix + "_"):
            lab = k[len(prefix) + 1:]
            # cell keys look like 'sparse|S=3|phi=0.2'; stratum keys do not.
            if "|" not in lab and "=" not in lab:
                out.append(lab)
    return sorted(out)


def show(path):
    storage = JournalStorage(JournalFileBackend(path))
    names = [s.study_name for s in optuna.get_all_study_summaries(storage=storage)]
    for name in names:
        st = optuna.load_study(study_name=name, storage=storage)
        done = [t for t in st.trials if t.state == optuna.trial.TrialState.COMPLETE]
        print(f"\n=== {name}: {len(done)}/{len(st.trials)} complete ===")
        for t in done:
            a = t.user_attrs
            pm = t.params.get("prior_mass", a.get("prior_mass", "n/a"))
            print(f"\n[trial {t.number}]  prior_mass = {pm}")
            feas = a.get("mass_feasible")
            print(f"  mass_feasible = {feas}"
                  + ("   (so it is ABSENT from --report's Pareto front)"
                     if feas is False else ""))

            labs = strata_of(a, "ap")
            if labs:
                print(f"  {'stratum':<12} {'AP':>7} {'FDR@.95':>9} "
                      f"{'ece_hi':>7} {'mass':>7} {'n@.95':>6}")
                for lab in labs:
                    print(f"  {lab:<12} {a.get('ap_'+lab, float('nan')):>7.3f} "
                          f"{a.get('fdr_at_95_'+lab, float('nan')):>9.3f} "
                          f"{a.get('ece_hi_'+lab, float('nan')):>7.3f} "
                          f"{a.get('total_mass_ratio_'+lab, float('nan')):>7.3f} "
                          f"{a.get('nsel_at_95_'+lab, float('nan')):>6}")

            print(f"  WORST  FDR@.95 {a.get('worst_fdr_at_95', float('nan')):.3f}"
                  f"  (stratum {a.get('worst_stratum', '?')})")
            print(f"  WORST CELL  AP {a.get('min_cell_ap', float('nan')):.3f}"
                  f"  FDR@.95 {a.get('worst_cell_fdr_at_95', float('nan')):.3f}"
                  f"  ({a.get('worst_cell', '?')})")
            print(f"  mean   AP {a.get('mean_ap', float('nan')):.3f}"
                  f"  mass {a.get('mean_total_mass_ratio', float('nan')):.3f}"
                  f"  reliab {a.get('mean_hi_pip_reliab', float('nan')):.3f}"
                  f"   <- do NOT compare candidates on these: averaging across"
                  f" strata is what broke v4")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for p in sys.argv[1:]:
        show(p)
