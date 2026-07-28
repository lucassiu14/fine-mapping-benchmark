#!/usr/bin/env python
"""
scripts/tuning/optuna_fb.py

Optuna hyperparameter search for Functional BEATRICE, designed for a PBS cluster:

  RESUMABLE   The study lives in an on-disk Optuna storage. Every worker opens it
              with load_if_exists=True, so re-submitting the job simply continues
              the SAME study - completed trials are remembered, the sampler keeps
              its history, nothing is recomputed. Killed-mid-trial jobs leave a
              stale RUNNING trial, which TPE ignores (only COMPLETE trials inform
              the sampler), so no cleanup is required.

  WALLTIME    --timeout seconds is passed to study.optimize(); each worker stops
              starting new trials once it is reached and exits cleanly, so the
              job ends under its walltime instead of being killed mid-write.
              Set it to (walltime - one typical trial) - the submit script does.

  PARALLEL    Many workers share one study through the storage file. We use
              JournalStorage, which is the backend Optuna recommends for
              networked filesystems (NFS/GPFS/RDS) - plain SQLite is prone to
              lock errors there.

  EFFICIENT   Each trial is scored on a small FIXED tuning simulation, one
              scenario at a time. In single-objective mode the running mean is
              reported to Optuna after each scenario and a MedianPruner abandons
              configurations that are already clearly worse than the median -
              usually after the first scenario, cutting most of the cost of a
              bad trial.

The heavy lifting stays in R (scripts/tuning/fb_objective.R) so the objective is
computed by the package's own evaluate_methods() with the audited metric
formulas, rather than a Python re-implementation that could drift.

Example (one worker, 8h budget):
  python scripts/tuning/optuna_fb.py \
      --sim tuning/sim_binary_ref500.rds --study fb_binary_ref500 \
      --storage tuning/fb_binary_ref500.journal --timeout 28000
"""
import argparse
import os
import subprocess
import sys
import time

import optuna


# --------------------------------------------------------------------------
# storage: JournalStorage is the NFS-safe backend. Its class moved between
# Optuna versions, so try the current location first and fall back.
# --------------------------------------------------------------------------
def make_storage(path):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    from optuna.storages import JournalStorage
    backend = None
    try:                                     # Optuna >= 4.0
        from optuna.storages.journal import JournalFileBackend
        backend = JournalFileBackend(path)
    except Exception:
        pass
    if backend is None:
        try:                                 # Optuna 3.4 - 3.x
            from optuna.storages import JournalFileBackend
            backend = JournalFileBackend(path)
        except Exception:
            pass
    if backend is None:                      # older Optuna
        from optuna.storages import JournalFileStorage
        backend = JournalFileStorage(path)
    return JournalStorage(backend)


def parse_kv(text):
    """Parse the `key=value` lines printed by fb_objective.R."""
    out = {}
    for line in text.splitlines():
        line = line.strip()
        if "=" in line and not line.startswith("#"):
            k, _, v = line.partition("=")
            out[k.strip()] = v.strip()
    return out


def run_one_scenario(args, params, scenario):
    """Invoke the R objective for a single scenario; return a dict of metrics."""
    cmd = [args.rscript, args.objective_script,
           "--sim", args.sim,
           "--scenario", str(scenario),
           "--regions", str(args.regions),
           "--beatrice_dir", args.beatrice_dir,
           "--python", args.python,
           "--max_iter", str(args.max_iter)]
    for k, v in params.items():
        cmd += ["--" + k, str(v)]

    proc = subprocess.run(cmd, capture_output=True, text=True,
                          timeout=args.trial_timeout or None)
    kv = parse_kv(proc.stdout)
    if kv.get("status") != "ok":
        msg = kv.get("message", (proc.stderr or "").strip()[-400:])
        raise RuntimeError(f"objective failed (scenario {scenario}): {msg}")

    def num(key):
        v = kv.get(key, "NA")
        try:
            return float(v)
        except ValueError:
            return float("nan")

    return {k: num(k) for k in
            ("ap", "max_fdr_violation_n20", "total_mass_ratio", "hi_pip_reliab", "seconds")}


def suggest_params(trial, args):
    """The search space: every knob beatrice_annot.py actually consumes.

    max_iter is deliberately NOT tuned - it trades compute for convergence, so
    including it would let the search buy objective value with runtime and make
    trials incomparable. It is fixed at the production value.
    """
    p = {
        "lambda_l1":            trial.suggest_float("lambda_l1", 1e-4, 1.0, log=True),
        "prior_regularisation": trial.suggest_float("prior_regularisation", 0.1, 50.0, log=True),
        "sigma_sq":             trial.suggest_float("sigma_sq", 5e-3, 0.5, log=True),
        "hierarchy_M":          trial.suggest_float("hierarchy_M", 1.0, 100.0, log=True),
        "n_caus":               trial.suggest_int("n_caus", 1, 10),
        "sparse_concrete":      trial.suggest_int("sparse_concrete", 10, 300, step=10),
    }
    if args.tune_gamma_coverage:
        p["gamma_coverage"] = trial.suggest_float("gamma_coverage", 0.8, 0.99)
    return p


def build_objective(args):
    n_scen = args.scenarios

    def objective(trial):
        params = suggest_params(trial, args)
        aps, viols, masses, reliabs = [], [], [], []

        for s in range(1, n_scen + 1):
            m = run_one_scenario(args, params, s)
            aps.append(m["ap"])
            if m["max_fdr_violation_n20"] == m["max_fdr_violation_n20"]:   # not NaN
                viols.append(m["max_fdr_violation_n20"])
            if m["total_mass_ratio"] == m["total_mass_ratio"]:
                masses.append(m["total_mass_ratio"])
            if m["hi_pip_reliab"] == m["hi_pip_reliab"]:
                reliabs.append(m["hi_pip_reliab"])

            # Pruning only applies to single-objective studies; Optuna does not
            # support pruners for multi-objective optimisation.
            if args.objective == "scalar":
                running = sum(aps) / len(aps)
                if viols:
                    running -= args.penalty * (sum(viols) / len(viols))
                trial.report(running, s)
                if trial.should_prune():
                    raise optuna.TrialPruned()

        mean_ap    = sum(aps) / len(aps) if aps else float("nan")
        mean_viol  = sum(viols) / len(viols) if viols else 1.0
        mean_mass  = sum(masses) / len(masses) if masses else float("nan")
        mean_relb  = sum(reliabs) / len(reliabs) if reliabs else float("nan")

        # Record everything, so any objective can be re-derived from the study
        # later without re-running trials.
        trial.set_user_attr("mean_ap", mean_ap)
        trial.set_user_attr("mean_fdr_violation_n20", mean_viol)
        trial.set_user_attr("mean_total_mass_ratio", mean_mass)
        trial.set_user_attr("mean_hi_pip_reliab", mean_relb)

        if args.objective == "multi":
            return mean_ap, mean_viol            # maximise AP, minimise violation
        return mean_ap - args.penalty * mean_viol

    return objective


def main():
    ap = argparse.ArgumentParser(description="Optuna tuning for Functional BEATRICE")
    ap.add_argument("--sim", required=True, help="fixed tuning sim .rds (make_tuning_sim.R)")
    ap.add_argument("--study", required=True, help="study name (one per stratum)")
    ap.add_argument("--storage", required=True, help="journal file path (persists the study)")
    ap.add_argument("--objective", choices=["multi", "scalar"], default="multi",
                    help="multi: Pareto front of (AP, FDR violation), no arbitrary "
                         "weighting but no pruning. scalar: AP - penalty*violation, "
                         "enables pruning. Default: multi.")
    ap.add_argument("--penalty", type=float, default=0.5,
                    help="scalar mode only: weight on the FDR violation")
    ap.add_argument("--scenarios", type=int, default=4, help="scenarios per trial")
    ap.add_argument("--regions", type=int, default=4, help="regions per scenario (0 = all)")
    ap.add_argument("--max_iter", type=int, default=2000, help="FB max_iter (fixed, not tuned)")
    ap.add_argument("--tune-gamma-coverage", action="store_true",
                    help="also tune the credible-set coverage threshold")
    ap.add_argument("--n-trials", type=int, default=None, help="stop after this many trials")
    ap.add_argument("--timeout", type=float, default=None,
                    help="seconds; stop launching trials after this (walltime guard)")
    ap.add_argument("--trial-timeout", type=float, default=None,
                    help="seconds; hard cap on a single R invocation")
    ap.add_argument("--rscript", default="Rscript")
    ap.add_argument("--objective-script", default="scripts/tuning/fb_objective.R")
    ap.add_argument("--beatrice_dir", default="BEATRICE_annot_sparse")
    ap.add_argument("--python", default=os.path.expanduser("~/tools/py-venv-runner.sh"))
    ap.add_argument("--startup-trials", type=int, default=20,
                    help="random trials before TPE starts modelling")
    ap.add_argument("--report", action="store_true",
                    help="print the study's best trials and exit (no optimisation)")
    args = ap.parse_args()

    storage = make_storage(args.storage)
    sampler = (optuna.samplers.NSGAIISampler(seed=None) if args.objective == "multi"
               else optuna.samplers.TPESampler(n_startup_trials=args.startup_trials))
    pruner = (optuna.pruners.MedianPruner(n_startup_trials=args.startup_trials,
                                          n_warmup_steps=1)
              if args.objective == "scalar" else optuna.pruners.NopPruner())

    kwargs = dict(study_name=args.study, storage=storage,
                  load_if_exists=True, sampler=sampler)
    if args.objective == "multi":
        kwargs["directions"] = ["maximize", "minimize"]
    else:
        kwargs["direction"] = "maximize"
        kwargs["pruner"] = pruner
    study = optuna.create_study(**kwargs)

    if args.report:
        report(study, args)
        return

    done = len([t for t in study.trials
                if t.state == optuna.trial.TrialState.COMPLETE])
    print(f"[optuna] study '{args.study}' opened: {done} completed trial(s) so far",
          flush=True)
    print(f"[optuna] objective={args.objective}  scenarios/trial={args.scenarios}  "
          f"regions/scenario={args.regions}  timeout={args.timeout}s", flush=True)

    t0 = time.time()
    study.optimize(build_objective(args),
                   n_trials=args.n_trials,
                   timeout=args.timeout,
                   catch=(RuntimeError, subprocess.TimeoutExpired),
                   gc_after_trial=True)

    done_after = len([t for t in study.trials
                      if t.state == optuna.trial.TrialState.COMPLETE])
    print(f"[optuna] worker finished after {time.time()-t0:.0f}s; "
          f"study now has {done_after} completed trial(s)", flush=True)
    report(study, args)


def report(study, args):
    complete = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
    print(f"\n=== study '{study.study_name}': {len(complete)} completed trials ===")
    if not complete:
        return
    if args.objective == "multi":
        print("Pareto front (AP up, FDR violation down):")
        for t in sorted(study.best_trials, key=lambda x: -x.values[0]):
            print(f"  AP={t.values[0]:.4f}  viol={t.values[1]:.4f}  "
                  f"mass={t.user_attrs.get('mean_total_mass_ratio', float('nan')):.2f}  "
                  f"reliab={t.user_attrs.get('mean_hi_pip_reliab', float('nan')):.2f}  "
                  f"#{t.number}  {t.params}")
    else:
        b = study.best_trial
        print(f"best score={b.value:.4f} (trial #{b.number})")
        print(f"  AP={b.user_attrs.get('mean_ap')}  "
              f"viol={b.user_attrs.get('mean_fdr_violation_n20')}  "
              f"mass={b.user_attrs.get('mean_total_mass_ratio')}  "
              f"reliab={b.user_attrs.get('mean_hi_pip_reliab')}")
        print(f"  params={b.params}")


if __name__ == "__main__":
    sys.exit(main())
