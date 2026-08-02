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
    if args.rotate_regions:
        cmd += ["--rotate_regions", "TRUE"]
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
            ("ap", "max_fdr_violation_n20", "fdr_at_90", "fdr_at_95",
             "total_mass_ratio", "hi_pip_reliab", "ece_hi", "seconds")}


def suggest_params(trial, args):
    """The search space: every knob beatrice_annot.py actually consumes.

    max_iter is deliberately NOT tuned - it trades compute for convergence, so
    including it would let the search buy objective value with runtime and make
    trials incomparable. It is fixed at the production value.
    """
    # BOUNDS REVISION (2026-07-30), after study fb_binary_ref500 v1 (1586 trials).
    # The v1 optimum sat ON two of these bounds, which means the bound - not the
    # objective - was choosing the value:
    #   sigma_sq        0.005323 against a 5e-3 floor  (1.4% into a 2-decade range)
    #   sparse_concrete 10       against a 10 floor    (exactly on it)
    # Both floors are now opened up. sparse_concrete goes to the package's own
    # declared lower_bound=1 and switches to a log-uniform integer, which puts the
    # sampling density where the optimum actually lives instead of spreading it
    # evenly over 10..300. (BEATRICE's README examples all use 10, so v1 landed on
    # the authors' value - but standing on a floor is not evidence it is optimal.)
    # Kept as-is: lambda_l1, prior_regularisation and n_caus all found interior
    # optima in v1, so their ranges are doing their job.
    fixed = args.fix or {}

    def sfloat(name, lo, hi, log=True):
        return fixed[name] if name in fixed else trial.suggest_float(name, lo, hi, log=log)

    def sint(name, lo, hi, log=False):
        return fixed[name] if name in fixed else trial.suggest_int(name, lo, hi, log=log)

    p = {
        "lambda_l1":            sfloat("lambda_l1", 1e-4, 1.0),
        "prior_regularisation": sfloat("prior_regularisation", 0.1, 50.0),
        "sigma_sq":             sfloat("sigma_sq", 1e-4, 0.5),
        "hierarchy_M":          sfloat("hierarchy_M", 1.0, 100.0),
        "n_caus":               sint("n_caus", 1, 10),
        "sparse_concrete":      sint("sparse_concrete", 1, 300, log=True),
    }
    if args.tune_gamma_coverage:
        p["gamma_coverage"] = sfloat("gamma_coverage", 0.8, 0.99, log=False)
    return p



def _abort_thresholds(study, n_scen, args):
    """Per-stage AP cut-offs, from trials that have already completed.

    Optuna supports pruners only for single-objective studies, and Trial.report()
    raises for multi-objective ones - but with a 20-cell grid a trial is 60
    BEATRICE fits, so running every hopeless configuration to completion is the
    dominant cost. We keep our own per-stage record (ap_after_<k> user attrs) and
    abort a trial whose running AP is below the `--abort-quantile` of what
    completed trials had at the SAME stage.

    Deliberately conservative: no aborting until `--abort-min-trials` trials have
    finished, and never before `--abort-after` scenarios, so a configuration that
    starts slow on the easy cells still gets a hearing.
    """
    if args.abort_quantile <= 0:
        return {}
    done = [t for t in study.get_trials(deepcopy=False)
            if t.state == optuna.trial.TrialState.COMPLETE]
    if len(done) < args.abort_min_trials:
        return {}
    cuts = {}
    for k in range(args.abort_after, n_scen):
        vals = sorted(v for v in (t.user_attrs.get(f"ap_after_{k}") for t in done)
                      if isinstance(v, float) and v == v)
        if len(vals) >= args.abort_min_trials:
            idx = int(args.abort_quantile * (len(vals) - 1))
            cuts[k] = vals[idx]
    return cuts


def build_objective(args):
    n_scen = args.scenarios

    def objective(trial):
        params = suggest_params(trial, args)
        cuts = _abort_thresholds(trial.study, n_scen, args)
        t_trial = time.time()
        aps, viols, masses, reliabs, eces, fdr95s = [], [], [], [], [], []

        for s in range(1, n_scen + 1):
            m = run_one_scenario(args, params, s)
            aps.append(m["ap"])
            if m["max_fdr_violation_n20"] == m["max_fdr_violation_n20"]:   # not NaN
                viols.append(m["max_fdr_violation_n20"])
            if m["total_mass_ratio"] == m["total_mass_ratio"]:
                masses.append(m["total_mass_ratio"])
            if m["hi_pip_reliab"] == m["hi_pip_reliab"]:
                reliabs.append(m["hi_pip_reliab"])
            if m["ece_hi"] == m["ece_hi"]:
                eces.append(m["ece_hi"])
            if m["fdr_at_95"] == m["fdr_at_95"]:
                fdr95s.append(m["fdr_at_95"])

            # Record the running AP so later trials can be compared against this
            # one at the same stage, then abort if we are already clearly behind.
            # Emit timing to OUR stdout (the PBS log). fb_objective.R's own
            # `seconds=` goes to a captured pipe, so without this line a trial
            # killed at walltime leaves no timing evidence whatsoever - which is
            # exactly the case a short canary run is trying to measure.
            print(f"[trial {trial.number}] scenario {s}/{n_scen} "
                  f"seconds={m['seconds']:.1f} ap={m['ap']:.4f} "
                  f"cum_min={(time.time() - t_trial) / 60:.1f}", flush=True)

            running_ap = sum(aps) / len(aps)
            trial.set_user_attr(f"ap_after_{s}", running_ap)
            if s in cuts and running_ap < cuts[s]:
                raise optuna.TrialPruned()

            # Pruning only applies to single-objective studies; Optuna does not
            # support pruners for multi-objective optimisation.
            if args.objective == "scalar":
                # Must use the SAME formula as the final score below, or the
                # pruner ranks trials on a different quantity than the one the
                # study is optimising.
                running = sum(aps) / len(aps)
                if viols:
                    running -= args.penalty * (sum(viols) / len(viols))
                if args.mass_penalty and masses:
                    running -= args.mass_penalty * abs(sum(masses) / len(masses) - 1.0)
                trial.report(running, s)
                if trial.should_prune():
                    raise optuna.TrialPruned()

        mean_ap    = sum(aps) / len(aps) if aps else float("nan")
        mean_viol  = sum(viols) / len(viols) if viols else 1.0
        mean_mass  = sum(masses) / len(masses) if masses else float("nan")
        mean_relb  = sum(reliabs) / len(reliabs) if reliabs else float("nan")
        mean_ece   = sum(eces) / len(eces) if eces else float("nan")
        mean_fdr95 = sum(fdr95s) / len(fdr95s) if fdr95s else float("nan")

        # Record everything, so any objective can be re-derived from the study
        # later without re-running trials.
        trial.set_user_attr("mean_ap", mean_ap)
        trial.set_user_attr("mean_fdr_violation_n20", mean_viol)
        trial.set_user_attr("mean_total_mass_ratio", mean_mass)
        trial.set_user_attr("mean_hi_pip_reliab", mean_relb)
        trial.set_user_attr("mean_ece_hi", mean_ece)
        trial.set_user_attr("mean_fdr_at_95", mean_fdr95)

        if args.objective == "triple":
            # The three quantities that jointly define a usable fine-mapper:
            # ranking (AP, up), selection safety (FDR violation, down) and
            # honesty of the reported probabilities (ece_hi, down). NaN scores
            # as clearly bad so a broken trial is dominated, not unsortable.
            # Excess FDR over the 0.05 budget promised at PIP >= 0.95. This
            # replaces max_fdr_violation_n20 as the safety objective: the
            # guarded max is a mid-range statistic (see fb_objective.R).
            fdr_excess = (max(0.0, mean_fdr95 - 0.05)
                          if mean_fdr95 == mean_fdr95 else 1.0)
            return (mean_ap,
                    fdr_excess,
                    mean_ece if mean_ece == mean_ece else 1.0)

        if args.objective == "multi":
            # Second objective, always minimised. 'violation' produced a
            # single-point front on fb_binary_ref500 because it never conflicts
            # with AP. 'mass' targets the trade-off that IS real once AP
            # saturates, and needs no weight chosen up front.
            # NaN would make the trial unsortable against the front, so a
            # mass ratio that could not be computed scores as clearly bad
            # rather than poisoning the study.
            if args.multi_second == "mass":
                second = abs(mean_mass - 1.0) if mean_mass == mean_mass else 10.0
            elif args.multi_second == "reliab":
                second = (1.0 - mean_relb) if mean_relb == mean_relb else 1.0
            else:
                second = mean_viol
            return mean_ap, second

        # Optional calibration term. Study fb_binary_ref500 (3102 trials) found
        # AP saturated at 0.809 while the mass ratio sat at ~2.0 in EVERY top
        # trial - i.e. the only metric with room left to improve is the one no
        # objective could see. |mass - 1| is symmetric because both directions
        # are miscalibration: >1 is over-confident (Iteration 002's 9-11x
        # failure mode), <1 leaves detectable mass on the table.
        score = mean_ap - args.penalty * mean_viol
        if args.mass_penalty and mean_mass == mean_mass:
            score -= args.mass_penalty * abs(mean_mass - 1.0)
        return score

    return objective


def main():
    ap = argparse.ArgumentParser(description="Optuna tuning for Functional BEATRICE")
    ap.add_argument("--sim", required=True, help="fixed tuning sim .rds (make_tuning_sim.R)")
    ap.add_argument("--study", required=True, help="study name (one per stratum)")
    ap.add_argument("--storage", required=True, help="journal file path (persists the study)")
    ap.add_argument("--objective", choices=["multi", "scalar", "triple"], default="multi",
                    help="multi: Pareto front of (AP, FDR violation), no arbitrary "
                         "weighting but no pruning. scalar: AP - penalty*violation, "
                         "enables pruning. Default: multi.")
    ap.add_argument("--rotate-regions", action="store_true",
                    help="score each scenario on ONE region, rotating which region "
                         "(and so which region size) across scenarios. Cost is "
                         "quadratic in p, so this is what makes a grid spanning "
                         "p=2000 affordable.")
    ap.add_argument("--abort-quantile", type=float, default=0.0,
                    help="abort a trial whose running AP falls below this quantile "
                         "of completed trials at the same scenario index. 0 (default) "
                         "disables. 0.25 aborts the worst quarter. Works in every "
                         "objective mode, including multi-objective, where Optuna "
                         "has no pruner of its own.")
    ap.add_argument("--abort-after", type=int, default=3,
                    help="never abort before this many scenarios have been scored.")
    ap.add_argument("--abort-min-trials", type=int, default=30,
                    help="never abort until this many trials have completed.")
    ap.add_argument("--multi-second", choices=["violation", "mass", "reliab"],
                    default="violation",
                    help="multi mode only: the second (minimised) objective. "
                         "'violation' = max_fdr_violation_n20 (the original); "
                         "'mass' = |total_mass_ratio - 1|; "
                         "'reliab' = 1 - hi_pip_reliab, i.e. the error rate among "
                         "PIP>=0.9 calls. 'reliab' is the one that matches the "
                         "measured Iteration 003 failure - but only if the tuning "
                         "sim spans phi=0.4, where that failure lives.")
    ap.add_argument("--mass-penalty", type=float, default=0.0,
                    help="scalar mode only: subtract this x |total_mass_ratio - 1| "
                         "from the score. 0 (default) reproduces the old objective. "
                         "Use when AP has saturated and calibration is the "
                         "remaining headroom.")
    ap.add_argument("--fix", nargs="*", default=None, metavar="NAME=VALUE",
                    help="remove parameters from the search space and hold them at "
                         "VALUE, e.g. --fix hierarchy_M=10.15 lambda_l1=0.1223. Use "
                         "for knobs an importance analysis shows are inert: tuning "
                         "them spends budget without moving the objective.")
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

    # --fix NAME=VALUE ... -> {name: value}. n_caus and sparse_concrete are
    # integer flags in beatrice_annot.py, so they must not be coerced to float.
    if args.fix:
        int_params = {"n_caus", "sparse_concrete"}
        fixed = {}
        for item in args.fix:
            if "=" not in item:
                ap.error(f"--fix expects NAME=VALUE, got '{item}'")
            k, _, v = item.partition("=")
            k = k.strip()
            try:
                fixed[k] = int(v) if k in int_params else float(v)
            except ValueError:
                ap.error(f"--fix {k}: '{v}' is not a number")
        args.fix = fixed
        print(f"[fix] held out of the search: {fixed}")

    storage = make_storage(args.storage)
    multi_obj = args.objective in ("multi", "triple")
    sampler = (optuna.samplers.NSGAIISampler(seed=None) if multi_obj
               else optuna.samplers.TPESampler(n_startup_trials=args.startup_trials))
    pruner = (optuna.pruners.MedianPruner(n_startup_trials=args.startup_trials,
                                          n_warmup_steps=1)
              if args.objective == "scalar" else optuna.pruners.NopPruner())

    kwargs = dict(study_name=args.study, storage=storage,
                  load_if_exists=True, sampler=sampler)
    if args.objective == "triple":
        kwargs["directions"] = ["maximize", "minimize", "minimize"]
    elif args.objective == "multi":
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


def _trial_ap(t):
    """Mean AP of a completed trial, however the study was scored.

    'multi' studies return (AP, violation) so AP is values[0]; 'scalar' studies
    return a penalised scalar, so AP has to come from the user attribute. Both
    modes set mean_ap, so prefer it and fall back to values[0].
    """
    v = t.user_attrs.get("mean_ap")
    if isinstance(v, float) and v == v:
        return v
    return t.values[0] if t.values else float("nan")


def report(study, args):
    complete = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
    by_state = {}
    for t in study.trials:
        by_state[t.state.name] = by_state.get(t.state.name, 0) + 1
    print(f"\n=== study '{study.study_name}': {len(complete)} completed trials ===")
    print("  states: " + ", ".join(f"{k}={v}" for k, v in sorted(by_state.items())))
    if not complete:
        return

    durs = sorted((t.datetime_complete - t.datetime_start).total_seconds() / 60.0
                  for t in complete
                  if t.datetime_start is not None and t.datetime_complete is not None)
    if durs:
        med = durs[len(durs) // 2]
        print(f"  trial duration (min): median {med:.1f}  min {durs[0]:.1f}  "
              f"max {durs[-1]:.1f}")
        print(f"  => {1440.0 / med * 20 / 24:.0f} trials per 72h x 20-worker run "
              f"(at {med:.0f} min/trial)")

    # Is the FDR-violation objective actually doing anything? If (nearly) every
    # trial achieves violation 0 it exerts no trade-off pressure, the Pareto
    # front degenerates to argmax(AP), and 'multi' is costing us the pruner
    # (Optuna supports pruners only for single-objective studies) for nothing.
    viols = [t.user_attrs.get("mean_fdr_violation_n20") for t in complete]
    viols = [v for v in viols if isinstance(v, float) and v == v]
    if viols:
        n_zero = sum(1 for v in viols if v <= 1e-12)
        print(f"  FDR-violation objective: {n_zero}/{len(viols)} trials at exactly 0 "
              f"({100.0 * n_zero / len(viols):.0f}%)"
              + ("  <- degenerate; prefer --objective scalar" if n_zero > 0.9 * len(viols) else ""))

    # Mass ratio is recorded but NOT optimised, so the search is free to inflate
    # it. Show its spread among the strongest trials so that stays visible.
    ranked = sorted(complete, key=lambda t: -_trial_ap(t))
    print("  top 5 by AP (mass ratio is NOT an objective - watch it):")
    for t in ranked[:5]:
        print(f"    AP={_trial_ap(t):.4f}  viol={t.user_attrs.get('mean_fdr_violation_n20', float('nan')):.4f}"
              f"  mass={t.user_attrs.get('mean_total_mass_ratio', float('nan')):.2f}"
              f"  reliab={t.user_attrs.get('mean_hi_pip_reliab', float('nan')):.2f}"
              f"  ece_hi={t.user_attrs.get('mean_ece_hi', float('nan')):.3f}  #{t.number}")

    # Which knobs actually move AP? A parameter whose importance is ~0 is worth
    # FIXING rather than tuning - it spends search budget for no return.
    try:
        imp = optuna.importance.get_param_importances(study, target=_trial_ap)
        print("  param importance for AP (fANOVA):")
        for k, v in imp.items():
            print(f"    {v:6.3f}  {k}")
    except Exception as e:                       # too few trials, missing sklearn, ...
        print(f"  (param importance unavailable: {e})")

    if args.objective == "triple":
        front = sorted(study.best_trials, key=lambda x: -x.values[0])
        print(f"Pareto front (AP up, FDR-excess@0.95 down, ece_hi down)  "
              f"[{len(front)} points]:")
        for t in front:
            print(f"  AP={t.values[0]:.4f}  fdr_xs={t.values[1]:.4f}  ece_hi={t.values[2]:.4f}"
                  f"  mass={t.user_attrs.get('mean_total_mass_ratio', float('nan')):.2f}"
                  f"  reliab={t.user_attrs.get('mean_hi_pip_reliab', float('nan')):.2f}"
                  f"  #{t.number}  {t.params}")
        return

    if args.objective == "multi":
        second = getattr(args, "multi_second", "violation")
        label = {"mass": "|mass-1|", "reliab": "1-hi_pip_reliab"}.get(second, "FDR violation")
        print(f"Pareto front (AP up, {label} down)  [{len(study.best_trials)} points]:")
        for t in sorted(study.best_trials, key=lambda x: -x.values[0]):
            print(f"  AP={t.values[0]:.4f}  obj2={t.values[1]:.4f}  "
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
