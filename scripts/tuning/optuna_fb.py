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


def parse_scenarios(args, n_sims):
    """--scenarios: one integer for every sim, or one per sim in order."""
    parts = [q.strip() for q in str(args.scenarios).split(",") if q.strip()]
    vals = [int(q) for q in parts]
    if len(vals) == 1:
        vals = vals * n_sims
    if len(vals) != n_sims:
        raise SystemExit(f"--scenarios gave {len(vals)} value(s) for {n_sims} --sim")
    if any(v < 1 for v in vals):
        raise SystemExit("--scenarios must be >= 1")
    return vals


def parse_sims(args):
    """--sim is repeatable, as LABEL=PATH, one per STRATUM.

    A stratum is a generative model x annotation type. They are never pooled:
    each contributes its own AP and its own FDR, and the objective aggregates
    those per-stratum numbers. One --sim without a label keeps the old
    single-stratum behaviour.
    """
    out = []
    for i, spec in enumerate(args.sim):
        if "=" in spec:
            label, _, path = spec.partition("=")
        else:
            label, path = (f"stratum{i+1}" if len(args.sim) > 1 else "default"), spec
        out.append((label, path))
    if any(not l for l, _ in out):
        raise SystemExit("--sim labels must be non-empty")
    # Labels MAY repeat: several sims can belong to one stratum, and then their
    # counts pool into it. That is how the rho sweep enters - sparse_inf is one
    # stratum whose rho levels the benchmark pools over, so rho=0.2 and rho=0.6
    # are two sims both labelled sparse_inf. Give them distinct labels instead
    # if you want the objective to target the worst rho rather than their pool.
    return out


def run_one_scenario(args, params, sim_path, stratum, scenario):
    """Invoke the R objective for a single scenario; return a dict of metrics."""
    cmd = [args.rscript, args.objective_script,
           "--sim", sim_path,
           "--stratum", stratum,
           "--method", args.method,
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
             "tp_at_95", "fp_at_95",
             "total_mass_ratio", "hi_pip_reliab", "ece_hi", "seconds",
             # scen_S is the scenario's TRUE causal count, reported by
             # fb_objective.R from simulation$scenarios[[k]]$S rather than
             # re-derived here from the expand.grid(S, phi) ordering. It keys
             # the per-S cells the `cells` objective worst-cases over.
             "scen_S", "scen_phi")}


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
    # Both floors were opened up. sparse_concrete went to lower_bound=1 and a
    # log-uniform integer, on the stated grounds of "the package's own declared
    # lower_bound=1".
    #
    # CORRECTED 2026-11-13. That bound is OUR FORK's. The published BEATRICE
    # (github.com/sayangsep/Beatrice-Finemapping, beatrice.py) declares
    #   DEFINE_integer('sparse_concrete', 50, ..., lower_bound=10)
    # and abseil refuses to start below 10, so v1 landing on 10 was landing on
    # the published MINIMUM, not on an arbitrary floor. sparse_concrete is the
    # top-K cap on how many variants may be simultaneously causal in one
    # configuration (joint_trainer.py:190 topk, original trainer.py:211), so a
    # value below the true causal count makes the model structurally unable to
    # represent the truth. Study v3 duly walked to sparse_concrete=2, which
    # scored a 5x better FDR by becoming a two-variant detector: measured
    # against the shipped config it gained 26% AP at S=1 and lost 48% at S=10.
    # The floor is now an argument, defaulting to the published 10.
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
        "sparse_concrete":      sint("sparse_concrete", args.sc_min, args.sc_max,
                                     log=True),
        # The prior's SCALE, unidentified until anchor_prior_mass existed.
        # Published BEATRICE does not learn the prior at all (upstream
        # scripts/trainer.py:386 sets p_0 = 1/m, so sum(p_0) = 1 exactly);
        # LassoNetPrior fixes the prior's SHAPE from annotations and left its
        # total free, and because `abf` uses a Bernoulli configuration prior,
        # sum_j PIP_j is governed by sum_j p0_j. Mass was therefore the one
        # constrained quantity that NO search parameter controlled directly -
        # gamma and sigma_sq could only push it sideways, at a cost in AP, which
        # is exactly what v4's lowest-mass trials did (gamma pinned at its 0.8
        # ceiling, min-cell AP 0.154 against a 0.225 baseline). Measured
        # locally, prior_mass=1.0 took sparse_inf S=10 mass from 2.00 to 1.00.
        # --pm-max <= 0 keeps the anchor out of the space entirely, which is
        # every run before this.
        "prior_mass":           (0.0 if args.pm_max <= 0 else
                                 sfloat("prior_mass", max(args.pm_min, 1e-3),
                                        args.pm_max, log=True)),
        # ADDED 2026-11-10. Both act BEFORE pip.csv is written, so both move the
        # PIPs and therefore AP and FDR; neither had ever been exposed by the R
        # wrappers, so every run so far used the Python defaults (0.1, 0.2).
        # gamma builds the reduced space of binary vectors during training.
        # gamma_key seeds the stepwise search whose abf() calls enlarge the memo
        # before the PIPs are recomputed - and reformat_memo's own docstring ties
        # it to the saturation bug ("no PIP ever reached gamma_key"), so its
        # default was implicitly calibrated against flattened PIPs.
        # gamma_coverage, gamma_selection and purity are deliberately NOT here:
        # generate_credible_sets.py reads them only after pip.csv is written, so
        # they cannot change AP or FDR at all.
        # CEILING RAISED 0.5 -> 0.8 (2026-11-13, user). With sparse_concrete
        # floored at its published 10, gamma is the only remaining lever on
        # configuration size: _abf keeps the top-K entries with cc > gamma, so
        # raising it shrinks configurations from INSIDE the top-K instead of
        # truncating it. The original declares gamma with no bounds at all, and
        # v3's lowest-mass sparse_concrete=7 configuration sat at gamma 0.478,
        # 4% under the old ceiling - the direction v4 needs. 0.8 matches
        # gamma_key's existing ceiling. Above ~0.8 few concrete values clear the
        # threshold, the training memo can come back near-empty and the PIPs
        # fall out of the stepwise search rather than the variational fit.
        "gamma":                sfloat("gamma", 0.01, 0.8, log=True),
        "gamma_key":            sfloat("gamma_key", 0.02, 0.8, log=True),
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
    sims   = parse_sims(args)
    counts = parse_scenarios(args, len(sims))
    # The flat list of work units a trial does: (stratum, sim path, scenario).
    units  = [(lab, path, sc)
              for (lab, path), n in zip(sims, counts)
              for sc in range(1, n + 1)]
    n_scen = len(units)

    def objective(trial):
        params = suggest_params(trial, args)
        cuts = _abort_thresholds(trial.study, n_scen, args)
        t_trial = time.time()
        aps = []
        # Per-stratum accumulators. Nothing is ever pooled ACROSS strata: AP and
        # ece are averaged within a stratum and then over strata (equal weight
        # per stratum), and a stratum's FDR is pooled from its own COUNTS.
        per = {lab: {"ap": [], "ece": [], "relb": [], "mass": [],
                     "viol": [], "tp": 0.0, "fp": 0.0} for lab, _ in sims}
        # Cells are (stratum, true causal count). Study v3 worst-cased FDR over
        # strata only, which left the S axis poolable - and a config that buys
        # S=1 with S=10 then averages into looking fine. Keyed lazily because S
        # arrives with the scenario, not with the --sim flag.
        cells = {}

        for s, (lab, sim_path, sc) in enumerate(units, start=1):
            m = run_one_scenario(args, params, sim_path, lab, sc)
            d = per[lab]
            aps.append(m["ap"])
            d["ap"].append(m["ap"])
            for key, dst in (("ece_hi", "ece"), ("hi_pip_reliab", "relb"),
                             ("total_mass_ratio", "mass"),
                             ("max_fdr_violation_n20", "viol")):
                if m[key] == m[key]:                      # not NaN
                    d[dst].append(m[key])
            # Counts, not the rate: mean(fp/n) != sum(fp)/sum(n).
            if m["tp_at_95"] == m["tp_at_95"] and m["fp_at_95"] == m["fp_at_95"]:
                d["tp"] += m["tp_at_95"]; d["fp"] += m["fp_at_95"]

            sS = m["scen_S"]
            key = f"{lab}|S={int(sS)}" if sS == sS else f"{lab}|S=?"
            cd = cells.setdefault(key, {"ap": [], "ece": [], "tp": 0.0, "fp": 0.0})
            cd["ap"].append(m["ap"])
            if m["ece_hi"] == m["ece_hi"]:
                cd["ece"].append(m["ece_hi"])
            if m["tp_at_95"] == m["tp_at_95"] and m["fp_at_95"] == m["fp_at_95"]:
                cd["tp"] += m["tp_at_95"]; cd["fp"] += m["fp_at_95"]

            print(f"[trial {trial.number}] {lab} scenario {sc} "
                  f"({s}/{n_scen}) seconds={m['seconds']:.1f} ap={m['ap']:.4f} "
                  f"cum_min={(time.time() - t_trial) / 60:.1f}", flush=True)

            running_ap = sum(aps) / len(aps)
            trial.set_user_attr(f"ap_after_{s}", running_ap)
            if s in cuts and running_ap < cuts[s]:
                raise optuna.TrialPruned()

            if args.objective == "scalar":
                running = running_ap
                vs = [v for dd in per.values() for v in dd["viol"]]
                ms = [v for dd in per.values() for v in dd["mass"]]
                if vs:
                    running -= args.penalty * (sum(vs) / len(vs))
                if args.mass_penalty and ms:
                    running -= args.mass_penalty * abs(sum(ms) / len(ms) - 1.0)
                trial.report(running, s)
                if trial.should_prune():
                    raise optuna.TrialPruned()

        def mean(xs, default=float("nan")):
            return sum(xs) / len(xs) if xs else default

        # --- per-stratum, then across strata --------------------------------
        strat_ap, strat_ece, strat_fdr = {}, {}, {}
        for lab, _ in sims:
            d = per[lab]
            strat_ap[lab]  = mean(d["ap"])
            strat_ece[lab] = mean(d["ece"])
            n = d["tp"] + d["fp"]
            strat_fdr[lab] = (d["fp"] / n) if n > 0 else float("nan")
            trial.set_user_attr(f"ap_{lab}", strat_ap[lab])
            trial.set_user_attr(f"fdr_at_95_{lab}", strat_fdr[lab])
            trial.set_user_attr(f"ece_hi_{lab}", strat_ece[lab])
            trial.set_user_attr(f"nsel_at_95_{lab}", n)

        mean_ap    = mean([v for v in strat_ap.values() if v == v])
        mean_ece   = mean([v for v in strat_ece.values() if v == v])
        mean_relb  = mean([v for dd in per.values() for v in dd["relb"]])
        # Mass is reported PER STRATUM and constrained on the WORST one.
        # v4 pooled it - the mean over every scenario of every stratum - which
        # is the one metric here that broke the never-pool rule, and since it
        # was the binding constraint it decided the whole search. Attainable
        # mass differs by generative model (measured locally at the shipped
        # config: 1.04 under sparse on real LD, 1.8-6.4 under sparse_inf), so a
        # window on their average can only be met by trading one stratum off
        # against the other - and with [0.5, 2.0] it was met by neither:
        # 0/76 trials feasible, Pareto front empty. Worst-stratum is the same
        # choice already made for FDR below, for the same reason.
        strat_mass = {}
        for lab, _ in sims:
            mm = [v for v in per[lab]["mass"] if v == v]
            strat_mass[lab] = mean(mm) if mm else float("nan")
            trial.set_user_attr(f"total_mass_ratio_{lab}", strat_mass[lab])
        mean_mass  = mean([v for v in strat_mass.values() if v == v])
        mean_viol  = mean([v for dd in per.values() for v in dd["viol"]], 1.0)
        # The reported FDR is the WORST stratum, not the average. Averaged, a
        # configuration that is honest under `sparse` and badly overconfident
        # under `sparse_inf` looks mid-range, and the search can buy objective
        # value by fixing only the easy stratum. The max cannot be satisfied
        # that way. A stratum that selected nothing at 0.95 contributes no rate.
        fdrs = [v for v in strat_fdr.values() if v == v]
        worst_fdr = max(fdrs) if fdrs else float("nan")

        trial.set_user_attr("mean_ap", mean_ap)
        trial.set_user_attr("mean_fdr_violation_n20", mean_viol)
        trial.set_user_attr("mean_total_mass_ratio", mean_mass)
        trial.set_user_attr("mean_hi_pip_reliab", mean_relb)
        trial.set_user_attr("mean_ece_hi", mean_ece)
        trial.set_user_attr("worst_fdr_at_95", worst_fdr)
        trial.set_user_attr("worst_stratum",
                            max(strat_fdr, key=lambda k: (strat_fdr[k] == strat_fdr[k],
                                                          strat_fdr[k]))
                            if fdrs else "")

        # --- per-cell summaries, for the `cells` objective ------------------
        cell_ap, cell_fdr = {}, {}
        for key, cd in cells.items():
            cell_ap[key] = mean(cd["ap"])
            n = cd["tp"] + cd["fp"]
            cell_fdr[key] = (cd["fp"] / n) if n > 0 else float("nan")
            trial.set_user_attr(f"ap_{key}", cell_ap[key])
            trial.set_user_attr(f"fdr_at_95_{key}", cell_fdr[key])
            trial.set_user_attr(f"nsel_at_95_{key}", n)
        aps_c = [v for v in cell_ap.values() if v == v]
        fdr_c = [v for v in cell_fdr.values() if v == v]
        min_cell_ap   = min(aps_c) if aps_c else float("nan")
        worst_cell_fdr = max(fdr_c) if fdr_c else float("nan")
        trial.set_user_attr("min_cell_ap", min_cell_ap)
        trial.set_user_attr("worst_cell_fdr_at_95", worst_cell_fdr)
        trial.set_user_attr("worst_cell",
                            max(cell_fdr, key=lambda k: (cell_fdr[k] == cell_fdr[k],
                                                         cell_fdr[k]))
                            if fdr_c else "")

        # Posterior mass as a HARD constraint, not a watched number. v3 printed
        # it as a warning and the search walked straight past it: every point at
        # the low-FDR end reached mass ~1 by capping sparse_concrete below the
        # true causal count. NSGA-II reads this through constraints_func;
        # <= 0 is feasible.
        _mv = [max(args.mass_lo - v, v - args.mass_hi)
               for v in strat_mass.values() if v == v]
        mass_c = (max(_mv) if _mv and len(_mv) == len(strat_mass) else 1e3)
        trial.set_user_attr("constraint", (float(mass_c),))
        trial.set_user_attr("mass_feasible", bool(mass_c <= 0))

        if args.objective == "cells":
            # As `triple`, but both AP and FDR are worst-cased over
            # (stratum x S) cells rather than over strata alone, so no
            # configuration can pay for one causal count with another.
            fdr_excess = (max(0.0, worst_cell_fdr - 0.05)
                          if worst_cell_fdr == worst_cell_fdr else 1.0)
            return (min_cell_ap if min_cell_ap == min_cell_ap else 0.0,
                    fdr_excess,
                    mean_ece if mean_ece == mean_ece else 1.0)

        if args.objective == "triple":
            # Ranking (AP, up), selection safety (worst-stratum FDR excess over
            # the 0.05 budget promised at PIP >= 0.95, down) and honesty of the
            # reported probabilities (ece_hi, down). NaN scores as clearly bad so
            # a broken trial is dominated rather than unsortable.
            fdr_excess = (max(0.0, worst_fdr - 0.05)
                          if worst_fdr == worst_fdr else 1.0)
            return (mean_ap,
                    fdr_excess,
                    mean_ece if mean_ece == mean_ece else 1.0)

        if args.objective == "multi":
            if args.multi_second == "mass":
                second = abs(mean_mass - 1.0) if mean_mass == mean_mass else 10.0
            elif args.multi_second == "reliab":
                second = (1.0 - mean_relb) if mean_relb == mean_relb else 1.0
            elif args.multi_second == "fdr":
                second = (max(0.0, worst_fdr - 0.05)
                          if worst_fdr == worst_fdr else 1.0)
            else:
                second = mean_viol
            return mean_ap, second

        score = mean_ap - args.penalty * mean_viol
        if args.mass_penalty and mean_mass == mean_mass:
            score -= args.mass_penalty * abs(mean_mass - 1.0)
        return score

    return objective


def main():
    ap = argparse.ArgumentParser(description="Optuna tuning for Functional BEATRICE")
    ap.add_argument("--sim", required=True, action="append", metavar="LABEL=PATH",
                    help="fixed tuning sim .rds (make_tuning_sim.R). Repeatable, "
                         "once per STRATUM, as LABEL=PATH. Strata are never "
                         "pooled: each yields its own AP and its own FDR.")
    ap.add_argument("--method", default="fb_xregion",
                    help="which FB to tune. Default fb_xregion, the cross-region "
                         "method the benchmark ships. 'functional_beatrice' is the "
                         "per-region method and tunes a DIFFERENT model.")
    ap.add_argument("--study", required=True, help="study name (one per stratum)")
    ap.add_argument("--storage", required=True, help="journal file path (persists the study)")
    ap.add_argument("--objective",
                    choices=["multi", "scalar", "triple", "cells"], default="multi",
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
    ap.add_argument("--multi-second", choices=["violation", "mass", "reliab", "fdr"],
                    default="violation",
                    help="multi mode only: the second (minimised) objective. "
                         "'violation' = max_fdr_violation_n20 (the original); "
                         "'mass' = |total_mass_ratio - 1|; "
                         "'reliab' = 1 - hi_pip_reliab, i.e. the error rate among "
                         "PIP>=0.9 calls. 'reliab' is the one that matches the "
                         "measured Iteration 003 failure - but only if the tuning "
                         "sim spans phi=0.4, where that failure lives.")
    # The published BEATRICE declares sparse_concrete with lower_bound=10 and
    # refuses to start below it. Defaulting the floor there keeps the search
    # inside the method as published; raise sc_max, never lower sc_min, unless
    # you mean to leave the published parameter space.
    ap.add_argument("--sc-min", type=int, default=10,
                    help="sparse_concrete floor (published lower_bound is 10)")
    ap.add_argument("--sc-max", type=int, default=300,
                    help="sparse_concrete ceiling (the original declares none)")
    ap.add_argument("--pm-min", type=float, default=0.5,
                    help="Lower bound for prior_mass, the per-region anchor on "
                         "sum(p_0). Only used when --pm-max > 0.")
    ap.add_argument("--pm-max", type=float, default=0.0,
                    help="Upper bound for prior_mass. <= 0 (the default) keeps "
                         "the anchor off and out of the search space, i.e. the "
                         "unanchored learned scale every run so far used. A "
                         "sensible range is 0.5-20: 1.0 is published BEATRICE's "
                         "scale, and the prior expected causal count should be "
                         "able to reach S.")
    ap.add_argument("--mass-lo", type=float, default=0.5,
                    help="feasible sum(PIP)/#causal floor (objective 'cells')")
    ap.add_argument("--mass-hi", type=float, default=2.0,
                    help="feasible sum(PIP)/#causal ceiling (objective 'cells')")
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
    ap.add_argument("--scenarios", default="4",
                    help="scenarios per trial: one integer for every sim, or a "
                         "comma-separated list, one per --sim, in order.")
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
    multi_obj = args.objective in ("multi", "triple", "cells")
    # `cells` carries a hard feasibility constraint (posterior mass), which
    # NSGA-II honours through constrained domination: an infeasible trial loses
    # to every feasible one but still informs the search, which is why this is
    # a constraint rather than a pruned trial. A trial that died before setting
    # the attr counts as infeasible.
    def _constraints(frozen):
        c = frozen.user_attrs.get("constraint")
        try:
            return tuple(float(v) for v in c)
        except (TypeError, ValueError):
            return (1e3,)
    if args.objective == "cells":
        sampler = optuna.samplers.NSGAIISampler(seed=None,
                                                constraints_func=_constraints)
    elif multi_obj:
        sampler = optuna.samplers.NSGAIISampler(seed=None)
    else:
        sampler = optuna.samplers.TPESampler(n_startup_trials=args.startup_trials)
    pruner = (optuna.pruners.MedianPruner(n_startup_trials=args.startup_trials,
                                          n_warmup_steps=1)
              if args.objective == "scalar" else optuna.pruners.NopPruner())

    kwargs = dict(study_name=args.study, storage=storage,
                  load_if_exists=True, sampler=sampler)
    if args.objective in ("triple", "cells"):
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
    print("[optuna] strata: " + ", ".join(f"{l}={pth}" for l, pth in parse_sims(args))
          + f"   method={args.method}")
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



def _rank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    for pos, i in enumerate(order):
        r[i] = float(pos)
    return r


def _corr(a, b):
    n = len(a)
    if n < 3:
        return float("nan")
    ma, mb = sum(a) / n, sum(b) / n
    va = sum((x - ma) ** 2 for x in a) ** 0.5
    vb = sum((y - mb) ** 2 for y in b) ** 0.5
    if va == 0 or vb == 0:
        return float("nan")
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (va * vb)


def _compute_confound(complete):
    """Rank correlation of each hyperparameter with trial WALL-CLOCK time.

    max_iter is excluded from the search space because it trades compute for
    convergence: tuning it lets the search buy objective value with runtime and
    makes trials incomparable. Any parameter that strongly predicts runtime is
    open to the same objection, so it should be visible rather than inferred.
    """
    rows = []
    durs, names = [], set()
    recs = []
    for t in complete:
        if t.datetime_start is None or t.datetime_complete is None:
            continue
        d = (t.datetime_complete - t.datetime_start).total_seconds()
        recs.append((t.params, d, _trial_ap(t)))
        names.update(t.params.keys())
    if len(recs) < 10:
        return rows
    for nm in sorted(names):
        pv, dv, av = [], [], []
        for params, d, ap in recs:
            if nm in params and isinstance(params[nm], (int, float)):
                pv.append(float(params[nm])); dv.append(d); av.append(ap)
        if len(pv) < 10:
            continue
        rows.append((nm,
                     _corr(_rank(pv), _rank(dv)),
                     _corr(_rank(pv), _rank(av))))
    return rows


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
        # 72h x 20 workers = 86400 core-MINUTES. (An earlier version divided
        # 1440 core-hours by a duration in minutes, mixing units and
        # understating throughput by ~70x.)
        core_minutes = 72 * 60 * 20
        print(f"  => ~{core_minutes / med:.0f} trials per 72h x 20-worker run "
              f"(at {med:.0f} min/trial, before any early-abort saving)")

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

    conf = _compute_confound(complete)
    if conf:
        print("  rank correlation of each parameter with (runtime, AP):")
        for nm, cd, ca in sorted(conf, key=lambda r: -abs(r[1])):
            flag = "   <- also predicts RUNTIME: compute confound?" if abs(cd) > 0.5 else ""
            print(f"    {nm:<22} runtime {cd:+.2f}   AP {ca:+.2f}{flag}")

    # Which knobs actually move AP? A parameter whose importance is ~0 is worth
    # FIXING rather than tuning - it spends search budget for no return.
    try:
        imp = optuna.importance.get_param_importances(study, target=_trial_ap)
        caveat = ("   <- TOO FEW TRIALS TO INTERPRET"
                  if len(complete) < 50 else "")
        print(f"  param importance for AP (fANOVA):{caveat}")
        for k, v in imp.items():
            print(f"    {v:6.3f}  {k}")
    except Exception as e:                       # too few trials, missing sklearn, ...
        print(f"  (param importance unavailable: {e})")

    if args.objective == "cells":
        feas = [t for t in complete if t.user_attrs.get("mass_feasible")]
        print(f"  mass-feasible trials: {len(feas)}/{len(complete)} "
              f"(sum(PIP)/#causal within [{args.mass_lo}, {args.mass_hi}])")
        if not feas:
            print("  NO feasible trial. Either the window is too tight or the")
            print("  post-fix model cannot reach calibrated mass inside the")
            print("  published sparse_concrete range - which is a finding, not")
            print("  a config. Check the sparse_concrete spread above.")
        front = sorted(study.best_trials, key=lambda x: -x.values[0])
        print(f"Pareto front (WORST-CELL AP up, WORST-CELL FDR-excess@0.95 down, "
              f"ece_hi down)  [{len(front)} points]:")
        for t in front:
            ua = t.user_attrs
            flag = "" if ua.get("mass_feasible") else "  [INFEASIBLE mass]"
            print(f"  minAP={t.values[0]:.4f}  worst_fdr_xs={t.values[1]:.4f}  "
                  f"ece_hi={t.values[2]:.4f}  meanAP={ua.get('mean_ap', float('nan')):.4f}  "
                  f"mass={ua.get('mean_total_mass_ratio', float('nan')):.2f}  "
                  f"worst={ua.get('worst_cell', '')}  #{t.number}{flag}")
            cells = sorted(k[3:] for k in ua if k.startswith("ap_") and "|S=" in k)
            if cells:
                print("      " + "  ".join(
                    f"{c}:AP={ua.get('ap_' + c, float('nan')):.3f}/"
                    f"FDR={ua.get('fdr_at_95_' + c, float('nan')):.3f}"
                    for c in cells))
            print(f"      {t.params}")
        return

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
