#!/usr/bin/env python3
# =============================================================================
# scripts/tuning/test_optuna_fb.py
#
# Unit tests for the `cells` objective in optuna_fb.py, with run_one_scenario
# stubbed out so no trainer runs. Checks the two properties the objective exists
# to enforce:
#   * nothing is pooled ACROSS strata - AP, FDR, ECE and (since 2026-11-14)
#     posterior mass are each kept per stratum, and the mass constraint is the
#     WORST stratum's violation. v4 pooled mass, which let one stratum's
#     calibration pay for the other's and left 0/76 trials feasible.
#   * FDR and AP are worst-cased over (stratum x S) cells, not averaged, so a
#     configuration that collapses at S=10 cannot hide behind S=1.
#
# Run from the project root:  python3 scripts/tuning/test_optuna_fb.py
# =============================================================================
import sys, types, importlib.util as iu
sys.argv = ["x"]
spec = iu.spec_from_file_location("ofb", "scripts/tuning/optuna_fb.py")
ofb = iu.module_from_spec(spec); spec.loader.exec_module(ofb)

FIX = dict(lambda_l1=0.01, prior_regularisation=1.0, sigma_sq=0.05,
           hierarchy_M=10.0, n_caus=5, sparse_concrete=10, gamma=0.1,
           gamma_key=0.2)

def mkargs(objective="cells", mass_lo=0.5, mass_hi=2.0):
    return types.SimpleNamespace(
        sim=["sparse=/s.rds", "sparse_inf=/i1.rds", "sparse_inf=/i2.rds"],
        scenarios="3", regions=10, objective=objective,
        mass_lo=mass_lo, mass_hi=mass_hi, abort_quantile=0.0, abort_after=3,
        abort_min_trials=30, penalty=0.0, mass_penalty=0.0, fix=dict(FIX),
        multi_second="fdr", tune_gamma_coverage=False, sc_min=10, sc_max=300,
        pm_min=0.5, pm_max=0.0)

class Study:
    def get_trials(self, deepcopy=False): return []
class Trial:
    def __init__(self): self.number = 0; self.attrs = {}; self.study = Study()
    def set_user_attr(self, k, v): self.attrs[k] = v
    @property
    def user_attrs(self): return self.attrs

# scenario k -> S (expand.grid varies S fastest): 1->1, 2->3, 3->10
S_OF = {1: 1, 2: 3, 3: 10}

def runner(ap_by_S, mass, fdr_by_S):
    def f(args, params, sim_path, stratum, scenario):
        S = S_OF[scenario]
        fdr = fdr_by_S[S]
        tp, fp = 100.0 * (1 - fdr), 100.0 * fdr
        return {"ap": ap_by_S[S], "ece_hi": 0.2, "hi_pip_reliab": 0.8,
                "total_mass_ratio": mass, "max_fdr_violation_n20": 0.0,
                "tp_at_95": tp, "fp_at_95": fp, "seconds": 1.0,
                "fdr_at_90": fdr, "fdr_at_95": fdr,
                "scen_S": float(S), "scen_phi": 0.2}
    return f

def run(objective, ap_by_S, mass, fdr_by_S):
    ofb.run_one_scenario = runner(ap_by_S, mass, fdr_by_S)
    a = mkargs(objective)
    t = Trial()
    vals = ofb.build_objective(a)(t)
    return vals, t.attrs

ok = True
def chk(name, got, want, tol=1e-9):
    global ok
    good = abs(got - want) < tol if isinstance(want, float) else got == want
    ok &= good
    print(f"  [{'PASS' if good else 'FAIL'}] {name}: got {got!r} want {want!r}")

# --- A: the #227 pathology - a two-variant detector -------------------------
A_ap  = {1: 0.9, 3: 0.5, 10: 0.15}
A_fdr = {1: 0.05, 3: 0.08, 10: 0.10}
print("A: two-variant detector (high S=1, collapsed S=10), mass 0.94")
(vA, aA) = run("cells", A_ap, 0.94, A_fdr)
chk("mean AP (pooled over S)", round(aA["mean_ap"], 4), 0.5167, 1e-3)
chk("min cell AP == the S=10 cell", round(vA[0], 4), 0.15, 1e-9)
chk("6 cells (2 strata x 3 S)", len([k for k in aA if k.startswith("ap_") and "|S=" in k]), 6)
chk("worst cell is an S=10 cell", "S=10" in aA["worst_cell"], True)
chk("mass feasible", aA["mass_feasible"], True)
chk("constraint <= 0 when feasible", aA["constraint"][0] <= 0, True)

# --- B: balanced, LOWER mean AP but HIGHER worst cell -----------------------
B_ap  = {1: 0.55, 3: 0.5, 10: 0.45}
B_fdr = {1: 0.05, 3: 0.08, 10: 0.10}
print("B: balanced (lower mean AP, much higher worst cell), mass 1.1")
(vB, aB) = run("cells", B_ap, 1.1, B_fdr)
chk("mean AP lower than A", aB["mean_ap"] < aA["mean_ap"], True)
chk("min cell AP higher than A", vB[0] > vA[0], True)
print("  -> under 'cells' B beats A on objective 1; under 'triple' A would win on AP")
(tA, _) = run("triple", A_ap, 0.94, A_fdr)
(tB, _) = run("triple", B_ap, 1.1, B_fdr)
chk("triple would have ranked A above B", tA[0] > tB[0], True)

# --- C: mass inflation is INFEASIBLE ----------------------------------------
print("C: same as B but mass 4.7")
(vC, aC) = run("cells", B_ap, 4.7, B_fdr)
chk("mass infeasible", aC["mass_feasible"], False)
chk("constraint > 0", aC["constraint"][0] > 0, True)
chk("constraint == mass - hi", round(aC["constraint"][0], 4), 2.7, 1e-9)

# --- D: FDR is worst CELL, pooled from counts, not mean of rates ------------
D_fdr = {1: 0.02, 3: 0.02, 10: 0.60}
print("D: FDR 0.02 at S=1,3 and 0.60 at S=10")
(vD, aD) = run("cells", B_ap, 1.0, D_fdr)
chk("worst-cell FDR == 0.60 (not the 0.213 mean)", round(aD["worst_cell_fdr_at_95"], 4), 0.6, 1e-9)
chk("objective 2 == worst excess over 0.05", round(vD[1], 4), 0.55, 1e-9)
chk("an S=10 cell pools BOTH sparse_inf sims", aD["nsel_at_95_sparse_inf|S=10"], 200.0)
chk("the sparse S=10 cell has only one sim", aD["nsel_at_95_sparse|S=10"], 100.0)


# --- E: mass is constrained on the WORST stratum, never pooled --------------
# sparse sits at 1.0 (inside [0.5, 2.0]); sparse_inf sits at 6.0 (outside).
# Pooling the three sims would give mean (1.0 + 6.0 + 6.0)/3 = 4.333 and a
# violation of 2.333. Worst-stratum gives 6.0 and a violation of 4.0. This is
# the v4 defect: the pooled number let one stratum's calibration pay for the
# other's, and with the real spread no point satisfied the window at all.
def runner_by_stratum(mass_by_stratum):
    def f(args, params, sim_path, stratum, scenario):
        S = S_OF[scenario]
        return {"ap": 0.5, "ece_hi": 0.2, "hi_pip_reliab": 0.8,
                "total_mass_ratio": mass_by_stratum[stratum],
                "max_fdr_violation_n20": 0.0,
                "tp_at_95": 95.0, "fp_at_95": 5.0, "seconds": 1.0,
                "fdr_at_90": 0.05, "fdr_at_95": 0.05,
                "scen_S": float(S), "scen_phi": 0.2}
    return f

print("E: sparse mass 1.0, sparse_inf mass 6.0 (pooled mean would be 4.333)")
ofb.run_one_scenario = runner_by_stratum({"sparse": 1.0, "sparse_inf": 6.0})
aE = mkargs("cells"); tE = Trial()
ofb.build_objective(aE)(tE)
chk("per-stratum mass recorded: sparse", tE.attrs["total_mass_ratio_sparse"], 1.0)
chk("per-stratum mass recorded: sparse_inf",
    tE.attrs["total_mass_ratio_sparse_inf"], 6.0)
chk("constraint is the WORST stratum (6.0 - 2.0), not the pooled 4.333 - 2.0",
    tE.attrs["constraint"][0], 4.0)
chk("infeasible on the worst stratum", tE.attrs["mass_feasible"], False)

# and the converse: both strata inside the window -> feasible
print("E2: both strata at 1.0 -> feasible, constraint <= 0")
ofb.run_one_scenario = runner_by_stratum({"sparse": 1.0, "sparse_inf": 1.0})
aE2 = mkargs("cells"); tE2 = Trial()
ofb.build_objective(aE2)(tE2)
chk("feasible when every stratum is inside", tE2.attrs["mass_feasible"], True)
chk("constraint <= 0", tE2.attrs["constraint"][0] <= 0, True)

print("\nALL PASS" if ok else "\nFAILURES ABOVE")
sys.exit(0 if ok else 1)
