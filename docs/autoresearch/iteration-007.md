# Iteration 007 — LD misspecification: reference panels of four sizes

**Status:** ready to submit.
**Scope:** exploratory and TEMPORARY, like Iterations 005 and 006. Nothing in `R/`
changes and Iteration 004's grid is untouched; see
[iteration-007-REVERT.md](iteration-007-REVERT.md). Iteration 004 remains the
standard. **The design below (n = 1000, small regions, reference-panel LD) is for
Iteration 007 only.**

> **Not comparable to Iteration 004**: n, region sizes, the S and φ levels, the
> enrichment and the LD all differ. Nor is it comparable to Iterations 005 or 006.
> Within this iteration, each method's in-sample rows are its control.

---

## 1. The question

Every result in the LSR uses in-sample LD, the most favourable case. How does
fine-mapping change when a method is given the LD of an independent reference
panel instead of the GWAS sample's own, and how does that depend on the size of
the panel?

## 2. Design (user-specified 2026-09-14)

**20 rows** = 5 LD levels × 2 models × 2 annotation types.

| factor | levels |
|---|---|
| LD | in-sample; independent reference panel of 500, 750, 1500 or 2000 individuals |
| model | sparse; sparse_inf with p_causal = 0.6 |
| annotations | binary, continuous; 10 tracks, the first 5 at fold 5.4 and the rest inert |

Within each row:

| | |
|---|---|
| S | 1, 3, 5 |
| φ | 0.1, 0.4 |
| iterations | 10 |
| regions | 10: two each of 100, 150, 200, 300 and 400 variants |
| n (GWAS) | 1000 |

**60 scenarios/row × 20 rows = 1,200 scenarios.**

The region sizes are Iteration 004's divided by 5, so p/n ≤ 0.4, as in Iteration
004. Every region is smaller than every sample, GWAS or panel.

**Unpaired.** Each row is a separate `run_simulation()` call with seed 1000 + row,
which is how `n_ref` has always worked. The LD levels therefore do not share
genotypes, causal variants or z-scores, and comparisons between levels include
draw-to-draw variation. Pairing them would need simulator changes, for two reasons
(`R/simulate_genotypes.R:413-449`):
- the panel draw consumes random numbers before the phenotypes are simulated;
- variants that are monomorphic in the panel are dropped from both samples.

**Methods (12)**, the user's list, selected with `FMB_ITER007_METHODS=1`:
- `susie`, `susie_inf`
- `finemap`, `finemap_inf`
- `finimom`, `sparsepro`
- `funmap`, `paintor`
- `polyfun_oracle`, `polyfun_ldsc`
- `beatrice`, `fb_xregion`

`fb_xregion` runs the current implementation. The lambda sweep showed it
reproduces Iteration 004's average precision on average.

**Caveat: tuned at n = 5000.** The BEATRICE-family hyperparameters (σ² = 0.0899 and
the rest) were tuned on a grid at n = 5000. The LSR's hyperparameter section notes
that the appropriate σ² depends on the GWAS sample size. How each method changes
across LD levels is unaffected. Comparisons between methods at n = 1000 carry this
caveat.

## 3. What each method receives

`run_simulation()` keeps `LD_true = cor(X)` and hands methods `LD = cor(X_ref)`
whenever a panel is drawn (`R/run_simulation.R:330-337`). The z-scores always come
from the GWAS sample: `R/simulate_phenotypes.R:1023-1029` computes them from `X`,
and that file never references `X_ref`.

Checked wrapper by wrapper:

- **`susie`, `susie_inf`:** z plus LD. The individual-level path needs
  `use_individual = TRUE` (`R/wrapper_susie.R:295-303`,
  `R/wrapper_susie_inf.R:326-334`), and `METHOD_ARGS` never sets it.
- **`finemap`, `finemap_inf`, `funmap`, `paintor`, `beatrice`, `sparsepro`:** `region_geno$LD`.
- **`polyfun_ldsc`, `polyfun_oracle`:** `susie_rss` with `region_geno$LD`.
  `polyfun_ldsc`'s LD scores come from the same matrix
  (`R/wrapper_polyfun_ldsc.R:267-271`).
- **`finimom`:** `region_geno$LD`. Its `insampleLD` flag is set to FALSE
  automatically when a panel exists (`R/wrapper_finimom.R:208-209`).
- **`fb_xregion`:** `genotypes[[i]]$LD` for every region of the joint fit
  (`R/wrapper_fb_joint.R:152`).

## 4. Checks

- **Grid:** `generate_params_grid_iter007.R` writes 20 rows. They pass
  `check_grid_columns.R` against the standing generator, and the labels keep
  Iteration 004's format.
- **Method set:** the 12 names match the user's list exactly, and all are in the
  standard set and the method registry. The worker parses, and the submitter
  passes `bash -n`.
- **Smoke test:** `run_simulation()` with the worker's arguments at n = 1000 and
  regions of 100 and 400 variants, using the bundled VCFs:

  | LD | LD given to methods | mean (LD − LD_true)² | NaN |
  |---|---|---|---|
  | in-sample | identical to LD_true | 0 | 0 |
  | panel 500 | cor(X_ref), 500 × p | 3.4e-3 | 0 |
  | panel 2000 | cor(X_ref), 2000 × p | 1.6e-3 | 0 |

- **Rank:** LD matrices are singular through perfectly correlated variant pairs,
  not through sample size, and at the same rate as in Iteration 004:

  | case | eigenvalues < 1e-8 |
  |---|---|
  | p = 400: in-sample, panel 500, panel 2000 | 15–17% |
  | Iteration 004's scale (n = 5000, p = 500), in-sample | 6–15% |

## 5. Running it

No `R/` file changed since Iteration 006's install (5590b95), so no reinstall is
needed. `check_pkg_current.R` stops the submitter if the installed package is
older anyway.

```bash
cd ~/fine-mapping-benchmark && git pull
module load R/4.5.2-gfbf-2025b GSL/2.8-GCC-14.3.0
bash scripts/hpc/check_toolchain.sh        # must print RESULT: PASS

FMB_GRID_GENERATOR=$PWD/scripts/hpc/generate_params_grid_iter007.R \
FMB_ITER007_METHODS=1 FMB_SCENARIOS_PER_TASK=5 \
PBS_QUEUE=v1_small72a PBS_WALLTIME=72:00:00 \
FMB_SCRATCH=$EPHEMERAL/fmbench_iter007 \
  bash scripts/hpc/submit_benchmark_pbs.sh
```

Check the echoed design before the array queues. It must show:

```
DESIGN: n=1000   regions={100,100,150,150,200,200,300,300,400,400}
        S={1,3,5}   phi={0.1,0.4}
        LD: in-sample, ref1500, ref2000, ref500, ref750
        p_causal: 0.6   annotations: binary, continuous
Grid: 20 rows x 60 scenarios; chunk=5 scenario(s)/task
      -> 12 tasks/row x 20 rows = 240 array tasks
```

The job logs print `[iter007] reduced method set: 12 methods: ...`.

**Cost.** In Iteration 004 these twelve methods took 1.5–1.8 CPU-h per scenario, at
regions five times larger and n five times larger. The expectation is well below
that here, so a 5-scenario task should finish far inside 72 h and the run should
cost under about 2,000 CPU-h. A mismatched LD can slow some methods' convergence,
so check the runtimes of the first tasks to finish.

## 6. Afterwards — before ephemeral deletes anything

Results land on `$EPHEMERAL`, which deletes files after 30 days. The completion
checker and the analysis are written when the jobs finish. The analysis must
stratify by LD level as well as by model × annotation type.
