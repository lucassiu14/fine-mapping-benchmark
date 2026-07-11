# Iteration 001 — Phase 1 HPC grid

**Date:** 2026-07-05
**Scope:** Phase 1 setup (§1.1, §1.3) — extend the SLURM grid + worker to
cover the plan's Phase 1 axes without launching the array yet.
**Outcome:** Grid, worker, and smoke tests all pass locally. Ready to
submit on the cluster once the operator confirms VCF availability.

## What landed

| # | File | Change |
|---|---|---|
| 1 | `scripts/hpc/generate_params_grid.R` | Rewritten for §1.1 / §1.3 axes |
| 2 | `scripts/hpc/run_benchmark_job.R` | Consumes new columns (`p_values`, `enrichment_values`, `annotation_correlation`, `p_causal`) and registers `polyfun_ldsc` + `sbayesrc` |
| 3 | `scripts/hpc/smoke_test.R` | Runs all 9 Tier-1 methods incl. new ones (touched last iteration; verified this iteration) |
| 4 | `scripts/hpc/params_grid.csv` | Regenerated: 25 tasks × 125 scenarios = 3125 scenarios (matches §1.3) |

## The Phase 1 grid (§1.1)

Fixed constants baked into every row:
- Per-region `p` = `c(100,100,100,100,200,200,200,200,400,400,400,400,500,500,500,500,1000,1000,1000,1000)` (length 20)
- `n_regions` = 20 (matches length of `p`)
- `n_annotations` = 20
- `enrichment_values` = `c(7.4, 7.4, 2.7, 2.7, rep(1, 16))` — two enriched
  pairs plus 16 null decoys (§1.1)
- Within-job sweep: `S ∈ {1,2,3,5,10}`, `phi ∈ {0.0075, 0.05, 0.1, 0.2, 0.4}`,
  `n_iter = 5`, `n = 1000` → 125 scenarios per task

Swept across tasks:
- `model ∈ {sparse, sparse_inf}`
- `annotation_regime`: `{none} ∪ {binary × ac ∈ {0, 0.25, 0.5, 0.75}}` (§1.3)
- `p_causal ∈ {0.5, 0.7, 0.9, 1.0}` for `sparse_inf` only

Row counts: `sparse × 5 + sparse_inf × 5 × 4 = 5 + 20 = 25` tasks.
Method-fits at full array: `25 × 125 × 20 regions × 9 methods = 562,500`.

## Worker changes (`run_benchmark_job.R`)

1. Row-parsing block now decodes four pipe-separated fields:
   `S_values`, `phi_values`, `p_values`, `enrichment_values`.
2. Echo block prints per-region `p`, per-annotation `enrichment`, and
   `annotation_correlation` (formatted `NA` for the `none` regime).
3. `p_causal` is only forwarded when the model is `sparse_inf` **and**
   the CSV cell is non-`NA` — `run_simulation()` uses its own default
   for the `sparse` model.
4. `METHODS` gained `polyfun_ldsc` and `sbayesrc`; `METHOD_ARGS` gained
   entries for both (`sbayesrc = list(n_iter=300, burn_in=150,
   gamma_update_every=10)`).
5. Backward-compat: the legacy `n_ref` column is still forwarded if it
   appears in a row, so older grids keep working.

## Smoke verification (local, laptop)

Two staged tests were run before flagging Phase 1 code as
submission-ready:

### 1. Dry-run harness (`scratchpad/dryrun_worker.R`)

Sources `run_benchmark_job.R` inside an env that shadows
`run_simulation`/`run_methods`/`evaluate_methods` and just captures
the args they would receive. Sampled `job_id ∈ {1, 2, 6, 10, 25}` (the
structurally distinct rows: sparse-none, sparse-binary, sparse_inf-none,
sparse_inf-binary-highcorr, and the last row).

Assertions that passed for every sample:
- `length(p) == 20`, `n_regions == 20`
- `length(S) == 5`, `length(phi) == 5`
- `length(methods) == 9`; `polyfun_ldsc` and `sbayesrc` are both present
- `sparse_inf` rows carry a valid `p_causal ∈ (0, 1]`
- `binary` rows carry length-20 `enrichment` and a non-`NA`
  `annotation_correlation`
- `sparse` + `none` rows do **not** pass `p_causal` to `run_simulation`

### 2. Miniature end-to-end (`scratchpad/mini_end_to_end.R`)

Three shape-representative rows (sparse+none, sparse+binary+ac0.5,
sparse_inf+binary+ac0.25+pc0.7) run through the full
`run_simulation → run_methods → evaluate_methods` chain with tiny
sizes (`p = c(60,60,100,100)`, `n = 300`, `n_iter = 1`, single
scenario). All 9 methods completed 16/16 fits with 0 failures.

## Findings from the smoke tests (recorded, not blockers)

1. **VCF panel does not have 20 regions.** With
   `VCF_DIR = "data/vcf"` and the current 2-file panel,
   `run_simulation()` errors out at
   `vcf_dir contains 2 VCF file(s) but n_regions = 20`.
   Before submitting on HPC, either:
   - Populate `data/vcf/` with ≥ 20 VCFs and rerun
     `inst/scripts/prepare_vcfs.R`, **or**
   - Leave `VCF_DIR` set but ship the extra VCFs via the cluster's data
     staging step (see `scripts/hpc/README.md`), **or**
   - Fall back to synthetic genotypes by setting `VCF_DIR = NULL` in the
     worker (worse LD realism; only appropriate for a Phase 1 shakedown).

   Decision: defer to the operator (user) — the smoke test does not
   need real LD to exercise the code paths.

2. **SBayesRC AUPRC is much lower than the other methods** on the
   miniature run (0.073 vs. 0.25–0.52 for the rest). Likely causes:
   - `n_iter = 100`, `burn_in = 50` in the smoke test — much shorter
     than the HPC default (300/150).
   - The four-slab default `sigma2_scale = c(0.05, 0.005, 5e-4, 5e-5)`
     was tuned against SBayesRC's paper defaults and has never been
     evaluated against our simulator's effect distribution.
   - Only `n_annotations = 6` with a very small sample size (`n = 300`)
     — pooled annotation regression is data-starved.

   All three are Phase 2 tuning questions. Recorded here so the auto-
   research loop can flag SBayesRC calibration in the first iteration.

## Functional BEATRICE — reviewer feedback state

External review of `BEATRICE_annot_sparse/` flagged four items.
Status after the Phase 1 audit + fix pass:

1. **Annotation-drop in the genome-wide path** — already fixed in
   Phase 0 (PR #19). `wrapper_functional_beatrice.R` prefers
   `region_geno$annotations_matrix` and falls back to
   `region_pheno$annotations_matrix`; regression-tested.
2. **Variational objective is only an approximation** (Bernoulli KL
   only over selected top-K; `∑ p_{0j}^2` regularisation shrinks all
   prior probabilities). Legitimate theoretical critique of the
   BEATRICE design, not a bug we can quick-fix. Logged as a Phase 2
   novel-method candidate: "Cardinality-consistent BEATRICE" (proper
   full-KL + Beta-Binomial K, or SuSiE-like SER decomposition of the
   prior).
3. **LassoNet identifiability** — fixed in the Phase-1 numpy-2 branch:
   `feature_importance` now uses the identifiable logit contrast
   `|theta[:, 1] - theta[:, 0]|`, not `||theta_j||_2`. The full
   reparameterisation to a single Bernoulli logit is deferred as a
   follow-up; also flagged that the manual proximal step differs from
   the LassoNet paper's projected-proximal-gradient path — same
   category of concern, kept as a documented limitation for now.
4. **numpy-2.x crash in `calculate_pip` + missing `return` in
   `reformat_memo`** (found during this iteration's audit, not in the
   reviewer note but same code region) — fixed in the same PR. Vanilla
   `beatrice` now runs the FB fork with `--annot` omitted, giving
   BEATRICE semantics without the upstream late-training crash.

## What's next

The grid + worker are ready. Next steps in order:
1. Operator: stage ≥ 20 VCFs onto the HPC and confirm `data/vcf/`
   contains them (see finding 1 above).
2. Submit the array: `sbatch scripts/hpc/submit_benchmark.sh`
   with `--array=1-25`.
3. When it returns, aggregate results and start Iteration 002
   (Phase 2 — auto-research loop, calibration gate, SBayesRC tuning).

## Reproducibility

- Grid CSV: `scripts/hpc/params_grid.csv` (checksum below).
- Smoke scripts (throwaway, not in git): `scratchpad/dryrun_worker.R`,
  `scratchpad/mini_end_to_end.R`.
- Seeds: worker uses `1000 + job_id` so every row is reproducible.
