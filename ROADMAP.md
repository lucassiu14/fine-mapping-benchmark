# Roadmap

Outstanding gaps in the fine-mapping-benchmark package, identified during the
May 2026 audit. Items are grouped by reviewer-importance for a publication and
roughly ordered by impact-per-effort within each tier.

**Publication framing (clarified May 2026):** the eventual paper is a
**fine-mapping benchmarking paper**, not a method paper for Functional
BEATRICE. This shifts priorities toward **comprehensiveness** (more methods,
more experimental axes, explicit external cross-checks) rather than
focused-contrast. Some items in the "Deferred entirely" list at the bottom of
this document deserve a second look under that framing — see the
"Re-prioritisation note for benchmark-paper framing" section before
implementing anything in Phase 4 or 5.

---

## Tier 1 — Reviewer-defensible benchmark

These are the gaps a referee or grant panel would most likely raise.

### 1. Baselines

Currently no naive comparator, so reported AUPRC numbers aren't anchored.

**Scoped decision (May 2026):** add a single `marginal_z` baseline. Random PIPs
and a top-pvalue / lead-SNP baseline were considered and deferred — they're
sanity-check floors, not informative comparators. **For a benchmark paper,
random / top-pvalue baselines deserve a second look** — they're cheap to add
and benchmark papers conventionally include them.

**Design:**

- Method name: `"marginal_z"`
- Normalisation: `pip_j = |z_j| / sum_k(|z_k|)` (sums to 1; single-causal-style).
  Avoids redundancy with ABF, which uses `z²/sum(z²)` with shrinkage.
- Credible set: greedy by PIP until cumulative PIP ≥ `coverage`, identical to
  ABF / PAINTOR. One CS returned per region.
- No model fitting; runtime_seconds ≈ 0.
- `params` records only `coverage` (default 0.95).

**Implementation checklist (estimated 1 hour):**

1. New file `R/wrappers/marginal_z.R` with `run_marginal_z(z, coverage = 0.95)`
   and `run_marginal_z_region(region_geno, region_pheno, ...)`. Pattern matches
   `R/wrappers/abf.R` exactly (input: z; output: standard list).
2. Register in `R/run_methods.R` `.FM_REGISTRY`:
   `marginal_z = "run_marginal_z_region"`.
3. Add a colour to `.FM_COLORS` in `R/plot_results.R`.
4. Add a row to the methods table in `README.md`.
5. Add a section to `docs/methods.md` describing the baseline and its purpose.
6. Add to `scripts/test_evaluate.R` and / or `scripts/test_comprehensive.R`.
7. Add to the HPC `METHODS` vector in `scripts/hpc/run_benchmark_job.R` and
   `scripts/hpc/smoke_test.R`.

**Acceptance: it shows up in PR curves and PIP-calibration plots, gets an AUPRC
column in `evaluation_summary.csv`, and no fine-mapping method should be
appreciably worse than it on any benchmark setting.** If a method does worse,
either the method is mis-tuned or the metric is misreporting.

### 2. LD-mismatch robustness

Every method currently gets `cor(X)` as LD — i.e. in-sample LD. Real GWAS uses
external reference panels (gnomAD, UK Biobank, 1000G). This is the open problem
behind CARMA's outlier detection and motivates much of the recent literature.

**Scoped decision (May 2026):**

- **Axis: sample-size mismatch only** (Axis A). Same population, smaller ref
  panel. Ancestry mismatch (Axis B) deferred — same machinery, follow-up
  data prep.
- **Strategy: independent draw.** Run `simulate_genotypes()` twice per region:
  once for the GWAS sample (size `n`), once for the reference panel (size
  `n_ref`). Phenotype + z-scores use the GWAS sample; LD passed to methods
  uses the reference panel. ~2× genotype generation time but
  methodologically clean.
- **Backwards compatible.** `n_ref = NULL` keeps the current in-sample
  behaviour exactly.
- **Default experimental grid:** `n_ref / n_gwas ∈ {0.1, 0.25, 0.5, 1.0}`.
  Proportional sweep scales with whatever `n_gwas` is set.

**Interface:**

```r
run_simulation(
  ...                       # existing args unchanged
  n             = 5000,     # GWAS sample size
  n_ref         = NULL,     # NEW: ref panel size. NULL = in-sample (current)
  ref_vcf_dir   = NULL,     # NEW: ref panel VCFs. NULL = same VCFs as GWAS
  ref_seed_offset = 1L      # NEW: offset for independent ref-panel draw
)
```

Same arguments added to `simulate_gwfm_data()`.

**Output structure changes (per region):**

- `genotypes[[i]]$LD` — what methods receive (ref-panel-derived when n_ref set,
  else in-sample as today)
- `genotypes[[i]]$LD_true` — NEW: always the in-sample LD, for diagnostics
- `genotypes[[i]]$n_ref` — NEW: ref panel size used (NA = in-sample)

`LD_true` lets evaluation compute LD-mismatch diagnostics (e.g.
`mean((LD - LD_true)^2)`) at analysis time without re-simulating.

**Implementation checklist:**

1. Extend `simulate_genotypes()` to accept and propagate an optional ref-panel
   draw. Re-uses `simulate_single_region()` internally.
2. Update `run_simulation()` and `simulate_gwfm_data()` signatures + plumbing.
3. Add `LD_true` storage in both pipelines.
4. Update test scripts to cover the new code path (`n_ref < n_gwas` case).
5. Add `n_ref` as a recognised stratification variable in `evaluate.R`
   (alongside `S`, `phi`, `p_causal`) so by-mismatch plots come for free.
6. Document in `docs/methods.md` and the per-pipeline docs.
7. Add a sweep config for the HPC grid (likely a new param column in
   `params_grid.csv`).

**Expected findings (to motivate the experiment):**

- SuSiE / SuSiE-inf should degrade monotonically as `n_ref` drops.
- CARMA should degrade less due to outlier detection.
- ABF should be relatively insensitive (no LD-inversion).
- BEATRICE / Functional BEATRICE plausibly handle mismatch well via the
  regularised variational prior — this is a methodological observation worth
  testing.

**Follow-up (Axis B, deferred):** ancestry mismatch is the same code path,
different inputs. Preprocess 1000G VCFs by super-population (AFR / AMR / EAS /
EUR / SAS) and pass via `ref_vcf_dir`. Add only after Axis A lands.

### 3. Multi-ancestry simulation

**Scoped decision (May 2026): deferred entirely** for the initial paper.
**Reconsider for benchmark-paper framing** — a benchmarking paper without
multi-ancestry in 2026 is incomplete and reviewers will likely demand it.
Currently parked but flagged as the most important deferred item to revisit.

When picked up later, the design direction is fixed:

- **VCF strategy: filter on-the-fly** inside `simulate_single_region()` via a
  `populations` argument (e.g. `populations = "EUR"` or
  `populations = c("AFR", "AMR")`). Lazier than pre-filtering, no disk
  pre-staging.
- **Reuses LD-mismatch machinery** from #2 — ancestry mismatch is Axis B of
  the same `n_ref / ref_vcf_dir` parameter framework.
- **Three experiments to design when revisited:**
  - Q1: within-ancestry power (does EUR-only vs AFR-only matter?)
  - Q2: cross-ancestry LD transfer (EUR causals, AFR ref panel)
  - Q3: trans-ethnic meta-fine-mapping (would need SuSiEx, MsCAVIAR wrappers)

Until revisited: the simulator continues to draw from "all populations"
unfiltered, which is what every existing run uses.

### 4. Binary / case-control traits

**Scoped decision (May 2026): deferred entirely** for the initial pass.
**Reconsider for benchmark-paper framing** — most large GWAS are
case-control; a benchmark that only covers continuous traits is hard to
defend. This is the second-most-important deferred item to revisit.

When picked up later, the design direction is fixed:

- **Liability-threshold model:** simulate continuous liability (as today),
  threshold at `(1 - K)`-quantile where K is the prevalence parameter, then
  sample cases and controls.
- **New parameters:** `trait_type = c("continuous", "binary")`,
  `prevalence` (default 0.10), `case_control_ratio` (default 1:1).
- **Z-score computation:** glm() per variant — slow but matches real GWAS;
  acceptable since summary stats are computed once per simulation.
- **Heritability scale:** report both observed-scale and liability-scale
  PVE so downstream comparisons are unambiguous.

Until revisited: only continuous traits.

### 5. More realistic genetic architecture

**Scoped decision (May 2026): deferred entirely** for the initial pass. Keep
just `"normal"` and `"equal"`. Methodologically nice-to-have but not
essential — major fine-mapping methods aren't very sensitive to the exact
tail distribution. **For a benchmark paper, adding Laplace is cheap (~30
lines) and pre-empts the obvious reviewer question.** Worth reconsidering
alongside the other deferred items.

When picked up later, the obvious additions are Laplace (~30 lines) and
mixture-of-normals (~60 lines). Power-law and empirical distributions
are lower priority.

---

## Tier 2 — Methodological additions

### 6. Missing methods

**Scoped decision (May 2026):** PolyFun is the priority. Ship alongside the
baselines (i.e. before LD-mismatch), as the foundational annotation-aware
comparator — PolyFun is the most-cited method in this category and every
serious benchmark includes it.

#### 6a. PolyFun-style methods — designed

Canonical PolyFun depends on pre-computed UK Biobank annotation priors
(~25 GB) and S-LDSC infrastructure calibrated on UKB. Neither applies to
simulated data with our 3-column synthetic annotations. We implement **two
defensible variants** instead, clearly labelled:

**`polyfun_oracle` — methodological upper bound (~50 lines)**

- Reconstructs the true per-SNP causal probability
  `π_j ∝ exp(A_j' log γ)` from `region_geno$annotations_matrix` and
  `region_pheno$truth$enrichment` (both already stored by the simulator).
- Feeds those as `prior_weights` to `susieR::susie_rss(...)`.
- Documented as "PolyFun with perfect annotation priors — ceiling for any
  annotation-aware method."

**`polyfun_est` — fair comparator (~150 lines)**

- Estimates per-annotation contributions `τ_k` by linear regression of
  `z_j² - 1` on `A_{j,·}` (LDSC-lite; no external reference LD scores
  needed since we're working per-region).
- Computes per-SNP prior `σ²(j) ∝ max(Σ_k τ_k · A_{j,k}, ε)`, normalises,
  feeds to `susieR::susie_rss` as `prior_weights`.
- Uses the **same** annotation matrix as the simulator generated — apples
  to apples with Funmap / Functional BEATRICE.
- Pure R, no external dependencies.

**Implementation checklist:**

1. New file `R/wrappers/polyfun_oracle.R`:
   - `run_polyfun_oracle(z, LD, n, annotations, enrichment, ...)` — explicit form
   - `run_polyfun_oracle_region(region_geno, region_pheno, ...)` — pulls
     `region_geno$annotations_matrix` and
     `region_pheno$truth$enrichment`
   - Calls `susieR::susie_rss(z, R, n, prior_weights = pi_vec, ...)`
2. New file `R/wrappers/polyfun_est.R`:
   - `run_polyfun_est(z, LD, n, annotations, ...)`
   - Internal `.estimate_per_snp_priors(z, A)` helper
   - Calls `susieR::susie_rss(..., prior_weights = sigma2_vec)`
3. Register both in `R/run_methods.R` `.FM_REGISTRY`.
4. Add colours to `.FM_COLORS` in `R/plot_results.R` (distinct from
   funmap / paintor / functional_beatrice).
5. Methods table in `README.md` and `docs/methods.md`.
6. Test scripts include both wrappers.
7. Add to HPC `METHODS` vector.

**Acceptance:**
- `polyfun_oracle` should be the **best or equal-best** annotation-aware
  method on any setting where the annotation signal is non-trivial. If
  not, the prior-injection plumbing into SuSiE is broken.
- `polyfun_est` should beat plain SuSiE on annotation-positive settings
  and roughly match plain SuSiE on annotation-null settings.
- Both should be appreciably worse than `polyfun_oracle` (the gap is the
  "cost of estimation").

**Resulting comparator landscape:**

| Method | Annotation prior source |
|---|---|
| polyfun_oracle | Truth (oracle ceiling) |
| polyfun_est | LDSC-lite regression of z² on A |
| Funmap | Joint random-effects fit |
| PAINTOR | EM enrichment weights |
| Functional BEATRICE | LassoNet prior network |

Five annotation-aware methods with distinct prior-estimation strategies.

#### 6b. Other comparators

**Scoped decision (May 2026):**

- **SparsePro added in Phase 4.** Tier 3 (Python) wrapper mirroring
  `R/wrappers/beatrice.R` — SparsePro is **CLI-only** (the upstream repo
  ships `sparsepro_zld.py` as a script, not a pip-installable module), so
  the wrapper uses `system2()` to invoke the script rather than
  `reticulate::import()`. Pure Python (numpy + scipy + pandas, no
  PyTorch), so the dependency footprint adds little on top of the
  existing BEATRICE/Funmap Python env. Adds a modern variational method
  to the comparator set, rounding out the catalogue.
- **DAP-G deferred** — eQTL-specific; reconsider for benchmark-paper
  framing if eQTL fine-mapping is in scope.
- **SBayesRC deferred** (Zheng 2024) — Bayesian fine-mapping with
  functional annotation categories via a BayesR-family prior plus
  per-annotation effect-size shrinkage. Distributed as part of GCTB
  (`gctb --sbayes-rc`); native binary, no R interface. A wrapper would
  follow the FINEMAP / PAINTOR pattern (system2 call to the binary).
  Worth adding once SparsePro is in, as another modern
  annotation-aware comparator distinct from the PolyFun / Funmap
  family.
- **MsCAVIAR / SuSiEx deferred** — couples with #3 multi-ancestry, also
  deferred.
- **eCAVIAR deferred** — colocalisation, not pure fine-mapping; out of
  scope for the first paper.

**SparsePro implementation:**

- `R/wrappers/sparsepro.R` with `setup_sparsepro()`, `run_sparsepro()`,
  `run_sparsepro_region()`. Pattern matches `R/wrappers/beatrice.R`
  (CLI invocation via `system2()`).
- Setup: user clones `https://github.com/zhwm/SparsePro`,
  `pip install -r requirements.txt` into the existing conda env, then
  passes `sparsepro_dir = "/path/to/SparsePro"` (and optionally
  `python = ...`) via `method_args`. No auto-clone, mirrors the
  BEATRICE pattern.
- Per-region call writes three temp files (zscore table, LD matrix,
  --zld summary file), invokes `python sparsepro_zld.py ...`, then
  parses the `.pip` (variant-level PIPs) and `.cs` (credible sets)
  output files.
- Register in `.FM_REGISTRY`, add to `.FM_COLORS`, etc.
- Effort: ~1 day.

### 7. MAF-stratified evaluation

**Scoped decision (May 2026):**

- **Bins (3):** `(0, 0.01]` rare, `(0.01, 0.05]` low-frequency, `(0.05, 0.50]`
  common — standard human-genetics partition.
- **Metric stratified:** AUPRC by causal-variant MAF only. (PIP calibration
  and CS size deferred — fewer cells per bin, less reliable.)
- **No new simulation code:** `genotypes[[i]]$maf` already exists.
- Add as `by_causal_maf` alongside the existing `by_S` / `by_phi` /
  `by_p_causal` axes in `evaluate.R`.

**Implementation:**

- In `.annotate_fits_with_truth()`, attach
  `f$causal_maf = region_geno$maf[truth$causal_indices]` to each fit.
- In `.stratify_metrics()` or a new sibling helper, bin by the **first**
  causal-variant MAF (or by the *minimum* — pick one and document; minimum
  is more conservative).
- Decision: bin by **minimum** causal MAF per region — captures "the
  rarest causal is the bottleneck for this region's fine-mapping".
- Output: `eval_out[[method]]$by_causal_maf` with keys `"rare"`, `"low"`,
  `"common"`.
- Add a corresponding `.print_section()` call in `plot_results()`.

**Acceptance:** for any method, AUPRC in the `rare` bin should be lower
than in the `common` bin under realistic phi. If they're equal, something
is wrong with the binning.

### 8. Model-misspecification evaluation

**Scoped decision (May 2026):** focus on the **annotation-misspecification**
case: test all annotation-aware methods (Funmap, PAINTOR,
functional_beatrice, polyfun_est) on `annotations = "none"` simulations.
This is the most realistic misspecification — practitioners often run
annotation-aware methods without strong functional priors.

SuSiE/SuSiE-inf cross-pairings (the classic 2×2) are deferred — they're
useful but the annotation-null case is more story-rich for the Functional
BEATRICE paper.

**No new simulation code, no new wrappers.** Pure experiment design + an
extra stratification axis.

**Implementation:**

- Expand the HPC params grid: add `annot_name = "none"` rows for both
  `sparse` and `sparse_inf` models (already partly present — confirm).
- Add `by_true_annotation_type` stratification in `evaluate.R`: keys
  `"none"`, `"binary"`, `"continuous"`. Pulled from
  `simulation$params$annotation_type`.
- Annotation-aware methods should remain runnable on `annotations = "none"`
  simulations (they should just run with uniform priors or refuse
  gracefully — Funmap currently refuses; PAINTOR / Functional BEATRICE
  fall back to uniform prior; polyfun_est should fall back to no
  prior boost). Audit each wrapper for this behaviour before running.

**Acceptance:** annotation-aware methods should match plain SuSiE
performance (not exceed it) when annotations are null. If they exceed it,
they're using information from elsewhere (suspicious). If they're worse,
they're over-fitting noise (informative — that's the misspecification cost).

### 7 + 8 timing

Bundle both with #2 LD mismatch as a single Phase 2 PR. All three are
analysis-side additions to `evaluate.R` + a few HPC grid changes. One
coherent update.

### 9. Anchor trait architectures to literature

`simulate_gwfm_data` sweeps π × h² freely. Anchor to known trait-level values:
- height-like: π ≈ 1e-3, h² ≈ 0.5
- T2D-like: π ≈ 5e-4, h² ≈ 0.4
- schizophrenia-like: π ≈ 1e-3, h² ≈ 0.3

Add to docs / smoke test, with citations.

---

## Tier 3 — Engineering / distribution

### 10. Convert to a proper R package

**Scoped decision (May 2026):**

- **Target: GitHub-installable only**, not CRAN. Most genomics packages live
  on GitHub or Bioconductor; CRAN's restrictions on `download.file()` would
  break the existing `setup_finemap()` auto-download pattern. Users will
  install via `remotes::install_github("lucassiu14/fine-mapping-benchmark")`.
- **Bioconductor not pursued** — useful in principle but the 6-monthly release
  schedule is too rigid for a research-stage package.
- **Tests: convert to testthat.** R-ecosystem standard, plays with
  `devtools::check()` and CI.

**Implementation steps:**

1. `usethis::create_package(".")` — creates `DESCRIPTION`, `NAMESPACE`,
   `.Rbuildignore`. Non-destructive on existing code.
2. Move `scripts/test_*.R` → `tests/testthat/test-*.R`, rewriting test
   bodies as `test_that()` blocks. The current custom assertion macros
   map cleanly to `expect_*()` functions.
3. Move `scripts/prepare_vcfs.R`, `prepare_gwfm_vcfs.R`,
   `download_ldetect_regions.R`, and `test_pipeline.R` → `inst/scripts/`
   (ships with install but stays out of the package namespace).
4. Keep `scripts/hpc/` at repo root — these are orchestration tooling
   that lives next to the package, not inside it.
5. Run `devtools::document()` to generate `man/*.Rd` from existing
   roxygen comments (most wrappers already have them).
6. Fill in `DESCRIPTION`:
   - Title, Description, License (MIT — see #11), Authors@R
   - Imports: susieR, sim1000G, ggplot2, gridExtra, grid, stats, utils
   - Suggests: CARMA (GitHub), reticulate (Python optional), testthat
   - SystemRequirements: htslib (tabix, bgzip)
7. Set `@export` tags correctly across all wrapper files and pipeline
   functions (already mostly done).
8. Run `devtools::check()` and address NOTES / WARNINGS.
9. `pkgdown::build_site()` → GitHub Pages.

**Effort:** ~1–2 days.

**Acceptance:**
- `remotes::install_github("lucassiu14/fine-mapping-benchmark")` works
- `library(finemapbenchmark); ?run_simulation` shows docs
- `devtools::check()` is clean (or has only documented NOTEs)
- pkgdown site builds and deploys

### 11. License

**Scoped decision: MIT.** Permissive, what susieR / CARMA / ggplot2 use,
maximum adoption friendliness. None of our dependencies require copyleft.

Add `LICENSE` at repo root containing the MIT text + author + year, and
declare `License: MIT + file LICENSE` in `DESCRIPTION` once the package is
created.

### 12. CITATION

Standard R package mechanism: `inst/CITATION` with `citEntry()` calls.

**Plan:**

1. Add `inst/CITATION` once the package is created (#10). Initial entry
   cites the GitHub repo + version.
2. Tag a v0.1.0 release on GitHub once Phase 1 lands.
3. Archive each tagged release to Zenodo for a DOI. Add the DOI to
   `inst/CITATION` and the README badge.
4. When the benchmarking paper is on bioRxiv / accepted, replace the
   manual entry with a `bibentry("Article", ...)` for the paper.

Effort: ~30 min for the initial CITATION; DOI archival is automated once
Zenodo + GitHub are linked.

### 13. GitHub Actions CI

Run `scripts/test_evaluate.R` on every push / PR. `r-lib/actions/setup-renv`
makes it straightforward. Catch regressions automatically.

### 14. Container

**Scoped decision (May 2026):**

- **Contents: deps-only (~2-3 GB).** R + susieR + CARMA + ggplot2 + sim1000G +
  reticulate + conda env (PyTorch, numpy, scipy, pandas, funmap) + FINEMAP
  binary + PAINTOR binary + htslib. **No reference VCFs or pre-baked
  simulations** — users mount their work and the simulator downloads
  what it needs.
- **Format: Docker source, Apptainer artefact.** Write a Dockerfile (single
  source of truth), build the Docker image for development, then convert
  to Apptainer (`apptainer build image.sif docker://...`) for HPC use.
- **Timing: before paper submission**, not now. Container is the
  reproducibility artefact that ships with the paper. Doing it during
  active wrapper development causes too many rebuilds.

**Implementation:**

1. `Dockerfile` at repo root, multi-stage:
   - Base: `rocker/r-ver:4.4.0` (or matching version)
   - Install system deps: htslib, conda
   - Install R packages via renv::restore() against the lockfile
   - Install PAINTOR (`conda install -c bioconda paintor`)
   - Auto-download FINEMAP at build time
   - Set up the BEATRICE / Funmap conda env
   - Clone Beatrice-Finemapping into a fixed location
2. `apptainer.def` *or* convert from Docker at build time
3. GitHub Actions workflow: build + push to GitHub Container Registry on
   every tagged release
4. README: add a "Container" section with both `docker run` and
   `apptainer exec` invocations

**Acceptance:** a reviewer can `apptainer pull ghcr.io/lucassiu14/finemapbenchmark:v0.1.0` and immediately run
`apptainer exec image.sif Rscript scripts/hpc/smoke_test.R` and have it pass.

### 15. Reference results / regression test fixture

A tiny canned simulation + expected evaluation output stored in `tests/fixtures/`.
Anyone can verify their installation produces bit-identical results.

---

## Tier 4 — Documentation / reporting

### 16. Worked vignette

**Scoped decision (May 2026):**

- **Format: R Markdown (`.Rmd`).** Stable, widely known, plays with pkgdown.
- **Scope: minimal getting-started.** Single file that walks through:
  small per-locus simulation → run two methods → evaluate → plot. Mirrors
  what `scripts/test_pipeline.R` does today but with prose commentary and
  inline output.
- **Location: `vignettes/getting-started.Rmd`** (after the R package
  conversion in #10).
- **Effort: ~4–6 hours.**

Don't try to demonstrate every feature — that's what `docs/methods.md` and
`docs/gw_simulation_documentation.md` are for. The vignette's job is "get a
new user productive in 10 minutes."

### 17. Method-selection guide

**Scoped decision (May 2026):**

- **Format: empirical table only.** Rows = simulation scenarios (e.g.
  "annotations + LD mismatch + S>1"), columns = recommended methods with
  brief justification from the benchmark numbers.
- **No methodological flowchart** — that's deferred. Empirical findings
  are stronger than theory-based recommendations and require less
  hedging.
- **Cannot be written until Phases 1–4 produce results.** This is a
  post-experiment deliverable.

When the time comes, write it as a section in the eventual paper and a
table in `docs/method_selection.md` so both audiences can reference it.

### 18. Cross-check against published benchmarks

**Scoped decision (May 2026): deferred entirely.** Useful for credibility but
not blocking. Park; revisit if reviewers raise the question or if a
clear-cut comparison opportunity appears (e.g. a published benchmark
re-using the same simulation framework we use).

### 19. Restructure README into install tiers

**Scoped decision (May 2026):** the current README presents installation as one
intimidating block. Restructure into three explicit tiers so users can see how
far they can get with minimal install. Costs ~30 minutes; pays off immediately.

**The tiers:**

- **Tier 1 — pure R (~15 min, no external tools):** `renv::restore()` + GitHub
  PAT. Covers SuSiE, SuSiE-inf, ABF, CARMA, marginal_z, polyfun_oracle,
  polyfun_est, plus all evaluation and plotting. **7 of 9 methods, all
  evaluation, all plotting.**
- **Tier 2 — adds binaries (~30 min extra):** FINEMAP (auto-downloads via
  `setup_finemap()`), PAINTOR (conda or compile from source).
- **Tier 3 — adds Python (~1 hour extra):** Conda env for BEATRICE, Functional
  BEATRICE, and Funmap. Clone Beatrice-Finemapping repo.

**Implementation:**

- Rewrite README "Installation" section with three clearly labelled subsections
- Add a "What works at each tier?" table mapping methods → tier
- Move the htslib / VCF setup into a separate "Reference data" section
  (it's not tied to a method tier — needed only for simulating from real LD)
- Add an "I just want to try it" path: clone repo → `renv::restore()` →
  `Rscript scripts/test_pipeline.R` with bundled VCF (Tier 1 only)

**Acceptance:** a new user can read the README, decide which tier they need,
and ignore the rest. Apple Silicon / cluster users can pick Tier 1 only
without ever touching PAINTOR or conda.

---

## Implementation order (decided May 2026)

**Phase 0 — Free wins (do anytime, no dependencies):**
0a. **#19 Restructure README into install tiers** — ~30 min; makes the package
    look 10× less intimidating to new users.

**Phase 1 — Foundational comparators (ship together):**
1. **#1 marginal_z baseline** — anchors AUPRC numbers
2. **#6a polyfun_oracle + polyfun_est** — most-cited annotation-aware
   method; benchmark papers without it are incomplete. Both land in Tier
   1 (pure R) so no install pain added.

These three wrappers share a code pattern and should be a single PR.

**Phase 2 — Experimental axes:**

3. **#2 LD mismatch** — sample-size mismatch (Axis A). Independent ref-panel
   draw; `n_ref / n_gwas ∈ {0.1, 0.25, 0.5, 1.0}`.
4. **#7 + #8 MAF stratification + misspecification** — analysis-only;
   no new simulation code, just new stratification axes in evaluate.R.

**Phase 3 — Engineering / distribution:**

5. **#10–12** R package conversion + LICENSE + CITATION.
6. **#14** Apptainer container.

**Phase 4 — Comparator additions:**

7. **#6b SparsePro wrapper** — Tier 3 (Python); strengthens the paper against
   the "you only beat old methods" critique. ~1 day.

**Phase 5 — Reporting (post-experiment):**

8. **#16 Vignette** — R Markdown `vignettes/getting-started.Rmd`, minimal
   scope. ~4–6 hours.
9. **#17 Method-selection table** — empirical findings only, written after
   Phases 1–4 produce real numbers.

**Deferred entirely (parked, not in current plan):**

- **#3 Multi-ancestry** — too large for the first paper; reuses #2's
  machinery when picked up.
- **#4 Binary / case-control traits** — liability-threshold model is the
  intended approach when picked up.
- **#5 Realistic effect-size distributions** — keep `normal` and `equal`;
  document as a known limitation.
- **#6b DAP-G / MsCAVIAR / SuSiEx / eCAVIAR** — audience-specific or coupled
  to deferred features.
- **#18 Cross-check against published benchmarks** — useful for credibility
  but not blocking; revisit if reviewers ask.

---

## Re-prioritisation note for benchmark-paper framing

The deferral decisions above were made under an implicit method-paper framing
("ship the minimum needed to demonstrate Functional BEATRICE"). The actual
paper is a **fine-mapping benchmark paper**, which calls for
**comprehensiveness over focused contrast**. Before starting Phase 4, this
deferral list should be revisited. Under benchmark-paper priorities:

| Item | Method-paper priority | Benchmark-paper priority |
|---|---|---|
| #3 Multi-ancestry | Deferred ("too big") | **Likely required** — 2026 benchmarks without it are hard to defend |
| #4 Binary traits | Deferred ("not blocking") | **Likely required** — half of GWAS is case-control |
| #5 Effect distributions | Deferred ("not essential") | **Worth adding Laplace** (~30 lines) — pre-empts reviewer questions |
| #6b DAP-G | Deferred ("audience-specific") | **Reconsider** — depends on whether eQTL fine-mapping is in scope |
| #18 Cross-check | Deferred ("not blocking") | **Should be done** — benchmark papers need external validation |

#17 Method-selection guide also gets promoted under benchmark-paper framing:
it's the headline deliverable of a benchmark paper, not an optional extra.

**Decision: hold this re-prioritisation as a checkpoint** — after Phase 3 lands
and before committing to Phase 4 work, revisit which deferred items to
unblock based on time / writing energy available.

**Not yet sequenced (but designed):**

- **#13 GitHub Actions CI** — fits naturally with Phase 3 (alongside R
  package conversion). Run `tests/testthat/` on push/PR.
- **#15 Reference results / regression test fixture** — Phase 3-ish, small
  test fixture as part of the testthat conversion.
