# UK Biobank cholesterol application

Genome-wide fine-mapping of UK Biobank total cholesterol over LDetect EUR
blocks, comparing **fb_xregion** (annotation selection) against **Funmap**
(annotation shrinkage) and **PolyFun + SuSiE**.

Everything below is reproducible from public data. No UK Biobank application
or credentials are needed: all three inputs are openly hosted.

---

## 1. Data

### 1.1 GWAS summary statistics

```
https://storage.googleapis.com/seqnn-share/sniff/_20260810_tmp_files/fixed_n_337k/sample_munged_gwas/biochemistry_Cholesterol.sumstats_munged.parquet   (526 MB)
```

UK Biobank total cholesterol, 15,221,393 variants on all 22 autosomes,
max |Z| = 82.9, 31,962 variants at p < 5e-8, median MAF 0.033.

**Provenance.** The `biochemistry_Cholesterol` naming, and the sibling
phenotypes in the same bucket (`body_HEIGHTz`, `bmd_HEEL_TSCOREz`,
`bp_DIASTOLICadjMEDz`, `pigment_HAIR`, `repro_MENARCHE_AGE`, ...), identify
this as the **Alkes Price group's 15-trait UK Biobank set** used in PolyFun
(Weissbrod et al. 2020) — the Loh et al. BOLT-LMM summary statistics computed
on UK Biobank **white-British unrelated** individuals.

**Ancestry consistency (checked).** The GWAS cohort (UKB white British), the
LD panel (§1.2, UKB white British unrelated) and the LDetect blocks (§1.3,
1000G EUR) are all European. GWAS and LD panel are the *same* cohort, so the
LD is in-sample rather than merely ancestry-matched.

**Sample size.** The `N` column reads **373,220** and that is what PolyFun
itself used on this data (§1.4). We use `n = 373220`. Note the LD panel holds
~337,000 individuals, so the panel is a subset of the GWAS sample; the
directory name `fixed_n_337k` refers to the panel, not the GWAS. The earlier
pilot used 337,000 — the difference is 10.7% in `n` and has only a small
effect on the Bayes factors, but 373,220 is the GWAS size and the value the
comparator used.

### 1.2 LD

```
s3://broad-alkesgroup-ukbb-ld/UKBB_LD/{window}.npz   (~0.94 GB each)
s3://broad-alkesgroup-ukbb-ld/UKBB_LD/{window}.gz    (variant table)
aws s3 cp --no-sign-request ...      # public, no credentials, no requester-pays
```

UK Biobank in-sample LD from the Alkes Price group
(<https://registry.opendata.aws/ukbb-ld/>), computed on ~337,000 unrelated
white-British individuals. Windows are `chr{c}_{k*1e6+1}_{(k+3)*1e6+1}` —
3 Mb wide, stepping 1 Mb, so consecutive windows overlap by 2 Mb. `.npz` is a
COO float32 **lower triangle**; the loader symmetrises it and sets the
diagonal to 1.

This is the same LD resource Funmap's own paper uses for its cholesterol
example, so all three methods share one LD definition.

**262 windows, ~246 GB.**

### 1.3 Loci

LDetect EUR approximately-independent LD blocks (Berisa & Pickrell 2016),
shipped with this package at `inst/extdata/gwfm_regions_ldetect_EUR.csv`
(1,701 blocks).

### 1.4 Annotations

```
https://storage.googleapis.com/seqnn-share/sniff/processed_borzoi_parquets/borzoi.{CHR}.annot.parquet   (0.2-1.3 GB each, ~15 GB total)
```

**Borzoi-102**: 102 principal components of Borzoi variant-effect predictions
over 6,385 datasets — ATAC 20, DNASE 20, RNA 20, CHIP 20, RNAGTEX 13, CAGE 9.
Public GCS, no authentication.

### 1.5 PolyFun priors (pre-computed by the real implementation)

```
.../sniff/sniff_results_ukb_15/sniff/biochemistry_Cholesterol.{CHR}.txt.gz      # Borzoi-102 priors
.../sniff/sniff_results_ukb_15/baseline/biochemistry_Cholesterol.{CHR}.txt.gz   # baselineLF priors
```

Columns: `CHR SNP BP A1 A2 SNPVAR MAF N Z P PIP BETA_MEAN BETA_SD CREDIBLE_SET`.

`SNPVAR` is PolyFun's **stage-1 output**: the per-SNP prior causal variance
from its S-LDSC + ridge model, computed by the authors' own code on exactly
this trait and these annotations. Using it lets us run the real
`finemapper.py` without re-running S-LDSC genome-wide for 102 annotations.
The two sets are genuinely different (corr(SNPVAR) = 0.247, and 3 vs 1 SNPs
at PIP > 0.5 on chr22), confirming `sniff/` is the Borzoi-based prior.

Note their `PIP`/`CREDIBLE_SET` columns fine-map over **3 Mb windows**, not
LDetect blocks, so they are not directly comparable to our loci — we use
`SNPVAR` only, and refit on our blocks.

### 1.6 Software (all upstream, unmodified)

| | |
|---|---|
| Functional BEATRICE | this project, branch **`fix/memo-normalisation`** or a descendant — **not `main`**, which lacks the `reformat_memo` max-offset fix and silently produces diluted PIPs. `03_fb_xregion.pbs` checks for it and refuses to run otherwise. |
| Funmap | <https://github.com/LeeHITsz/Funmap> |
| PolyFun | <https://github.com/omerwe/polyfun> |
| BEATRICE (if run) | <https://github.com/sayangsep/Beatrice-Finemapping> |

**Do not** use the copies under `Beatrice_annot/alt_methods/` — the Funmap
copy there has its input readers rewritten (originals commented out). The
algorithm is untouched, but comparators must be the authors' code.

---

## 2. Preprocessing, exactly as performed

### Step A — score every LDetect block (`00a_score_blocks.py`)
For each of the 1,701 LDetect EUR blocks, take `max(Z^2)` over the summary
statistics inside it and convert to a two-sided p. Writes `blocks_scored.csv`.

### Step B — select loci and assign LD windows (`00_make_blocks.py`)
Keep blocks with `min_p < 1e-6`: **265 blocks**.

For each, choose the 3 Mb UKBB window that (a) contains the block's **lead
variant** and (b) retains the most of the block, then trim the block to that
window.

> A block up to 2 Mb always fits some window; wider ones may not, because
> window boundaries sit on integer Mb. **37 of the 265 need trimming: 18.81 Mb
> total, 3.89% of the span, with zero lead variants lost.**
>
> An earlier version of this pipeline silently dropped all 58 blocks wider
> than 3 Mb. **13 of them pass p < 1e-6**, including `ldEUR_11_0034`
> (chr11:58.8–62.2 Mb, lead |Z| = 19.1 — the **FADS1/FADS2** cluster) at
> p = 2.4e-81 and `ldEUR_17_0039` (lead |Z| = 14.3, the **ABCA6/8/9/10**
> cluster) at p = 1.9e-46. Dropping them would have removed two of the
> strongest cholesterol loci in the genome. They are now trimmed like any
> other wide block and retained.
>
> The alternative — stitching two overlapping windows and zeroing pairs more
> than 2 Mb apart — was rejected: it keeps more variants but makes the LD
> inexact for exactly the long-range pairs a 3 Mb locus exists to model.

Writes `blocks.csv` (with `orig_start`, `orig_end`, `trimmed_bp`, `lead_bp`,
`lead_abs_z` per block) and `ld_windows.txt`.

### Step C — build per-block inputs (`prep_block.py`)
Per block:

1. Read the LD panel's variant table; restrict to the block and to
   **biallelic SNVs**; drop duplicate positions entirely.
2. Match to the summary statistics on **position + unordered allele pair**,
   and flip `Z` wherever the sumstats' A1 differs from the LD panel's A1, so
   `Z` is harmonised to the LD allele coding.
3. Join the Borzoi annotations on the **ordered** allele pair, as coded by the
   LD panel.
4. Take the LD submatrix **by the row index carried through every merge**,
   symmetrise, set the diagonal to 1, store float32.
5. Standardise annotations per block (mean 0, sd 1) and clip to ±10.

Each of those rules exists because of a bug that actually occurred:

- **Row index, never position.** Position is not unique. A position-keyed LD
  lookup hands every allele at a site the *same* LD row, producing r = 1 pairs
  carrying different z — logically impossible input, and the worst case for an
  RSS model. It passed every global diagnostic (`estimate_s_rss` s = 0.01–0.03)
  while capturing up to 100% of the posterior mass.
- **Ordered allele pair for annotations.** Borzoi lists some indels in both
  allele orders with different values; an unordered key silently duplicates
  rows.
- **Biallelic SNVs, no duplicate positions.** What Funmap's own pipeline does.
  Removes the problem at source, at a cost of ~10% of variants.

`prep_block.py` asserts no two retained variants share an LD row. Verified on
the pilot blocks: zero duplicate positions, and max |dZ| = 0.000 within every
genuine r = 1 group.

### Resulting scale
**265 blocks, 464 Mb, ~2.36M variants**, median p 8,997, max p 15,219.

---

## 3. Running

Three Python environments, because the pinned dependencies conflict:

| env | for | pins |
|---|---|---|
| `APP_VENV` | prep, fb_xregion | numpy, pandas, pyarrow, scipy, torch |
| `FUNMAP_VENV` | Funmap | numpy 1.26.2, scipy 1.11.4, pandas 2.1.3, matplotlib 3.8.2 |
| `POLYFUN_ENV` | PolyFun | see `polyfun.yml` upstream; needs rpy2 + susieR |

```bash
export APP_ROOT=/rds/.../chol_app
export N_LIKELIHOOD=373220

python 00a_score_blocks.py                              # -> blocks_scored.csv
python 00_make_blocks.py                                # -> blocks.csv, ld_windows.txt
qsub -v APP_ROOT 00b_download_annot.pbs                 # 22 tasks,  ~15 GB
qsub -v APP_ROOT 01_download_ld.pbs                     # 262 tasks, ~246 GB
qsub -v APP_ROOT,APP_VENV 02_prep.pbs                   # 265 tasks
qsub -v APP_ROOT,APP_VENV,FB_DIR 03_fb_xregion.pbs      # ONE job, 160 GB
qsub -v APP_ROOT,FUNMAP_VENV,FUNMAP_DIR 04_funmap.pbs   # 265 tasks
qsub -v APP_ROOT,POLYFUN_ENV,POLYFUN_DIR 05_polyfun.pbs # 265 tasks
python 06_collect.py
```

Steps 03–05 are independent and can run concurrently. Every step is
idempotent: a block with a `DONE` marker or an existing output is skipped, so
any array can be resubmitted after a partial failure. Update the array sizes
in the `.pbs` headers if you change the p-value threshold.

---

## 4. Resources

**Measured** (5 pilot blocks, 35,621 variants, 2000 joint iterations):
366 s wall with `.npy` LD against 493 s with text; peak RSS 3.25 GB of which
1.06 GB is the retained float32 LD.

**Projected** for 265 blocks: LD 94 GB + annotations ~1 GB + ~2 GB fixed
overhead. **Request 160 GB for `03_fb_xregion`**; mini-batching is not needed.

**Measured per-method, on one block (`ldEUR_20_0038`, chr20, p = 4,238,
8-core laptop):**

| method | block | p | wall time |
|---|---|---|---|
| fb_xregion | 5 blocks **jointly**, 2000 iterations | 35,621 total | 366 s |
| Funmap | ldEUR_22_0024_part | 7,777 | **702 s (11.7 min)** |
| Funmap | ldEUR_20_0038 | 4,238 | 3,107 s — *see below* |
| PolyFun `finemapper.py` | ldEUR_20_0038 | 4,238 | 46 s (6–8 s in SuSiE) |
| PolyFun `finemapper.py` | ldEUR_22_0024_part | 7,777 | 29 s |

**The 3,107 s figure is not a clean measurement** — that run spent most of its
life competing with a stray process for the CPU. The clean number is chr22:
**11.7 minutes at p = 7,777**, despite that block being larger *and* taking
more stage-2 iterations (19 vs 13). Funmap's cost tracks iteration count to
convergence as much as p, so it varies block to block.

Budget on the clean figure: roughly **15–25 min per block** at the median
p ≈ 9,000, more for the largest. A 12 h walltime on `04_funmap.pbs` is ample,
and it parallelises trivially as an array job.

**Still to measure: fb_xregion total runtime.** `regularize_ld` always falls
through to a full eigendecomposition on real float32 UKBB LD (min eigenvalue
≈ −0.029, so the Cholesky fast path never succeeds), which is O(p³) *per
block*. At p ≈ 9,000–15,000 over 265 blocks this dominates setup. Run `03`
against a 10-block manifest and extrapolate before trusting the walltime.

### Verified end to end

All three methods were run on `ldEUR_20_0038` with identical inputs:

**ldEUR_20_0038** (chr20:62.19–62.97 Mb, p = 4,238, max |Z| = 6.05)

| method | sum PIP | max PIP | top variant | \|Z\| | rank |
|---|---|---|---|---|---|
| fb_xregion | 1.72 | 0.728 | 62,896,110 | 5.98 | 27 |
| funmap | 9.99 | **1.000** | 62,370,732 | 6.00 | 13 |
| polyfun_susie | 9.04 | 0.407 | 62,370,732 | 6.00 | 13 |

**ldEUR_22_0024_part** (chr22:50.00–51.24 Mb, p = 7,777, max |Z| = 5.52)

| method | sum PIP | max PIP | top variant | \|Z\| | rank |
|---|---|---|---|---|---|
| fb_xregion | 1.00 | 0.033 | 50,858,813 | 5.46 | 2 |
| funmap | 9.99 | **1.000** | 50,858,813 | 5.46 | 2 |
| polyfun_susie | 9.97 | 0.133 | 50,838,241 | 5.42 | 11 |

Spearman between PIP vectors on chr20: funmap–polyfun 0.722,
fb_xregion–polyfun 0.533, fb_xregion–funmap 0.518.

Three things to look at before reading anything into the full run:

1. **Funmap returns PIP of exactly 1.000**, on the lead variant in both blocks
   and on two variants in chr20. At |Z| ≈ 6 inside a region of several thousand
   variants in LD, a posterior of exactly 1 is not credible; it suggests
   numerical saturation rather than certainty.
2. **fb_xregion is far more conservative** — sum(PIP) 1.0–1.7 against ~10 for
   the other two, and on chr22 its maximum is 0.033, i.e. it declines to
   localise at all where funmap is certain. Whether that is appropriate
   caution or a failure to commit is the thing to establish.
3. **Funmap and PolyFun both put sum(PIP) ≈ 9.99**, their L = 10 prior, in
   every block. That is SuSiE-like behaviour and means the posterior expected
   number of causal variants is essentially the prior.

A real 50 kb gap at chr22:50,364,571–50,414,983 is visible in the z-scores. It
is genuine: zero variants there in both the summary statistics and the LD
panel, not an artefact of the filtering.

### Building the PolyFun environment

`polyfun.yml` (python 3.8, `r-susier==0.11.92`) solves on x86_64 Linux, which
is what the cluster runs. It does **not** solve on osx-arm64. A working pip
environment needs four specific pins, each found the hard way:

| pin | why |
|---|---|
| `pyarrow<=17` | newer pyarrow requires numpy ≥ 2, which conflicts with the rest |
| `rpy2==3.5.x` | 3.6 turned `numpy2ri.activate()` into a raised `DeprecationWarning`, and finemapper.py calls it |
| `susieR 0.11.92` | 0.12+ moved `bhat/shat/R/var_y` out of `susie_suff_stat` into `susie_rss`; finemapper.py calls the old signature. Install with `remotes::install_github("stephenslab/susieR@v0.11.92")` — note the `v` prefix on the tag |
| `RPY2_CFFI_MODE=ABI` | only if the rpy2 wheel was built against a different R minor version than the one installed |

Set `R_LIBS` so the 0.11.92 library comes first, and `R_HOME` to the R
installation rpy2 should embed.

---

## 5. Open items

- **Pre-registered prediction.** `RNAGTEX-PC5` (0-based) separates liver at
  3.5 SD, so a lipid trait should select it. Check
  `results/annotation_weights.csv` against Funmap's shrunk weights.
- **Correctness anchor.** Funmap's paper fine-maps cholesterol at
  chr8:5–8 Mb using this same LD bucket. Reproducing that locus would validate
  the whole chain end to end.
