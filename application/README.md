# UK Biobank cholesterol application

Genome-wide fine-mapping of UK Biobank total cholesterol over LDetect EUR
blocks, comparing four methods on identical loci:

| method | annotations | how it uses them |
|---|---|---|
| `fb_xregion` | Borzoi-102 | LassoNet — sparse **selection** |
| `funmap` | Borzoi-102 | random effects — **shrinkage**, no selection |
| `polyfun_sniff` | Borzoi-102 | S-LDSC prior on the same 102 |
| `polyfun_combined` | Borzoi-102 + baselineLF (289) | S-LDSC prior, standard specification |

The two PolyFun variants separate the method from the annotation set:
`polyfun_sniff` sees exactly what the other two see, `polyfun_combined` is
PolyFun at its strongest.

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
Public GCS, no authentication. Used by all four methods.

```
.../sniff/processed_baselineLF_parquets/baselineLF2.2.UKB.{CHR}.annot.parquet   (~5 GB total)
```

**baselineLF v2.2**: 187 annotations, used only by `polyfun_combined`. It
shares the Borzoi files' variant set and row order exactly — verified on
chr22: identical SNP, BP and A1 sequences, 271,699 rows each — so `00b2`
combines them by column concatenation into 289 annotations, asserting the
alignment rather than assuming it. baselineLF has no tissue-specific
annotations.

### 1.5 PolyFun prior — computed here, not downloaded

The prior is computed from our own annotation set by `00d`–`00g`, using the
authors' `polyfun.py`. It is **not** downloaded.

Precomputed priors do exist in the same bucket
(`sniff_results_ukb_15/{sniff,baseline}/biochemistry_Cholesterol.{CHR}.txt.gz`)
and they are for the **same GWAS** — verified on chr22: N, Z and MAF agree to
within text rounding across 209,974 clean 1:1 matches on position+alleles,
with zero allele flips. But **which annotation set produced them cannot be
established**: the bucket holds only stage-2 intermediates (30 `snpvar_bin`
columns), not the annotation LD-scores, and a scalar `SNPVAR` cannot be
reverse-engineered. `00c` can still fetch them as a cross-check against our
own, which is the only honest use for them.

Reference data for computing it ourselves (0.30 GB, from the Zenodo record
"S-LDSC reference files", 10515792 — the old
`alkesgroup.broadinstitute.org/LDSCORE` URLs 404 after redirect):

```
1000G_Phase3_plinkfiles.tgz            0.29 GB  -> 1000G_EUR_Phase3_plink/
1000G_Phase3_weights_hm3_no_MHC.tgz    0.01 GB
```

**Two deliberate choices.**

*LD-scores from the 1000G EUR panel, not the UKB LD.* Computing them from the
UKB LD would mean all 2,763 regions, **3.06 TB** of transfer (mean 1.11 GB
each). LD-scores are a genome-wide summary rather than a per-locus quantity,
and a reference panel is the conventional S-LDSC choice. The UKB LD is still
used for the fine-mapping itself, where it matters. EUR is the right panel:
the GWAS is UK Biobank white-British and the loci are LDetect EUR blocks.

*Two annotation sets, run separately.* `ANNOT_SET=sniff` gives PolyFun the
102 Borzoi PCs only — non-standard for S-LDSC, since enrichment estimates may
be biased without baseline controls, but it means PolyFun sees exactly what
`fb_xregion` and `funmap` see. `ANNOT_SET=combined` adds baselineLF v2.2 for
289 annotations, the standard specification and PolyFun at its strongest.
Running both removes the objection either way: if `polyfun_combined` wins, the
289-annotation advantage is visible rather than hidden; if `polyfun_sniff`
already wins, the extra annotations were not what mattered.

Each set gets its own LD-scores (`ldscores_<set>/`), prior (`pf_prior_<set>/`)
and fine-mapping output (`polyfun_out_<set>/`).

S-LDSC is a genome-wide regression, so this runs over all 22 chromosomes.
Restricting it to the 265 significant blocks would be selection on the outcome.

### 1.6 Software (all upstream, unmodified)

| | |
|---|---|
| Functional BEATRICE | this project, `main` at **`7a46a13`** or later. Anything older lacks the `reformat_memo` max-offset fix and silently produces diluted PIPs, so `03_fb_xregion.pbs` tests the checkout for the fix by content and refuses to run without it. |
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

### Bootstrap once

```bash
export APP_ROOT=/rds/general/user/<you>/home/chol_app
bash setup.sh
```

`setup.sh` creates the directory tree, downloads the summary statistics,
clones upstream Funmap and PolyFun, builds the three Python environments,
installs susieR 0.11.92 for PolyFun, and writes `$APP_ROOT/env.sh`. It is
safe to re-run: every step skips work already done.

Three environments are needed because the pins genuinely conflict:

| env | for | pins |
|---|---|---|
| `APP_VENV` | prep, fb_xregion | numpy, pandas, pyarrow, scipy, torch |
| `FUNMAP_VENV` | Funmap | numpy 1.26.2, scipy 1.11.4, pandas 2.1.3, matplotlib 3.8.2 |
| `POLYFUN_ENV` | PolyFun | numpy<2, pyarrow 17, rpy2 3.5.17, plus susieR 0.11.92 in `POLYFUN_RLIB` |

Then put functional-beatrice at `$FB_DIR` — **`main` at `d0fb022` or later**;
`03_fb_xregion.pbs` verifies this by content and refuses to run otherwise.

### Submit

```bash
source "$APP_ROOT/env.sh"
bash submit_all.sh --dry-run     # inspect the qsub commands and dependencies
bash submit_all.sh
```

`submit_all.sh` chains the stages with PBS dependencies, so nothing starts
before its inputs exist:

```
  00b annotations ─┬──────────────────────> 02 prep ─┬─> 03 fb_xregion (ONE job, 160 GB)
  01  LD          ─┘                                  ├─> 04 funmap     (265 tasks)
       │                                              │
       └─> 00b2 combine (289)                         │
                                                      │
  00d LDSC ref ─> 00e munge ─┬─> 00f ldscores[sniff]    ─> 00g prior[sniff]    ─┼─> 05 polyfun[sniff]
                             └─> 00f ldscores[combined] ─> 00g prior[combined] ─┴─> 05 polyfun[combined]
```

`00d`–`00g` compute each PolyFun prior from our own annotations and run in
parallel with the LD download. `00f`, `00g` and `05` are each submitted twice,
once per `ANNOT_SET`; `submit_all.sh` does this for you.

The three downloads run concurrently. 03, 04 and 05 are independent of one
another. `--from 03` resumes at a stage with dependencies dropped.

### Check before committing CPU

```bash
python preflight.py data      # before prep:    downloads complete and non-truncated?
python preflight.py prep      # before methods: all 265 blocks built?
python preflight.py methods   # progress per method
```

`preflight.py` exits non-zero and names what is missing. The point is to fail
in seconds rather than a third of the way through a 265-task array.

### Collect

```bash
python 06_collect.py && python preflight.py collect
python 07_plot_blocks.py      # z-scores and PIPs per block
```

`06_collect.py` is deliberately **not** chained, so a partial failure in one
method does not quietly produce a half-empty results table.

Every stage is idempotent — a block with a `DONE` marker or an existing output
is skipped — so any array can be resubmitted after a partial failure. If you
change the p-value threshold, re-run `00a`/`00` and update the `#PBS -J`
array sizes to match the new block and window counts.

---

## 4. Diagnostics

Every stage sources `lib.sh`, which gives it:

- a **banner** recording job id, host, cpus, cwd, free disk, every relevant
  environment variable, and — after the venv is activated — the python that
  actually ran and its package versions. Most failures are an environment
  pointing somewhere unexpected, and this makes that visible immediately.
- an **ERR trap** naming the failing line number, the command, the working
  directory and the elapsed time. No more silent `set -e` exits.
- **preconditions** (`require_file`, `require_dir`, `require_env`) checked
  *before* expensive work. `require_file` takes a minimum size, because the
  failure that actually happens is a truncated download that is present but
  unusable.
- **output validation** (`check_output`) after it, so a stage cannot report
  success having written a zero-row file. `01` opens each `.npz` to prove it
  loads; `00f` checks the LD-score file has rows and annotation columns; `00g`
  asserts all 22 chromosomes exist and that SNPVAR is finite and positive;
  `04` checks funmap wrote one PIP per input variant.

Abnormal lines are prefixed `FATAL` or `WARN`, so an entire array triages with:

```bash
bash triage.sh            # all stages
bash triage.sh 02_prep    # one stage
```

It separates completed tasks, tasks that failed with a `FATAL` (showing the
message), and tasks with neither — which means the scheduler killed them, so
look at walltime and memory. It also groups the distinct failure modes, which
is what tells you whether 40 failures are one problem or forty.

A worked example of what a failure looks like:

```
[14:02:11] chol_prep[57]  ===== start =====
[14:02:11] chol_prep[57]  host     : cx3-6-21.cx3.hpc.ic.ac.uk  (4 cpus)
[14:02:11] chol_prep[57]  APP_ROOT : /rds/general/user/xx/home/chol_app
[14:02:11] chol_prep[57]  disk free: 812 GB
[14:02:12] chol_prep[57]  python   : .../appenv/bin/python (Python 3.11.6)
[14:02:12] chol_prep[57]  block    : ldEUR_04_0048  chr4  window chr4_68000001_71000001
[14:02:12] chol_prep[57]  FATAL missing file: ld/chr4_68000001_71000001.npz  (stage 01 incomplete for this window)
```

---

## 5. Resources

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
| fb_xregion | 1.55 | 0.561 | 62,362,116 | 6.05 | 1 |
| funmap | 9.99 | **1.000** | 62,370,732 | 6.00 | 13 |
| polyfun_susie | 9.04 | 0.407 | 62,370,732 | 6.00 | 13 |

This block holds **two independent signals** of near-identical strength 534 kb
apart — peak A at 62.36 Mb (max |Z| 6.05) and peak B at 62.91 Mb (max |Z|
6.04), r = 0.097 between their leads. All three methods now settle on peak A.
In blocks with multiple peaks, report per-signal rank rather than global rank:
a correct call on a second signal scores badly on a global |Z| rank.

### The LD regularisation defect (fixed)

Before `d0fb022`, `regularize_ld` added `(1e-3 − min(eig))` to the LD
**diagonal**, scaling every implied correlation by `1/(1+c)`. On ldEUR_01_0067
(min eig −0.069) that turned r = 0.996 into 0.931 — and the approximate Bayes
factor has a crossover near r = 0.955, above which co-inclusion of collinear
variants is penalised (−11.1 nats at 0.996) and below which it is **favoured**
(+8.2 nats at 0.931). 100% of sampled pairs with raw |r| > 0.955 were pushed
across it.

The effect was real and visible: at SORT1, the single highest-weighted
configuration in the memo was `{rs12740374, rs7528419}` (r = 0.996), 71% of
posterior mass sat on configurations containing both, and PIP summed to 1.97
across the four variants at r > 0.99.

Clipping the eigenvalues and renormalising the diagonal fixes it:

| | sum over r>0.99 | rs12740374 | rs7528419 | rs646776 | rs629301 |
|---|---|---|---|---|---|
| diagonal shift | 1.97 | 0.856 | 0.772 | 0.184 | 0.160 |
| **eigenvalue clip** | **1.00** | **0.942** | 0.059 | 0.002 | 0.000 |
| PolyFun (published) | 1.02 | 0.917 | 0.064 | 0.022 | 0.016 |

rs12740374 is the functionally validated causal variant for this locus
(Musunuru et al., *Nature* 2010). The four variants are in r = 0.996–1.000 LD
with |Z| spanning only 35.04–35.14, so the marginal statistics cannot separate
them — and fb_xregion and PolyFun, independently, both put ~0.92–0.94 on the
validated one.

**The benchmark is unaffected.** The shift only fired when the LD was not
positive definite. Simulated LD is PD (min eig +0.042) and takes the Cholesky
fast path; `examples/data` reproduces bit-identically after the change.

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
   the other two. Across five prepared blocks its confidence tracks signal
   strength: max PIP 0.86 at |Z| 35, 0.94 at |Z| 10, 0.73 at |Z| 6, but only
   0.16 at |Z| 5.0 and 0.03 at |Z| 5.5. Being unconfident at |Z| 5.5 among
   7,777 variants in LD is defensible; funmap returning exactly 1.000 on the
   same variant is the harder number to justify. Which of the two is better
   calibrated is an open question, and not one to settle by assuming the
   confident method is right.
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

## 6. Open items

- **Pre-registered prediction.** `RNAGTEX-PC5` (0-based) separates liver at
  3.5 SD, so a lipid trait should select it. Check
  `results/annotation_weights.csv` against Funmap's shrunk weights.
- **Correctness anchor.** Funmap's paper fine-maps cholesterol at
  chr8:5–8 Mb using this same LD bucket. Reproducing that locus would validate
  the whole chain end to end.
