# Genome-Wide Fine-Mapping Simulation: Technical Documentation

This document describes in full mathematical detail how `simulate_gwfm_data()` simulates genome-wide fine-mapping data. It covers the origin of the genotype data, every step of the phenotype simulation, the parameter space, the annotation model, and the precise distributions sampled at each stage. All notation is defined in the table below.

See also Section 14 for input data format requirements, genome coverage options, and instructions for downloading region sets and VCF files.

---

## 1. Notation

| Symbol | Definition |
|---|---|
| $K$ | Number of genomic regions |
| $n$ | Number of individuals |
| $p_i$ | Number of SNPs in region $i$ |
| $p_{\text{total}} = \sum_i p_i$ | Total number of SNPs across all regions |
| $\mathbf{X}_i \in \mathbb{R}^{n \times p_i}$ | Standardised genotype matrix for region $i$ |
| $\boldsymbol{\beta}_i \in \mathbb{R}^{p_i}$ | True (scaled) effect size vector for region $i$ |
| $\mathbf{y} \in \mathbb{R}^n$ | Phenotype vector, shared across all regions |
| $\pi$ | Genome-wide polygenicity: probability that any given variant is causal |
| $h^2$ | Target total SNP heritability |
| $\rho$ | `p_causal`: proportion of $h^2$ explained by the sparse component (sparse_inf model only) |
| $\sigma^2_\beta$ | `effect_variance`: prior variance of causal effect sizes under the normal model |
| $m$ | Number of functional annotation columns |
| $\gamma_k$ | Enrichment fold-factor for annotation $k$ |
| $\mathbf{A}_i \in \mathbb{R}^{p_i \times m}$ | Annotation matrix for region $i$ |
| $S_{\text{total}}$ | Realised total number of causal variants across all regions |
| $S_i$ | Realised number of causal variants in region $i$ (may be 0) |

---

## 2. Overview of the Simulation Pipeline

The simulation proceeds in four stages, each performed once per call to `simulate_gwfm_data()`:

1. **Genotype simulation** — genotype matrices are drawn from real 1000 Genomes haplotypes and standardised. LD matrices are pre-computed. This is done **once** and shared across all scenarios.
2. **Annotation simulation** — functional annotation matrices are drawn. Also done **once** and shared across all scenarios.
3. **Scenario loop** — for each combination of $(\pi, h^2, \text{iter})$, causal variants are assigned genome-wide, effect sizes drawn, the shared phenotype simulated, and per-region summary statistics computed from that shared phenotype.
4. **Output packaging** — results are stored in a format compatible with `run_gwfm_methods()` and `evaluate_methods()`.

The key structural difference from the locus-based pipeline (`run_simulation`) is that **a single phenotype vector $\mathbf{y}$ is shared across all $K$ regions within each scenario**, rather than generating independent phenotypes per region. This reflects how a real GWAS works: one phenotype is measured in $n$ individuals, and marginal associations are computed for every variant in the genome against that same $\mathbf{y}$.

---

## 3. Genomic Regions

### 3.1 Region set

The default region set (`regions = "representative"`) consists of **128 pre-defined 300 kb windows** spread across all 22 autosomes, bundled with the package as `inst/extdata/gwfm_regions.csv`. The number of regions per chromosome is approximately proportional to chromosome length (roughly 3–10 per chromosome).

Regions were selected to satisfy the following criteria:

- **LD independence**: regions are separated by at least 5 Mb, which is well beyond the typical range of linkage disequilibrium in human populations (LD decays to near-zero within ~500 kb in outbred populations).
- **Centromere avoidance**: regions avoid ±5 Mb around known centromere positions (GRCh37/hg19 coordinates).
- **MHC avoidance**: the major histocompatibility complex (chr6:28–34 Mb) is excluded due to exceptionally long-range and complex LD.
- **Acrocentric chromosome handling**: chromosomes 13, 14, 15, 21, and 22 have large heterochromatic short arms; all regions start at least 16 Mb from the chromosome start.
- **Telomere avoidance**: regions end at least 3 Mb from the chromosome end.

Users may supply their own regions as a data frame with columns `region_id`, `chrom`, `start`, `end`, or restrict to a single chromosome via `regions = "1"` (useful for small-scale testing).

### 3.2 VCF data source

Genotype data are derived from the **1000 Genomes Project Phase 3** reference panel (GRCh37 / hg19 build), which comprises 2,504 unrelated individuals from 26 global populations. The remote VCF files are hosted at:

```
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
ALL.chr{CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
```

The script `inst/scripts/prepare_gwfm_vcfs.R` uses `tabix` to stream precisely the 300 kb window for each region from the remote VCF — no whole-chromosome download is required. Each extracted file is stored as a bgzipped, tabix-indexed VCF in `data/gwfm_vcf/<region_id>.vcf.gz`. The total size across all 128 regions is approximately 400 MB.

---

## 4. Genotype Simulation

Genotype simulation uses the **sim1000G** R package, which generates new individuals by recombining haplotypes from the 1000 Genomes reference panel. This preserves realistic LD patterns within each region.

For each region $i$:

### 4.1 VCF loading

The VCF for region $i$ is read via `sim1000G::readVCF()`. To allow for MAF filtering and monomorphic removal, up to $1.5 \times p_i$ variants are initially read.

### 4.2 MAF filtering

Only variants with minor allele frequency $\geq$ `min_maf` (default 0.01) are retained. This removes rare variants that carry very little information and can cause numerical issues in LD estimation.

### 4.3 Recombination-based individual simulation

The genetic map for the chromosome is loaded from the HapMap GRCh37 recombination rate files (downloaded automatically by sim1000G and cached in `data/genetic_maps/`). sim1000G uses this map to simulate breakpoint positions during meiosis, then assembles $n$ simulated individuals by independently drawing two haplotypes from the 1000 Genomes reference panel and recombining them.

The raw genotype matrix $\mathbf{G}_i \in \{0, 1, 2\}^{n \times p_i}$ contains dosage counts (0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate).

### 4.4 Monomorphic removal

Any column of $\mathbf{G}_i$ that is monomorphic in the simulated sample (sample MAF = 0) is removed. This can occur even after VCF-level MAF filtering because the simulated sample is a random draw.

**Realised vs requested $p_i$.** The requested per-region target `p` is treated as an upper bound: if the VCF window contains fewer variants than `p`, or if extra variants are lost to MAF filtering and monomorphic removal, then $p_i$ < requested. The realised $p_i$ for every region is recorded in `params$p_actual` (a length-$K$ integer vector); the requested values are kept in `params$p_requested`. Downstream evaluation always uses $p_i$ from `genotypes[[i]]$p`, but stratifying analyses by `p` should use `p_actual` rather than the requested value.

### 4.5 Standardisation

Each column $j$ of $\mathbf{G}_i$ is standardised to mean 0 and variance 1:

$$X_{ij} = \frac{G_{ij} - \bar{G}_j}{s_j}$$

where $\bar{G}_j = \frac{1}{n}\sum_{l=1}^n G_{lj}$ is the sample mean and $s_j = \sqrt{\frac{1}{n-1}\sum_{l=1}^n (G_{lj} - \bar{G}_j)^2}$ is the sample standard deviation. The resulting matrix $\mathbf{X}_i \in \mathbb{R}^{n \times p_i}$ has columns with mean 0 and variance 1.

Standardisation to unit variance ensures that the scale of effect sizes $\boldsymbol{\beta}_i$ is interpretable per-standard-deviation of genotype, and that the genetic variance $\text{Var}(\mathbf{X}_i \boldsymbol{\beta}_i)$ does not depend on allele frequencies.

### 4.6 LD matrix

The LD matrix for region $i$ is the **sample Pearson correlation matrix** of $\mathbf{X}_i$:

$$\boldsymbol{\Sigma}_i = \frac{1}{n-1} \mathbf{X}_i^\top \mathbf{X}_i \in \mathbb{R}^{p_i \times p_i}$$

Since $\mathbf{X}_i$ is already column-standardised, this is equivalent to `cor(X_i)` in R. The diagonal entries are 1 and off-diagonal entries $\Sigma_{i,jk}$ are the sample correlations between variants $j$ and $k$ within region $i$.

LD matrices are pre-computed once after genotype simulation and reused across all $(\pi, h^2, \text{iter})$ scenarios.

### 4.7 Optional reference-panel LD (sample-size mismatch)

To simulate the practical case in which fine-mapping methods receive an
LD matrix computed from a **smaller, independent reference panel** rather
than the GWAS sample itself (e.g. 1000 Genomes EUR with $n \approx 2{,}500$
used to fine-map a $n = 500{,}000$ UK Biobank GWAS), the `n_ref` argument
to `simulate_genotypes()`, `run_simulation()`, and `simulate_gwfm_data()`
triggers an additional independent draw of size $n_{\text{ref}}$ from the
same VCF, in the same sim1000G session.

When `n_ref` is set:

- A second sample $\mathbf{G}_i^{\text{ref}} \in \{0, 1, 2\}^{n_{\text{ref}} \times p_i}$
  is drawn alongside the GWAS sample. Both draws use the same VCF load,
  the same genetic map, and the same MAF-filtered variant set.
- The polymorphic filter (§4.4) is applied to the **intersection** of
  variants polymorphic in both samples — variants monomorphic in either
  the GWAS or the reference are dropped from both, mirroring how a real
  reference-panel workflow handles missing or absent variants.
- $\mathbf{X}_i^{\text{ref}}$ is standardised by its **own** column means
  and SDs, not the GWAS sample's. This mirrors how a real reference panel
  produces an LD matrix on its own scale.
- The reference-panel LD is
  $\boldsymbol{\Sigma}_i^{\text{ref}} = \mathrm{cor}(\mathbf{X}_i^{\text{ref}})$
  and is stored as `genotypes[[i]]$LD`. The in-sample LD
  $\boldsymbol{\Sigma}_i^{\text{true}} = \mathrm{cor}(\mathbf{X}_i)$ is
  also stored as `genotypes[[i]]$LD_true` for diagnostics.
- Fine-mapping methods receive `LD` (ref-panel-derived). `LD_true` lets
  evaluation compute mismatch diagnostics — e.g.
  $\frac{1}{p_i^2}\|\boldsymbol{\Sigma}_i^{\text{ref}} - \boldsymbol{\Sigma}_i^{\text{true}}\|_F^2$
  — without re-simulating.

When `n_ref = NULL` (the default), no second sample is drawn. `LD == LD_true`
exactly, and `LD_true` is the only correlation matrix stored. This keeps
the behaviour backwards-compatible with pre-LD-mismatch simulations.

`params$n_ref` records the chosen reference-panel size (or `NULL`) for
reproducibility. The two LD matrices roughly double the per-region memory
footprint when `n_ref` is set; the LD memory guard accounts for this.

---

## 5. Functional Annotation Simulation

Annotation matrices, if requested, are also simulated once and shared across all scenarios.

### 5.1 Binary annotations (`annotations = "binary"`)

For each region $i$ and annotation $k = 1, \ldots, m$, a binary indicator vector is drawn:

$$A_{i,jk} \sim \text{Bernoulli}(q_k), \quad j = 1, \ldots, p_i$$

independently for all variants $j$. The proportion $q_k$ is either:

- **User-specified** via `annotation_proportions[k]`
- **Randomly drawn** from $q_k \sim \text{Uniform}(0.01, 0.30)$ if `annotation_proportions = NULL`

The same $q_k$ is used for all regions (annotation $k$ has the same marginal density genome-wide), but the actual indicator values differ across variants and regions due to independent sampling.

### 5.2 Continuous annotations (`annotations = "continuous"`)

Each entry of the annotation matrix is drawn independently from a standard normal:

$$A_{i,jk} \sim \mathcal{N}(0, 1), \quad j = 1, \ldots, p_i, \quad k = 1, \ldots, m$$

### 5.3 User-supplied annotations (`annotations = <matrix>`)

The user may supply a single $p_{\text{total}} \times m$ matrix directly. This is split row-wise by region in order (rows $1, \ldots, p_1$ to region 1, rows $p_1+1, \ldots, p_1+p_2$ to region 2, etc.).

---

## 6. Causal Variant Assignment

For each scenario, causal variants are assigned **genome-wide** rather than per-region. This is the key distributional difference from the locus-based pipeline.

### 6.1 Annotation enrichment weights

If annotations are present, per-variant enrichment weights are first computed. A single set of enrichment values $\gamma_k > 0$ is shared across all regions. Each $\gamma_k$ is either:

- **User-specified** via `enrichment[k]`
- **Randomly drawn** from $\gamma_k \sim \text{Uniform}(2, 10)$ if `enrichment = NULL`

For variant $j$ in region $i$, the unnormalised log-weight is:

$$\log w_{ij} = \sum_{k=1}^{m} A_{i,jk} \log \gamma_k$$

and the weight is $w_{ij} = \exp(\log w_{ij} - \max_{j'} \log w_{ij'})$ (the subtraction of the maximum is for numerical stability and does not affect the relative weights). The per-variant causal probability is then:

$$\pi_{ij} = \pi \cdot \frac{w_{ij}}{\bar{w}_i}, \quad \text{where } \bar{w}_i = \frac{1}{p_i} \sum_{j=1}^{p_i} w_{ij}$$

This normalisation ensures that $\frac{1}{p_i}\sum_j \pi_{ij} = \pi$ within each region, so the **expected number of causal variants per region is $\pi \cdot p_i$**, and the **expected total number of causal variants is $\pi \cdot p_{\text{total}}$**, regardless of annotation structure.

**Numerical clamping caveat.** Before sampling, $\pi_{ij}$ is clamped to $[10^{-10}, 1-10^{-10}]$ for safety. When `enrichment` is large and an annotation is sparse, a small number of variants can have $\pi_{ij}$ that would exceed 1 under the rescaling; the clamp shaves off that excess and breaks the $\bar{\pi}_i = \pi$ guarantee. The realised expected number of causals can then fall slightly below $\pi \cdot p_{\text{total}}$. The realised total $S_{\text{total}}$ is always recorded in `scenarios[[s]]$S_total`, so downstream analysis should use that rather than $\pi \cdot p_{\text{total}}$ when comparing.

When no annotations are present, $\pi_{ij} = \pi$ for all variants.

### 6.2 Bernoulli assignment

Each variant $j$ in region $i$ is independently assigned as causal:

$$z_{ij} \sim \text{Bernoulli}(\pi_{ij}), \quad j = 1, \ldots, p_i, \quad i = 1, \ldots, K$$

The causal index set for region $i$ is $\mathcal{C}_i = \{j : z_{ij} = 1\}$, with realised size $S_i = |\mathcal{C}_i|$.

**Important**: unlike the locus-based pipeline where every region is guaranteed exactly $S$ causal variants, here $S_i$ is a random variable. The marginal distribution of $S_i$ (without annotations) is $\text{Binomial}(p_i, \pi)$, with:

$$\mathbb{E}[S_i] = \pi p_i, \qquad \text{Var}(S_i) = \pi(1-\pi)p_i$$

In particular, **$S_i = 0$ is possible** (a null region), especially when $\pi$ is small. Null regions are included in evaluation and contribute true negatives (all PIPs should be low).

The total realised causal count is $S_{\text{total}} = \sum_i S_i$.

---

## 7. Effect Size Simulation

### 7.1 Raw effect sizes

For each causal variant $j \in \mathcal{C}_i$, a raw effect size $\tilde{\beta}_{ij}$ is drawn independently. Non-causal variants have $\tilde{\beta}_{ij} = 0$. Two distributions are available:

**Normal** (`effect_distribution = "normal"`):
$$\tilde{\beta}_{ij} \mid z_{ij} = 1 \sim \mathcal{N}(0,\ \sigma^2_\beta)$$
where $\sigma^2_\beta$ = `effect_variance` (default 0.36, i.e. standard deviation 0.6).

**Equal magnitudes** (`effect_distribution = "equal"`):
$$\tilde{\beta}_{ij} \mid z_{ij} = 1 = \pm 1 \quad \text{with equal probability}$$
Signs are drawn independently as $\text{Bernoulli}(0.5)$ mapped to $\{-1, +1\}$.

### 7.2 Genome-wide genetic signal

The raw genome-wide genetic signal is:

$$\tilde{\mathbf{g}} = \sum_{i=1}^{K} \mathbf{X}_i \tilde{\boldsymbol{\beta}}_i \in \mathbb{R}^n$$

with sample variance $\widetilde{V}_g = \widehat{\text{Var}}(\tilde{\mathbf{g}})$ (computed as the sample variance across $n$ individuals).

---

## 8. Phenotype Simulation

### 8.1 Sparse model (`model = "sparse"`)

**Scaling**: The raw genetic signal is scaled so that its sample variance equals the target heritability $h^2$:

$$c = \sqrt{\frac{h^2}{\widetilde{V}_g}}, \qquad \mathbf{g} = c \cdot \tilde{\mathbf{g}}, \qquad \boldsymbol{\beta}_i = c \cdot \tilde{\boldsymbol{\beta}}_i$$

By construction, $\widehat{\text{Var}}(\mathbf{g}) = h^2$.

**Residual noise**: The residual variance is fixed at $\sigma^2_\varepsilon = 1 - h^2$, so that $\widehat{\text{Var}}(\mathbf{y}) \approx h^2 + (1 - h^2) = 1$ in expectation. The noise vector is:

$$\boldsymbol{\varepsilon} \sim \mathcal{N}(\mathbf{0},\ (1 - h^2)\mathbf{I}_n)$$

**Phenotype**:

$$\mathbf{y} = \mathbf{g} + \boldsymbol{\varepsilon} = \sum_{i=1}^K \mathbf{X}_i \boldsymbol{\beta}_i + \boldsymbol{\varepsilon}$$

**Realised heritability**: After sampling, the realised $h^2$ is reported as:

$$\hat{h}^2 = \frac{\widehat{\text{Var}}(\mathbf{g})}{\widehat{\text{Var}}(\mathbf{y})}$$

This differs from the target $h^2$ by a finite-sample fluctuation that shrinks as $n \to \infty$.

### 8.2 Sparse + infinitesimal model (`model = "sparse_inf"`)

The phenotype is decomposed into three independent components:

$$\mathbf{y} = \mathbf{g}_{\text{sparse}} + \mathbf{g}_{\text{inf}} + \boldsymbol{\varepsilon}$$

where the target variances are:

$$\widehat{\text{Var}}(\mathbf{g}_{\text{sparse}}) = \rho \cdot h^2, \qquad \widehat{\text{Var}}(\mathbf{g}_{\text{inf}}) = (1-\rho) \cdot h^2, \qquad \sigma^2_\varepsilon = 1 - h^2$$

and $\rho$ = `p_causal` $\in (0, 1]$.

**Sparse component**: drawn and scaled exactly as in Section 8.1, but targeting $\rho \cdot h^2$:

$$c_{\text{sparse}} = \sqrt{\frac{\rho \cdot h^2}{\widetilde{V}_g}}, \qquad \mathbf{g}_{\text{sparse}} = c_{\text{sparse}} \cdot \tilde{\mathbf{g}}$$

**Infinitesimal component** (see Section 8.3 for derivation):

$$c_{\text{inf}} = \sqrt{\frac{(1-\rho) \cdot h^2}{\widetilde{V}_{\text{inf}}}}, \qquad \mathbf{g}_{\text{inf}} = c_{\text{inf}} \cdot \tilde{\mathbf{g}}_{\text{inf}}$$

**Residual**:

$$\boldsymbol{\varepsilon} \sim \mathcal{N}(\mathbf{0},\ (1 - h^2)\mathbf{I}_n)$$

### 8.3 Infinitesimal component

Two formulations are implemented, corresponding to the BEATRICE and SuSiE-inf models.

#### BEATRICE formulation (`inf_model = "beatrice"`)

Infinitesimal effects are placed **only on non-causal variants**. Let $\mathcal{C}_i^c = \{1,\ldots,p_i\} \setminus \mathcal{C}_i$ be the non-causal index set in region $i$, with $|\mathcal{C}_i^c| = p_i - S_i$ non-causal variants. Let $\mathbf{X}_{i,\mathcal{C}_i^c}$ denote the submatrix of $\mathbf{X}_i$ restricted to non-causal columns.

For each non-causal variant $j \in \mathcal{C}_i^c$ across all regions, an infinitesimal effect is drawn:

$$\alpha_j \sim \mathcal{N}\!\left(0,\ \frac{1}{p_{\text{total}}}\right)$$

independently. The raw infinitesimal signal is:

$$\tilde{\mathbf{g}}_{\text{inf}} = \sum_{i=1}^{K} \mathbf{X}_{i,\mathcal{C}_i^c}\, \boldsymbol{\alpha}_{\mathcal{C}_i^c}$$

The choice of variance $1/p_{\text{total}}$ follows BEATRICE (Kadie & Heckerman, 2021). Under this parameterisation, with standardised genotypes and independent effects, $\mathbb{E}[\widehat{\text{Var}}(\tilde{\mathbf{g}}_{\text{inf}})] \approx (p_{\text{total}} - S_{\text{total}}) / p_{\text{total}} \approx 1 - \pi$. The signal is then rescaled to $(1-\rho) h^2$ using the sample-based scale factor $c_{\text{inf}}$ above.

#### SuSiE-inf formulation (`inf_model = "susie_inf"`)

Infinitesimal effects are placed on **all variants**, including causal ones. For each variant $j$ in region $i$:

$$\alpha_j \sim \mathcal{N}\!\left(0,\ \frac{1}{p_{\text{total}}}\right)$$

The raw infinitesimal signal is:

$$\tilde{\mathbf{g}}_{\text{inf}} = \sum_{i=1}^{K} \mathbf{X}_i\, \boldsymbol{\alpha}_i$$

This follows the SuSiE-inf formulation (Karber et al., 2024), in which every variant has a small polygenicity background in addition to any sparse effect.

---

## 9. Summary Statistics

After simulating $\mathbf{y}$, marginal summary statistics are computed **per region** from the **shared phenotype**. For each region $i$ and variant $j \in \{1, \ldots, p_i\}$:

Let $\mathbf{x}_{ij}$ denote column $j$ of $\mathbf{X}_i$ and $\mathbf{y}_c = \mathbf{y} - \bar{y}$ the mean-centred phenotype. The marginal OLS estimator (univariate regression of $\mathbf{y}_c$ on $\mathbf{x}_{ij}$ without intercept, valid because $\mathbf{x}_{ij}$ is mean-zero by standardisation) is:

$$\hat{\beta}_{ij} = \frac{\mathbf{x}_{ij}^\top \mathbf{y}_c}{\mathbf{x}_{ij}^\top \mathbf{x}_{ij}}$$

The residual vector is $\mathbf{r}_{ij} = \mathbf{y}_c - \mathbf{x}_{ij} \hat{\beta}_{ij}$, and the residual variance estimate (with $n - 2$ degrees of freedom for the implicit intercept) is:

$$\hat{\sigma}^2_{ij} = \frac{\|\mathbf{r}_{ij}\|^2}{n - 2}$$

The standard error and $z$-score are:

$$\widehat{\text{se}}_{ij} = \sqrt{\frac{\hat{\sigma}^2_{ij}}{\mathbf{x}_{ij}^\top \mathbf{x}_{ij}}}, \qquad z_{ij} = \frac{\hat{\beta}_{ij}}{\widehat{\text{se}}_{ij}}$$

These are computed independently for all $j$ (no joint model). The critical point is that $\mathbf{y}$ is the same vector for all regions: the $z$-score of a variant in region $i$ is influenced by the genetic signal from all other regions through the shared residual.

In contrast, the locus-based pipeline computes $\hat{\beta}_{ij}$ by regressing each region's own independently-simulated phenotype $\mathbf{y}_i = \mathbf{X}_i \boldsymbol{\beta}_i + \boldsymbol{\varepsilon}_i$ against $\mathbf{x}_{ij}$, so the summary statistics for region $i$ carry no information about other regions.

---

## 10. Parameter Space and Scenario Grid

The full parameter grid is constructed by `expand.grid` over all combinations of the swept parameters:

**Sparse model**:

$$\{(\pi, h^2, t) : \pi \in \Pi,\ h^2 \in \mathcal{H},\ t \in \{1,\ldots,T\}\}$$

giving $|\Pi| \times |\mathcal{H}| \times T$ scenarios.

**Sparse + infinitesimal model** adds `p_causal` ($\rho$) as a swept parameter:

$$\{(\pi, h^2, \rho, t) : \pi \in \Pi,\ h^2 \in \mathcal{H},\ \rho \in \mathcal{P},\ t \in \{1,\ldots,T\}\}$$

giving $|\Pi| \times |\mathcal{H}| \times |\mathcal{P}| \times T$ scenarios.

Within each scenario:
- Genotypes $\{\mathbf{X}_i\}$ and annotations $\{\mathbf{A}_i\}$ are **fixed** (shared from the single pre-simulation).
- Causal indices $\{\mathcal{C}_i\}$, effect sizes $\{\boldsymbol{\beta}_i\}$, and phenotype $\mathbf{y}$ are **freshly drawn** from their respective distributions.

The iteration index $t$ thus provides $T$ independent realisations of the genetic architecture and phenotype for each $(\pi, h^2)$ parameter combination, enabling stable estimation of mean metrics and their standard errors across replicates.

---

## 11. Ground Truth Storage

For each region $i$ in each scenario, the following ground truth quantities are stored alongside the summary statistics:

| Field | Type | Definition |
|---|---|---|
| `causal_indices` | integer vector | Indices $\mathcal{C}_i$ of causal variants within region $i$ (may be empty) |
| `causal_effects` | numeric vector | Scaled effect sizes $\boldsymbol{\beta}_i[\mathcal{C}_i]$ |
| `beta_true` | numeric vector (length $p_i$) | Full effect vector: $\boldsymbol{\beta}_i$ (0 for non-causal) |
| `pve_region` | numeric | $\widehat{\text{Var}}(\mathbf{X}_i\boldsymbol{\beta}_i) / \widehat{\text{Var}}(\mathbf{y})$: region's contribution to total phenotypic variance |
| `S_realized` | integer | $S_i = |\mathcal{C}_i|$ (may be 0 for null regions) |
| `pi` | numeric | Genome-wide $\pi$ used in this scenario |
| `h2` | numeric | Target $h^2$ used in this scenario |

---

## 12. Key Differences from the Locus-Based Pipeline

| Aspect | Locus-based (`run_simulation`) | Genome-wide (`simulate_gwfm_data`) |
|---|---|---|
| **Causal assignment** | Fixed $S$ causal variants per region | Bernoulli($\pi$) per variant; $S_i$ is random |
| **Null regions** | None — every region has $S \geq 1$ | $S_i = 0$ possible, especially for small $\pi$ |
| **Phenotype** | Independent $\mathbf{y}_i = \mathbf{X}_i\boldsymbol{\beta}_i + \boldsymbol{\varepsilon}_i$ per region | Shared $\mathbf{y} = \sum_i \mathbf{X}_i\boldsymbol{\beta}_i + \boldsymbol{\varepsilon}$ |
| **Heritability** | $\phi$ = per-region PVE; each region explains $\phi$ of its own variance | $h^2$ = total genome-wide SNP heritability |
| **Summary stats** | Regressed against per-region $\mathbf{y}_i$ | Regressed against shared genome-wide $\mathbf{y}$ |
| **Cross-region confounding** | None by construction | Present: $z_{ij}$ is influenced by signal in all other regions through shared $\mathbf{y}$ |
| **Parameter swept** | $(S, \phi, t)$ grid | $(\pi, h^2, t)$ grid |
| **Enrichment normalisation** | Per-region: weights normalised within region | Per-region: weights normalised within region; $\bar{\pi}_{ij} = \pi$ within each region |

---

## 13. Function Arguments Reference

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `n` | integer | — | Number of individuals |
| `n_iter` | integer | 5 | Number of replicates $T$ per $(\pi, h^2)$ combination |
| `pi` | numeric vector | `c(1e-4, 1e-3)` | Genome-wide polygenicity values $\Pi$ to sweep |
| `h2` | numeric vector | `c(0.1, 0.3)` | Total heritability values $\mathcal{H}$ to sweep; must be in $(0,1)$ |
| `model` | character | `"sparse"` | `"sparse"` or `"sparse_inf"` |
| `p_causal` | numeric | 0.5 | $\rho$: sparse fraction of $h^2$; only for `sparse_inf` |
| `inf_model` | character | `"beatrice"` | `"beatrice"` (non-causal only) or `"susie_inf"` (all variants) |
| `effect_distribution` | character | `"normal"` | `"normal"` ($\mathcal{N}(0,\sigma^2_\beta)$) or `"equal"` ($\pm 1$) |
| `effect_variance` | numeric | 0.36 | $\sigma^2_\beta$; only for `"normal"` distribution |
| `annotations` | character or matrix | `"none"` | Annotation type: `"none"`, `"binary"`, `"continuous"`, or matrix |
| `n_annotations` | integer | 3 | Number of annotation columns $m$ |
| `annotation_proportions` | numeric vector or NULL | NULL | $q_k$ values for binary annotations; NULL = random from $U(0.01, 0.30)$ |
| `enrichment` | numeric vector or NULL | NULL | $\gamma_k$ values; NULL = random from $U(2, 10)$ |
| `regions` | character or data frame | `"representative"` | Region set: `"representative"`, `"<chr>"`, or custom data frame (see Section 14) |
| `coverage` | numeric in $(0,1]$ | `1` | Fraction of the loaded region set to use; stratified by chromosome (see Section 14) |
| `p` | integer | 200 | Target SNPs per region (after MAF filtering) |
| `min_maf` | numeric | 0.01 | Minimum MAF filter applied to VCF variants |
| `vcf_dir` | character | `"data/gwfm_vcf"` | Directory containing per-region VCF files |
| `genetic_map_dir` | character | `"data/genetic_maps"` | Directory for cached HapMap recombination maps |
| `seed` | integer or NULL | NULL | Master random seed (set before genotype simulation) |
| `save` | logical | FALSE | If TRUE, save result as `.rds` in `output_dir` |
| `output_dir` | character | `"results"` | Directory for saved output |
| `verbose` | logical | TRUE | Print progress messages |

---

## 14. Input Data Requirements, Coverage Options, and Download Instructions

### 14.1 Required input: VCF files

`simulate_gwfm_data()` requires one bgzipped, tabix-indexed VCF file per region. These are **not downloaded automatically** — they must be prepared in advance by running one of the preparation scripts below. The function will error clearly if any expected VCF file is missing.

**Format requirements for each VCF file:**

| Property | Requirement |
|---|---|
| Genome build | GRCh37 / hg19 (must match 1000G Phase 3 coordinates) |
| Compression | bgzipped (`.vcf.gz`) — use `bgzip`, not `gzip` |
| Index | tabix index (`.vcf.gz.tbi`) — produced by `tabix -p vcf file.vcf.gz` |
| Naming | `<region_id>.vcf.gz`, where `region_id` matches the `region_id` column in the region CSV |
| Location | All files in a single directory, passed as `vcf_dir` |
| Content | Standard VCF v4.x with GT fields; multi-allelic sites are handled by sim1000G |

A VCF file for a 300 kb region from 1000 Genomes Phase 3 is typically **2–5 MB** on disk.

### 14.2 Required input: region CSV file

The region set is specified by a CSV file with the following columns:

| Column | Type | Description |
|---|---|---|
| `region_id` | character | Unique identifier; must match the VCF filename (e.g. `gw001` → `gw001.vcf.gz`) |
| `chrom` | integer | Autosomal chromosome number (1–22); no `chr` prefix |
| `start` | integer | Window start position (GRCh37, 0-based or 1-based — consistent with the VCF) |
| `end` | integer | Window end position |
| `notes` | character | Optional description; ignored by the function |

Additional columns are permitted and ignored.

### 14.3 Genome coverage options

Three region sets are available, covering different fractions of the autosome:

| Region set | Regions | LD blocks covered | SNPs (at p=200) | VCF download | Genotype sim. time\* |
|---|---|---|---|---|---|
| Bundled (`"representative"`) | 128 | ~7.5% of LDetect EUR | ~25,600 | ~400 MB | ~15 min |
| LDetect EUR (full) | ~1,703 | ~100% | ~340,000 | ~5 GB | ~3 hr |
| LDetect AFR (full) | ~1,445 | ~100% | ~290,000 | ~4.5 GB | ~2.5 hr |
| LDetect ASN (full) | ~1,647 | ~100% | ~330,000 | ~5 GB | ~3 hr |

\* Approximate, for $n = 1{,}000$ on a modern laptop. Time scales linearly with $n$.

The `coverage` argument controls what fraction of the **loaded** region set is used:

- `coverage = 1` (default): use all regions in the loaded set
- `coverage = 0.5`: use approximately half, sampled proportionally within each chromosome
- `coverage = 0.1`: use approximately 10%

Subsampling is **stratified by chromosome**: for chromosome $c$ with $n_c$ regions, $\text{round}(\texttt{coverage} \times n_c)$ regions are sampled without replacement. Chromosomes where this rounds to zero are omitted. This preserves genome-wide spread at all coverage values.

**Combining `regions` and `coverage`:** use `regions` to select the upper bound (which region set to load) and `coverage` to select the density within that set. For example, to use 20% of the full LDetect EUR partition:

```r
ldetect <- read.csv(system.file("extdata", "gwfm_regions_ldetect_EUR.csv", package = "fmbenchmark"))
sim <- simulate_gwfm_data(
  n       = 2000,
  regions = ldetect,
  coverage = 0.20,          # ~340 of 1703 regions
  vcf_dir = "data/gwfm_vcf_ldetect_EUR",
  ...
)
```

### 14.4 Downloading the bundled region set VCFs (128 regions, ~400 MB)

Run once from the project root (source checkout). Requires `tabix` and `bgzip` (install via `brew install htslib` on macOS or `conda install -c bioconda htslib`):

```bash
Rscript inst/scripts/prepare_gwfm_vcfs.R
```

This streams 128 × 300 kb windows from the 1000 Genomes Phase 3 FTP server using tabix, compresses them with bgzip, and saves them to `data/gwfm_vcf/`. No full chromosome is downloaded.

### 14.5 Downloading a full LDetect partition (~1,703 blocks, ~5 GB VCF)

**LDetect** (Berisa & Pickrell, 2016, *Bioinformatics* 32:283–285; doi: [10.1093/bioinformatics/btv546](https://doi.org/10.1093/bioinformatics/btv546)) provides approximately independent LD block partitions for three ancestry groups, estimated from 1000 Genomes Phase 3 haplotypes using a Fourier approach.

Data source: [https://bitbucket.org/nygcresearch/ldetect-data](https://bitbucket.org/nygcresearch/ldetect-data)

Available partitions:

| File (on Bitbucket) | Population | Blocks |
|---|---|---|
| `EUR/fourier_ls-all.bed` | European | ~1,703 |
| `AFR/fourier_ls-all.bed` | African | ~1,445 |
| `ASN/fourier_ls-all.bed` | East Asian | ~1,647 |

**Step 1 — Download and convert the LDetect block coordinates** (fast, ~seconds, no large files):

```bash
# Edit POPULATION <- "EUR" (or "AFR"/"ASN") in the script first
Rscript inst/scripts/download_ldetect_regions.R
```

This produces `inst/extdata/gwfm_regions_ldetect_EUR.csv` (same format as `inst/extdata/gwfm_regions.csv`). Each block's central 300 kb window is extracted as the simulation region.

**Step 2 — Download VCF files** (slow, ~5 GB, run only when needed):

Either set `DOWNLOAD_VCFS <- TRUE` inside `download_ldetect_regions.R` and re-run, or adapt `prepare_gwfm_vcfs.R` to point at the LDetect CSV:

```r
# In inst/scripts/prepare_gwfm_vcfs.R, change:
REGIONS <- find_extdata("gwfm_regions_ldetect_EUR.csv")
VCF_DIR <- "data/gwfm_vcf_ldetect_EUR"
```

Then:

```bash
Rscript inst/scripts/prepare_gwfm_vcfs.R
```

The VCF download can be interrupted and resumed — already-downloaded files are skipped unless `OVERWRITE <- TRUE`.

### 14.6 Alternative region sources

Users may supply any data frame of regions meeting the format in Section 14.2. Other published LD block resources that are compatible with this format include:

- **GCTB LD blocks** — used by SBayesRC and the GWFM paper (Wu et al., 2026). Available alongside the GCTB software at [https://cnsgenomics.com/software/gctb](https://cnsgenomics.com/software/gctb). Coordinates must be converted from their format to our CSV format.
- **1000 Genomes LD blocks** — any custom partition of the 1000G reference panel is valid, provided VCF data are available for each window.
- **Custom loci** — for targeted simulation around specific GWAS loci, supply a hand-curated data frame with the relevant windows.
