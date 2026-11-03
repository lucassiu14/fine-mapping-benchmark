#!/bin/bash
# Build application/README.pdf from README.md.
#
# This TinyTeX has no fvextra, and `tlmgr install fvextra` is refused because the
# local installation is TeX Live 2024 against a 2026 repository. Verbatim therefore
# cannot wrap, so three long lines and one long inline path are rewritten in a
# throwaway copy of the README; the committed README is never touched.
set -euo pipefail
APP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SC="$(mktemp -d)"
trap 'rm -rf "$SC"' EXIT
python3 - "$APP/README.md" "$SC/README_pdf.md" <<'PY'
import sys
src = open(sys.argv[1]).read()
pairs = [
# The two long download URLs leave the code fence so pandoc emits \url{}, which
# breaks at slashes and stays copy-pastable.
("""```
https://storage.googleapis.com/seqnn-share/sniff/_20260810_tmp_files/fixed_n_337k/sample_munged_gwas/biochemistry_Cholesterol.sumstats_munged.parquet   (526 MB)
```""",
 """<https://storage.googleapis.com/seqnn-share/sniff/_20260810_tmp_files/fixed_n_337k/sample_munged_gwas/biochemistry_Cholesterol.sumstats_munged.parquet> (526 MB)"""),
("""```
https://storage.googleapis.com/seqnn-share/sniff/processed_borzoi_parquets/borzoi.{CHR}.annot.parquet   (0.2-1.3 GB each, ~15 GB total)
```""",
 """<https://storage.googleapis.com/seqnn-share/sniff/processed_borzoi_parquets/borzoi.%7BCHR%7D.annot.parquet> (0.2-1.3 GB each, ~15 GB total)"""),
# A 74-character path inside running prose: \texttt does not break, so it ran 2.75 cm
# into the margin. Displayed instead, where it fits.
("""Precomputed priors do exist in the same bucket
(`sniff_results_ukb_15/{sniff,baseline}/biochemistry_Cholesterol.{CHR}.txt.gz`)
and they are for the **same GWAS**""",
 """Precomputed priors do exist in the same bucket

```
sniff_results_ukb_15/{sniff,baseline}/biochemistry_Cholesterol.{CHR}.txt.gz
```

and they are for the **same GWAS**"""),
# One log line 10 characters too wide; wrapped where the log itself would wrap.
("[14:02:12] chol_prep[57]  FATAL missing file: ld/chr4_68000001_71000001.npz  (stage 01 incomplete for this window)",
 "[14:02:12] chol_prep[57]  FATAL missing file: ld/chr4_68000001_71000001.npz\n[14:02:12] chol_prep[57]        (stage 01 incomplete for this window)"),
]
for a, b in pairs:
    assert a in src, "pattern no longer in README.md: " + a[:70]
    src = src.replace(a, b)
# The three figures live in application/figures/ (gitignored, rebuilt by
# 12_figures.py) and are not referenced inline by the README. Appended here so
# the PDF is a self-contained read; the README body is unchanged.
src += """

------------------------------------------------------------------------

## Appendix: figures

Produced by `12_figures.py` from the export bundle, not from the cluster.
Colour marks the four methods; grey marks the three ablations of our own code,
which are states, not methods (§6.1).

![Posterior inclusion probabilities across all 254 LDetect blocks,
UK Biobank total cholesterol. Four methods above the rule, three ablations below
it. Each panel’s right-hand note counts the variants that panel calls at
PIP ≥ 0.95 and the blocks holding them (§6.3).](figures/fig1_pips_genome.pdf){width=82%}

![Method behaviour. (a) the |Z| carried by each region’s own
top-PIP variant (§6.4); (b) every variant called at PIP ≥ 0.95, against
its |Z| (§6.3); (c) sum(PIP) per block against each model’s own cap, the
saturation check (§6.6); (d) 95% credible-set size.](figures/fig2_method_behaviour.pdf)

![The seven positive controls (§6.11). Methods left of the
rule, ablations right. The last column gives k, the number of variants at
r² ≥ 0.99 with the control: where k > 1 the association data alone cannot
put more than about 1/k on any one of them, so a PIP near 1 there is the prior
talking.](figures/fig3_positive_controls.pdf)
"""
open(sys.argv[2], "w").write(src)
PY
cat > $SC/pdfhead.tex <<'EOF'
\setlength{\tabcolsep}{3.5pt}
\usepackage{etoolbox}
\AtBeginEnvironment{longtable}{\footnotesize}
EOF
# Lucida Grande, not Helvetica Neue: it is the only font on this machine with
# U+2192, and a missing glyph would silently delete the arrow from every
# "653 -> 136" in section 6.
pandoc "$SC/README_pdf.md" -o "$APP/README.pdf" \
  --pdf-engine=xelatex --toc --toc-depth=3 \
  -f markdown-tex_math_dollars-raw_tex \
  -H "$SC/pdfhead.tex" --resource-path="$APP" \
  -V documentclass=article -V classoption=10pt \
  -V geometry:margin=1.7cm -V colorlinks=true -V linkcolor=RoyalBlue -V toccolor=black \
  -V mainfont="Lucida Grande" -V monofont="Menlo" -V monofontoptions="Scale=0.80" \
  -V title="UK Biobank cholesterol application" \
  -V subtitle="fine-mapping-benchmark, branch application/ukb-cholesterol" \
  -V date="generated $(date '+%-d %B %Y') from application/README.md at $(git -C "$APP" rev-parse --short HEAD)"
echo "built: $(ls -la "$APP/README.pdf" | awk '{print $5" bytes"}')"
