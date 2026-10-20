#!/bin/bash
# Submit the whole pipeline with PBS dependencies, so nothing starts before
# its inputs exist. Run `source env.sh` first.
#
#   bash submit_all.sh              # submit everything
#   bash submit_all.sh --dry-run    # print the qsub commands and stop
#   bash submit_all.sh --from 03    # resume from a stage (deps are dropped)
#
# Every stage is idempotent, so resubmitting after a partial failure re-runs
# only what is missing.
set -euo pipefail
cd "${APP_ROOT:?source env.sh first}"

DRY=0; FROM=""
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) DRY=1 ;;
    --from) FROM="$2"; shift ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac; shift
done

PASS=APP_ROOT,N_LIKELIHOOD,APP_VENV,FB_DIR,FUNMAP_VENV,FUNMAP_DIR,POLYFUN_ENV,POLYFUN_DIR,POLYFUN_RLIB,R_HOME,RPY2_CFFI_MODE

sub () {  # sub <stage> <script> [dependency-expression] [VAR=value ...]
  local stage="$1" script="$2" dep="${3:-}"; shift 3 2>/dev/null || shift $#
  local extra=""
  for kv in "$@"; do extra="${extra},${kv}"; done
  if [ -n "$FROM" ] && [ "$stage" \< "$FROM" ]; then echo "skip $script"; return; fi
  if [ -n "$FROM" ] && [ "$stage" = "$FROM" ]; then dep=""; fi   # resuming: no deps
  local cmd=(qsub -v "${PASS}${extra}")
  [ -n "$dep" ] && cmd+=(-W "depend=$dep")
  cmd+=("$script")
  if [ "$DRY" = 1 ]; then
    # the command goes to stderr; ONLY the job id goes to stdout, because the
    # caller captures stdout into a variable it then builds dependencies from
    { echo; echo "  \$ ${cmd[*]}"; } >&2
    echo "dry-$stage"
    return
  fi
  "${cmd[@]}"
}

echo "== downloads (independent, run concurrently) =="
J_ANNOT=$(sub 00b 00b_download_annot.pbs)   ; echo "  annotations : $J_ANNOT  (Borzoi + baselineLF)"
J_LD=$(   sub 01  01_download_ld.pbs)       ; echo "  LD          : $J_LD"
J_REF=$(  sub 00d 00d_ldsc_reference.pbs)   ; echo "  LDSC ref    : $J_REF"

echo "== annotation sets =="
J_COMB=$(sub 00b2 00b2_combine_annot.pbs "afterokarray:${J_ANNOT}[]")
echo "  combined    : $J_COMB  (Borzoi-102 + baselineLF = 289)"

echo "== PolyFun prior, computed from OUR annotations, once per set =="
J_MUNGE=$(sub 00e 00e_munge.pbs "afterok:${J_REF}")
echo "  munge       : $J_MUNGE"
J_LDSC_S=$(sub 00f 00f_ldscores.pbs "afterok:${J_MUNGE},afterokarray:${J_ANNOT}[]" ANNOT_SET=sniff)
echo "  LD-scores   : $J_LDSC_S  (sniff)"
J_LDSC_C=$(sub 00f 00f_ldscores.pbs "afterok:${J_MUNGE},afterokarray:${J_COMB}[]"  ANNOT_SET=combined)
echo "  LD-scores   : $J_LDSC_C  (combined)"
J_PRIOR_S=$(sub 00g 00g_polyfun_prior.pbs "afterokarray:${J_LDSC_S}[]" ANNOT_SET=sniff)
echo "  prior       : $J_PRIOR_S  (sniff)"
J_PRIOR_C=$(sub 00g 00g_polyfun_prior.pbs "afterokarray:${J_LDSC_C}[]" ANNOT_SET=combined)
echo "  prior       : $J_PRIOR_C  (combined)"

echo "== prep (needs annotations + LD) =="
J_PREP=$(sub 02 02_prep.pbs "afterokarray:${J_ANNOT}[],afterokarray:${J_LD}[]")
echo "  prep        : $J_PREP"

echo "== methods (independent of each other) =="
J_FB=$(sub 03 03_fb_xregion.pbs "afterokarray:${J_PREP}[]")
echo "  fb_xregion  : $J_FB"
J_FM=$(sub 04 04_funmap.pbs     "afterokarray:${J_PREP}[]")
echo "  funmap      : $J_FM"
J_PF_S=$(sub 05 05_polyfun.pbs "afterokarray:${J_PREP}[],afterok:${J_PRIOR_S}" ANNOT_SET=sniff)
echo "  polyfun     : $J_PF_S  (sniff, 102 annotations)"
J_PF_C=$(sub 05 05_polyfun.pbs "afterokarray:${J_PREP}[],afterok:${J_PRIOR_C}" ANNOT_SET=combined)
echo "  polyfun     : $J_PF_C  (combined, 289 annotations)"

cat <<NOTE

Submitted. Useful checks while it runs:

  python preflight.py data      # before prep: downloads complete?
  python preflight.py prep      # before the methods: all 265 blocks built?
  python preflight.py methods   # progress per method
  qstat -tu \$USER

06_collect.py is deliberately NOT chained - run it once the methods finish,
so a partial failure in one of them does not produce a half-empty table:

  python 06_collect.py && python preflight.py collect
NOTE
