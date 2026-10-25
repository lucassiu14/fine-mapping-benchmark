#!/bin/bash
# Submit the whole pipeline with PBS dependencies, so nothing starts before
# its inputs exist. Run `source env.sh` first.
#
#   bash submit_all.sh              # submit everything
#   bash submit_all.sh --dry-run    # print the qsub commands and stop
#   bash submit_all.sh --from 03    # resume from a stage (deps are dropped)
#   bash submit_all.sh --prior-only # ONLY the PolyFun-prior chain
#   bash submit_all.sh --prep-only  # everything up to and including 02 prep
#
# --prior-only submits 00b -> 00b2 -> 00d -> 00h1 -> 00h2 -> 00e -> 00f -> 00g
# and stops: the
# S-LDSC chain that produces the prior, which is the part of the pipeline that
# has never been exercised end to end. It deliberately leaves out 01 (the LD
# download, ~429 GB peak), 02 and the methods, so the chain can be proven for
# the cost of the ~20 GB annotation download alone. Re-run without the flag
# afterwards: every stage is idempotent and the finished ones are skipped.
#
# Every stage is idempotent, so resubmitting after a partial failure re-runs
# only what is missing.
set -euo pipefail
cd "${APP_ROOT:?source env.sh first}"

DRY=0; FROM=""; PRIOR_ONLY=0; PREP_ONLY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) DRY=1 ;;
    --prior-only) PRIOR_ONLY=1 ;;
    --prep-only) PREP_ONLY=1 ;;
    --from) FROM="$2"; shift ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac; shift
done

# Variables the jobs need. Build the -v list from the ones that are ACTUALLY
# SET: naming an unset variable makes qsub refuse the whole submission with
#   qsub: cannot send environment with the job
# which is what happened when PF_NUM_BINS joined this list - it is optional
# and only set when pinning the bin count, so an ordinary run has it unset.
PASS_VARS="APP_ROOT APP_MODULES AWS_MODULE N_LIKELIHOOD APP_VENV FB_DIR
           FUNMAP_VENV FUNMAP_DIR POLYFUN_ENV POLYFUN_DIR POLYFUN_RLIB
           R_HOME R_LIBS RPY2_CFFI_MODE PF_NUM_BINS"
PASS=""; PASS_SKIPPED=""
for v in $PASS_VARS; do
  if [ -n "${!v:-}" ]; then PASS="${PASS:+$PASS,}$v"
  else PASS_SKIPPED="${PASS_SKIPPED:+$PASS_SKIPPED }$v"; fi
done
[ -n "$PASS" ] || { echo "none of $PASS_VARS is set - source env.sh first" >&2; exit 2; }
[ -n "$PASS_SKIPPED" ] && echo "note: not set, so not passed: $PASS_SKIPPED"

dep_arr () {
  # Dependency on an ARRAY job. MEASURED on this scheduler with
  # probe_depend.sh (PBS Pro pbs_version = 2024.1.1.20240608133507):
  #
  #   afterokarray:4424651[].pbs-7     rejected: illegal -W value
  #   afterokarray:4424651.pbs-7[]     rejected
  #   afterokarray:4424651[]           rejected
  #   afterokarray:4424651.pbs-7       rejected
  #   afterokarray:4424651             rejected
  #   afterok:4424651[].pbs-7          ACCEPTED
  #
  # Every afterokarray spelling is refused, so this build does not accept
  # that dependency type; afterok with the id exactly as qsub returns it
  # does. Do NOT change this back to afterokarray without re-running
  # probe_depend.sh on the scheduler you are submitting to.
  #
  # Caveat, not measured: whether afterok on an array parent requires every
  # subjob to have succeeded, or only the parent to have finished. The
  # stages do not rely on it either way - each one verifies its own inputs
  # (00f checks its annotation file exists, 00g requires all 22
  # snpvar_constrained files) and dies rather than running on partial data.
  printf 'afterok:%s' "$1"
}

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
    # Mimic the id shape PBS Pro actually returns, brackets and all, so a
    # dry run exercises the dependency expressions the same way a real run does.
    if grep -qE '^#PBS +-J ' "$script"; then echo "dry-$stage[].pbs-0"
    else echo "dry-$stage.pbs-0"; fi
    return
  fi
  "${cmd[@]}"
}

echo "== downloads (independent, run concurrently) =="
J_ANNOT=$(sub 00b 00b_download_annot.pbs)   ; echo "  annotations : $J_ANNOT  (Borzoi + baselineLF)"
if [ "$PRIOR_ONLY" = 1 ]; then
  echo "  LD          : SKIPPED (--prior-only)"
else
  J_LD=$(sub 01 01_download_ld.pbs)         ; echo "  LD          : $J_LD"
fi
J_REF=$(  sub 00d 00d_ldsc_reference.pbs)   ; echo "  LDSC ref    : $J_REF"

echo "== annotation sets =="
J_COMB=$(sub 00b2 00b2_combine_annot.pbs "$(dep_arr "$J_ANNOT")")
echo "  combined    : $J_COMB  (Borzoi-102 + baselineLF = 289)"

echo "== PolyFun annotations: base + quintile indicators per Borzoi PC =="
# The raw PCs are signed, so their column sums (= M, which S-LDSC reads as a
# per-annotation SNP count) came out negative, and neither source set carries
# a base column. 00h1 computes genome-wide quintile breakpoints once; 00h2
# builds base + 4 indicators per PC per chromosome. 00f then computes the
# LD-scores FROM those built files.
J_BP=$(sub 00h1 00h1_annot_breakpoints.pbs "$(dep_arr "$J_ANNOT")")
echo "  breakpoints : $J_BP  (genome-wide, shared by both sets)"
J_BLD_S=$(sub 00h2 00h2_build_annot.pbs "afterok:${J_BP},$(dep_arr "$J_ANNOT")" ANNOT_SET=sniff)
echo "  build annot : $J_BLD_S  (sniff, 409)"
J_BLD_C=$(sub 00h2 00h2_build_annot.pbs "afterok:${J_BP},$(dep_arr "$J_COMB")"  ANNOT_SET=combined)
echo "  build annot : $J_BLD_C  (combined, 596)"

echo "== PolyFun prior, computed from OUR annotations, once per set =="
J_MUNGE=$(sub 00e 00e_munge.pbs "afterok:${J_REF}")
echo "  munge       : $J_MUNGE"
J_LDSC_S=$(sub 00f 00f_ldscores.pbs "afterok:${J_MUNGE},$(dep_arr "$J_BLD_S")" ANNOT_SET=sniff)
echo "  LD-scores   : $J_LDSC_S  (sniff)"
J_LDSC_C=$(sub 00f 00f_ldscores.pbs "afterok:${J_MUNGE},$(dep_arr "$J_BLD_C")" ANNOT_SET=combined)
echo "  LD-scores   : $J_LDSC_C  (combined)"
J_PRIOR_S=$(sub 00g 00g_polyfun_prior.pbs "$(dep_arr "$J_LDSC_S")" ANNOT_SET=sniff)
echo "  prior       : $J_PRIOR_S  (sniff)"
J_PRIOR_C=$(sub 00g 00g_polyfun_prior.pbs "$(dep_arr "$J_LDSC_C")" ANNOT_SET=combined)
echo "  prior       : $J_PRIOR_C  (combined)"

if [ "$PRIOR_ONLY" = 1 ]; then
  cat <<'PNOTE'

--prior-only: stopping after the PolyFun-prior chain
(00b -> 00b2 -> 00d -> 00h1 -> 00h2 -> 00e -> 00f -> 00g). Not submitted:
01 (LD download), 02 (prep), 03/04/05 (methods).

Watch it with:

  qstat -tu $USER

PBS names each log by JOB ID, not job name - logs/<jobid>.pbs-N.OU - so a
glob on the job name matches nothing. To follow the newest logs:

  tail -f $(ls -t $APP_ROOT/logs/*.OU | head -2)

or, for one job you have the id of:

  tail -f $APP_ROOT/logs/<jobid>.*.OU

The prior chain is proven when BOTH sets have 22 chromosomes:

  for S in sniff combined; do
    echo "$S: $(ls $APP_ROOT/pf_prior_$S/chol.*.snpvar_constrained.gz 2>/dev/null | wc -l) of 22"
  done

Then submit the rest - finished stages are skipped:

  bash submit_all.sh
PNOTE
  exit 0
fi

echo "== prep (needs annotations + LD) =="
J_PREP=$(sub 02 02_prep.pbs "$(dep_arr "$J_ANNOT"),$(dep_arr "$J_LD")")
echo "  prep        : $J_PREP"

if [ "$PREP_ONLY" = 1 ]; then
  cat <<'PNOTE'

--prep-only: stopping after 02. Not submitted: 03/04/05 (the methods).

03_fb_xregion is a SINGLE job over all blocks (the prior head is shared, so
it cannot be an array) with a 48 h request that has never been measured. Time
it on a stratified subset first, which takes about 15 minutes:

  qsub -v APP_ROOT,APP_MODULES,APP_VENV,FB_DIR,N_LIKELIHOOD 03t_fb_timing.pbs

It reports two extrapolations to the full block set and a suggested walltime.
Then submit the methods:

  bash submit_all.sh
PNOTE
  exit 0
fi

echo "== methods (independent of each other) =="
J_FB=$(sub 03 03_fb_xregion.pbs "$(dep_arr "$J_PREP")")
echo "  fb_xregion  : $J_FB"
J_FM=$(sub 04 04_funmap.pbs     "$(dep_arr "$J_PREP")")
echo "  funmap      : $J_FM"
J_PF_S=$(sub 05 05_polyfun.pbs "$(dep_arr "$J_PREP"),afterok:${J_PRIOR_S}" ANNOT_SET=sniff)
echo "  polyfun     : $J_PF_S  (sniff, 102 annotations)"
J_PF_C=$(sub 05 05_polyfun.pbs "$(dep_arr "$J_PREP"),afterok:${J_PRIOR_C}" ANNOT_SET=combined)
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
