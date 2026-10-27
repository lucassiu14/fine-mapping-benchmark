# Shared diagnostics for every stage. Sourced, not executed.
#
# Gives each job: a timestamped log, an environment banner, an ERR trap that
# names the failing line and command, precondition helpers that fail loudly
# before expensive work, and output validation after it. Every abnormal line
# is prefixed FATAL or WARN so a whole array can be triaged with
#   grep -l FATAL logs/*
set -Eeuo pipefail

_T0=$(date +%s)
_STAGE="${PBS_JOBNAME:-$(basename "${BASH_SOURCE[1]:-stage}")}"
_TASK="${PBS_ARRAY_INDEX:-}"
_TAG="${_STAGE}${_TASK:+[$_TASK]}"

log()     { printf '[%s] %s  %s\n' "$(date +%H:%M:%S)" "$_TAG" "$*"; }
warn()    { printf '[%s] %s  WARN  %s\n' "$(date +%H:%M:%S)" "$_TAG" "$*" >&2; }
die()     { printf '[%s] %s  FATAL %s\n' "$(date +%H:%M:%S)" "$_TAG" "$*" >&2; exit 1; }
section() { printf '\n[%s] %s  ===== %s =====\n' "$(date +%H:%M:%S)" "$_TAG" "$*"; }

on_err() {
  local rc=$? line="$1" cmd="$2"
  printf '\n[%s] %s  FATAL exit %d at line %s\n' "$(date +%H:%M:%S)" "$_TAG" "$rc" "$line" >&2
  printf '[%s] %s  FATAL command: %s\n' "$(date +%H:%M:%S)" "$_TAG" "$cmd" >&2
  printf '[%s] %s  FATAL cwd: %s\n' "$(date +%H:%M:%S)" "$_TAG" "$PWD" >&2
  printf '[%s] %s  FATAL ran for %ds\n' "$(date +%H:%M:%S)" "$_TAG" "$(( $(date +%s) - _T0 ))" >&2
  exit "$rc"
}
trap 'on_err "$LINENO" "$BASH_COMMAND"' ERR

# Load cluster modules before anything else. The venvs are built on the login
# node against a particular python; if that python comes from a module, the
# compute node needs the same module loaded or the venv will not start. Set
# APP_MODULES in env.sh, e.g.
#   export APP_MODULES="tools/prod Python/3.11.2 R/4.3.2 aws-cli"
# Left unset, nothing is loaded and the system defaults are used.
load_modules() {
  [ -n "${APP_MODULES:-}" ] || return 0
  command -v module >/dev/null 2>&1 || { warn "APP_MODULES is set but 'module' is unavailable"; return 0; }
  for m in $APP_MODULES; do try_module "$m"; done
}

# Load one module, for stages that need something the others do not. Keeping
# such a module out of APP_MODULES matters when it is built against a
# different toolchain than the rest: EasyBuild will often refuse to load
# awscli/...-GCCcore-13.2.0 alongside Python/...-GCCcore-12.3.0, and loading
# it only where it is needed sidesteps the conflict entirely.
try_module() {
  command -v module >/dev/null 2>&1 || { warn "cannot load $1: no 'module' command"; return 0; }
  if module load "$1" 2>/dev/null; then log "module   : loaded $1"; else warn "module load $1 failed"; fi
}

banner() {
  section "start"
  load_modules
  log "job      : ${PBS_JOBID:-interactive}${_TASK:+  task $_TASK}"
  log "host     : $(hostname)  ($(nproc 2>/dev/null || echo '?') cpus)"
  log "cwd      : $PWD"
  for v in APP_ROOT N_LIKELIHOOD ANNOT_SET APP_MODULES FB_DIR FUNMAP_DIR POLYFUN_DIR \
           APP_VENV FUNMAP_VENV POLYFUN_ENV POLYFUN_RLIB R_HOME; do
    [ -n "${!v:-}" ] && log "$(printf '%-9s' "$v"): ${!v}"
  done
  local avail
  avail=$(df -Pk . 2>/dev/null | awk 'NR==2{printf "%.0f GB", $4/1048576}')
  log "disk free: ${avail:-?}"
}

# Report python/R provenance AFTER the environment is activated, so the log
# records what actually ran rather than what was meant to.
report_python() {
  log "python   : $(command -v python) ($(python -V 2>&1))"
  python - <<'PY' 2>/dev/null || warn "could not report package versions"
import importlib
out=[]
for m in ("numpy","scipy","pandas","pyarrow","torch","bitarray","rpy2","sklearn"):
    try: out.append(f"{m} {importlib.import_module(m).__version__}")
    except Exception: pass
print("            packages : " + ", ".join(out))
PY
}

require_cmd()  { command -v "$1" >/dev/null || die "command not found: $1"; }
require_dir()  { [ -d "$1" ] || die "missing directory: $1${2:+  ($2)}"; }
require_file() {
  [ -e "$1" ] || die "missing file: $1${3:+  ($3)}"
  local min="${2:-1}" sz
  sz=$(wc -c < "$1" 2>/dev/null || echo 0)
  [ "$sz" -ge "$min" ] || die "file too small (${sz}B < ${min}B), likely truncated: $1"
}
require_env()  { [ -n "${!1:-}" ] || die "environment variable $1 is not set${2:+  ($2)}"; }

# Validate an output before declaring success, so a stage cannot report OK
# while having written a zero-row or truncated file.
check_output() {
  local f="$1" min="${2:-1000}"
  [ -e "$f" ] || die "stage produced no output: $f"
  local sz; sz=$(wc -c < "$f")
  [ "$sz" -ge "$min" ] || die "output suspiciously small (${sz}B): $f"
  log "wrote    : $f ($(du -h "$f" | cut -f1))"
}

require_fb_fixes() {
  # --- guard: refuse to run against a checkout missing the memo fix ----------
  # These fixes landed on main in 7a46a13. Any checkout older than that silently
  # produces diluted PIPs: reformat_memo would subtract the MEAN
  # log-BF instead of the max, saturating 27-53% of the memo at the 1e15 clamp
  # and turning calculate_pip into a count of visited configurations. Max PIP
  # collapses to ~0.02-0.05 and no credible set is ever formed. It fails
  # quietly, looks like a modelling result, and is the single most expensive
  # mistake available here - so check for it rather than trusting the checkout.
  : "${FB_DIR:?set FB_DIR to the functional-beatrice checkout}"
  grep -q 'm0 = max(' "${FB_DIR}/scripts/trainer_annot.py" || {
    echo "FATAL: ${FB_DIR} is missing the reformat_memo max-offset fix." >&2
    echo "       Update the checkout to main at d0fb022 or later." >&2
    exit 1
  }
  # Search the whole scripts dir, not trainer_annot.py: regularize_ld existed
  # as two byte-identical copies (trainer_annot and generate_credible_sets) and
  # now lives once in scripts/ld_regularisation.py. Pinning this check to one
  # file made it fail on a checkout that HAS the fix.
  grep -rq 'nearest PD correlation matrix' "${FB_DIR}/scripts/" || {
    echo "FATAL: ${FB_DIR} is missing the LD regularisation fix." >&2
    echo "       The diagonal shift scales every implied correlation by 1/(1+c)," >&2
    echo "       which makes the Bayes factor FAVOUR co-including collinear variants." >&2
    echo "       Update the checkout to main at d0fb022 or later." >&2
    exit 1
  }
  # The adaptive eigenvalue floor. A fixed 1e-3 clip lets near-singular
  # directions amplify z-noise ~1000x: on ldEUR_08_0045 a 36-variant
  # configuration scored a quadratic term of 745.88 against 85.12 for the true
  # lead at |Z|=13.05, and took PIP 1.0 on a variant at |Z|=0.04.
  grep -q 'eig_floor' "${FB_DIR}/scripts/ld_regularisation.py" 2>/dev/null || {
    echo "FATAL: ${FB_DIR} is missing the adaptive LD eigenvalue floor." >&2
    echo "       A fixed 1e-3 clip manufactures evidence in near-singular" >&2
    echo "       directions. Update to fc82a9c or later." >&2
    exit 1
  }
  # The learned prior must reach the credible-set sweep. prior_location is ''
  # on every path, so generate_credible_sets' bare except silently substituted
  # a UNIFORM prior and the annotation prior never reached a reported PIP.
  grep -q 'p0_t' "${FB_DIR}/scripts/generate_credible_sets.py" \
    && grep -q "'p0_t'" "${FB_DIR}/scripts/joint_trainer.py" || {
    echo "FATAL: ${FB_DIR} does not pass the learned prior to the sweep." >&2
    echo "       Without it the reported PIPs are computed under a UNIFORM" >&2
    echo "       prior and the annotations do not enter them." >&2
    echo "       Update to fc82a9c or later." >&2
    exit 1
  }
  grep -q 'def read_ld' "${FB_DIR}/scripts/joint_trainer.py" || {
    echo "FATAL: ${FB_DIR} is missing .npy LD support; the manifest points at .npy files." >&2
    exit 1
  }
  echo "[guard] fb checkout has: memo fix, LD regularisation + adaptive eigenvalue floor, learned prior passed to the sweep, .npy LD support"
}

rpy2_env() {
  # rpy2 embeds R IN-PROCESS, and that embedded R does not inherit a library
  # the R module did not install. Measured: its .libPaths() is only
  #   /rds/easybuild/rome/apps/software/R/4.3.2-gfbf-2023a/lib64/R/library
  # so anything installed into POLYFUN_RLIB is invisible to it unless R_LIBS
  # says otherwise. This is easy to get wrong because an
  #   Rscript -e ".libPaths(c('$POLYFUN_RLIB', .libPaths())); library(X)"
  # test PASSES while rpy2 still fails - the two do not see the same
  # libraries. 00g reached PolyFun's binning step after 20 minutes of ridge
  # regression and died with
  #   PackageNotInstalledError: The R package "Ckmeans.1d.dp" is not installed
  # while that very package loaded fine under Rscript.
  export R_HOME="${R_HOME:?set R_HOME to the R installation rpy2 should embed}"
  # Prepend only if absent: env.sh exports R_LIBS too, and prepending
  # unconditionally yields "<rlib>:<rlib>" in the log. R dedupes it, so this
  # is cosmetic - but a doubled path in a log is a false lead waiting to be
  # chased.
  : "${POLYFUN_RLIB:?set POLYFUN_RLIB}"
  case ":${R_LIBS:-}:" in
    *":$POLYFUN_RLIB:"*) ;;
    *) export R_LIBS="$POLYFUN_RLIB${R_LIBS:+:$R_LIBS}" ;;
  esac
  # Only needed when the rpy2 wheel was built against a different R minor
  # version than the one installed; harmless otherwise.
  export RPY2_CFFI_MODE="${RPY2_CFFI_MODE:-ABI}"
  log "R_HOME   : $R_HOME"
  log "R_LIBS   : $R_LIBS"
}

require_r_package() {
  # Verify through rpy2's own importr, exactly as PolyFun calls it. Checking
  # with Rscript instead is what let this through the first time.
  python - "$1" <<'RPY' || die "the embedded R cannot load '$1' - check R_LIBS above"
import sys
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
name = sys.argv[1]
print("            embedded R .libPaths():")
for p in ro.r('.libPaths()'):
    print(f"              {p}")
importr(name)
print(f"            '{name}' loads under rpy2", flush=True)
RPY
}

finish() {
  local el=$(( $(date +%s) - _T0 ))
  section "done in ${el}s ($((el/60)) min)"
}
