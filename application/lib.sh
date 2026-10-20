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

finish() {
  local el=$(( $(date +%s) - _T0 ))
  section "done in ${el}s ($((el/60)) min)"
}
