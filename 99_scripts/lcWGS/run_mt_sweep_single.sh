#!/usr/bin/env bash
set -euo pipefail

PIPELINE="/BIGDATA1/r_projects/mattia.baricordi2/lcWGS/lcWGS_mt_sweep.py"

R1=""; R2=""; SAMPLE=""; ROOT=""
THREADS=16
SWEEP_READS="50000,100000,250000,500000,1000000,2000000"
DEPTH_TARGET=10
ERASE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -1|--r1) R1="$2"; shift 2;;
    -2|--r2) R2="$2"; shift 2;;
    -s|--sample) SAMPLE="$2"; shift 2;;
    -r|--root) ROOT="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    --reads) SWEEP_READS="$2"; shift 2;;
    --depth) DEPTH_TARGET="$2"; shift 2;;
    --erase) ERASE=1; shift 1;;
    *) echo "Unknown arg: $1" >&2; exit 2;;
  esac
done

[[ -n "$R1" && -n "$R2" && -n "$SAMPLE" && -n "$ROOT" ]] || {
  echo "Usage: $0 -1 R1 -2 R2 -s SAMPLE -r ROOT [-t 16] [--reads list] [--depth 10] [--erase]" >&2
  exit 2
}

RUNID="$(date +%Y%m%d_%H%M%S_%N)"
OUT="$ROOT/runs/$SAMPLE/$RUNID"
LOGDIR="$ROOT/logs/$SAMPLE"
LOG="$LOGDIR/$RUNID.log"

mkdir -p "$ROOT/runs/$SAMPLE" "$LOGDIR" "$ROOT/reports"

cmd=(python -u "$PIPELINE"
  -1 "$R1" -2 "$R2" -o "$OUT" -t "$THREADS"
  --mt_sweep_reads "$SWEEP_READS"
  --mt_depth_target "$DEPTH_TARGET"
  --stop_after_pass
  --mtonly
  --verbose
)
if [[ $ERASE -eq 1 ]]; then cmd+=(--erase); fi

echo "OUT=$OUT"
echo "LOG=$LOG"
printf 'CMD:'; printf ' %q' "${cmd[@]}"; echo
"${cmd[@]}" |& tee "$LOG"
