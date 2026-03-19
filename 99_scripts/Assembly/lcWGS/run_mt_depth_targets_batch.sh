#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage:
  $(basename "$0") -i <fastq_dir> -o <out_root>
                   [-p <lcWGS_mt_sweep.py>] [-t <threads>]
                   [--reads 50000,100000,...] [--depths 5,10,20]
                   [--skip 16-CI1f,OTHER] [--dry-run] [--force]

Defaults:
  -p ./lcWGS_mt_sweep.py
  -t 16
  --reads 50000,100000,250000,500000,1000000,2000000
  --depths 5,10,20,30
  --skip 16-CI1f
EOF
}

INDIR=""
OUTROOT=""
PIPELINE="./lcWGS_mt_sweep.py"
THREADS=16
MT_SWEEP_READS="50000,100000,250000,500000,1000000,2000000"
DEPTHS_CSV="5,10,20,30"
SKIP_CSV="16-CI1f"
DRYRUN=0
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)    INDIR="$2"; shift 2;;
    -o|--out)      OUTROOT="$2"; shift 2;;
    -p|--pipeline) PIPELINE="$2"; shift 2;;
    -t|--threads)  THREADS="$2"; shift 2;;
    --reads)       MT_SWEEP_READS="$2"; shift 2;;
    --depths)      DEPTHS_CSV="$2"; shift 2;;
    --skip)        SKIP_CSV="$2"; shift 2;;
    --dry-run)     DRYRUN=1; shift;;
    --force)       FORCE=1; shift;;
    -h|--help)     usage; exit 0;;
    *) echo "ERROR: unknown option: $1" >&2; usage; exit 2;;
  esac
done

[[ -n "$INDIR" && -n "$OUTROOT" ]] || { echo "ERROR: -i and -o are required" >&2; usage; exit 2; }
[[ -f "$PIPELINE" ]] || { echo "ERROR: pipeline not found: $PIPELINE" >&2; exit 2; }

IFS=',' read -r -a DEPTHS <<< "$DEPTHS_CSV"
IFS=',' read -r -a SKIPS  <<< "$SKIP_CSV"

is_skipped() {
  local s="$1"
  for x in "${SKIPS[@]}"; do [[ -n "$x" && "$s" == "$x" ]] && return 0; done
  return 1
}

mkdir -p "$OUTROOT"/{runs,logs,reports}

ALL_TSV="$OUTROOT/reports/all_runs.tsv"
PASS_TSV="$OUTROOT/reports/selected_min_pass.tsv"
FAIL_TSV="$OUTROOT/reports/failed_runs.tsv"

init_reports() {
  [[ -f "$ALL_TSV"  ]] || echo -e "sample\tdepth_target\trun_dir\tstatus\tselected_sweep_reads\tmt_primary_len\tmean_depth\tpct_sites_ge_target\tlog" > "$ALL_TSV"
  [[ -f "$PASS_TSV" ]] || echo -e "sample\tdepth_target\trun_dir\tselected_sweep_reads\tmt_primary_len\tmean_depth\tpct_sites_ge_target\tlog" > "$PASS_TSV"
  [[ -f "$FAIL_TSV" ]] || echo -e "sample\tdepth_target\trun_dir\texit_code\tlog" > "$FAIL_TSV"
}
init_reports

# Find R1 files: accept *_1.fastq.gz (your case). (easy to extend later)
shopt -s nullglob
R1S=( "$INDIR"/*_1.fastq.gz )
if ((${#R1S[@]}==0)); then
  echo "ERROR: no R1 files found like *_1.fastq.gz in: $INDIR" >&2
  exit 2
fi

for r1 in "${R1S[@]}"; do
  base=$(basename "$r1")
  sample="${base%_1.fastq.gz}"
  r2="$INDIR/${sample}_2.fastq.gz"
  [[ -f "$r2" ]] || { echo "WARN: missing R2 for $sample ($r2); skipping" >&2; continue; }

  if is_skipped "$sample"; then
    echo "SKIP sample=$sample (in skip list)"
    continue
  fi

  for depth in "${DEPTHS[@]}"; do
    tag="mtDepth_${depth}x"
    sample_root="$OUTROOT/runs/$sample/$tag"
    doneflag="$sample_root/DONE.ok"
    logdir="$OUTROOT/logs/$sample/$tag"
    mkdir -p "$sample_root" "$logdir"

    if [[ $FORCE -eq 0 && -f "$doneflag" ]]; then
      echo "SKIP sample=$sample depth=${depth}x (already DONE: $doneflag)"
      continue
    fi

    stamp="$(date +%Y%m%d_%H%M%S_%N)"
    outdir="$sample_root/run_${stamp}"
    logfile="$logdir/run_${stamp}.log"

    cmd=( python -u "$PIPELINE"
          -1 "$r1" -2 "$r2"
          -o "$outdir"
          -t "$THREADS"
          --mt_sweep_reads "$MT_SWEEP_READS"
          --mt_depth_target "$depth"
          --stop_after_pass
          --mtonly
          --verbose )

    echo "==> [$sample] depth=${depth}x"
    echo "    OUT: $outdir"
    echo "    LOG: $logfile"
    echo "    CMD: ${cmd[*]}"

    if [[ $DRYRUN -eq 1 ]]; then
      continue
    fi

    set +e
    "${cmd[@]}" |& tee "$logfile"
    rc=${PIPESTATUS[0]}
    set -e

    if [[ $rc -ne 0 ]]; then
      echo -e "${sample}\t${depth}\t${outdir}\t${rc}\t${logfile}" >> "$FAIL_TSV"
      echo "FAIL sample=$sample depth=${depth}x (exit=$rc). See: $logfile" >&2
      continue
    fi

    # mark done + store which run dir is the "done" one
    echo "$outdir" > "$doneflag"

    # parse selected min pass line from mt_sweep_summary.tsv (if present)
    sel_reads="NA"; mtlen="NA"; meand="NA"; pct="NA"; status="OK"
    summ="$outdir/mt_sweep_summary.tsv"
    if [[ -f "$summ" ]]; then
      read -r sel_reads mtlen meand pct < <(
        awk -F'\t' 'NR==1{next} $4=="YES"{print $1, $6, $9, $10; exit}' "$summ"
      )
      [[ -n "${sel_reads:-}" ]] || sel_reads="NA"
      [[ -n "${mtlen:-}" ]]     || mtlen="NA"
      [[ -n "${meand:-}" ]]     || meand="NA"
      [[ -n "${pct:-}" ]]       || pct="NA"
    else
      status="OK_NO_SUMMARY"
    fi

    echo -e "${sample}\t${depth}\t${outdir}\t${status}\t${sel_reads}\t${mtlen}\t${meand}\t${pct}\t${logfile}" >> "$ALL_TSV"
    echo -e "${sample}\t${depth}\t${outdir}\t${sel_reads}\t${mtlen}\t${meand}\t${pct}\t${logfile}" >> "$PASS_TSV"

    echo "DONE sample=$sample depth=${depth}x -> $outdir"
  done
done

echo "Batch reports:"
echo "  ALL:  $ALL_TSV"
echo "  PASS: $PASS_TSV"
echo "  FAIL: $FAIL_TSV"
