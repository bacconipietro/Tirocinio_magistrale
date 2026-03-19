#!/usr/bin/env python3
"""
lcWGS Pipeline (MT Sweep Edition)

Implements:
- randomized downsampling per sweep level (seqtk sample, deterministic via seed)
- MitoZ validation
- Bowtie2 -> samtools depth mean depth calculation
- automatic selection of the minimum sweep level that yields >= target mean depth (default 10x)
- TSV summary of all sweep runs

Notes:
- For paired-end, --mt_sweep_reads values are interpreted as READ PAIRS (records sampled from R1 and R2).
- Depth is computed on the LONGEST MitoZ-validated contig ("primary" MT contig).
- Uses SPAdes header k-mer coverage only as an initial filter; acceptance uses samtools depth (mean depth).
- BUSCO default lineage is insecta_odb10 (changeable via --busco_lineage).

Requires in PATH:
  trimmomatic, fastqc, spades.py, mitoz, bowtie2-build, bowtie2, samtools, seqtk
"""

import os
import re
import csv
import shutil
import datetime
import argparse
import subprocess
from os import path
from typing import List, Tuple, Optional, Dict

from Bio import SeqIO


# ------------------------- Helpers -------------------------

def ts() -> str:
    return datetime.datetime.now().strftime("%H:%M:%S")


def ensure_dir(p: str) -> None:
    if p and (not path.isdir(p)):
        os.makedirs(p)


def file_nonempty(p: str) -> bool:
    return path.exists(p) and path.getsize(p) > 0


def which_or_die(exe: str) -> str:
    pth = shutil.which(exe)
    if not pth:
        raise SystemExit(f"\nERROR: Required executable not found in PATH: {exe}\n")
    return pth


def run_logged(cmd: List[str], log_path: str, cwd: Optional[str] = None) -> None:
    """Run a command, redirecting stdout+stderr to a log file. Fail fast on error."""
    ensure_dir(path.dirname(log_path))
    with open(log_path, "w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=log, cwd=cwd)
    if proc.returncode != 0:
        raise RuntimeError( 
            f"\nERROR: Command failed (exit code {proc.returncode}).\n"
            f"See log: {log_path}\n"
            f"CMD: {' '.join(cmd)}\n"
        )


def run_logged_stdout_to_file(cmd: List[str], out_path: str, log_path: str, cwd: Optional[str] = None) -> None:
    """Run a command with stdout to a file and stderr to a log file."""
    ensure_dir(path.dirname(out_path))
    ensure_dir(path.dirname(log_path))
    with open(out_path, "w") as out_f, open(log_path, "w") as log:
        proc = subprocess.run(cmd, stdout=out_f, stderr=log, cwd=cwd)
    if proc.returncode != 0:
        raise RuntimeError(
            f"\nERROR: Command failed (exit code {proc.returncode}).\n"
            f"See log: {log_path}\n"
            f"CMD: {' '.join(cmd)}\n"
        )


def safe_link_or_copy(src: str, dst: str) -> None:
    """Prefer hardlink, then symlink, else copy."""
    if path.exists(dst):
        os.remove(dst)
    try:
        os.link(src, dst)
        return
    except Exception:
        pass
    try:
        os.symlink(src, dst)
        return
    except Exception:
        pass
    shutil.copy2(src, dst)


def parse_spades_id(record_id: str) -> Tuple[Optional[int], Optional[float]]:
    """Parse SPAdes header: NODE_1_length_12345_cov_56.7"""
    m = re.search(r"length_(\d+)_cov_([0-9.]+)", record_id)
    if not m:
        return None, None
    return int(m.group(1)), float(m.group(2))


def find_trimmomatic_adapter(filename: str) -> Optional[str]:
    """Auto-detect adapter path in Conda environment."""
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        return None
    share_dir = path.join(conda_prefix, "share")
    if not path.isdir(share_dir):
        return None

    for entry in os.listdir(share_dir):
        if entry.startswith("trimmomatic"):
            pth = path.join(share_dir, entry, "adapters", filename)
            if path.exists(pth):
                return pth
    return None


def parse_int_list(csv_str: str) -> List[int]:
    """Parse comma-separated integers (allow whitespace)."""
    vals = []
    for x in csv_str.split(","):
        x = x.strip()
        if not x:
            continue
        try:
            vals.append(int(x))
        except ValueError:
            raise SystemExit(f"\nERROR: Invalid integer in list: '{x}'\n")
    if not vals:
        raise SystemExit("\nERROR: No values parsed for --mt_sweep_reads.\n")
    return vals


def mitoz_valid_ids_from_summary(summary_path: str, candidate_ids: set) -> List[str]:
    """
    Conservative parser:
    - collects any first token that matches a candidate contig id
    This works across multiple MitoZ summary formats without relying on exact table markers.
    """
    valid = []
    seen = set()
    with open(summary_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tok = line.split()[0]
            if tok in candidate_ids and tok not in seen:
                valid.append(tok)
                seen.add(tok)
    return valid


def longest_contig_fasta(in_fa: str, out_fa: str) -> Tuple[str, int]:
    """Write the longest contig from in_fa to out_fa; return (id, length)."""
    longest = None
    longest_len = -1
    for rec in SeqIO.parse(in_fa, "fasta"):
        L = len(rec.seq)
        if L > longest_len:
            longest = rec
            longest_len = L
    if not longest:
        raise SystemExit("ERROR: No records found in validated MT fasta.")
    with open(out_fa, "w") as out_f:
        SeqIO.write(longest, out_f, "fasta")
    return longest.id, longest_len


def compute_depth_stats(depth_txt: str, ref_len: int, depth_target: float) -> Tuple[float, float]:
    """
    Compute mean depth and percent positions >= depth_target.
    Assumes samtools depth was run with -a; if missing sites occur, denominator remains ref_len.
    """
    if ref_len <= 0:
        return 0.0, 0.0
    if not path.exists(depth_txt) or path.getsize(depth_txt) == 0:
        return 0.0, 0.0

    depth_sum = 0
    ge_target = 0

    with open(depth_txt, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                d = int(parts[2])
            except ValueError:
                continue
            depth_sum += d
            if d >= depth_target:
                ge_target += 1

    mean_depth = float(depth_sum) / float(ref_len)
    pct_ge = 100.0 * float(ge_target) / float(ref_len)
    return mean_depth, pct_ge


# ------------------------- MT Sweep -------------------------

def run_mt_sweep_level(
    sweep_dir: str,
    paired_mode: bool,
    trim_r1: Optional[str],
    trim_r2: Optional[str],
    trim_se: Optional[str],
    reads_n: int,
    seed: int,
    threads: int,
    mtlen: float,
    mtcov_kmer: float,
    depth_target: float,
    sample_name: str,
) -> Dict[str, object]:
    """
    One sweep level:
    - downsample (seqtk sample)
    - SPAdes assemble
    - filter contigs by SPAdes header len/kmer-cov
    - MitoZ annotate + validate contig IDs
    - select longest validated contig (primary)
    - Bowtie2 map reads -> BAM -> samtools depth -> mean depth stats
    """
    ensure_dir(sweep_dir)
    logs_dir = path.join(sweep_dir, "logs")
    ensure_dir(logs_dir)

    result: Dict[str, object] = {
        "sweep_reads": reads_n,
        "status": "UNKNOWN",
        "mt_primary_id": "",
        "mt_primary_len": "",
        "mt_total_valid_len": "",
        "mt_num_valid_contigs": "",
        "mean_depth": "",
        "pct_sites_ge_target": "",
        "notes": "",
    }

    # ---- Downsample reads ----
    ds_r1 = path.join(sweep_dir, "ds_R1.fq")
    ds_r2 = path.join(sweep_dir, "ds_R2.fq")
    ds_se = path.join(sweep_dir, "ds_SE.fq")

    try:
        if reads_n == 0:
            # all reads
            if paired_mode:
                safe_link_or_copy(trim_r1, ds_r1)  # type: ignore
                safe_link_or_copy(trim_r2, ds_r2)  # type: ignore
            else:
                safe_link_or_copy(trim_se, ds_se)  # type: ignore
        else:
            which_or_die("seqtk")
            if paired_mode:
                # same seed + same N on matching-ordered R1/R2 preserves pairing in practice
                run_logged_stdout_to_file(
                    ["seqtk", "sample", "-s", str(seed), trim_r1, str(reads_n)],  # type: ignore
                    ds_r1,
                    path.join(logs_dir, "log_seqtk_R1.txt"),
                )
                run_logged_stdout_to_file(
                    ["seqtk", "sample", "-s", str(seed), trim_r2, str(reads_n)],  # type: ignore
                    ds_r2,
                    path.join(logs_dir, "log_seqtk_R2.txt"),
                )
            else:
                run_logged_stdout_to_file(
                    ["seqtk", "sample", "-s", str(seed), trim_se, str(reads_n)],  # type: ignore
                    ds_se,
                    path.join(logs_dir, "log_seqtk_SE.txt"),
                )
    except Exception as e:
        result["status"] = "FAILED_DOWNSAMPLE"
        result["notes"] = str(e)
        return result

    # ---- SPAdes MT assembly ----
    spades_out = path.join(sweep_dir, "spades_mt")
    try:
        spades_cmd = ["spades.py", "-o", spades_out, "-t", str(threads), "--cov-cutoff", "off"]
        if paired_mode:
            spades_cmd += ["-1", ds_r1, "-2", ds_r2]
        else:
            spades_cmd += ["-s", ds_se]
        run_logged(spades_cmd, path.join(logs_dir, "log_spades_mt.txt"))
    except Exception as e:
        result["status"] = "FAILED_SPADES"
        result["notes"] = str(e)
        return result

    spades_contigs = path.join(spades_out, "contigs.fasta")
    if not file_nonempty(spades_contigs):
        result["status"] = "FAILED_SPADES_NO_CONTIGS"
        result["notes"] = "SPAdes produced no contigs."
        return result

    # ---- Filter contigs by SPAdes header length/cov ----
    mt_candidate = path.join(sweep_dir, "mtcandidate.fa")
    kept = 0
    candidate_ids = set()
    try:
        with open(mt_candidate, "w") as out_f:
            for rec in SeqIO.parse(spades_contigs, "fasta"):
                L, C = parse_spades_id(rec.id)
                if L is None or C is None:
                    continue
                if float(L) > float(mtlen) and float(C) > float(mtcov_kmer):
                    rec.id = rec.id.split("_length")[0]  # NODE_#
                    rec.description = rec.id
                    SeqIO.write(rec, out_f, "fasta")
                    kept += 1
                    candidate_ids.add(rec.id)
    except Exception as e:
        result["status"] = "FAILED_FILTER"
        result["notes"] = str(e)
        return result

    if kept == 0:
        result["status"] = "FAILED_FILTER_NONE_PASSED"
        result["notes"] = f"No contigs passed filters (Len>{mtlen}, kmerCov>{mtcov_kmer})."
        return result

    # ---- MitoZ annotate + validate ----
    try:
        run_logged(
            [
                "mitoz", "annotate",
                "--fastafiles", mt_candidate,
                "--outprefix", "annot",
                "--thread_number", str(threads),
                "--species_name", sample_name,
            ],
            path.join(logs_dir, "log_mitoz.txt"),
            cwd=sweep_dir
        )
    except Exception as e:
        result["status"] = "FAILED_MITOZ"
        result["notes"] = str(e)
        return result

    mitoz_res_dir = path.join(sweep_dir, "annot.mtcandidate.fa.result")
    mitoz_summary = path.join(mitoz_res_dir, "summary.txt")
    if not path.exists(mitoz_summary):
        result["status"] = "FAILED_MITOZ_NO_SUMMARY"
        result["notes"] = "MitoZ summary.txt missing."
        return result

    valid_ids = mitoz_valid_ids_from_summary(mitoz_summary, candidate_ids)
    if not valid_ids:
        result["status"] = "FAILED_VALIDATION_NO_HITS"
        result["notes"] = "No validated contigs detected in MitoZ summary."
        return result

    mt_validated_all = path.join(sweep_dir, "mt_validated_all.fa")
    n_valid = 0
    valid_len_total = 0
    try:
        with open(mt_validated_all, "w") as out_f:
            for rec in SeqIO.parse(mt_candidate, "fasta"):
                if rec.id in valid_ids:
                    SeqIO.write(rec, out_f, "fasta")
                    n_valid += 1
                    valid_len_total += len(rec.seq)
    except Exception as e:
        result["status"] = "FAILED_VALIDATION_WRITE"
        result["notes"] = str(e)
        return result

    if n_valid == 0:
        result["status"] = "FAILED_VALIDATION_EMPTY"
        result["notes"] = "Validated contig list was non-empty but produced empty fasta."
        return result

    result["mt_num_valid_contigs"] = n_valid
    result["mt_total_valid_len"] = valid_len_total

    # ---- Primary = longest validated contig ----
    mt_primary = path.join(sweep_dir, "mt_primary.fa")
    try:
        primary_id, primary_len = longest_contig_fasta(mt_validated_all, mt_primary)
        result["mt_primary_id"] = primary_id
        result["mt_primary_len"] = primary_len
    except Exception as e:
        result["status"] = "FAILED_PRIMARY_SELECT"
        result["notes"] = str(e)
        return result

    # ---- Bowtie2 -> samtools depth (on primary) ----
    idx_prefix = path.join(sweep_dir, "mt_primary_idx")
    sam_path = path.join(sweep_dir, "mt_primary.sam")
    bam_path = path.join(sweep_dir, "mt_primary.bam")
    bam_sorted = path.join(sweep_dir, "mt_primary.sorted.bam")
    depth_txt = path.join(sweep_dir, "mt_primary.depth.txt")

    try:
        run_logged(["bowtie2-build", mt_primary, idx_prefix], path.join(logs_dir, "log_bowtie2_build.txt"))

        bt_cmd = ["bowtie2", "-x", idx_prefix, "-p", str(threads)]
        if paired_mode:
            bt_cmd += ["-1", ds_r1, "-2", ds_r2]
        else:
            bt_cmd += ["-U", ds_se]

        # bowtie2: alignments to stdout; summary to stderr
        with open(sam_path, "w") as sam_f, open(path.join(logs_dir, "log_bowtie2_map.txt"), "w") as log_f:
            proc = subprocess.run(bt_cmd, stdout=sam_f, stderr=log_f)
        if proc.returncode != 0:
            raise RuntimeError(f"bowtie2 failed (exit {proc.returncode}).")

        run_logged(["samtools", "view", "-bS", sam_path, "-o", bam_path], path.join(logs_dir, "log_samtools_view.txt"))
        run_logged(["samtools", "sort", "-o", bam_sorted, bam_path], path.join(logs_dir, "log_samtools_sort.txt"))
        run_logged(["samtools", "index", bam_sorted], path.join(logs_dir, "log_samtools_index.txt"))

        with open(depth_txt, "w") as out_f, open(path.join(logs_dir, "log_samtools_depth.txt"), "w") as log_f:
            proc = subprocess.run(["samtools", "depth", "-a", bam_sorted], stdout=out_f, stderr=log_f)
        if proc.returncode != 0:
            raise RuntimeError(f"samtools depth failed (exit {proc.returncode}).")

        mean_depth, pct_ge = compute_depth_stats(depth_txt, int(primary_len), float(depth_target))
        result["mean_depth"] = f"{mean_depth:.4f}"
        result["pct_sites_ge_target"] = f"{pct_ge:.2f}"
        result["status"] = "OK"

    except Exception as e:
        result["status"] = "FAILED_MAPPING_DEPTH"
        result["notes"] = str(e)
        return result

    return result


# ------------------------- Main -------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="lcWGS MT Sweep: trim/QC -> MT sweep -> select minimum >= target mean depth -> optional nuclear steps"
    )

    # Inputs
    parser.add_argument("-1", "--one", metavar="FILE", help="Paired-end Forward (R1)")
    parser.add_argument("-2", "--two", metavar="FILE", help="Paired-end Reverse (R2)")
    parser.add_argument("-u", "--unpaired", metavar="FILE", help="Single-end reads")

    # Outputs & control
    parser.add_argument("-o", "--output", required=True, metavar="DIR", help="Output folder")
    parser.add_argument("-e", "--erase", action="store_true", help="Erase output folder if it exists")
    parser.add_argument("-v", "--verbose", action="store_true", help="Keep temporary and intermediate files")
    parser.add_argument("-t", "--threads", type=int, required=True, metavar="INT", help="Number of threads")
    parser.add_argument("-m", "--mtonly", action="store_true", help="Stop after mitochondrial steps")

    # MT sweep parameters
    parser.add_argument("--mt_sweep_reads", type=str, default="100000,250000,500000,1000000",
                        help="Comma-separated sweep of read counts. Use 0 to represent ALL trimmed reads.")
    parser.add_argument("--seed", type=int, default=1,
                        help="Random seed base for downsampling (seed increments per sweep level).")
    parser.add_argument("--mt_depth_target", type=float, default=10.0,
                        help="Target mean depth for MT (default 10.0).")
    parser.add_argument("--mtlen", type=float, default=1000.0,
                        help="Min contig length filter (bp) on SPAdes headers (default 1000).")
    parser.add_argument("--mtcov_kmer", type=float, default=10.0,
                        help="Min SPAdes k-mer coverage filter from contig headers (default 10.0).")
    parser.add_argument("--stop_after_pass", action="store_true",
                        help="Stop sweep once the first passing level is found (requires ascending sweep).")

    # Optional nuclear inputs
    parser.add_argument("-d", "--db_uc", metavar="FILE", help="UCE probes FASTA (optional).")
    parser.add_argument("--busco_lineage", default="insecta_odb10",
                        help="BUSCO lineage dataset (default insecta_odb10).")

    # Adapters (optional override)
    parser.add_argument("--pe_adapters", help="Path to TruSeq3-PE.fa (Optional)")
    parser.add_argument("--se_adapters", help="Path to TruSeq3-SE.fa (Optional)")

    args = parser.parse_args()

    # ---- Validate mode ----
    paired_mode = args.unpaired is None
    if paired_mode and (not args.one or not args.two):
        raise SystemExit("\nERROR: Paired mode requires -1 and -2.\n")
    if (not paired_mode) and (args.one or args.two):
        raise SystemExit("\nERROR: Single-end mode (-u) cannot use -1/-2.\n")

    # ---- Resolve paths ----
    output_dir = path.abspath(args.output)
    in_r1 = path.abspath(args.one) if args.one else None
    in_r2 = path.abspath(args.two) if args.two else None
    in_se = path.abspath(args.unpaired) if args.unpaired else None
    uce_db = path.abspath(args.db_uc) if args.db_uc else None

    # ---- Required executables ----
    for exe in ["trimmomatic", "fastqc", "spades.py", "mitoz", "bowtie2-build", "bowtie2", "samtools", "seqtk"]:
        which_or_die(exe)

    # ---- Setup output ----
    print(f"\n[{ts()}] Pipeline Started")

    if path.exists(output_dir):
        if args.erase:
            shutil.rmtree(output_dir)
        else:
            raise SystemExit(f"\nERROR: Output '{output_dir}' exists. Use --erase to overwrite.\n")
    ensure_dir(output_dir)

    tmp_dir = path.join(output_dir, "tmp")
    ensure_dir(tmp_dir)

    # ---- Adapter detection ----
    pe_adap = args.pe_adapters or find_trimmomatic_adapter("TruSeq3-PE.fa")
    se_adap = args.se_adapters or find_trimmomatic_adapter("TruSeq3-SE.fa")
    if paired_mode and not pe_adap:
        raise SystemExit("ERROR: TruSeq3-PE.fa not found. Use --pe_adapters or ensure trimmomatic adapters exist in env.")
    if (not paired_mode) and not se_adap:
        raise SystemExit("ERROR: TruSeq3-SE.fa not found. Use --se_adapters or ensure trimmomatic adapters exist in env.")

    ADAPTER_STR = f"ILLUMINACLIP:{pe_adap}:2:30:10" if paired_mode else f"ILLUMINACLIP:{se_adap}:2:30:10"

    # ---------------- 1) TRIMMING ----------------
    print(f"[{ts()}] Reads Trimming")
    trim_log = path.join(tmp_dir, "log_trimmomatic.txt")
    if paired_mode:
        trim_cmd = [
            "trimmomatic", "PE", in_r1, in_r2,  # type: ignore
            "-baseout", "trimmomatic",
            ADAPTER_STR,
            "LEADING:20", "TRAILING:20", "SLIDINGWINDOW:5:20", "MINLEN:95"
        ]
    else:
        trim_cmd = [
            "trimmomatic", "SE", in_se,  # type: ignore
            "trimmomatic_U",
            ADAPTER_STR,
            "LEADING:20", "TRAILING:20", "SLIDINGWINDOW:5:20", "MINLEN:95"
        ]
    run_logged(trim_cmd, trim_log, cwd=tmp_dir)

    # Standardize filenames
    trim_r1 = trim_r2 = trim_se = None
    if paired_mode:
        os.rename(path.join(tmp_dir, "trimmomatic_1P"), path.join(tmp_dir, "trimmomatic_1P.fq"))
        os.rename(path.join(tmp_dir, "trimmomatic_2P"), path.join(tmp_dir, "trimmomatic_2P.fq"))
        trim_r1 = path.join(tmp_dir, "trimmomatic_1P.fq")
        trim_r2 = path.join(tmp_dir, "trimmomatic_2P.fq")
    else:
        os.rename(path.join(tmp_dir, "trimmomatic_U"), path.join(tmp_dir, "trimmomatic_U.fq"))
        trim_se = path.join(tmp_dir, "trimmomatic_U.fq")

    # ---------------- 2) QC ----------------
    print(f"[{ts()}] Reads QC")
    fastqc_dir = path.join(output_dir, "fastqc")
    ensure_dir(fastqc_dir)
    qc_log = path.join(tmp_dir, "log_fastqc.txt")
    if paired_mode:
        run_logged(["fastqc", "-o", fastqc_dir, trim_r1, trim_r2], qc_log)
    else:
        run_logged(["fastqc", "-o", fastqc_dir, trim_se], qc_log)

    # ---------------- 3) MT SWEEP ----------------
    print(f"[{ts()}] MT Sweep (random downsampling -> SPAdes -> MitoZ -> depth)")

    sweep_reads = sorted(parse_int_list(args.mt_sweep_reads))
    mt_sweep_root = path.join(output_dir, "mt_sweep")
    ensure_dir(mt_sweep_root)

    rows: List[Dict[str, object]] = []
    pass_levels: List[Tuple[int, float]] = []

    sample_name = path.basename(output_dir.rstrip("/"))

    for i, reads_n in enumerate(sweep_reads):
        level_name = "ALL" if reads_n == 0 else str(reads_n)
        sweep_dir = path.join(mt_sweep_root, f"reads_{level_name}")
        seed_i = int(args.seed) + int(i)

        print(f"[{ts()}]  Sweep level: reads={level_name} (seed={seed_i})")

        res = run_mt_sweep_level(
            sweep_dir=sweep_dir,
            paired_mode=paired_mode,
            trim_r1=trim_r1,
            trim_r2=trim_r2,
            trim_se=trim_se,
            reads_n=reads_n,
            seed=seed_i,
            threads=args.threads,
            mtlen=args.mtlen,
            mtcov_kmer=args.mtcov_kmer,
            depth_target=args.mt_depth_target,
            sample_name=sample_name,
        )

        passed = False
        mean_depth_val = None
        if res.get("status") == "OK":
            try:
                mean_depth_val = float(res["mean_depth"])  # type: ignore
                if mean_depth_val >= float(args.mt_depth_target):
                    passed = True
                    pass_levels.append((reads_n, mean_depth_val))
            except Exception:
                passed = False

        res["passed_target"] = "YES" if passed else "NO"
        res["selected_min_pass"] = "NO"
        rows.append(res)

        if args.stop_after_pass and passed:
            print(f"[{ts()}]  Stop-after-pass enabled; minimum passing level found at reads={level_name}.")
            break

    # Write TSV summary (always)
    tsv_path = path.join(output_dir, "mt_sweep_summary.tsv")
    with open(tsv_path, "w", newline="") as out_tsv:
        writer = csv.DictWriter(
            out_tsv,
            fieldnames=[
                "sweep_reads", "status", "passed_target", "selected_min_pass",
                "mt_primary_id", "mt_primary_len", "mt_num_valid_contigs", "mt_total_valid_len",
                "mean_depth", "pct_sites_ge_target", "notes"
            ],
            delimiter="\t"
        )
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    if not pass_levels:
        raise SystemExit(
            f"\nERROR: No sweep level achieved mean depth >= {args.mt_depth_target}x.\n"
            f"See: {tsv_path}\n"
        )

    # Select minimum passing reads level (prefer non-zero minima if any pass)
    pass_nonzero = [(r, d) for (r, d) in pass_levels if r != 0]
    if pass_nonzero:
        selected_reads, selected_depth = sorted(pass_nonzero, key=lambda x: x[0])[0]
    else:
        selected_reads, selected_depth = pass_levels[0]

    for r in rows:
        if int(r.get("sweep_reads")) == int(selected_reads):
            r["selected_min_pass"] = "YES"

    # Re-write TSV with selection marked
    with open(tsv_path, "w", newline="") as out_tsv:
        writer = csv.DictWriter(
            out_tsv,
            fieldnames=[
                "sweep_reads", "status", "passed_target", "selected_min_pass",
                "mt_primary_id", "mt_primary_len", "mt_num_valid_contigs", "mt_total_valid_len",
                "mean_depth", "pct_sites_ge_target", "notes"
            ],
            delimiter="\t"
        )
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    # Copy selected MT outputs to top-level output
    sel_level_name = "ALL" if selected_reads == 0 else str(selected_reads)
    selected_dir = path.join(mt_sweep_root, f"reads_{sel_level_name}")

    sel_mt_primary = path.join(selected_dir, "mt_primary.fa")
    sel_mitoz_summary = path.join(selected_dir, "annot.mtcandidate.fa.result", "summary.txt")
    sel_gbf = path.join(selected_dir, "annot.mtcandidate.fa.result", "annot_mtcandidate.fa_mitoscaf.fa.gbf")

    if file_nonempty(sel_mt_primary):
        shutil.copy2(sel_mt_primary, path.join(output_dir, "mt.fasta"))
    if path.exists(sel_mitoz_summary):
        shutil.copy2(sel_mitoz_summary, path.join(output_dir, "mt_summary.txt"))
    if path.exists(sel_gbf):
        shutil.copy2(sel_gbf, path.join(output_dir, "mt.gbf"))

    print(f"[{ts()}] Selected minimum passing sweep: reads={sel_level_name}, mean_depth≈{selected_depth:.3f}x")
    print(f"[{ts()}] MT sweep summary TSV: {tsv_path}")

    if args.mtonly:
        if not args.verbose:
            shutil.rmtree(tmp_dir)
        print(f"[{ts()}] Finished (MT Only)")
        return

    # ---------------- 4) REMOVE MT READS (map all trimmed reads to selected MT reference) ----------------
    print(f"[{ts()}] Removing MT Reads (map all trimmed reads to selected mt.fasta)")

    mt_ref = path.join(output_dir, "mt.fasta")
    if not file_nonempty(mt_ref):
        raise SystemExit("\nERROR: Selected mt.fasta missing; cannot proceed to nuclear steps.\n")

    mt_idx_prefix = path.join(tmp_dir, "mtref_idx")
    run_logged(["bowtie2-build", mt_ref, mt_idx_prefix], path.join(tmp_dir, "log_bowtie2_build_mtref.txt"))

    if paired_mode:
        bt_cmd = [
            "bowtie2", "-x", mt_idx_prefix,
            "--un-conc", path.join(tmp_dir, "mtless"),
            "-p", str(args.threads),
            "-1", trim_r1, "-2", trim_r2
        ]
        run_logged(bt_cmd, path.join(tmp_dir, "log_bowtie2_align_mtref.txt"))
        os.rename(path.join(tmp_dir, "mtless.1"), path.join(tmp_dir, "mtless.1.fq"))
        os.rename(path.join(tmp_dir, "mtless.2"), path.join(tmp_dir, "mtless.2.fq"))
        nc_in = ["-1", path.join(tmp_dir, "mtless.1.fq"), "-2", path.join(tmp_dir, "mtless.2.fq")]
    else:
        bt_cmd = [
            "bowtie2", "-x", mt_idx_prefix,
            "--un", path.join(tmp_dir, "mtless.U.fq"),
            "-p", str(args.threads),
            "-U", trim_se
        ]
        run_logged(bt_cmd, path.join(tmp_dir, "log_bowtie2_align_mtref.txt"))
        nc_in = ["-s", path.join(tmp_dir, "mtless.U.fq")]

    # ---------------- 5) NUCLEAR ASSEMBLY ----------------
    print(f"[{ts()}] Nuclear Assembly (MT-depleted reads)")
    nc_spades_out = path.join(tmp_dir, "spades_nc")
    run_logged(
        ["spades.py"] + nc_in + ["-o", nc_spades_out, "-t", str(args.threads), "--cov-cutoff", "off"],
        path.join(tmp_dir, "log_spades_nc.txt")
    )
    nc_contigs = path.join(nc_spades_out, "contigs.fasta")
    if not file_nonempty(nc_contigs):
        raise SystemExit("ERROR: Nuclear SPAdes produced no contigs.")
    shutil.copy2(nc_contigs, path.join(output_dir, "nc.fasta"))

    # Optional UCE match (Phyluce)
        # Optional UCE match (Phyluce)
    if uce_db:
        # Only require Phyluce if user requested UCE matching
        if not shutil.which("phyluce_assembly_match_contigs_to_probes"):
            raise SystemExit(
                "\nERROR: UCE probes provided (-d), but Phyluce is not available in the current environment.\n"
                "Activate your Phyluce environment and re-run, e.g.:\n"
                "  mamba activate phyluce_env\n"
            )

        print(f"[{ts()}] UCE Annotation (Phyluce)")
        phyluce_out = path.join(tmp_dir, "folder_phyluce")
        ensure_dir(phyluce_out)
        run_logged(
            [
                "phyluce_assembly_match_contigs_to_probes",
                "--contigs", nc_contigs,
                "--probes", uce_db,
                "--output", phyluce_out,
                "--log-path", phyluce_out
            ],
            path.join(tmp_dir, "log_phyluce.txt"),
            cwd=tmp_dir
        )
        uce_log = path.join(phyluce_out, "phyluce_assembly_match_contigs_to_probes.log")
        if path.exists(uce_log):
            shutil.copy2(uce_log, path.join(output_dir, "summary_uces.txt"))
    else:
        print(f"[{ts()}] Skipping Phyluce (no -d/--db_uc provided).")


    # BUSCO (default insecta_odb10)
    print(f"[{ts()}] BUSCO Annotation (lineage={args.busco_lineage})")
    busco_dirname = "folder_busco"
    run_logged(
        [
            "busco",
            "-i", nc_contigs,
            "-o", busco_dirname,
            "--mode", "geno",
            "--cpu", str(args.threads),
            "-l", args.busco_lineage,
            "-f"
        ],
        path.join(tmp_dir, "log_busco.txt"),
        cwd=tmp_dir
    )

    # Copy BUSCO short summary if present (BUSCO v5 commonly writes short_summary*.txt)
    busco_run_dir = path.join(tmp_dir, busco_dirname, f"run_{args.busco_lineage}")
    if path.isdir(busco_run_dir):
        for fname in os.listdir(busco_run_dir):
            if fname.startswith("short_summary") and fname.endswith(".txt"):
                shutil.copy2(path.join(busco_run_dir, fname), path.join(output_dir, "summary_busco.txt"))
                break

    if not args.verbose:
        shutil.rmtree(tmp_dir)

    print(f"[{ts()}] Pipeline Completed Successfully")


if __name__ == "__main__":
    try:
        main()
    except RuntimeError as e:
        raise SystemExit(str(e))
