#!/usr/bin/env python3
"""
lcWGS MT Post-processing (corrected parsing + safe best-contig selection)

For each successful run (from <batch_root>/reports/selected_min_pass.tsv):
  - locate the selected sweep level directory (reads_<selected>)
  - parse MitoZ GenBank output (.gbf) to summarize gene content per contig
      * robust PCG normalization (COX/COI/ND/NADH dehydrogenase product strings)
      * robust rRNA detection (rrnL/rrnS from 16S/12S/product strings)
      * robust tRNA counting (tRNA features OR gene features whose qualifiers indicate tRNA/trnX)
  - compute per-contig mean depth + covered fraction by mapping selected sweep reads to mt_validated_all.fa
  - score and select the best contig per run using weights:
        score = w_gene * gene_completeness + w_depth * depth_norm + w_len * len_norm
    with a guard:
        if ANY contig length >= --min_len_best, selection is restricted to those contigs only
        (prevents tiny fragments being selected)

Outputs (under <batch_root>/reports/mt_post/):
  - contigs.tsv
  - runs.tsv
  - best_per_sample.tsv
  - best_candidates.fasta
  - blast_hits.tsv  (if BLAST+ available and not --skip_blast)

Dependencies in PATH:
  bowtie2, bowtie2-build, samtools
  (optional) makeblastdb, blastn
Python deps:
  biopython
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO


# ------------------------- Canonical gene sets -------------------------

EXPECTED_PCGS = [
    "atp6", "atp8", "cob",
    "cox1", "cox2", "cox3",
    "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6",
]

PCG_SYNONYMS = {
    "cytb": "cob",
    "cyt b": "cob",
    "coi": "cox1",
    "coii": "cox2",
    "coiii": "cox3",
    "cox i": "cox1",
    "cox ii": "cox2",
    "cox iii": "cox3",
    "nd1": "nad1",
    "nd2": "nad2",
    "nd3": "nad3",
    "nd4": "nad4",
    "nd4l": "nad4l",
    "nd5": "nad5",
    "nd6": "nad6",
}

ROMAN_TO_ARABIC = {"i": "1", "ii": "2", "iii": "3"}


def canonical_pcg(raw: str) -> str:
    """
    Map a variety of gene/product strings to canonical PCG codes:
      atp6 atp8 cob cox1 cox2 cox3 nad1 nad2 nad3 nad4 nad4l nad5 nad6
    """
    if not raw:
        return ""
    s = raw.strip().lower()
    s = re.sub(r"[\(\)\[\],;:/]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()

    # direct map
    if s in PCG_SYNONYMS:
        return PCG_SYNONYMS[s]

    # cytochrome b
    if "cytochrome b" in s:
        return "cob"

    # ATP synthase
    if ("atp synthase" in s) or ("atpase" in s):
        if "subunit 6" in s:
            return "atp6"
        if "subunit 8" in s:
            return "atp8"

    # COX
    if ("cytochrome c oxidase" in s) or ("cytochrome oxidase" in s) or s.startswith("cox"):
        m = re.search(r"subunit\s+([0-9]+|i{1,3})\b", s)
        if m:
            x = ROMAN_TO_ARABIC.get(m.group(1), m.group(1))
            if x in {"1", "2", "3"}:
                return f"cox{x}"
        m2 = re.fullmatch(r"cox\s*([123])", s)
        if m2:
            return f"cox{m2.group(1)}"

    # ND / NAD forms
    m = re.fullmatch(r"nd\s*([0-9]+)\s*(l)?", s)
    if m:
        n = m.group(1)
        suf = m.group(2) or ""
        return f"nad{n}{'l' if suf else ''}"

    m = re.fullmatch(r"nad\s*([0-9]+)\s*(l)?", s)
    if m:
        n = m.group(1)
        suf = m.group(2) or ""
        return f"nad{n}{'l' if suf else ''}"

    # product: "NADH dehydrogenase subunit X"
    if ("nadh dehydrogenase" in s) and ("subunit" in s):
        m = re.search(r"subunit\s*([0-9]+)\s*(l)?", s)
        if m:
            n = m.group(1)
            suf = m.group(2) or ""
            return f"nad{n}{'l' if suf else ''}"

    return ""


def canonical_rrna(gene_q: str, prod_q: str) -> str:
    """
    Return rrnl/rrns if 16S/12S patterns are detected in gene/product text.
    """
    text = f"{gene_q} {prod_q}".strip().lower()
    text = re.sub(r"\s+", " ", text)

    # common patterns
    if "16s" in text or "rrnl" in text or "large subunit ribosomal" in text:
        return "rrnl"
    if "12s" in text or "rrns" in text or "small subunit ribosomal" in text:
        return "rrns"

    # sometimes "ribosomal rna" without 12S/16S label
    if "ribosomal rna" in text:
        # ambiguous; return empty rather than guessing
        return ""

    return ""


def is_trna_feature(feat_type: str, gene_q: str, prod_q: str) -> bool:
    """
    Identify tRNAs even if encoded as type=gene with product/gene qualifiers.
    """
    text = f"{gene_q} {prod_q}".strip().lower()
    text = re.sub(r"\s+", " ", text)

    if feat_type == "tRNA":
        return True
    # qualifier patterns
    if "trna" in text:
        return True
    # gene names often like trnL, trnS, etc.
    if text.startswith("trn"):
        return True
    return False


# Your sample->species mapping (used mainly for reporting/BLAST labels)
DEFAULT_SAMPLE_SPECIES = {
    "43_assoi": "Ameles assoi",
    "16-CI1f": "Pseudoyersinia betancuriae",
    "41-K26": "Ameles picteti",
    "59-T3": "Ameles dumonti",
    "23_spal": "Ameles spallanzania",
    "ADEC02": "Ameles decolor",
    "RBE05": "Riventina baetica",
    "AAP07": "Apteromantis aptera",
}


# ------------------------- Utilities -------------------------

def die(msg: str) -> None:
    raise SystemExit(f"ERROR: {msg}")


def which_or_none(exe: str) -> Optional[str]:
    return shutil.which(exe)


def run(cmd: List[str], log: Optional[Path] = None, cwd: Optional[Path] = None) -> None:
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        with open(log, "w") as fh:
            p = subprocess.run(cmd, stdout=fh, stderr=fh, cwd=str(cwd) if cwd else None)
    else:
        p = subprocess.run(cmd)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed (exit {p.returncode}): {' '.join(cmd)}")


def read_tsv(path: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            rows.append({k: (v if v is not None else "") for k, v in r.items()})
    return rows


def write_tsv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


# ------------------------- Data models -------------------------

@dataclass
class GeneStats:
    pcgs_present: List[str]
    rrnas_present: List[str]
    trna_count: int
    total_features: int

    @property
    def pcg_count(self) -> int:
        return len(set(self.pcgs_present))

    @property
    def rrna_count(self) -> int:
        return len(set(self.rrnas_present))

    @property
    def gene_completeness(self) -> float:
        """
        Completeness in [0,1] across 37 expected mitochondrial genes:
          13 PCGs + 2 rRNAs + 22 tRNAs.
        """
        pcg = min(self.pcg_count, 13)
        rrna = min(self.rrna_count, 2)
        trna = min(self.trna_count, 22)
        return float(pcg + rrna + trna) / 37.0

    @property
    def missing_pcgs(self) -> List[str]:
        have = set(self.pcgs_present)
        return [g for g in EXPECTED_PCGS if g not in have]


@dataclass
class ContigMetrics:
    contig_id: str
    length: int
    mean_depth: float
    covered_frac: float
    gene_stats: Optional[GeneStats]

    @property
    def gene_completeness(self) -> float:
        return self.gene_stats.gene_completeness if self.gene_stats else 0.0

    @property
    def pcg_count(self) -> int:
        return self.gene_stats.pcg_count if self.gene_stats else 0

    @property
    def rrna_count(self) -> int:
        return self.gene_stats.rrna_count if self.gene_stats else 0

    @property
    def trna_count(self) -> int:
        return self.gene_stats.trna_count if self.gene_stats else 0


# ------------------------- Parsing -------------------------

def parse_genbank_by_contig(gbf: Path) -> Dict[str, GeneStats]:
    """
    Parse multi-record GenBank produced by MitoZ and return gene content per record.
    """
    by: Dict[str, GeneStats] = {}
    if not gbf.exists() or gbf.stat().st_size == 0:
        return by

    for rec in SeqIO.parse(str(gbf), "genbank"):
        pcgs: List[str] = []
        rrnas: List[str] = []
        trna_names: List[str] = []
        feat_n = 0

        for feat in rec.features:
            # include 'gene' because tRNAs (and sometimes rRNAs) may be encoded as gene features
            if feat.type not in {"CDS", "rRNA", "tRNA", "gene"}:
                continue
            feat_n += 1

            q = feat.qualifiers
            gene_q = (q.get("gene", [""])[0] or "").strip()
            prod_q = (q.get("product", [""])[0] or "").strip()

            # ---- PCGs ----
            if feat.type == "CDS":
                g = canonical_pcg(gene_q) or canonical_pcg(prod_q)
                if g:
                    pcgs.append(g)
                continue

            # ---- rRNAs ----
            # rRNAs may appear as type=rRNA or as type=gene with product containing "16S/12S"
            rr = ""
            if feat.type == "rRNA" or feat.type == "gene":
                rr = canonical_rrna(gene_q, prod_q)
            if rr:
                rrnas.append(rr)
                continue

            # ---- tRNAs ----
            if is_trna_feature(feat.type, gene_q, prod_q):
                trna_names.append(gene_q or prod_q or "tRNA")
                continue

        trna_count = len([x for x in trna_names if x])
        by[rec.id] = GeneStats(
            pcgs_present=pcgs,
            rrnas_present=rrnas,
            trna_count=trna_count,
            total_features=feat_n
        )

    return by


def fasta_lengths(fa: Path) -> Dict[str, int]:
    d: Dict[str, int] = {}
    for rec in SeqIO.parse(str(fa), "fasta"):
        d[rec.id] = len(rec.seq)
    return d


def parse_selected_level(run_dir: Path) -> Tuple[str, Path]:
    """
    Determine selected sweep reads from mt_sweep_summary.tsv (selected_min_pass==YES)
    and return (level_name, sweep_dir).
    level_name is 'ALL' or the integer reads (as string).
    """
    summ = run_dir / "mt_sweep_summary.tsv"
    if not summ.exists():
        die(f"Missing mt_sweep_summary.tsv in run: {run_dir}")

    rows = read_tsv(summ)
    sel = [r for r in rows if r.get("selected_min_pass", "") == "YES"]
    if not sel:
        sel = [r for r in rows if r.get("passed_target", "") == "YES"]
    if not sel:
        die(f"No selected_min_pass or passed_target found in: {summ}")

    sweep_reads = sel[0].get("sweep_reads", "").strip()
    level = "ALL" if sweep_reads == "0" else sweep_reads

    sweep_dir = run_dir / "mt_sweep" / f"reads_{level}"
    if not sweep_dir.exists():
        die(f"Selected sweep dir not found: {sweep_dir}")

    return level, sweep_dir


# ------------------------- Coverage computation -------------------------

def map_reads_and_coverage(
    sweep_dir: Path,
    mt_validated_fa: Path,
    threads: int,
    force: bool,
) -> Dict[str, Tuple[float, float]]:
    """
    Map downsampled reads (ds_R1/ds_R2 or ds_SE) to mt_validated_all.fa and compute per-contig:
      - mean_depth
      - covered_frac (fraction of reference bases with depth > 0)
    Returns dict: contig -> (mean_depth, covered_frac)

    Writes/uses cached files in:
      sweep_dir / candidate_eval/
    """
    outdir = sweep_dir / "candidate_eval"
    outdir.mkdir(parents=True, exist_ok=True)

    cov_tsv = outdir / "per_contig_coverage.tsv"
    if cov_tsv.exists() and not force:
        out: Dict[str, Tuple[float, float]] = {}
        with open(cov_tsv, "r") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(hdr)}
            for line in fh:
                p = line.rstrip("\n").split("\t")
                rname = p[idx["rname"]]
                md = float(p[idx["mean_depth"]])
                cf = float(p[idx["covered_frac"]])
                out[rname] = (md, cf)
        return out

    ds_r1 = sweep_dir / "ds_R1.fq"
    ds_r2 = sweep_dir / "ds_R2.fq"
    ds_se = sweep_dir / "ds_SE.fq"

    paired = ds_r1.exists() and ds_r2.exists()
    single = ds_se.exists()

    if not (paired or single):
        die(f"Downsampled reads not found in sweep dir: {sweep_dir}")

    if not which_or_none("bowtie2-build") or not which_or_none("bowtie2") or not which_or_none("samtools"):
        die("bowtie2/bowtie2-build/samtools must be in PATH to compute per-contig coverage.")

    idx_prefix = outdir / "mt_validated_idx"
    bt2_any = idx_prefix.with_suffix(".1.bt2")
    if not bt2_any.exists() or force:
        run(["bowtie2-build", str(mt_validated_fa), str(idx_prefix)], log=outdir / "log_bowtie2_build.txt")

    bam_sorted = outdir / "aln.sorted.bam"
    if not bam_sorted.exists() or force:
        log_bt = open(outdir / "log_bowtie2_map.txt", "w")
        log_st = open(outdir / "log_samtools.txt", "w")

        if paired:
            bt_cmd = ["bowtie2", "-x", str(idx_prefix), "-p", str(threads), "-1", str(ds_r1), "-2", str(ds_r2)]
        else:
            bt_cmd = ["bowtie2", "-x", str(idx_prefix), "-p", str(threads), "-U", str(ds_se)]

        p1 = subprocess.Popen(bt_cmd, stdout=subprocess.PIPE, stderr=log_bt)
        p2 = subprocess.Popen(["samtools", "view", "-bS", "-"], stdin=p1.stdout, stdout=subprocess.PIPE, stderr=log_st)
        p3 = subprocess.Popen(["samtools", "sort", "-o", str(bam_sorted), "-"], stdin=p2.stdout, stderr=log_st)

        assert p1.stdout is not None
        assert p2.stdout is not None
        p1.stdout.close()
        p2.stdout.close()

        rc3 = p3.wait()
        rc2 = p2.wait()
        rc1 = p1.wait()
        log_bt.close()
        log_st.close()

        if rc1 != 0 or rc2 != 0 or rc3 != 0:
            raise RuntimeError(f"Mapping pipeline failed in {outdir}. See logs.")

        run(["samtools", "index", str(bam_sorted)], log=outdir / "log_samtools_index.txt")

    cov: Dict[str, Tuple[float, float]] = {}

    # Prefer samtools coverage when available
    cov_raw = outdir / "samtools_coverage.tsv"
    try:
        with open(cov_raw, "w") as out, open(outdir / "log_samtools_coverage.txt", "w") as log:
            p = subprocess.run(["samtools", "coverage", str(bam_sorted)], stdout=out, stderr=log)
        if p.returncode == 0 and cov_raw.stat().st_size > 0:
            with open(cov_raw, "r") as fh:
                header = fh.readline().rstrip("\n").split("\t")
                idx = {h: i for i, h in enumerate(header)}
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    rname = parts[idx.get("#rname", idx.get("rname", 0))]
                    meand = float(parts[idx.get("meandepth", idx.get("mean_depth", -1))])
                    covbases = float(parts[idx.get("covbases", -1)])
                    endpos = float(parts[idx.get("endpos", -1)])
                    startpos = float(parts[idx.get("startpos", 1)])
                    length = max(1.0, endpos - startpos + 1.0)
                    covered_frac = covbases / length
                    cov[rname] = (meand, covered_frac)
    except Exception:
        cov = {}

    # Fallback: samtools depth -a
    if not cov:
        depth_txt = outdir / "samtools_depth.txt"
        with open(depth_txt, "w") as out, open(outdir / "log_samtools_depth.txt", "w") as log:
            p = subprocess.run(["samtools", "depth", "-a", str(bam_sorted)], stdout=out, stderr=log)
        if p.returncode != 0:
            raise RuntimeError(f"samtools depth failed in {outdir}. See logs.")

        lengths = fasta_lengths(mt_validated_fa)
        sums: Dict[str, int] = {k: 0 for k in lengths}
        covpos: Dict[str, int] = {k: 0 for k in lengths}

        with open(depth_txt, "r") as fh:
            for line in fh:
                c, _pos, d = line.rstrip("\n").split("\t")
                di = int(d)
                sums[c] += di
                if di > 0:
                    covpos[c] += 1

        for c, L in lengths.items():
            md = float(sums.get(c, 0)) / float(L)
            cf = float(covpos.get(c, 0)) / float(L)
            cov[c] = (md, cf)

    with open(cov_tsv, "w") as fh:
        fh.write("rname\tmean_depth\tcovered_frac\n")
        for rname, (md, cf) in sorted(cov.items()):
            fh.write(f"{rname}\t{md:.6f}\t{cf:.6f}\n")

    return cov


# ------------------------- Best contig selection -------------------------

def select_best_contig(
    contigs: List[ContigMetrics],
    w_gene: float,
    w_depth: float,
    w_len: float,
    min_len_best: int,
) -> Tuple[ContigMetrics, float]:
    if not contigs:
        die("No contigs to select from.")

    # If there is at least one "long" contig, restrict selection to long contigs only.
    if any(c.length >= min_len_best for c in contigs):
        contigs = [c for c in contigs if c.length >= min_len_best]

    max_len = max(c.length for c in contigs) or 1
    max_depth = max(c.mean_depth for c in contigs) or 1.0

    best = contigs[0]
    best_score = -1e9

    for c in contigs:
        len_norm = float(c.length) / float(max_len)
        depth_norm = float(c.mean_depth) / float(max_depth) if max_depth > 0 else 0.0
        gene = c.gene_completeness

        score = (w_gene * gene) + (w_depth * depth_norm) + (w_len * len_norm)

        # deterministic tie-break
        if (score > best_score) or (abs(score - best_score) < 1e-12 and c.length > best.length):
            best, best_score = c, score

    return best, best_score


def write_best_contig_fasta(run_dir: Path, sweep_dir: Path, best_id: str, out_fa: Path) -> None:
    """
    Extract best contig sequence from mt_validated_all.fa or mt_primary.fa.
    """
    fa1 = sweep_dir / "mt_validated_all.fa"
    fa2 = sweep_dir / "mt_primary.fa"
    src = fa1 if fa1.exists() else fa2
    if not src.exists():
        die(f"No source fasta found for contig extraction in {sweep_dir}")

    for rec in SeqIO.parse(str(src), "fasta"):
        if rec.id == best_id:
            out_fa.parent.mkdir(parents=True, exist_ok=True)
            SeqIO.write(rec, str(out_fa), "fasta")
            return

    die(f"Best contig id {best_id} not found in {src}")


# ------------------------- BLAST helpers -------------------------

def make_blast_db(ref_fasta: Path, db_prefix: Path, force: bool) -> None:
    if not which_or_none("makeblastdb"):
        die("makeblastdb not found in PATH (install BLAST+ or use --skip_blast).")

    nin = db_prefix.with_suffix(".nin")
    if nin.exists() and not force:
        return

    run(
        ["makeblastdb", "-in", str(ref_fasta), "-dbtype", "nucl", "-out", str(db_prefix)],
        log=db_prefix.parent / "log_makeblastdb.txt",
    )


def blast_top_hit(query_fa: Path, db_prefix: Path, threads: int) -> Optional[Dict[str, object]]:
    if not which_or_none("blastn"):
        die("blastn not found in PATH (install BLAST+ or use --skip_blast).")

    out = query_fa.parent / f"{query_fa.stem}.blast.tsv"
    fmt = "6 qseqid sseqid pident length qlen slen evalue bitscore stitle"
    cmd = [
        "blastn", "-task", "megablast",
        "-query", str(query_fa),
        "-db", str(db_prefix),
        "-outfmt", fmt,
        "-max_target_seqs", "5",
        "-evalue", "1e-20",
        "-num_threads", str(threads),
    ]
    with open(out, "w") as fh, open(query_fa.parent / f"{query_fa.stem}.blast.log", "w") as log:
        p = subprocess.run(cmd, stdout=fh, stderr=log)
    if p.returncode != 0:
        return None

    with open(out, "r") as fh:
        line = fh.readline().rstrip("\n")
        if not line:
            return None
        parts = line.split("\t")
        return {
            "qseqid": parts[0],
            "sseqid": parts[1],
            "pident": float(parts[2]),
            "aln_len": int(parts[3]),
            "qlen": int(parts[4]),
            "slen": int(parts[5]),
            "evalue": parts[6],
            "bitscore": float(parts[7]),
            "stitle": parts[8] if len(parts) > 8 else "",
        }


def build_ref_fasta_from_samples(
    runs_rows: List[Dict[str, object]],
    ref_samples: List[str],
    out_fa: Path,
) -> None:
    """
    Build a FASTA containing best candidate contigs for specified samples from this batch.
    Sequence IDs are set to the sample name for easy BLAST attribution.
    """
    by_sample: Dict[str, Dict[str, object]] = {}
    for r in runs_rows:
        if r.get("status") != "OK":
            continue
        s = str(r["sample"])
        score = float(r["run_best_score"])
        if s not in by_sample or score > float(by_sample[s]["run_best_score"]):
            by_sample[s] = r

    out_fa.parent.mkdir(parents=True, exist_ok=True)
    records = []

    for s in ref_samples:
        if s not in by_sample:
            continue
        best_fa_path = Path(str(by_sample[s]["best_contig_fasta"]))
        if not best_fa_path.exists():
            continue
        rec = next(SeqIO.parse(str(best_fa_path), "fasta"))
        rec.id = s
        rec.name = s
        rec.description = DEFAULT_SAMPLE_SPECIES.get(s, "")
        records.append(rec)

    if not records:
        die("No reference records could be built from ref samples (check names and available runs).")

    SeqIO.write(records, str(out_fa), "fasta")


# ------------------------- Main -------------------------

def main() -> None:
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--batch_root",
        required=True,
        help="Root folder created by run_mt_depth_targets_batch.sh (has runs/, reports/).",
    )
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--force", action="store_true", help="Recompute cached coverage/BLAST outputs.")
    ap.add_argument(
        "--weights",
        default="0.60,0.30,0.10",
        help="Weights for gene,depth,length (comma-separated). Default 0.60,0.30,0.10",
    )
    ap.add_argument(
        "--min_len_best",
        type=int,
        default=8000,
        help="If any contig length >= this, restrict best-contig selection to contigs >= this (default 8000).",
    )

    # BLAST options
    ap.add_argument("--skip_blast", action="store_true")
    ap.add_argument("--blast_ref_fasta", default="", help="FASTA to build BLAST DB from (optional).")
    ap.add_argument(
        "--blast_ref_samples",
        default="43_assoi,41-K26,59-T3,23_spal,ADEC02,RBE05,AAP07",
        help="Comma list of samples to build BLAST DB from (best contigs from this batch). Ignored if --blast_ref_fasta is set.",
    )

    args = ap.parse_args()

    batch_root = Path(args.batch_root).resolve()
    pass_tsv = batch_root / "reports" / "selected_min_pass.tsv"
    if not pass_tsv.exists():
        die(f"Cannot find {pass_tsv}. Provide correct --batch_root.")

    w_gene, w_depth, w_len = [float(x) for x in args.weights.split(",")]

    out_root = batch_root / "reports" / "mt_post"
    out_root.mkdir(parents=True, exist_ok=True)

    pass_rows = read_tsv(pass_tsv)
    if not pass_rows:
        die(f"No rows found in {pass_tsv}")

    contigs_out_rows: List[Dict[str, object]] = []
    runs_out_rows: List[Dict[str, object]] = []
    best_fasta_records = []

    for r in pass_rows:
        sample = r["sample"]
        depth_target = float(r["depth_target"])
        run_dir = Path(r["run_dir"]).resolve()

        if not run_dir.exists():
            runs_out_rows.append({
                "sample": sample,
                "depth_target": depth_target,
                "run_dir": str(run_dir),
                "status": "MISSING_RUN_DIR",
            })
            continue

        try:
            selected_level, sweep_dir = parse_selected_level(run_dir)
        except Exception as e:
            runs_out_rows.append({
                "sample": sample,
                "depth_target": depth_target,
                "run_dir": str(run_dir),
                "status": f"FAILED_SELECT_LEVEL: {e}",
            })
            continue

        mt_validated = sweep_dir / "mt_validated_all.fa"
        if not mt_validated.exists():
            mt_validated = sweep_dir / "mt_primary.fa"
        if not mt_validated.exists():
            runs_out_rows.append({
                "sample": sample,
                "depth_target": depth_target,
                "run_dir": str(run_dir),
                "selected_sweep_reads": selected_level,
                "status": "NO_VALIDATED_FASTA",
            })
            continue

        gbf = sweep_dir / "annot.mtcandidate.fa.result" / "annot_mtcandidate.fa_mitoscaf.fa.gbf"
        genes_by_contig = parse_genbank_by_contig(gbf)

        # Per-contig coverage
        try:
            cov = map_reads_and_coverage(sweep_dir, mt_validated, threads=args.threads, force=args.force)
        except Exception:
            cov = {}

        lengths = fasta_lengths(mt_validated)

        contigs: List[ContigMetrics] = []
        for cid, L in lengths.items():
            md, cf = cov.get(cid, (0.0, 0.0))
            contigs.append(ContigMetrics(
                contig_id=cid,
                length=int(L),
                mean_depth=float(md),
                covered_frac=float(cf),
                gene_stats=genes_by_contig.get(cid),
            ))

        best_contig, best_score = select_best_contig(
            contigs,
            w_gene=w_gene,
            w_depth=w_depth,
            w_len=w_len,
            min_len_best=args.min_len_best,
        )

        best_fa = out_root / "best_contigs" / sample / f"depth_{int(depth_target)}x" / "best_contig.fa"
        write_best_contig_fasta(run_dir, sweep_dir, best_contig.contig_id, best_fa)

        # Combined FASTA
        rec = next(SeqIO.parse(str(best_fa), "fasta"))
        rec.id = f"{sample}|{int(depth_target)}x|{best_contig.contig_id}"
        rec.name = rec.id
        rec.description = DEFAULT_SAMPLE_SPECIES.get(sample, "")
        best_fasta_records.append(rec)

        # Per-contig rows
        for c in contigs:
            gs = c.gene_stats
            contigs_out_rows.append({
                "sample": sample,
                "depth_target": depth_target,
                "run_dir": str(run_dir),
                "selected_sweep_reads": selected_level,
                "contig_id": c.contig_id,
                "length": c.length,
                "mean_depth": f"{c.mean_depth:.6f}",
                "covered_frac": f"{c.covered_frac:.6f}",
                "gene_completeness": f"{c.gene_completeness:.6f}",
                "pcg_count": c.pcg_count,
                "rrna_count": c.rrna_count,
                "trna_count": c.trna_count,
                "missing_pcgs": ",".join(gs.missing_pcgs) if gs else ",".join(EXPECTED_PCGS),
            })

        runs_out_rows.append({
            "sample": sample,
            "species": DEFAULT_SAMPLE_SPECIES.get(sample, ""),
            "depth_target": depth_target,
            "run_dir": str(run_dir),
            "selected_sweep_reads": selected_level,
            "best_contig_id": best_contig.contig_id,
            "best_length": best_contig.length,
            "best_mean_depth": f"{best_contig.mean_depth:.6f}",
            "best_covered_frac": f"{best_contig.covered_frac:.6f}",
            "best_gene_completeness": f"{best_contig.gene_completeness:.6f}",
            "best_pcg_count": best_contig.pcg_count,
            "best_rrna_count": best_contig.rrna_count,
            "best_trna_count": best_contig.trna_count,
            "run_best_score": f"{best_score:.6f}",
            "status": "OK",
            "best_contig_fasta": str(best_fa),
        })

    # Write combined FASTA
    best_all_fa = out_root / "best_candidates.fasta"
    if best_fasta_records:
        SeqIO.write(best_fasta_records, str(best_all_fa), "fasta")

    # Write TSVs
    contigs_tsv = out_root / "contigs.tsv"
    runs_tsv = out_root / "runs.tsv"

    write_tsv(
        contigs_tsv,
        contigs_out_rows,
        fieldnames=[
            "sample", "depth_target", "run_dir", "selected_sweep_reads",
            "contig_id", "length", "mean_depth", "covered_frac",
            "gene_completeness", "pcg_count", "rrna_count", "trna_count",
            "missing_pcgs",
        ],
    )

    write_tsv(
        runs_tsv,
        runs_out_rows,
        fieldnames=[
            "sample", "species", "depth_target", "run_dir", "selected_sweep_reads",
            "best_contig_id", "best_length", "best_mean_depth", "best_covered_frac",
            "best_gene_completeness", "best_pcg_count", "best_rrna_count", "best_trna_count",
            "run_best_score", "status", "best_contig_fasta",
        ],
    )

    # Best run per sample across depth targets
    best_per_sample: Dict[str, Dict[str, object]] = {}
    for rr in runs_out_rows:
        if rr.get("status") != "OK":
            continue
        s = str(rr["sample"])
        sc = float(rr["run_best_score"])
        if s not in best_per_sample or sc > float(best_per_sample[s]["run_best_score"]):
            best_per_sample[s] = rr

    best_per_sample_tsv = out_root / "best_per_sample.tsv"
    write_tsv(
        best_per_sample_tsv,
        list(best_per_sample.values()),
        fieldnames=[
            "sample", "species", "depth_target", "run_dir", "selected_sweep_reads",
            "best_contig_id", "best_length", "best_mean_depth", "best_covered_frac",
            "best_gene_completeness", "best_pcg_count", "best_rrna_count", "best_trna_count",
            "run_best_score", "status", "best_contig_fasta",
        ],
    )

    # Optional BLAST
    blast_hits: List[Dict[str, object]] = []
    if not args.skip_blast:
        db_dir = out_root / "blast_db"
        db_dir.mkdir(parents=True, exist_ok=True)

        if args.blast_ref_fasta:
            ref_fa = Path(args.blast_ref_fasta).resolve()
            if not ref_fa.exists():
                die(f"--blast_ref_fasta not found: {ref_fa}")
        else:
            ref_samples = [x.strip() for x in args.blast_ref_samples.split(",") if x.strip()]
            ref_fa = db_dir / "refs_from_batch.fasta"
            build_ref_fasta_from_samples(runs_out_rows, ref_samples, ref_fa)

        db_prefix = db_dir / "mt_refs_db"
        make_blast_db(ref_fa, db_prefix, force=args.force)

        for rr in runs_out_rows:
            if rr.get("status") != "OK":
                continue
            qfa = Path(str(rr["best_contig_fasta"]))
            hit = blast_top_hit(qfa, db_prefix, threads=args.threads)
            if not hit:
                blast_hits.append({
                    "sample": rr["sample"],
                    "depth_target": rr["depth_target"],
                    "best_contig_id": rr["best_contig_id"],
                    "best_contig_fasta": rr["best_contig_fasta"],
                    "blast_status": "NO_HIT_OR_FAILED",
                })
                continue

            blast_hits.append({
                "sample": rr["sample"],
                "depth_target": rr["depth_target"],
                "best_contig_id": rr["best_contig_id"],
                "best_contig_fasta": rr["best_contig_fasta"],
                "blast_status": "OK",
                **hit,
            })

        blast_tsv = out_root / "blast_hits.tsv"
        fields = [
            "sample", "depth_target", "best_contig_id", "best_contig_fasta",
            "blast_status", "qseqid", "sseqid", "pident", "aln_len", "qlen", "slen", "evalue", "bitscore", "stitle"
        ]
        write_tsv(blast_tsv, blast_hits, fields)

    print("MT post-processing complete.")
    print(f"  contigs:         {contigs_tsv}")
    print(f"  runs:            {runs_tsv}")
    print(f"  best_per_sample: {best_per_sample_tsv}")
    print(f"  best_fasta:      {best_all_fa}")
    if not args.skip_blast:
        print(f"  blast_hits:      {out_root / 'blast_hits.tsv'}")


if __name__ == "__main__":
    main()
