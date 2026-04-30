#!/usr/bin/env python3
"""
tree + per-partition GC-content + root-to-tip lollipop plot.

This script can either:
  1. compute GC content directly from an alignment plus an IQ-TREE partition file, or
  2. read a precomputed partition-GC table.

Main outputs:
  - figure with tree | labels | GC panel | root-to-tip distance
  - wide GC-content table (one row per sequence/header, one column per partition)
  - long GC-content table
  - matched plot-data table
  - fastest-branches table

GC panel modes:
  - partitions : one box column per partition (full matrix)
  - overall    : one GC strip using GC_all (or a chosen overall column)
  - single     : one GC strip for a single selected partition/gene

This is useful for BUSCO-scale datasets, where plotting thousands of partition tiles
would be unreadable. In those cases, use:
  --gc-panel-mode overall
or
  --gc-panel-mode single --single-partition BUSCO_xxx
  
USAGE
python3 plot_tree_partition_gc_lollipop.py \
  -t input_tree \
  -a alignment.fasta \
  -p partitions.txt \
  --gc-panel-mode overall \
  --top-branches X \
  -o gc_overall
"""

import argparse
import re
import sys
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm, Normalize
from Bio import Phylo, SeqIO

STYLE = {
    "font_family": "sans-serif",
    "font_sans": ["Arial", "Helvetica", "DejaVu Sans"],
    "color_tree": "#596275",
    "color_guide": "#f1f2f6",
    "color_fast": "#e17055",
}

# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

def set_style() -> None:
    mpl.rcParams.update({
        "font.family": STYLE["font_family"],
        "font.sans-serif": STYLE["font_sans"],
        "axes.linewidth": 0.0,
        "xtick.major.width": 0.0,
        "ytick.major.width": 0.0,
        "pdf.fonttype": 42,
        "axes.grid": False,
    })


def clean_name(s) -> str:
    return str(s).strip().replace(" ", "_").replace("'", "").replace('"', "")


def parse_label_display(label: str) -> str:
    label = clean_name(label)
    parts = label.split("_")
    if len(parts) >= 3 and parts[0] in ("GCA", "GCF"):
        return " ".join(parts[2:])
    return label.replace("_", " ")


def extract_accession(name: str):
    match = re.search(r"(GC[AF]_[0-9]+\.[0-9]+)", str(name))
    return match.group(1) if match else None


def read_table_auto(path: Path, sep: str = "auto") -> pd.DataFrame:
    if sep == "\\t":
        sep = "\t"
    if sep != "auto":
        return pd.read_csv(path, sep=sep)

    suffix = path.suffix.lower()
    if suffix in {".tsv", ".tab"}:
        return pd.read_csv(path, sep="\t")
    if suffix == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep=None, engine="python")


def split_keywords(values) -> list[str]:
    if not values:
        return []
    out = []
    for value in values:
        if value is None:
            continue
        for x in str(value).split(","):
            x = x.strip()
            if x:
                out.append(x)
    return out


def load_keyword_file(path: Path | None) -> list[str]:
    if path is None:
        return []
    keywords = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            keywords.append(line)
    return keywords


def load_name_map(path: Path | None, key_col: str, value_col: str, sep: str) -> dict[str, str]:
    if path is None:
        return {}
    df = read_table_auto(path, sep=sep)
    missing = [c for c in (key_col, value_col) if c not in df.columns]
    if missing:
        raise ValueError(
            f"Mapping file {path} is missing column(s): {', '.join(missing)}. "
            f"Available columns: {', '.join(map(str, df.columns))}"
        )
    name_map = {}
    for _, row in df[[key_col, value_col]].dropna().iterrows():
        key = clean_name(row[key_col])
        val = clean_name(row[value_col])
        if key:
            name_map[key] = val
    return name_map


def read_newick(path: Path):
    try:
        return Phylo.read(str(path), "newick")
    except Exception:
        text = path.read_text()
        if ";" not in text:
            raise
        return Phylo.read(StringIO(text[: text.find(";") + 1]), "newick")


# -----------------------------------------------------------------------------
# Partition parsing and GC computation
# -----------------------------------------------------------------------------

def parse_partition_ranges(range_text: str) -> list[int]:
    positions = []
    for chunk in str(range_text).split(","):
        chunk = chunk.strip().rstrip(";")
        if not chunk:
            continue

        step = 1
        if "\\" in chunk:
            chunk, step_part = chunk.split("\\", 1)
            step = int(step_part)
        elif "/" in chunk:
            chunk, step_part = chunk.split("/", 1)
            step = int(step_part)

        if "-" in chunk:
            start, end = chunk.split("-", 1)
            start_i, end_i = int(start), int(end)
            positions.extend(list(range(start_i - 1, end_i, step)))
        else:
            positions.append(int(chunk) - 1)

    return sorted(set(positions))


def load_iqtree_partitions(path: Path) -> list[dict]:
    partitions = []
    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            low = line.lower()
            if low in {"begin sets;", "end;", "endblock;"}:
                continue
            if low.startswith("charset "):
                line = re.sub(r"^charset\s+", "", line, flags=re.IGNORECASE)

            if "=" not in line:
                continue

            left, right = line.split("=", 1)
            left = left.strip().rstrip(",")
            right = right.strip().rstrip(";")

            if "," in left:
                part_name = left.split(",")[-1].strip()
            else:
                part_name = left.strip()

            part_name = clean_name(part_name)
            positions = parse_partition_ranges(right)
            if not positions:
                continue

            partitions.append({
                "partition": part_name,
                "n_sites": len(positions),
                "positions": positions,
                "range_text": right,
            })

    if not partitions:
        raise ValueError(f"No partitions could be parsed from: {path}")

    return partitions


def compute_gc_for_sequence_substring(seq: str) -> tuple[float, int, int]:
    seq = str(seq).upper()
    gc_sites = 0
    valid_sites = 0
    for base in seq:
        if base in {"A", "C", "G", "T"}:
            valid_sites += 1
            if base in {"G", "C"}:
                gc_sites += 1
    gc = np.nan if valid_sites == 0 else gc_sites / valid_sites
    return gc, valid_sites, gc_sites


def compute_partition_gc_table(alignment_path: Path, alignment_format: str, partition_path: Path):
    records = list(SeqIO.parse(str(alignment_path), alignment_format))
    if not records:
        raise ValueError(f"No sequences read from alignment: {alignment_path}")

    lengths = {len(rec.seq) for rec in records}
    if len(lengths) != 1:
        raise ValueError("Alignment sequences do not all have the same length.")
    aln_len = list(lengths)[0]

    partitions = load_iqtree_partitions(partition_path)
    for part in partitions:
        max_pos = max(part["positions"])
        if max_pos >= aln_len:
            raise ValueError(
                f"Partition '{part['partition']}' uses site {max_pos + 1}, "
                f"but alignment length is only {aln_len}."
            )

    rows_wide = []
    rows_long = []
    for rec in records:
        seq = str(rec.seq)
        taxon = clean_name(rec.id)
        wide = {
            "taxon": taxon,
            "header": rec.id,
            "sequence_length": len(seq),
        }

        total_gc, total_valid, total_gc_sites = compute_gc_for_sequence_substring(seq)
        wide["GC_all"] = total_gc
        wide["valid_sites_all"] = total_valid
        wide["gc_sites_all"] = total_gc_sites

        for part in partitions:
            subseq = "".join(seq[i] for i in part["positions"])
            gc, valid_sites, gc_sites = compute_gc_for_sequence_substring(subseq)
            wide[part["partition"]] = gc
            rows_long.append({
                "taxon": taxon,
                "header": rec.id,
                "partition": part["partition"],
                "gc_content": gc,
                "valid_sites": valid_sites,
                "gc_sites": gc_sites,
                "n_sites_in_partition": part["n_sites"],
            })

        rows_wide.append(wide)

    wide_df = pd.DataFrame(rows_wide)
    long_df = pd.DataFrame(rows_long)
    part_df = pd.DataFrame([
        {"partition": p["partition"], "n_sites": p["n_sites"], "range_text": p["range_text"]}
        for p in partitions
    ])
    partition_order = [p["partition"] for p in partitions]
    return wide_df, long_df, part_df, partition_order


# -----------------------------------------------------------------------------
# Tree handling
# -----------------------------------------------------------------------------

def rename_tree_tips(tree, name_map: dict[str, str], use_accession_extraction: bool = True):
    if not name_map:
        return tree
    for tip in tree.get_terminals():
        original = clean_name(tip.name)
        if original in name_map:
            tip.name = name_map[original]
            continue
        if use_accession_extraction:
            acc = extract_accession(original)
            if acc and acc in name_map:
                tip.name = name_map[acc]
    return tree


def filter_tree_by_keywords(tree, keep_keywords: list[str]):
    if not keep_keywords:
        return tree
    keep_lower = [k.lower() for k in keep_keywords]
    terminals = list(tree.get_terminals())
    keep = []
    for tip in terminals:
        label = clean_name(tip.name).lower()
        if any(k in label for k in keep_lower):
            keep.append(tip)
    if not keep:
        raise ValueError("No tree tips matched --keep-keywords / --keep-keyword-file.")
    keep_ids = {id(t) for t in keep}
    for tip in terminals:
        if id(tip) not in keep_ids:
            tree.prune(tip)
    return tree


def get_layout_coords(tree):
    terminals = tree.get_terminals()
    y_map = {id(t): i for i, t in enumerate(terminals)}

    def calc_y(clade):
        if id(clade) in y_map:
            return
        for child in clade.clades:
            calc_y(child)
        y_map[id(clade)] = np.mean([y_map[id(child)] for child in clade.clades])

    calc_y(tree.root)

    x_map = {}
    def calc_x(clade, curr_x):
        x_map[id(clade)] = curr_x
        for child in clade.clades:
            bl = child.branch_length if child.branch_length else 0.0
            calc_x(child, curr_x + bl)

    calc_x(tree.root, 0.0)
    return x_map, y_map, terminals


def get_segments(tree, x_map, y_map):
    h_lines, v_lines, edge_data = [], [], []

    def walk(clade):
        x_start = x_map[id(clade)]
        if not clade.is_terminal():
            ys = [y_map[id(c)] for c in clade.clades]
            v_lines.append([(x_start, min(ys)), (x_start, max(ys))])
        for child in clade.clades:
            x_end, y_end = x_map[id(child)], y_map[id(child)]
            h_lines.append([(x_start, y_end), (x_end, y_end)])
            bl = child.branch_length if child.branch_length else 0.0
            edge_data.append({"clade": child, "x0": x_start, "x1": x_end, "y": y_end, "bl": bl})
            walk(child)

    walk(tree.root)
    return h_lines, v_lines, edge_data


# -----------------------------------------------------------------------------
# Matching GC table to tree tips
# -----------------------------------------------------------------------------

def build_reverse_map(name_map: dict[str, str]) -> dict[str, list[str]]:
    reverse = {}
    for source, target in name_map.items():
        reverse.setdefault(clean_name(target), []).append(clean_name(source))
    return reverse


def collapse_partition_rows(rows: pd.DataFrame, gc_columns: list[str], agg: str) -> dict[str, float]:
    out = {}
    for col in gc_columns:
        vals = pd.to_numeric(rows[col], errors="coerce").dropna()
        if vals.empty:
            out[col] = np.nan
            continue
        if agg == "first":
            out[col] = float(vals.iloc[0])
        elif agg == "mean":
            out[col] = float(vals.mean())
        elif agg == "median":
            out[col] = float(vals.median())
        elif agg == "min":
            out[col] = float(vals.min())
        elif agg == "max":
            out[col] = float(vals.max())
        else:
            raise ValueError(f"Unsupported --duplicate-agg value: {agg}")

        if pd.notna(out[col]) and out[col] > 1.5:
            out[col] = out[col] / 100.0
    return out


def find_gc_row_match(
    tip_name: str,
    gc_df: pd.DataFrame,
    gc_columns: list[str],
    name_map: dict[str, str],
    reverse_map: dict[str, list[str]],
    duplicate_agg: str,
    allow_substring_match: bool,
):
    current = clean_name(tip_name)
    candidates = [current]

    acc = extract_accession(current)
    if acc:
        candidates.append(acc)

    if current in reverse_map:
        candidates.extend(reverse_map[current])

    if current in name_map:
        candidates.append(clean_name(name_map[current]))

    seen = set()
    candidates = [x for x in candidates if not (x in seen or seen.add(x))]

    rows = gc_df[gc_df["clean_tax"].isin(candidates)]
    if not rows.empty:
        return collapse_partition_rows(rows, gc_columns, duplicate_agg), "exact_or_mapped"

    if allow_substring_match and len(current) > 5:
        rows = gc_df[gc_df["clean_tax"].str.contains(current, regex=False, na=False)]
        if not rows.empty:
            return collapse_partition_rows(rows, gc_columns, duplicate_agg), "substring"

    return {col: np.nan for col in gc_columns}, "unmatched"


def build_plot_data(
    terminals,
    y_map,
    x_map,
    gc_df: pd.DataFrame,
    gc_columns: list[str],
    name_map: dict[str, str],
    duplicate_agg: str,
    allow_substring_match: bool,
):
    reverse_map = build_reverse_map(name_map)
    r2t_source = {id(t): x_map[id(t)] for t in terminals}
    mean_r2t = np.mean(list(r2t_source.values())) if r2t_source else 1.0

    rows = []
    match_counts = {"exact_or_mapped": 0, "substring": 0, "unmatched": 0}

    for tip in terminals:
        gc_vals, match_type = find_gc_row_match(
            tip_name=tip.name,
            gc_df=gc_df,
            gc_columns=gc_columns,
            name_map=name_map,
            reverse_map=reverse_map,
            duplicate_agg=duplicate_agg,
            allow_substring_match=allow_substring_match,
        )
        match_counts[match_type] += 1

        row = {
            "id": clean_name(tip.name),
            "display_name": parse_label_display(clean_name(tip.name)),
            "y_coord": y_map[id(tip)],
            "r2t_val": r2t_source[id(tip)],
            "match_type": match_type,
        }
        row.update(gc_vals)
        rows.append(row)

    merged_df = pd.DataFrame(rows)
    return merged_df, mean_r2t, match_counts


# -----------------------------------------------------------------------------
# Plot functions
# -----------------------------------------------------------------------------

def draw_zebra_background(ax, n_tips, y_map, terminals):
    ax.set_ylim(-0.5, n_tips - 0.5)
    ax.invert_yaxis()
    for tip in terminals:
        y = y_map[id(tip)]
        if int(y) % 2 == 1:
            ax.axhspan(y - 0.5, y + 0.5, color=STYLE["color_guide"], alpha=0.6, zorder=-10)


def plot_tree(ax, h_lines, v_lines, n_tips, y_map, terminals):
    draw_zebra_background(ax, n_tips, y_map, terminals)
    ax.add_collection(LineCollection(v_lines, colors=STYLE["color_tree"], lw=1.0, alpha=0.9))
    ax.add_collection(LineCollection(h_lines, colors=STYLE["color_tree"], lw=1.2, alpha=0.9))
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.spines["bottom"].set_color("#bdc3c7")
    ax.tick_params(axis="x", labelsize=8, length=4, color="#bdc3c7", labelcolor="#2d3436")
    ax.set_xlabel("Substitutions / site", fontsize=9, fontweight="bold", color="#2d3436", labelpad=8)
    ax.axis("on")
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_yticks([])


def mark_fastest_branches(ax, edge_data, top_n):
    if top_n <= 0:
        return pd.DataFrame()
    valid_edges = [e for e in edge_data if e["bl"] > 0]
    if not valid_edges:
        return pd.DataFrame()
    mean_bl = np.mean([e["bl"] for e in valid_edges])
    top_edges = sorted(valid_edges, key=lambda x: x["bl"], reverse=True)[:top_n]

    table_rows = []
    for rank, edge in enumerate(top_edges, 1):
        ax.plot([edge["x0"], edge["x1"]], [edge["y"], edge["y"]], color=STYLE["color_fast"], lw=2.5, zorder=10, alpha=0.8)
        mid_x = (edge["x0"] + edge["x1"]) / 2
        ax.scatter(mid_x, edge["y"], s=100, facecolor="white", edgecolor=STYLE["color_fast"], zorder=11, lw=1.5)
        ax.text(mid_x, edge["y"], str(rank), color=STYLE["color_fast"], fontsize=7, fontweight="bold", ha="center", va="center", zorder=12)
        terms = edge["clade"].get_terminals()
        if len(terms) == 1:
            desc = f"Tip: {parse_label_display(clean_name(terms[0].name))}"
        else:
            desc = f"Node: Clade of {len(terms)} tips (e.g., {parse_label_display(clean_name(terms[0].name))})"
        table_rows.append({"Rank": rank, "Branch_Length": edge["bl"], "Relative_Rate": round(edge["bl"] / mean_bl, 2), "Description": desc})
    return pd.DataFrame(table_rows)


def plot_labels(ax, tips, y_map, fontsize):
    draw_zebra_background(ax, len(tips), y_map, tips)
    ax.axis("off")
    for tip in tips:
        y = y_map[id(tip)]
        name = parse_label_display(clean_name(tip.name))
        ax.text(0.0, y, name, va="center", fontsize=fontsize, ha="left", color="#2d3436")


def get_gc_norm_and_ticks(values: np.ndarray):
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return Normalize(vmin=0, vmax=1), [0.0, 0.5, 1.0]
    min_v, med_v, max_v = float(np.min(values)), float(np.median(values)), float(np.max(values))
    if min_v < med_v < max_v:
        norm = TwoSlopeNorm(vmin=min_v, vcenter=med_v, vmax=max_v)
    else:
        vmax = max_v if max_v > min_v else min_v + 0.01
        norm = Normalize(vmin=min_v, vmax=vmax)
    return norm, [min_v, med_v, max_v]


def plot_gc_strip(ax, merged_df, gc_col, y_map, terminals, gc_cmap_name: str, title: str):
    draw_zebra_background(ax, len(terminals), y_map, terminals)
    ax.set_xlim(0, 1)
    ax.axis("off")
    ax.set_title(title, fontsize=9, fontweight="bold", color="#2d3436", pad=10)

    vals = pd.to_numeric(merged_df[gc_col], errors="coerce").dropna().values
    norm, ticks = get_gc_norm_and_ticks(vals.astype(float) if len(vals) else np.array([]))
    cmap = plt.get_cmap(gc_cmap_name)

    for _, row in merged_df.iterrows():
        val = row.get(gc_col, np.nan)
        if pd.isna(val):
            continue
        y = row["y_coord"]
        color = cmap(norm(val))
        rect = plt.Rectangle((0.1, y - 0.4), 0.8, 0.8, facecolor=color, edgecolor='none')
        ax.add_patch(rect)

    return norm, cmap, ticks


def plot_partition_gc_boxes(ax, merged_df, partition_cols, y_map, terminals, gc_cmap_name: str):
    draw_zebra_background(ax, len(terminals), y_map, terminals)
    n_parts = len(partition_cols)
    ax.set_xlim(-0.5, n_parts - 0.5)
    ax.set_title("Partition GC", fontsize=9, fontweight="bold", color="#2d3436", pad=10)

    vals = merged_df[partition_cols].to_numpy(dtype=float).ravel()
    norm, ticks = get_gc_norm_and_ticks(vals)
    cmap = plt.get_cmap(gc_cmap_name)

    for _, row in merged_df.iterrows():
        y = row["y_coord"]
        for j, part in enumerate(partition_cols):
            val = row[part]
            if pd.isna(val):
                face = "#ecf0f1"
            else:
                face = cmap(norm(val))
            rect = plt.Rectangle((j - 0.45, y - 0.4), 0.9, 0.8, facecolor=face, edgecolor="white", linewidth=0.3)
            ax.add_patch(rect)

    ax.set_xticks(range(n_parts))
    font_size = 8 if n_parts <= 12 else 7 if n_parts <= 25 else 6
    ax.set_xticklabels(partition_cols, rotation=90, fontsize=font_size)
    ax.tick_params(axis="x", length=0)
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.spines["bottom"].set_color("#bdc3c7")

    return norm, cmap, ticks


def plot_rates(ax, merged_df, y_map, terminals, mean_r2t, rate_cmap_name: str):
    draw_zebra_background(ax, len(terminals), y_map, terminals)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.spines["bottom"].set_color("#bdc3c7")
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=8, length=4, color="#bdc3c7", labelcolor="#2d3436")
    ax.set_title("Root-to-Tip Distance\n(Dotted line = Mean)", fontsize=9, fontweight="bold", color="#2d3436", pad=10)

    rates = merged_df["r2t_val"].dropna().values
    if len(rates) == 0:
        return
    min_r, max_r = min(rates), max(rates)
    span = max_r - min_r if max_r > min_r else 0.1
    ax.set_xlim(min_r - (span * 0.05), max_r + (span * 0.05))
    norm = Normalize(vmin=min_r, vmax=max_r)
    cmap = plt.get_cmap(rate_cmap_name)
    ax.axvline(mean_r2t, color="#95a5a6", ls=":", lw=1.0, zorder=1)
    baseline = min_r - (span * 0.05)
    for _, row in merged_df.iterrows():
        y = row["y_coord"]
        val = row["r2t_val"]
        if pd.isna(val):
            continue
        col = cmap(norm(val))
        ax.hlines(y, baseline, val, color="#b2bec3", lw=1.0, zorder=2)
        ax.scatter(val, y, s=30, color=col, edgecolor="none", zorder=3, alpha=0.9)


def draw_footer(fig, gc_norm, gc_cmap, gc_ticks, fig_height_inches, gc_label):
    to_frac = lambda x: x / fig_height_inches

    rect_leg = [0.05, to_frac(0.2), 0.25, to_frac(0.5)]
    ax_leg = fig.add_axes(rect_leg)
    ax_leg.axis("off")
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="white", markeredgecolor=STYLE["color_fast"], markersize=8, label="Fastest Branch"),
        Line2D([0], [0], color=STYLE["color_fast"], lw=2.5, label="High Rate Branch"),
    ]
    ax_leg.legend(handles=legend_elements, loc="lower left", fontsize=8, frameon=False, ncol=1)

    rect_cb = [0.40, to_frac(0.35), 0.2, to_frac(0.12)]
    ax_cb = fig.add_axes(rect_cb)
    cb = plt.colorbar(mpl.cm.ScalarMappable(norm=gc_norm, cmap=gc_cmap), cax=ax_cb, orientation="horizontal")
    formatted_ticks = []
    for t in gc_ticks:
        if max(gc_ticks) <= 1.0:
            formatted_ticks.append(f"{t * 100:.0f}%")
        else:
            formatted_ticks.append(f"{t:.1f}")
    cb.set_ticks(gc_ticks)
    cb.set_ticklabels(formatted_ticks)
    cb.ax.tick_params(labelsize=7, length=0, color="#2d3436")
    cb.outline.set_visible(False)
    cb.set_label(gc_label, fontsize=8, labelpad=4, color="#2d3436")


# -----------------------------------------------------------------------------
# Input preparation
# -----------------------------------------------------------------------------

def prepare_precomputed_gc_table(path: Path, sep: str, taxon_col: str, partition_cols: list[str] | None):
    df = read_table_auto(path, sep=sep)
    if taxon_col not in df.columns:
        raise ValueError(f"Precomputed GC table is missing taxon column '{taxon_col}'. Available columns: {', '.join(map(str, df.columns))}")

    if partition_cols is None or len(partition_cols) == 0:
        partition_cols = [c for c in df.columns if c != taxon_col and c != 'GC_all' and not str(c).startswith("valid_") and not str(c).startswith("gc_")]
    missing = [c for c in partition_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Precomputed GC table is missing partition column(s): {', '.join(missing)}")

    out = df.copy()
    out["taxon"] = out[taxon_col].apply(clean_name)
    out["header"] = out[taxon_col]
    out["clean_tax"] = out["taxon"]

    if "GC_all" not in out.columns and partition_cols:
        out["GC_all"] = out[partition_cols].apply(pd.to_numeric, errors='coerce').mean(axis=1)

    numeric_cols = list(set(partition_cols + (["GC_all"] if "GC_all" in out.columns else [])))
    for col in numeric_cols:
        out[col] = pd.to_numeric(out[col], errors="coerce")
        mask = out[col] > 1.5
        out.loc[mask, col] = out.loc[mask, col] / 100.0
    return out, partition_cols


def resolve_gc_mode_columns(args, gc_wide_df, partition_cols):
    gc_mode = args.gc_panel_mode

    if gc_mode == 'partitions':
        if len(partition_cols) == 0:
            raise ValueError("No partition columns available for --gc-panel-mode partitions.")
        return partition_cols, 'partitions', 'Partition GC', max(1.6, min(8.0, 0.38 * len(partition_cols)))

    if gc_mode == 'overall':
        gc_col = args.overall_gc_col
        if gc_col not in gc_wide_df.columns:
            raise ValueError(f"Requested overall GC column '{gc_col}' not found. Available columns include: {', '.join(map(str, gc_wide_df.columns[:20]))}...")
        return [gc_col], 'single_strip', args.gc_panel_title or 'GC Content', 0.7

    if gc_mode == 'single':
        if args.single_partition is None:
            raise ValueError("--gc-panel-mode single requires --single-partition")
        gc_col = args.single_partition
        if gc_col not in gc_wide_df.columns:
            preview = ', '.join(map(str, partition_cols[:20]))
            raise ValueError(f"Requested single partition '{gc_col}' not found. First available partitions: {preview}")
        return [gc_col], 'single_strip', args.gc_panel_title or f'GC: {gc_col}', 0.7

    raise ValueError(f"Unsupported --gc-panel-mode: {gc_mode}")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Tree + GC panel + root-to-tip lollipop plot.")

    parser.add_argument("-t", "--tree", required=True, type=Path, help="Input Newick tree.")
    parser.add_argument("-o", "--outprefix", default="tree_partition_gc_lollipop", help="Output prefix.")

    # Mode 1: compute from alignment
    parser.add_argument("-a", "--alignment", type=Path, default=None, help="Input alignment used to compute GC.")
    parser.add_argument("-p", "--partition-file", type=Path, default=None, help="IQ-TREE partition file.")
    parser.add_argument("--alignment-format", default="fasta", help="Biopython SeqIO format for alignment. Default: fasta")

    # Mode 2: precomputed table
    parser.add_argument("--partition-gc-table", type=Path, default=None, help="Precomputed partition-GC table.")
    parser.add_argument("--gc-sep", default="auto", help="Delimiter for precomputed GC table. Default: auto")
    parser.add_argument("--taxon-col", default="taxon", help="Taxon column in the precomputed GC table. Default: taxon")
    parser.add_argument("--partition-cols", nargs="*", default=None, help="Partition columns in the precomputed GC table. If omitted, infer all non-taxon columns except GC_all.")

    # Mapping and filtering
    parser.add_argument("--map", type=Path, default=None, help="Optional accession/name mapping TSV/CSV.")
    parser.add_argument("--map-key-col", default="accession", help="Key column in --map. Default: accession")
    parser.add_argument("--map-value-col", default="species", help="Value column in --map. Default: species")
    parser.add_argument("--map-sep", default="auto", help="Delimiter for --map table. Default: auto")
    parser.add_argument("--no-accession-extraction", action="store_true", help="Do not extract GCA/GCF accessions from longer tip names during matching.")
    parser.add_argument("--keep-keywords", nargs="*", default=[], help="Optional comma-separated or space-separated keywords. If absent, all tips are kept.")
    parser.add_argument("--keep-keyword-file", type=Path, default=None, help="Optional file with one keep keyword per line.")
    parser.add_argument("--duplicate-agg", choices=["first", "mean", "median", "min", "max"], default="first", help="How to collapse duplicated alignment/GC rows per taxon.")
    parser.add_argument("--allow-substring-match", action="store_true", help="Enable fallback substring matching between tip labels and taxon labels.")

    # GC panel mode
    parser.add_argument("--gc-panel-mode", choices=["partitions", "overall", "single"], default="partitions", help="How to draw the GC panel. Default: partitions")
    parser.add_argument("--single-partition", default=None, help="Partition/gene name to plot when --gc-panel-mode single")
    parser.add_argument("--overall-gc-col", default="GC_all", help="Column to use when --gc-panel-mode overall. Default: GC_all")
    parser.add_argument("--gc-panel-title", default=None, help="Optional custom title for the GC panel")

    # Plot options
    parser.add_argument("--top-branches", type=int, default=10, help="Number of longest branches to mark. Use 0 to disable.")
    parser.add_argument("--title", default="Phylogenomic Profile (n={n})", help="Plot title. Supports {n}, {tree}, {n_partitions}, {gc_mode}.")
    parser.add_argument("--gc-cmap", default="BrBG", help="Matplotlib colormap for GC panel. Default: BrBG")
    parser.add_argument("--rate-cmap", default="plasma", help="Matplotlib colormap for root-to-tip points. Default: plasma")
    parser.add_argument("--fig-width", type=float, default=None, help="Figure width in inches. Default: auto")
    parser.add_argument("--row-height", type=float, default=0.20, help="Figure height per tip in inches.")
    parser.add_argument("--min-height", type=float, default=8.0, help="Minimum figure height in inches.")
    parser.add_argument("--label-fontsize", type=float, default=9.0, help="Tip label font size.")
    parser.add_argument("--formats", nargs="+", default=["png", "pdf", "svg"], choices=["png", "pdf", "svg"], help="Output figure formats.")
    parser.add_argument("--dpi", type=int, default=300, help="PNG dpi. Default: 300")

    return parser.parse_args()


def main():
    args = parse_args()
    set_style()

    if not args.tree.exists():
        raise FileNotFoundError(f"Tree file not found: {args.tree}")

    using_alignment_mode = args.alignment is not None or args.partition_file is not None
    using_table_mode = args.partition_gc_table is not None

    if using_alignment_mode and using_table_mode:
        raise ValueError("Use either alignment+partition-file OR --partition-gc-table, not both.")
    if not using_alignment_mode and not using_table_mode:
        raise ValueError("Provide either --alignment with --partition-file, or --partition-gc-table.")
    if using_alignment_mode and (args.alignment is None or args.partition_file is None):
        raise ValueError("When computing GC from alignment, both --alignment and --partition-file are required.")

    tree = read_newick(args.tree)
    name_map = load_name_map(args.map, args.map_key_col, args.map_value_col, args.map_sep)
    tree = rename_tree_tips(tree, name_map=name_map, use_accession_extraction=not args.no_accession_extraction)

    keep_keywords = split_keywords(args.keep_keywords) + load_keyword_file(args.keep_keyword_file)
    tree = filter_tree_by_keywords(tree, keep_keywords)

    x_map, y_map, terminals = get_layout_coords(tree)
    h_lines, v_lines, edge_data = get_segments(tree, x_map, y_map)

    if using_alignment_mode:
        if not args.alignment.exists():
            raise FileNotFoundError(f"Alignment file not found: {args.alignment}")
        if not args.partition_file.exists():
            raise FileNotFoundError(f"Partition file not found: {args.partition_file}")

        gc_wide_df, gc_long_df, partition_info_df, partition_cols = compute_partition_gc_table(
            alignment_path=args.alignment,
            alignment_format=args.alignment_format,
            partition_path=args.partition_file,
        )
        gc_wide_df["clean_tax"] = gc_wide_df["taxon"].apply(clean_name)
    else:
        if not args.partition_gc_table.exists():
            raise FileNotFoundError(f"Precomputed GC table not found: {args.partition_gc_table}")
        gc_wide_df, partition_cols = prepare_precomputed_gc_table(
            path=args.partition_gc_table,
            sep=args.gc_sep,
            taxon_col=args.taxon_col,
            partition_cols=args.partition_cols,
        )
        gc_long_df = gc_wide_df.melt(
            id_vars=[c for c in ["taxon", "header", "clean_tax"] if c in gc_wide_df.columns],
            value_vars=partition_cols,
            var_name="partition",
            value_name="gc_content",
        )
        partition_info_df = pd.DataFrame({"partition": partition_cols})

    gc_columns_for_matching = list(dict.fromkeys(["GC_all"] + partition_cols))
    merged_df, mean_r2t, match_counts = build_plot_data(
        terminals=terminals,
        y_map=y_map,
        x_map=x_map,
        gc_df=gc_wide_df,
        gc_columns=gc_columns_for_matching,
        name_map=name_map,
        duplicate_agg=args.duplicate_agg,
        allow_substring_match=args.allow_substring_match,
    )

    selected_gc_cols, gc_panel_kind, gc_panel_title, gc_panel_width = resolve_gc_mode_columns(args, gc_wide_df, partition_cols)

    n_tips = len(terminals)
    matched = n_tips - match_counts.get("unmatched", 0)
    print(f"Tips plotted: {n_tips}")
    print(
        "GC matches: "
        f"{matched}/{n_tips} "
        f"(exact/mapped={match_counts.get('exact_or_mapped', 0)}, "
        f"substring={match_counts.get('substring', 0)}, "
        f"unmatched={match_counts.get('unmatched', 0)})"
    )
    print(f"GC panel mode: {args.gc_panel_mode}")
    if args.gc_panel_mode == 'partitions':
        print(f"Partitions plotted: {len(partition_cols)}")
    elif args.gc_panel_mode == 'single':
        print(f"Single partition plotted: {args.single_partition}")
    else:
        print(f"Overall GC column plotted: {args.overall_gc_col}")

    header_in = 1.1
    footer_in = 1.6
    fig_h = max(args.min_height, (n_tips * args.row_height) + header_in + footer_in)
    fig_w = args.fig_width if args.fig_width is not None else 4 + 2 + gc_panel_width + 2.5 + 1.0
    top_frac = 1.0 - (header_in / fig_h)
    bottom_frac = footer_in / fig_h

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(
        1,
        4,
        width_ratios=[4, 2, gc_panel_width, 2.5],
        wspace=0.15,
        bottom=bottom_frac,
        top=top_frac,
        left=0.05,
        right=0.95,
    )

    ax_tree = fig.add_subplot(gs[0])
    ax_lbl = fig.add_subplot(gs[1])
    ax_gc = fig.add_subplot(gs[2])
    ax_rt = fig.add_subplot(gs[3])

    plot_tree(ax_tree, h_lines, v_lines, n_tips, y_map, terminals)
    top_df = mark_fastest_branches(ax_tree, edge_data, args.top_branches)
    plot_labels(ax_lbl, terminals, y_map, fontsize=args.label_fontsize)

    if gc_panel_kind == 'partitions':
        norm_gc, cmap_gc, ticks_gc = plot_partition_gc_boxes(ax_gc, merged_df, selected_gc_cols, y_map, terminals, args.gc_cmap)
        gc_legend_label = 'Partition GC content'
    else:
        norm_gc, cmap_gc, ticks_gc = plot_gc_strip(ax_gc, merged_df, selected_gc_cols[0], y_map, terminals, args.gc_cmap, gc_panel_title)
        gc_legend_label = gc_panel_title

    plot_rates(ax_rt, merged_df, y_map, terminals, mean_r2t, args.rate_cmap)
    draw_footer(fig, norm_gc, cmap_gc, ticks_gc, fig_h, gc_legend_label)

    title = args.title.format(
        n=n_tips,
        tree=args.tree.name,
        n_partitions=len(partition_cols),
        gc_mode=args.gc_panel_mode,
    )
    fig.suptitle(title, y=1.0 - (0.5 / fig_h), fontsize=14, fontweight="bold", color="#2c3e50")

    outprefix = str(args.outprefix)
    for fmt in args.formats:
        out = f"{outprefix}.{fmt}"
        print(f"Saving figure -> {out}")
        if fmt == "png":
            fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
        else:
            fig.savefig(out, bbox_inches="tight")

    fastest_out = f"{outprefix}.fastest_branches.tsv"
    plotdata_out = f"{outprefix}.plot_data.tsv"
    gc_wide_out = f"{outprefix}.partition_gc.tsv"
    gc_long_out = f"{outprefix}.partition_gc.long.tsv"
    partinfo_out = f"{outprefix}.partitions.tsv"

    print(f"Saving fastest branch table -> {fastest_out}")
    top_df.to_csv(fastest_out, sep="\t", index=False)
    print(f"Saving plot data -> {plotdata_out}")
    merged_df.to_csv(plotdata_out, sep="\t", index=False)
    print(f"Saving wide partition GC table -> {gc_wide_out}")
    gc_wide_df.to_csv(gc_wide_out, sep="\t", index=False)
    print(f"Saving long partition GC table -> {gc_long_out}")
    gc_long_df.to_csv(gc_long_out, sep="\t", index=False)
    print(f"Saving partition definitions -> {partinfo_out}")
    partition_info_df.to_csv(partinfo_out, sep="\t", index=False)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
