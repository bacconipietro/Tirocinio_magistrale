#!/usr/bin/env python3
"""
aa_freq_from_alignment.py
--------------------------
Three-pipeline tool for amino acid frequency analysis of AA alignments,
designed for invertebrate mitochondrial data.

─────────────────────────────────────────────────────────────────────────────
PIPELINE 1 — freq
  Read a FASTA AA alignment, compute per-species AA relative frequencies,
  and write a TSV table.

  Usage:
    python aa_freq_from_alignment.py freq -i alignment.fasta -o aa_freq.tsv

─────────────────────────────────────────────────────────────────────────────
PIPELINE 2 — heatmap  (no tree)
  Read the TSV from pipeline 1 and draw a standalone violet heatmap:
    Y axis → species  |  X axis → amino acids
    Color  → light violet (low) → dark violet (high)

  Usage:
    python aa_freq_from_alignment.py heatmap -i aa_freq.tsv -o heatmap.png

─────────────────────────────────────────────────────────────────────────────
PIPELINE 3 — treemap  (tree + labels + heatmap)
  Combine a Newick tree with the AA frequency TSV into one figure laid out as:

        [ tree topology ] | [ species labels ] | [ AA freq heatmap ]

  Tip order is driven by the tree (top to bottom).
  Tips in the TSV but absent from the tree are silently skipped.
  Tips in the tree but absent from the TSV are shown as blank rows.

  Usage:
    python aa_freq_from_alignment.py treemap \\
        -i  aa_freq.tsv \\
        -t  my_tree.nwk \\
        -o  treemap.png

─────────────────────────────────────────────────────────────────────────────
Invertebrate mitochondrial notes
  - W (Trp) is kept: UGA encodes Trp in invertebrate mito code (NCBI table 5).
  - Stop codons (*) are excluded from counts.
  - Ambiguous residues (X, B, Z, J) are excluded from the denominator by
    default; use --include-ambiguous to change this.
  - Gap characters (- and .) are always excluded.

Dependencies
  pip install biopython matplotlib pandas scipy
"""

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

# ─── Constants ────────────────────────────────────────────────────────────────

STANDARD_AA     = list("ACDEFGHIKLMNPQRSTVWY")
GAP_CHARS       = set("-.")
STOP_CHARS      = set("*")
AMBIGUOUS_CHARS = set("XBZJ")
EXCLUDE_AMBIGUOUS = True       # module-level default; overridden by CLI flag

# Dark-theme palette (shared by both heatmap and treemap)
BG_COLOR      = "#0D0D12"      # near-black canvas
TREE_COLOR    = "#8892A4"      # muted blue-grey for tree branches
LABEL_COLOR   = "#EDE0F5"      # light lavender for labels
ACCENT_COLOR  = "#C9A0DC"      # medium violet for axis titles / colorbar
STRIPE_COLOR  = "#1A1A2A"      # alternating row stripe (dark)
GRID_COLOR    = "#0D0D12"      # cell grid lines


# ─────────────────────────────────────────────────────────────────────────────
#  PIPELINE 1 — frequency computation
# ─────────────────────────────────────────────────────────────────────────────

def compute_aa_frequencies(sequence: str) -> dict:
    """
    Compute relative AA frequencies for one sequence.

    Returns a dict with:
      freq_<AA>        — relative frequency for each of the 20 standard AAs
      count_valid_AA   — denominator (residues used for freq calculation)
      count_gaps       — gap characters skipped
      count_stops      — stop codon symbols skipped
      count_ambiguous  — ambiguous residue symbols (X/B/Z/J) found
      alignment_length — total character length including gaps
    """
    seq_upper = sequence.upper()
    counts = Counter()
    n_gaps = n_stops = n_ambig = 0

    for char in seq_upper:
        if char in GAP_CHARS:
            n_gaps += 1
        elif char in STOP_CHARS:
            n_stops += 1
        elif char in AMBIGUOUS_CHARS:
            n_ambig += 1
            counts[char] += 1
        else:
            counts[char] += 1

    if EXCLUDE_AMBIGUOUS:
        denominator = sum(counts[aa] for aa in STANDARD_AA)
    else:
        denominator = sum(counts[aa] for aa in STANDARD_AA) + n_ambig

    result = {}
    for aa in STANDARD_AA:
        result[f"freq_{aa}"] = (counts[aa] / denominator) if denominator > 0 else 0.0

    result["count_valid_AA"]   = denominator
    result["count_gaps"]       = n_gaps
    result["count_stops"]      = n_stops
    result["count_ambiguous"]  = n_ambig
    result["alignment_length"] = len(seq_upper)
    return result


def parse_alignment(fasta_path: Path) -> list:
    """Parse a FASTA alignment; return one dict per sequence."""
    from Bio import SeqIO
    records = []
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        freqs = compute_aa_frequencies(str(record.seq))
        row = {"species": record.id}
        row.update(freqs)
        records.append(row)
    if not records:
        sys.exit(f"[ERROR] No sequences found in {fasta_path}.")
    return records


def write_tsv(records: list, out_path: Path):
    """Write per-species AA frequencies to a TSV file."""
    fieldnames = (
        ["species"]
        + [f"freq_{aa}" for aa in STANDARD_AA]
        + ["count_valid_AA", "count_gaps", "count_stops",
           "count_ambiguous", "alignment_length"]
    )
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    print(f"[OK] Written {len(records)} species → {out_path}")


def run_freq(args):
    """Entry point for the 'freq' subcommand."""
    global EXCLUDE_AMBIGUOUS
    if args.include_ambiguous:
        EXCLUDE_AMBIGUOUS = False
        print("[INFO] Ambiguous residues included in denominator.")

    fasta_path = Path(args.input)
    if not fasta_path.exists():
        sys.exit(f"[ERROR] File not found: {fasta_path}")

    print(f"[INFO] Parsing alignment: {fasta_path}")
    records = parse_alignment(fasta_path)
    print(f"[INFO] Found {len(records)} sequences.")
    write_tsv(records, Path(args.output))


# ─────────────────────────────────────────────────────────────────────────────
#  Shared heatmap helpers (used by both pipeline 2 and 3)
# ─────────────────────────────────────────────────────────────────────────────

def build_violet_cmap():
    """
    Light lavender → medium violet → deep indigo.
    Perceptually smooth light-to-dark violet gradient.
    """
    from matplotlib.colors import LinearSegmentedColormap
    stops = ["#EDE0F5", "#C9A0DC", "#7B3FA0", "#4B0082", "#1A003D"]
    return LinearSegmentedColormap.from_list("violet_freq", stops)


def load_freq_table(tsv_path: Path):
    """
    Read the TSV produced by pipeline 1.

    Returns
    -------
    species_labels : list[str]
    matrix         : numpy ndarray, shape (n_species, 20)
    aa_labels      : list[str] — the 20 single-letter AA codes
    freq_dict      : dict[species_id -> ndarray of shape (20,)]
    """
    import numpy as np
    import pandas as pd

    df = pd.read_csv(tsv_path, sep="\t")
    freq_cols = [f"freq_{aa}" for aa in STANDARD_AA]
    missing = [c for c in freq_cols if c not in df.columns]
    if missing:
        sys.exit(
            f"[ERROR] TSV is missing expected columns: {missing}\n"
            "Make sure you are passing the output of the 'freq' pipeline."
        )
    species_labels = df["species"].tolist()
    matrix         = df[freq_cols].values.astype(float)
    freq_dict      = {sp: matrix[i] for i, sp in enumerate(species_labels)}
    return species_labels, matrix, STANDARD_AA, freq_dict


def render_heatmap_axes(ax, matrix, species_labels, aa_labels,
                        show_values: bool, show_y_labels: bool,
                        y_fontsize: int):
    """
    Draw the AA-frequency heatmap onto *ax*.
    Returns (im, cbar-ready ScalarMappable).
    """
    import numpy as np
    import matplotlib.pyplot as plt

    cmap = build_violet_cmap()
    vmax = matrix.max() if matrix.size > 0 else 1.0
    im   = ax.imshow(
        matrix,
        aspect="auto",
        cmap=cmap,
        vmin=0.0,
        vmax=vmax,
        interpolation="nearest",
    )

    # Optional cell-value annotations
    if show_values:
        thresh = vmax * 0.60
        n_r, n_c = matrix.shape
        for r in range(n_r):
            for c in range(n_c):
                val   = matrix[r, c]
                color = LABEL_COLOR if val < thresh else BG_COLOR
                ax.text(c, r, f"{val:.3f}",
                        ha="center", va="center",
                        fontsize=5, color=color)

    # X axis — amino acids (top)
    n_aa = len(aa_labels)
    ax.set_xticks(range(n_aa))
    ax.set_xticklabels(aa_labels, fontsize=9, fontweight="bold",
                       color=LABEL_COLOR, fontfamily="monospace")
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.set_xlabel("Amino Acid", fontsize=10, color=ACCENT_COLOR,
                  labelpad=8, fontweight="bold")

    # Y axis — species (only when no separate label panel)
    n_sp = len(species_labels)
    ax.set_yticks(range(n_sp))
    if show_y_labels:
        ax.set_yticklabels(species_labels, fontsize=y_fontsize,
                           color=LABEL_COLOR, fontfamily="monospace")
        ax.set_ylabel("Species", fontsize=10, color=ACCENT_COLOR,
                      labelpad=8, fontweight="bold")
    else:
        ax.set_yticklabels([])
        ax.set_ylabel("")

    # Minor-grid between cells
    ax.set_xticks(np.arange(-0.5, n_aa, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_sp, 1), minor=True)
    ax.grid(which="minor", color=GRID_COLOR, linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False,
                   top=False, right=False)

    # Spines
    for spine in ax.spines.values():
        spine.set_edgecolor("#4B0082")
        spine.set_linewidth(1.0)

    return im


# ─────────────────────────────────────────────────────────────────────────────
#  PIPELINE 2 — standalone heatmap (no tree)
# ─────────────────────────────────────────────────────────────────────────────

def draw_heatmap(tsv_path: Path, out_path: Path,
                 fig_width: float, fig_height,
                 show_values: bool, cluster_species: bool, dpi: int):
    """Build and save the standalone violet heatmap."""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist

    species_labels, matrix, aa_labels, _ = load_freq_table(tsv_path)
    n_species = len(species_labels)

    # Optional hierarchical clustering
    if cluster_species and n_species > 1:
        dist   = pdist(matrix, metric="euclidean")
        Z      = linkage(dist, method="average")
        order  = leaves_list(Z)
        matrix = matrix[order]
        species_labels = [species_labels[i] for i in order]
        print("[INFO] Species rows reordered by hierarchical clustering.")

    auto_h  = max(6.0, n_species * 0.28 + 3.0)
    height  = fig_height if fig_height else auto_h
    y_fs    = max(4, min(9, 200 // n_species))

    fig, ax = plt.subplots(figsize=(fig_width, height))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor(BG_COLOR)

    im = render_heatmap_axes(ax, matrix, species_labels, aa_labels,
                             show_values=show_values,
                             show_y_labels=True,
                             y_fontsize=y_fs)

    ax.set_title("Amino Acid Frequency — Invertebrate Mitochondrial Alignment",
                 fontsize=12, color=LABEL_COLOR, pad=18, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, aspect=40)
    cbar.set_label("Relative AA frequency", fontsize=9,
                   color=ACCENT_COLOR, labelpad=8)
    cbar.ax.yaxis.set_tick_params(color=ACCENT_COLOR, labelcolor=ACCENT_COLOR,
                                  labelsize=7)
    cbar.outline.set_edgecolor("#4B0082")

    plt.tight_layout()
    fig.savefig(out_path, dpi=dpi,
                bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"[OK] Heatmap saved → {out_path}")
    print(f"     {n_species} species × {len(aa_labels)} amino acids")


def run_heatmap(args):
    tsv_path = Path(args.input)
    if not tsv_path.exists():
        sys.exit(f"[ERROR] File not found: {tsv_path}")
    draw_heatmap(
        tsv_path        = tsv_path,
        out_path        = Path(args.output),
        fig_width       = args.width,
        fig_height      = args.height,
        show_values     = args.show_values,
        cluster_species = args.cluster,
        dpi             = args.dpi,
    )


# ─────────────────────────────────────────────────────────────────────────────
#  PIPELINE 3 — tree + label strip + heatmap
# ─────────────────────────────────────────────────────────────────────────────

# ── Tree geometry ─────────────────────────────────────────────────────────────

def read_newick(path: Path):
    """Read a Newick tree with Bio.Phylo."""
    from Bio import Phylo
    from io import StringIO
    try:
        return Phylo.read(str(path), "newick")
    except Exception:
        text = path.read_text()
        if ";" not in text:
            raise
        return Phylo.read(StringIO(text[: text.rfind(";") + 1]), "newick")


def get_tree_layout(tree):
    """
    Compute (x, y) layout coordinates for every clade in the tree.

    Y coordinates: tips are assigned integer positions 0..n-1 in
    in-order (top-to-bottom) traversal; internal nodes get the mean
    of their children's Y.

    X coordinates: root at 0, each node at parent_x + branch_length.

    Returns
    -------
    x_map   : dict[id(clade) -> float]
    y_map   : dict[id(clade) -> float]
    terminals : list[Clade]  — tips in top-to-bottom display order
    """
    import numpy as np

    terminals = tree.get_terminals()           # in-order traversal order
    y_map = {id(t): float(i) for i, t in enumerate(terminals)}

    def _calc_y(clade):
        if id(clade) in y_map:
            return
        for child in clade.clades:
            _calc_y(child)
        y_map[id(clade)] = float(np.mean([y_map[id(c)] for c in clade.clades]))

    _calc_y(tree.root)

    x_map = {}

    def _calc_x(clade, curr_x):
        x_map[id(clade)] = curr_x
        for child in clade.clades:
            bl = child.branch_length if child.branch_length else 0.0
            _calc_x(child, curr_x + bl)

    _calc_x(tree.root, 0.0)
    return x_map, y_map, terminals


def get_tree_segments(tree, x_map, y_map):
    """
    Return two lists of line segments for drawing the tree:
      h_segs — horizontal branches (parent x → child x, at child y)
      v_segs — vertical connectors (min child y → max child y, at parent x)
    """
    h_segs, v_segs = [], []

    def _walk(clade):
        px = x_map[id(clade)]
        if not clade.is_terminal():
            ys = [y_map[id(c)] for c in clade.clades]
            v_segs.append([(px, min(ys)), (px, max(ys))])
        for child in clade.clades:
            cx = x_map[id(child)]
            cy = y_map[id(child)]
            h_segs.append([(px, cy), (cx, cy)])
            _walk(child)

    _walk(tree.root)
    return h_segs, v_segs


# ── Zebra-stripe background helper ───────────────────────────────────────────

def draw_zebra(ax, n_tips):
    """Paint alternating faint horizontal bands on *ax* for readability."""
    ax.set_ylim(-0.5, n_tips - 0.5)
    ax.invert_yaxis()
    for i in range(n_tips):
        if i % 2 == 1:
            ax.axhspan(i - 0.5, i + 0.5,
                       color=STRIPE_COLOR, alpha=1.0, zorder=0)


# ── Tree panel ────────────────────────────────────────────────────────────────

def draw_tree_panel(ax, h_segs, v_segs, n_tips, x_map, tree):
    """
    Draw the phylogenetic tree topology onto *ax*.
    The X axis shows branch-length scale; Y is species order (inverted top→bottom).
    """
    from matplotlib.collections import LineCollection
    import numpy as np

    draw_zebra(ax, n_tips)

    ax.add_collection(LineCollection(v_segs, colors=TREE_COLOR, lw=1.0, zorder=2))
    ax.add_collection(LineCollection(h_segs, colors=TREE_COLOR, lw=1.2, zorder=2))

    # Auto-scale X to the data
    all_x = list(x_map.values())
    x_min, x_max = min(all_x), max(all_x)
    x_span = x_max - x_min if x_max > x_min else 1.0
    ax.set_xlim(x_min - x_span * 0.02, x_max + x_span * 0.05)

    # Minimal axis decoration — only the bottom scale bar
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=7, length=3,
                   color=TREE_COLOR, labelcolor=LABEL_COLOR)
    ax.set_xlabel("Substitutions / site", fontsize=8,
                  color=ACCENT_COLOR, labelpad=6)
    for spine in ("top", "left", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_color(TREE_COLOR)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.set_facecolor(BG_COLOR)


# ── Label strip panel ────────────────────────────────────────────────────────

def draw_label_panel(ax, tip_order, n_tips, fontsize):
    """
    Draw species names as a plain text strip.
    *tip_order* is the list of tip names in top-to-bottom display order.
    """
    draw_zebra(ax, n_tips)
    ax.set_xlim(0, 1)
    ax.axis("off")
    ax.set_facecolor(BG_COLOR)

    for i, name in enumerate(tip_order):
        # Clean display name: strip accession prefix if present
        display = _clean_display_name(name)
        ax.text(0.02, i, display,
                va="center", ha="left",
                fontsize=fontsize,
                color=LABEL_COLOR,
                fontfamily="monospace")


def _clean_display_name(label: str) -> str:
    """
    Human-readable tip label.
    If the name looks like GCA_XXXX.1_Genus_species_..., return 'Genus species'.
    Otherwise return the raw label with underscores replaced by spaces.
    """
    parts = label.replace("'", "").replace('"', "").split("_")
    if len(parts) >= 3 and parts[0] in ("GCA", "GCF"):
        return " ".join(parts[2:])
    return label.replace("_", " ")


# ── Main treemap builder ──────────────────────────────────────────────────────

def draw_treemap(tsv_path: Path, tree_path: Path, out_path: Path,
                 fig_width, fig_height,
                 show_values: bool,
                 label_fontsize: float,
                 dpi: int):
    """
    Compose the three-panel figure:
        [ tree ] | [ species labels ] | [ AA heatmap ]

    Tip order is dictated by the tree's in-order traversal (top → bottom).
    """
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    # ── Load data ──────────────────────────────────────────────────────────
    _, _, aa_labels, freq_dict = load_freq_table(tsv_path)

    tree = read_newick(tree_path)
    x_map, y_map, terminals = get_tree_layout(tree)
    h_segs, v_segs          = get_tree_segments(tree, x_map, y_map)

    # Tip names in display order (top → bottom as laid out by tree)
    tip_order = [t.name.strip().replace("'", "").replace('"', "")
                 for t in terminals]
    n_tips    = len(tip_order)

    # ── Build matrix aligned to tree tip order ─────────────────────────────
    # Tips present in both tree and TSV → fill matrix row
    # Tips in tree but missing in TSV  → row of NaN (rendered as "no data")
    nan_row = np.full(len(aa_labels), np.nan)
    matrix  = np.vstack([
        freq_dict.get(name, nan_row) for name in tip_order
    ])

    n_matched = sum(1 for name in tip_order if name in freq_dict)
    n_missing = n_tips - n_matched
    if n_missing:
        print(f"[WARN] {n_missing} tree tips not found in TSV "
              "(shown as blank rows).")
    print(f"[INFO] Matched {n_matched}/{n_tips} tree tips to TSV rows.")

    # ── Figure geometry ────────────────────────────────────────────────────
    row_h      = 0.28                          # inches per tip row
    auto_h     = max(8.0, n_tips * row_h + 3.0)
    height     = fig_height if fig_height else auto_h
    width      = fig_width  if fig_width  else 20.0

    label_fs   = label_fontsize if label_fontsize else max(4, min(9, 200 // n_tips))

    # Column width ratios: tree | labels | heatmap
    # Tree panel: ~30 % of width; labels: ~20 %; heatmap: remaining ~50 %
    fig = plt.figure(figsize=(width, height))
    fig.patch.set_facecolor(BG_COLOR)

    gs = gridspec.GridSpec(
        1, 3,
        width_ratios=[3, 2, 5],
        wspace=0.04,
        left=0.04, right=0.96,
        bottom=0.06, top=0.94,
        figure=fig,
    )

    ax_tree  = fig.add_subplot(gs[0])
    ax_label = fig.add_subplot(gs[1])
    ax_heat  = fig.add_subplot(gs[2])

    # ── Panel 1 — tree ────────────────────────────────────────────────────
    draw_tree_panel(ax_tree, h_segs, v_segs, n_tips, x_map, tree)
    ax_tree.set_ylim(-0.5, n_tips - 0.5)
    ax_tree.invert_yaxis()

    # ── Panel 2 — label strip ─────────────────────────────────────────────
    draw_label_panel(ax_label, tip_order, n_tips, fontsize=label_fs)
    ax_label.set_ylim(-0.5, n_tips - 0.5)
    ax_label.invert_yaxis()

    # ── Panel 3 — heatmap ─────────────────────────────────────────────────
    ax_heat.set_facecolor(BG_COLOR)

    # Replace NaN rows with zeros for display (they'll look like "no data"
    # because 0 maps to the lightest violet)
    display_matrix = np.where(np.isnan(matrix), 0.0, matrix)

    im = render_heatmap_axes(
        ax_heat,
        display_matrix,
        tip_order,
        aa_labels,
        show_values  = show_values,
        show_y_labels= False,          # labels drawn in the dedicated panel
        y_fontsize   = label_fs,
    )

    # Synchronise Y limits across all panels (imshow sets its own)
    for ax in (ax_tree, ax_label):
        ax.set_ylim(ax_heat.get_ylim())

    # ── Colorbar ──────────────────────────────────────────────────────────
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.025, pad=0.02, aspect=40)
    cbar.set_label("Relative AA frequency", fontsize=9,
                   color=ACCENT_COLOR, labelpad=8)
    cbar.ax.yaxis.set_tick_params(color=ACCENT_COLOR, labelcolor=ACCENT_COLOR,
                                  labelsize=7)
    cbar.outline.set_edgecolor("#4B0082")

    # ── Title ─────────────────────────────────────────────────────────────
    fig.suptitle(
        f"Amino Acid Frequency — Invertebrate Mitochondrial Alignment  "
        f"(n = {n_tips})",
        fontsize=13, color=LABEL_COLOR, fontweight="bold", y=0.975,
    )

    # ── Save ──────────────────────────────────────────────────────────────
    fig.savefig(out_path, dpi=dpi,
                bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"[OK] Treemap saved → {out_path}")
    print(f"     Layout: tree | labels | {n_tips} species × {len(aa_labels)} AAs")


def run_treemap(args):
    """Entry point for the 'treemap' subcommand."""
    tsv_path  = Path(args.input)
    tree_path = Path(args.tree)

    for p, label in [(tsv_path, "TSV"), (tree_path, "tree")]:
        if not p.exists():
            sys.exit(f"[ERROR] {label} file not found: {p}")

    draw_treemap(
        tsv_path      = tsv_path,
        tree_path     = tree_path,
        out_path      = Path(args.output),
        fig_width     = args.width,
        fig_height    = args.height,
        show_values   = args.show_values,
        label_fontsize= args.label_fontsize,
        dpi           = args.dpi,
    )


# ─────────────────────────────────────────────────────────────────────────────
#  CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        prog="aa_freq_from_alignment.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ── subcommand: freq ──────────────────────────────────────────────────
    p_freq = sub.add_parser(
        "freq",
        help="Pipeline 1 — compute AA frequencies from a FASTA alignment",
    )
    p_freq.add_argument("-i", "--input",  required=True,
                        help="Input FASTA file (AA alignment)")
    p_freq.add_argument("-o", "--output", default="aa_frequencies.tsv",
                        help="Output TSV  [default: aa_frequencies.tsv]")
    p_freq.add_argument("--include-ambiguous", action="store_true",
                        help="Include ambiguous residues (X/B/Z/J) in denominator")

    # ── subcommand: heatmap ───────────────────────────────────────────────
    p_heat = sub.add_parser(
        "heatmap",
        help="Pipeline 2 — draw a violet heatmap from a frequency TSV",
    )
    p_heat.add_argument("-i", "--input",  required=True,
                        help="Input TSV (output of 'freq' pipeline)")
    p_heat.add_argument("-o", "--output", default="aa_heatmap.png",
                        help="Output image (.png/.pdf/.svg)  "
                             "[default: aa_heatmap.png]")
    p_heat.add_argument("--width",  type=float, default=12.0,
                        help="Figure width in inches  [default: 12]")
    p_heat.add_argument("--height", type=float, default=None,
                        help="Figure height in inches  [default: auto]")
    p_heat.add_argument("--dpi", type=int, default=300,
                        help="Resolution for PNG output  [default: 300]")
    p_heat.add_argument("--show-values", action="store_true",
                        help="Annotate each cell with its numeric value")
    p_heat.add_argument("--cluster", action="store_true",
                        help="Reorder species rows by hierarchical clustering")

    # ── subcommand: treemap ───────────────────────────────────────────────
    p_tree = sub.add_parser(
        "treemap",
        help="Pipeline 3 — tree topology + label strip + AA heatmap",
    )
    p_tree.add_argument("-i", "--input",  required=True,
                        help="Input TSV (output of 'freq' pipeline)")
    p_tree.add_argument("-t", "--tree",   required=True,
                        help="Input Newick tree file")
    p_tree.add_argument("-o", "--output", default="aa_treemap.png",
                        help="Output image (.png/.pdf/.svg)  "
                             "[default: aa_treemap.png]")
    p_tree.add_argument("--width",  type=float, default=None,
                        help="Figure width in inches  [default: 20]")
    p_tree.add_argument("--height", type=float, default=None,
                        help="Figure height in inches  [default: auto]")
    p_tree.add_argument("--dpi", type=int, default=300,
                        help="Resolution for PNG output  [default: 300]")
    p_tree.add_argument("--show-values", action="store_true",
                        help="Annotate each heatmap cell with its numeric value "
                             "(recommended only for small datasets)")
    p_tree.add_argument("--label-fontsize", type=float, default=None,
                        help="Species label font size  [default: auto]")

    args = parser.parse_args()

    if args.command == "freq":
        run_freq(args)
    elif args.command == "heatmap":
        run_heatmap(args)
    elif args.command == "treemap":
        run_treemap(args)


if __name__ == "__main__":
    main()
