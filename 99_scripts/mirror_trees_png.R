#conda install -n phylo_mirror -c conda-forge r-base r-ape r-phytools -y

## -------------------------
## SETTINGS
## -------------------------
tree1_path <- "/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/04_IQtree/01_runs/00_nucleotide/00_genafpair/mantids_nu_genaf_ML_MP.contree"
tree2_path <- "/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/04_IQtree/01_runs/00_nucleotide/00_genafpair/without_3rd/mantids_nu_genaf_ML_W3_MP.contree"

output_base <- "mirrored_phylo_comparison"

## custom titles
left_title  <- "cp123"
right_title <- "cp12"
main_title  <- "Mirrored comparison of phylogenetic trees"

## Plot style
plot_type <- "phylogram"      # "phylogram" or "cladogram"
rotate_to_match <- TRUE
italic_tip_labels <- FALSE
highlight_top_n <- 15

## Figure aesthetics
tip_cex <- 0.35
tree_lwd <- 1.2

base_edge_col <- "grey35"
base_link_col <- "grey20"
base_link_alpha <- 0.20
base_link_lwd <- 0.8

highlight_col <- "#B2182B"
highlight_link_lwd <- 2.4

## Output sizes
pdf_width  <- 16
pdf_height <- 22

png_width  <- 6000
png_height <- 8000
png_res    <- 300
## -------------------------


## -------------------------
## Package setup
## -------------------------
required_pkgs <- c("ape", "phytools")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install them with conda first:\n",
    "  conda install -c conda-forge r-base r-ape r-phytools",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

## -------------------------
## READ TREES
## -------------------------
tree1 <- ape::read.tree(tree1_path)
tree2 <- ape::read.tree(tree2_path)


## -------------------------
## ROOT TREES (unrooted .contree input)
## -------------------------

# Set your outgroup tip label exactly as it appears in the tree
# This should be the same taxon you use in iTOL/FigTree
outgroup_path <- "/home/STUDENTI/pietro.bacconi/Tirocinio_magistrale/00_data/03_outgroups/rerooting_outgroup.tsv"

# Read outgroup tips from file
if (!file.exists(outgroup_path)) {
  stop("Outgroup file not found: ", outgroup_path, call. = FALSE)
}

outgroup <- readLines(outgroup_path)
outgroup <- trimws(outgroup)               # removes accidental spaces/tabs
outgroup <- outgroup[nchar(outgroup) > 0]  # removes empty lines

message("Outgroup taxa loaded: ", length(outgroup), " tips")

## -------------------------
## DIAGNOSTIC - verify names and MRCA before rooting
## -------------------------
cat("\n--- Tips in your outgroup file ---\n")
print(outgroup)

cat("\n--- Matching tips found in tree1 ---\n")
print(tree1$tip.label[tree1$tip.label %in% outgroup])

cat("\n--- NOT matched in tree1 ---\n")
not_matched <- outgroup[!outgroup %in% tree1$tip.label]
if (length(not_matched) > 0) {
  print(not_matched)
  stop("Fix tip name mismatches above before rooting.", call. = FALSE)
} else {
  cat("All outgroup tips matched successfully.\n")
}

# Find MRCA node of the outgroup clade
mrca_node_t1 <- ape::getMRCA(tree1, tip = outgroup)
mrca_node_t2 <- ape::getMRCA(tree2, tip = outgroup)

cat("\n--- MRCA node in tree1:", mrca_node_t1, "---\n")
cat("--- MRCA node in tree2:", mrca_node_t2, "---\n")

# Verify descendants of MRCA - should be ONLY your outgroup tips
descendants_t1 <- ape::extract.clade(tree1, node = mrca_node_t1)$tip.label
descendants_t2 <- ape::extract.clade(tree2, node = mrca_node_t2)$tip.label

cat("\n--- Tips descending from MRCA in tree1 ---\n")
print(descendants_t1)

cat("\n--- Tips descending from MRCA in tree2 ---\n")
print(descendants_t2)

# Check for unexpected mantid tips inside the outgroup clade
unexpected_t1 <- descendants_t1[!descendants_t1 %in% outgroup]
unexpected_t2 <- descendants_t2[!descendants_t2 %in% outgroup]

if (length(unexpected_t1) > 0 | length(unexpected_t2) > 0) {
  cat("\nWARNING: Non-outgroup tips found inside outgroup clade!\n")
  cat("Tree1 unexpected tips:", unexpected_t1, "\n")
  cat("Tree2 unexpected tips:", unexpected_t2, "\n")
  cat("Outgroup may not be monophyletic - rooting may be unreliable.\n")
} else {
  cat("\nOutgroup clade is clean - rooting is safe to proceed.\n")
}

## ---------------------
## ROOT on outgroup MRCA
## ---------------------
tree1 <- ape::root(tree1, node = mrca_node_t1, resolve.root = TRUE)
tree2 <- ape::root(tree2, node = mrca_node_t2, resolve.root = TRUE)

message("Trees successfully rooted on ", length(outgroup), " outgroup taxa")


## -------------------------
## OPTIONAL: rotate right tree to better match left
## -------------------------
if (rotate_to_match) {
  tree2 <- ape::rotateConstr(tree2, tree1$tip.label)
}

## -------------------------
## IDENTIFY TIPS TO HIGHLIGHT
## -------------------------
tip_counts <- table(c(tree1$tip.label, tree2$tip.label))
top_tips   <- names(sort(tip_counts, decreasing = TRUE))[seq_len(highlight_top_n)]

## -------------------------
## PLOT FUNCTION (reused for PDF and PNG)
## -------------------------
draw_mirror <- function() {
  layout(matrix(c(1, 2), nrow = 1))
  par(oma = c(1, 0, 3, 0))

  # --- Left tree ---
  par(mar = c(1, 0, 2, 8), font = if (italic_tip_labels) 3 else 1)  # 3 = italic, 1 = regular
  phytools::plotTree(
    tree1,
    type       = plot_type,
    direction  = "rightwards",
    lwd        = tree_lwd,
    edge.color = base_edge_col,
    fsize      = tip_cex,
    add        = FALSE
  )
  if (!is.null(left_title)) title(main = left_title, cex.main = 0.9, line = 0.5)

  left_idx <- which(tree1$tip.label %in% top_tips)
  if (length(left_idx) > 0)
    tiplabels(tip = left_idx, pch = 21,
              bg = highlight_col, col = highlight_col, cex = 0.55)

  # --- Right tree (mirrored) ---
  par(mar = c(1, 8, 2, 0), font = if (italic_tip_labels) 3 else 1)
  phytools::plotTree(
    tree2,
    type       = plot_type,
    direction  = "leftwards",
    lwd        = tree_lwd,
    edge.color = base_edge_col,
    fsize      = tip_cex,
    add        = FALSE
  )
  if (!is.null(right_title)) title(main = right_title, cex.main = 0.9, line = 0.5)

  right_idx <- which(tree2$tip.label %in% top_tips)
  if (length(right_idx) > 0)
    tiplabels(tip = right_idx, pch = 21,
              bg = highlight_col, col = highlight_col, cex = 0.55)

  # --- Shared main title ---
  mtext(main_title, outer = TRUE, side = 3, cex = 1.1, font = 2, line = 1)
}


## -------------------------
## SAVE PDF
## -------------------------
pdf(paste0(output_base, ".pdf"), width = pdf_width, height = pdf_height)
draw_mirror()
dev.off()
message("Saved: ", output_base, ".pdf")

## -------------------------
## SAVE PNG
## -------------------------
png(paste0(output_base, ".png"),
    width = png_width, height = png_height, res = png_res)
draw_mirror()
dev.off()
message("Saved: ", output_base, ".png")


