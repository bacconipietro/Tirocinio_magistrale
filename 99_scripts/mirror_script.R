## -------------------------
## INSTALLING PACKAGES
## conda install r-base r-ape r-phytools 
## -------------------------
## -------------------------
## SETTINGS
## -------------------------
tree1_path <- "/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/04_IQtree/01_runs/00_nucleotide/00_genafpair/mantids_nu_genaf_ML_MP.contree"
tree2_path <- "/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/04_IQtree/01_runs/00_nucleotide/00_genafpair/without_3rd/mantids_nu_genaf_ML_W3_MP.contree"

output_base <- "mirrored_phylo_comparison"

## Custom titles
left_title  <- "cp123"
right_title <- "cp12"
main_title  <- "Mirrored comparison of phylogenetic trees"

## Plot style
plot_type        <- "phylogram"   # "phylogram" or "cladogram"
rotate_to_match  <- TRUE
italic_tip_labels <- FALSE
highlight_top_n  <- 15

## Figure aesthetics
tip_cex        <- 0.35
tree_lwd       <- 1.2

base_edge_col  <- "grey35"
base_link_col  <- "grey20"
base_link_alpha <- 0.20
base_link_lwd  <- 0.8

highlight_col     <- "#B2182B"
highlight_link_lwd <- 2.4

## Output sizes
pdf_width  <- 16
pdf_height <- 22

png_width  <- 6000
png_height <- 8000
png_res    <- 300
## -------------------------


## -------------------------
## PACKAGE SETUP
## -------------------------
required_pkgs <- c("ape", "phytools")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

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

message("Trees loaded successfully.")
message("Tree1 tips: ", length(tree1$tip.label))
message("Tree2 tips: ", length(tree2$tip.label))


## -------------------------
## OPTIONAL: rotate right tree to better match left
## -------------------------
if (rotate_to_match) {
  tree2 <- ape::rotateConstr(tree2, tree1$tip.label)
  message("Tree2 rotated to match Tree1 tip order.")
}


## -------------------------
## IDENTIFY TIPS TO HIGHLIGHT
## -------------------------
tip_counts <- table(c(tree1$tip.label, tree2$tip.label))
top_tips   <- names(sort(tip_counts, decreasing = TRUE))[seq_len(highlight_top_n)]


## -------------------------
## PLOT FUNCTION
## -------------------------
draw_mirror <- function() {
  layout(matrix(c(1, 2), nrow = 1))
  par(oma = c(1, 0, 3, 0))

  # --- Left tree ---
  par(mar = c(1, 0, 2, 8), font = if (italic_tip_labels) 3 else 1)
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
## SAVE PNG
## -------------------------
png(paste0(output_base, ".png"),
    width = png_width, height = png_height, res = png_res)
draw_mirror()
dev.off()
message("Saved: ", output_base, ".png")
