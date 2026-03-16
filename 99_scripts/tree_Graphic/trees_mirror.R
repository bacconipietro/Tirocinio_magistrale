#!/usr/bin/env Rscript

## -------------------------
## INSTALLING PACKAGES
## conda install -c conda-forge r-base r-ape r-phytools r-optparse
## -------------------------

## -------------------------
## PACKAGE SETUP
## -------------------------
required_pkgs <- c("ape", "phytools", "optparse")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install them with conda first:\n",
    "  conda install -c conda-forge r-base r-ape r-phytools r-optparse",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(optparse)
})

## -------------------------
## COMMAND LINE ARGUMENTS
## -------------------------
option_list <- list(
  make_option(c("-f", "--tree1"), type = "character", default = NULL,
              help = "Path to the first tree file (.nwk)", metavar = "FILE"),
  make_option(c("-s", "--tree2"), type = "character", default = NULL,
              help = "Path to the second tree file (.nwk)", metavar = "FILE"),
  make_option(c("-o", "--out"), type = "character", default = "mirrored_comparison.png",
              help = "Output file name (e.g., name.png) [default = %default]", metavar = "FILE"),
  make_option(c("-L", "--left_title"), type = "character", default = "Tree 1",
              help = "Title for the left tree [default = %default]", metavar = "STRING"),
  make_option(c("-R", "--right_title"), type = "character", default = "Tree 2",
              help = "Title for the right tree [default = %default]", metavar = "STRING"),
  make_option(c("-M", "--main_title"), type = "character", default = "Mirrored comparison of phylogenetic trees",
              help = "Main overarching title [default = %default]", metavar = "STRING")
)

opt_parser <- OptionParser(option_list = option_list, 
                           description = "Mirrored comparison of two phylogenetic trees.")
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$tree1) || is.null(opt$tree2)) {
  print_help(opt_parser)
  stop("ERROR: You must provide both input trees using the -1 and -2 flags.", call. = FALSE)
}

tree1_path <- opt$tree1
tree2_path <- opt$tree2
output_filename <- opt$out

# Ensure output has a .png extension
if (!grepl("\\.png$", output_filename, ignore.case = TRUE)) {
  output_filename <- paste0(output_filename, ".png")
}

## -------------------------
## SETTINGS (Aesthetics)
## -------------------------
## Plot style
plot_type        <- "phylogram"   # "phylogram" or "cladogram"
italic_tip_labels <- FALSE
highlight_top_n  <- 15

## Figure aesthetics
tip_cex        <- 0.35
tree_lwd       <- 1.2

base_edge_col  <- "grey35"
highlight_col  <- "#B2182B"

## Output sizes
png_width  <- 6000
png_height <- 8000
png_res    <- 300
## -------------------------

## -------------------------
## READ TREES
## -------------------------
tree1 <- ape::read.tree(tree1_path)
tree2 <- ape::read.tree(tree2_path)

message("Trees loaded successfully.")

## -------------------------
## OPTIMIZE ROTATION (MIRRORING)
## -------------------------
message("Calculating optimal node rotations to mirror the trees...")
comp_trees <- phytools::cophylo(tr1 = tree1, tr2 = tree2)

tree1_rotated <- comp_trees$trees[[1]]
tree2_rotated <- comp_trees$trees[[2]]
message("Trees successfully rotated to match topologies.")

## -------------------------
## IDENTIFY TIPS TO HIGHLIGHT
## -------------------------
tip_counts <- table(c(tree1_rotated$tip.label, tree2_rotated$tip.label))
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
    tree1_rotated,
    type       = plot_type,
    direction  = "rightwards",
    lwd        = tree_lwd,
    edge.color = base_edge_col,
    fsize      = tip_cex,
    add        = FALSE
  )
  if (!is.null(opt$left_title)) title(main = opt$left_title, cex.main = 0.9, line = 0.5)

  left_idx <- which(tree1_rotated$tip.label %in% top_tips)
  if (length(left_idx) > 0)
    tiplabels(tip = left_idx, pch = 21,
              bg = highlight_col, col = highlight_col, cex = 0.55)

  # --- Right tree (mirrored) ---
  par(mar = c(1, 8, 2, 0), font = if (italic_tip_labels) 3 else 1)
  phytools::plotTree(
    tree2_rotated,
    type       = plot_type,
    direction  = "leftwards",
    lwd        = tree_lwd,
    edge.color = base_edge_col,
    fsize      = tip_cex,
    add        = FALSE
  )
  if (!is.null(opt$right_title)) title(main = opt$right_title, cex.main = 0.9, line = 0.5)

  right_idx <- which(tree2_rotated$tip.label %in% top_tips)
  if (length(right_idx) > 0)
    tiplabels(tip = right_idx, pch = 21,
              bg = highlight_col, col = highlight_col, cex = 0.55)

  # --- Shared main title ---
  mtext(opt$main_title, outer = TRUE, side = 3, cex = 1.1, font = 2, line = 1)
}

## -------------------------
## SAVE PNG
## -------------------------
png(output_filename, width = png_width, height = png_height, res = png_res)
draw_mirror()
invisible(dev.off())
message("Saved: ", output_filename)
