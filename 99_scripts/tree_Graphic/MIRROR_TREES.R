#!/usr/bin/env Rscript

## -------------------------
## INSTALLING PACKAGES
## conda install -c conda-forge r-base r-ape r-phytools r-optparse
## -------------------------

required_pkgs <- c("ape", "phytools", "optparse")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("Missing R packages. Install with conda: conda install -c conda-forge r-base r-ape r-phytools r-optparse", call. = FALSE)
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
  make_option(c("-f", "--tree1"), type = "character", default = NULL),
  make_option(c("-s", "--tree2"), type = "character", default = NULL),
  make_option(c("-o", "--out"), type = "character", default = "Colored_Families.png"),
  make_option(c("-L", "--left_title"), type = "character", default = "Tree 1"),
  make_option(c("-R", "--right_title"), type = "character", default = "Tree 2"),
  make_option(c("-M", "--main_title"), type = "character", default = "Mirrored comparison of phylogenetic trees")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$tree1) || is.null(opt$tree2)) {
  stop("ERROR: You must provide both input trees using the -f and -s flags.", call. = FALSE)
}

output_filename <- opt$out
if (!grepl("\\.png$", output_filename, ignore.case = TRUE)) {
  output_filename <- paste0(output_filename, ".png")
}

## -------------------------
## SETTINGS & HARDCODED PALETTE
## -------------------------
plot_type        <- "phylogram"   
italic_tip_labels <- FALSE
tip_cex          <- 0.35
tree_lwd         <- 1.2
base_edge_col    <- "grey35"

png_width  <- 6000
png_height <- 8000
png_res    <- 300

family_palette <- c(
  "Acanthopidae"     = "#B22222", "Amelidae"         = "#0078D7",
  "Amorphoscelidae"  = "#556B2F", "Coptopterygidae"  = "#CD5C5C",
  "Deroplatyidae"    = "#0000CD", "Empusidae"        = "#8A2BE2",
  "Eremiaphilidae"   = "#008B8B", "Galinthiadidae"   = "#00008B",
  "Gonypetidae"      = "#808000", "Haaniidae"        = "#D2691E",
  "Hymenopodidae"    = "#C71585", "Leptomantellidae" = "#2E8B57",
  "Mantidae"         = "#4B0082", "Metallyticidae"   = "#8B0000",
  "Miomantidae"      = "#4169E1", "Nanomantidae"     = "#228B22",
  "Rivetinidae"      = "#008080", "Thespidae"        = "#B8860B",
  "Toxoderidae"      = "#4682B4", "Termitoidae"      = "#708090",
  "Blaberidae"       = "#708090", "Ectobiidae"       = "#708090",
  "Blattidae"        = "#708090", "Corydiidae"       = "#708090",
  "Cryptocercidae"   = "#708090"
)
# Sort to prevent "Mantidae" from overwriting "Nanomantidae"
family_palette <- family_palette[order(-nchar(names(family_palette)))]

## -------------------------
## READ TREES & MIRROR
## -------------------------
tree1 <- ape::read.tree(opt$tree1)
tree2 <- ape::read.tree(opt$tree2)
message("\nCalculating optimal node rotations...")
comp_trees <- phytools::cophylo(tr1 = tree1, tr2 = tree2)

tree1_rotated <- comp_trees$trees[[1]]
tree2_rotated <- comp_trees$trees[[2]]

## -------------------------
## COLOR MATCHING LOGIC
## -------------------------
search_tip_colors <- function(tip_labels) {
  tip_cols <- rep("black", length(tip_labels))
  for (i in seq_along(tip_labels)) {
    for (fam_name in names(family_palette)) {
      if (grepl(fam_name, tip_labels[i], ignore.case = TRUE)) {
        tip_cols[i] <- family_palette[[fam_name]]
        break
      }
    }
  }
  return(tip_cols)
}

left_colors  <- search_tip_colors(tree1_rotated$tip.label)
right_colors <- search_tip_colors(tree2_rotated$tip.label)

## -------------------------
## TERMINAL DIAGNOSTICS
## -------------------------
message("\n--- COLOR MATCHING DIAGNOSTICS ---")
message("Successfully matched ", sum(left_colors != "black"), " out of ", length(left_colors), " tips in Left Tree.")
message("Successfully matched ", sum(right_colors != "black"), " out of ", length(right_colors), " tips in Right Tree.")
message("----------------------------------\n")

## -------------------------
## PLOT FUNCTION (USING BASE APE)
## -------------------------
draw_mirror <- function() {
  layout(matrix(c(1, 2), nrow = 1))
  par(oma = c(1, 0, 4, 0))

  # --- LEFT TREE ---
  par(mar = c(1, 0, 2, 8), font = if (italic_tip_labels) 3 else 1)
  # Swapped plotTree for plot.phylo to bypass any phytools bugs
  ape::plot.phylo(
    tree1_rotated, type = plot_type, direction = "rightwards",
    edge.width = tree_lwd, edge.color = base_edge_col, cex = tip_cex, 
    tip.color = left_colors, no.margin = FALSE
  )
  mtext(opt$left_title, side = 3, line = 0, cex = 1.1, font = 2) 

  # --- RIGHT TREE ---
  par(mar = c(1, 8, 2, 0), font = if (italic_tip_labels) 3 else 1)
  ape::plot.phylo(
    tree2_rotated, type = plot_type, direction = "leftwards",
    edge.width = tree_lwd, edge.color = base_edge_col, cex = tip_cex, 
    tip.color = right_colors, no.margin = FALSE
  )
  mtext(opt$right_title, side = 3, line = 0, cex = 1.1, font = 2)

  # --- SHARED MAIN TITLE ---
  mtext(opt$main_title, outer = TRUE, side = 3, cex = 1.4, font = 2, line = 1)
}

## -------------------------
## SAVE PNG
## -------------------------
png(output_filename, width = png_width, height = png_height, res = png_res)
draw_mirror()
invisible(dev.off())
message("Success! Image saved as: ", output_filename)
