#!/usr/bin/env Rscript

## -------------------------
## INSTALLING/CHECKING PACKAGES
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
  make_option(c("-t", "--trees"), type = "character", default = NULL, 
              help = "Comma-separated list of .nwk files (e.g., 'tree1.nwk,tree2.nwk,tree3.nwk')", metavar = "character"),
  make_option(c("-n", "--titles"), type = "character", default = NULL, 
              help = "Comma-separated list of titles matching the trees (e.g., 'Title 1,Title 2,Title 3')", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = "combined_trees", 
              help = "Base name for output files (no extension)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$trees) || is.null(opt$titles)) {
  print_help(opt_parser)
  stop("You must provide comma-separated lists for both --trees and --titles.", call. = FALSE)
}

# Split the comma-separated strings into vectors and remove extra whitespace
tree_files <- trimws(unlist(strsplit(opt$trees, ",")))
tree_titles <- trimws(unlist(strsplit(opt$titles, ",")))

# Ensure the number of trees matches the number of titles
if (length(tree_files) != length(tree_titles)) {
  stop(sprintf("Mismatch error: You provided %d trees but %d titles. Please make sure they match.", 
               length(tree_files), length(tree_titles)), call. = FALSE)
}

num_trees <- length(tree_files)

## -------------------------
## TYPO AUTO-CORRECTION FUNCTION
## -------------------------
fix_tip_names <- function(tree) {
  tree$tip.label[tree$tip.label == "Haanidae"]  <- "Haaniidae"
  tree$tip.label[tree$tip.label == "Rivetidae"] <- "Rivetinidae"
  return(tree)
}

## -------------------------
## READING, LADDERIZING & FIXING TREES
## -------------------------
# lapply reads and processes all trees in the list automatically
trees_list <- lapply(tree_files, function(file) {
  if (!file.exists(file)) stop(paste("Cannot find file:", file), call. = FALSE)
  
  t <- read.tree(file)
  t <- ladderize(t, right = FALSE)
  t <- fix_tip_names(t)
  return(t)
})

## -------------------------
## COLOR PALETTE
## -------------------------
family_colors <- c(
  "Acanthopidae"      = "#B22222", "Amelidae"          = "#0078D7",
  "Amorphoscelidae"   = "#556B2F", "Coptopterygidae"   = "#CD5C5C",
  "Deroplatyidae"     = "#0000CD", "Empusidae"         = "#8A2BE2",
  "Eremiaphilidae"    = "#20B2AA", "Galinthiadidae"    = "#00008B",
  "Gonypetidae"       = "#808000", "Haaniidae"         = "#D2691E", 
  "Hymenopodidae"     = "#C71585", "Leptomantellidae"  = "#2E8B57",
  "Mantidae"          = "#4B0082", "Metallyticidae"    = "#8B0000",
  "Miomantidae"       = "#4169E1", "Nanomantidae"      = "#228B22",
  "Rivetinidae"       = "#008080", "Thespidae"         = "#B8860B",
  "Toxoderidae"       = "#4682B4", "Outgroup"          = "#708090",
  "Angelidae"         = "#A52A2A", "Chaeteessidae"     = "#2F4F4F",
  "Chroicopteridae"   = "#9400D3", "Dactylopterygidae" = "#BDB76B",
  "Epaphroditidae"    = "#483D8B", "Hoplocoryphidae"   = "#008B8B",
  "Liturgusidae"      = "#8B4513", "Majangidae"        = "#DB7093",
  "Mantoididae"       = "#6B8E23", "Photinaidae"       = "#CD853F"
)

get_tip_colors <- function(tree) {
  colors <- family_colors[tree$tip.label]
  colors[is.na(colors)] <- "#000000"
  return(colors)
}

## -------------------------
## DYNAMIC PLOTTING FUNCTION
## -------------------------
generate_plot <- function() {
  # Calculate rows and columns based on total number of trees
  cols <- ceiling(sqrt(num_trees))
  rows <- ceiling(num_trees / cols)
  
  par(mfrow = c(rows, cols), mar = c(1, 1, 3, 1))
  
  for (i in 1:num_trees) {
    plot(trees_list[[i]], main = tree_titles[i], cex = 0.8, tip.color = get_tip_colors(trees_list[[i]]))
  }
}

## -------------------------
## SAVING OUTPUTS
## -------------------------
output_base <- opt$out
pdf_file <- paste0(output_base, ".pdf")
png_file <- paste0(output_base, ".png")

# Dynamic scaling based on how many rows/cols we have
pdf(pdf_file, width = 6 * ceiling(sqrt(num_trees)), height = 5 * ceiling(num_trees / ceiling(sqrt(num_trees))))
generate_plot()
invisible(dev.off())
cat("Saved PDF version to:", pdf_file, "\n")

png(png_file, width = 600 * ceiling(sqrt(num_trees)), height = 500 * ceiling(num_trees / ceiling(sqrt(num_trees))), res = 100)
generate_plot()
invisible(dev.off())
cat("Saved PNG version to:", png_file, "\n")

cat("Done!\n")
