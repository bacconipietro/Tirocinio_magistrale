#!/usr/bin/env Rscript

## -------------------------
## INSTALLING PACKAGES
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
  make_option(c("-o", "--out"), type = "character", default = "Ultimate_Mirrored_Trees"),
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
output_base <- sub("\\.(png|pdf)$", "", output_filename, ignore.case = TRUE)

## -------------------------
## SETTINGS & AESTHETICS
## -------------------------
plot_type         <- "phylogram"   
italic_tip_labels <- FALSE

# Tree Lines & Colors
tip_cex           <- 0.35
tree_lwd          <- 1.2
base_edge_col     <- "grey35"

# Bootstraps
show_bootstraps   <- TRUE
boot_cex          <- 0.25      
boot_font         <- 1         
boot_col          <- "grey20"  

# Brackets
family_cex        <- 0.55      

# Output Dimensions
png_width  <- 6000
png_height <- 8000
png_res    <- 300

pdf_width  <- 16
pdf_height <- 22

family_palette <- c(
  "Acanthopidae"     = "#B22222", "Amelidae"         = "#0078D7",
  "Amorphoscelidae"  = "#556B2F", "Coptopterygidae"  = "#CD5C5C",
  "Deroplatyidae"    = "#0000CD", "Empusidae"        = "#8A2BE2",
  "Eremiaphilidae"   = "#20B2AA", "Galinthiadidae"   = "#00008B", 
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
family_palette <- family_palette[order(-nchar(names(family_palette)))]

## -------------------------
## PARSERS & HELPERS
## -------------------------
read_iqtree_nwk <- function(file_path) {
  raw_txt <- paste(readLines(file_path, warn = FALSE), collapse = "")
  mod_txt <- gsub("\\):([0-9\\.eE-]+)\\[([^]]+)\\]", ")\\2:\\1", raw_txt)
  return(ape::read.tree(text = mod_txt))
}

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

draw_family_brackets <- function(tree, side = "left") {
  env <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  N <- length(tree$tip.label)
  y_coords <- env$yy[1:N]
  tips <- tree$tip.label
  
  fams <- rep("Unknown", N)
  for(i in 1:N) {
    for(fam in names(family_palette)) {
      if(grepl(fam, tips[i], ignore.case=TRUE)) { fams[i] <- fam; break }
    }
  }
  
  df <- data.frame(tip = tips, y = y_coords, fam = fams, stringsAsFactors = FALSE)
  df <- df[order(df$y), ]
  
  runs <- rle(df$fam)
  end_idx <- cumsum(runs$lengths)
  start_idx <- c(1, end_idx[-length(end_idx)] + 1)
  
  usr <- par("usr")
  max_width <- max(strwidth(tips, cex = tip_cex))
  
  par(xpd = NA) 
  
  for(i in seq_along(runs$values)) {
    fam_name <- runs$values[i]
    if(fam_name == "Unknown") next
    
    y_start <- df$y[start_idx[i]]
    y_end <- df$y[end_idx[i]]
    y_mid <- (y_start + y_end) / 2
    col <- family_palette[[fam_name]]
    
    if (y_start == y_end) {
      y_start <- y_start - 0.3
      y_end <- y_end + 0.3
    }
    
    if(side == "left") {
      x_bar <- usr[2] + max_width + (usr[2] - usr[1]) * 0.01
      x_txt <- x_bar + (usr[2] - usr[1]) * 0.015
      
      segments(x0 = x_bar, y0 = y_start, x1 = x_bar, y1 = y_end, col = col, lwd = 2.5)
      segments(x0 = x_bar - (usr[2]-usr[1])*0.005, y0 = c(y_start, y_end), x1 = x_bar, y1 = c(y_start, y_end), col = col, lwd = 2.5)
      text(x = x_txt, y = y_mid, labels = fam_name, col = col, adj = 0, cex = family_cex, font = 2)
      
    } else {
      x_bar <- usr[1] - max_width - (usr[2] - usr[1]) * 0.01
      x_txt <- x_bar - (usr[2] - usr[1]) * 0.015
      
      segments(x0 = x_bar, y0 = y_start, x1 = x_bar, y1 = y_end, col = col, lwd = 2.5)
      segments(x0 = x_bar + (usr[2]-usr[1])*0.005, y0 = c(y_start, y_end), x1 = x_bar, y1 = c(y_start, y_end), col = col, lwd = 2.5)
      text(x = x_txt, y = y_mid, labels = fam_name, col = col, adj = 1, cex = family_cex, font = 2)
    }
  }
  par(xpd = FALSE)
}

## -------------------------
## READ & PROCESS TREES
## -------------------------
tree1 <- read_iqtree_nwk(opt$tree1)
tree2 <- read_iqtree_nwk(opt$tree2)

message("\nRotating trees to mirror...")
comp_trees <- phytools::cophylo(tr1 = tree1, tr2 = tree2)
tree1_rotated <- comp_trees$trees[[1]]
tree2_rotated <- comp_trees$trees[[2]]

left_tip_colors  <- search_tip_colors(tree1_rotated$tip.label)
right_tip_colors <- search_tip_colors(tree2_rotated$tip.label)

## -------------------------
## PLOT FUNCTION
## -------------------------
draw_mirror <- function() {
  layout(matrix(c(1, 2), nrow = 1))
  par(oma = c(1, 0, 4, 0))

  # --- LEFT TREE ---
  par(mar = c(1, 0, 2, 14), font = if (italic_tip_labels) 3 else 1)
  ape::plot.phylo(
    tree1_rotated, type = plot_type, direction = "rightwards",
    edge.width = tree_lwd, edge.color = base_edge_col, cex = tip_cex, 
    tip.color = left_tip_colors, no.margin = FALSE
  )
  mtext(opt$left_title, side = 3, line = 0, cex = 1.1, font = 2) 
  
  # Bootstraps hover slightly above the node
  if (show_bootstraps && !is.null(tree1_rotated$node.label)) {
    nl <- tree1_rotated$node.label
    nl[nl == ""] <- NA
    nodelabels(nl, frame = "none", cex = boot_cex, font = boot_font, col = boot_col, adj = c(1.1, -0.4))
  }
  
  draw_family_brackets(tree1_rotated, side = "left")

  # --- RIGHT TREE ---
  par(mar = c(1, 14, 2, 0), font = if (italic_tip_labels) 3 else 1)
  ape::plot.phylo(
    tree2_rotated, type = plot_type, direction = "leftwards",
    edge.width = tree_lwd, edge.color = base_edge_col, cex = tip_cex, 
    tip.color = right_tip_colors, no.margin = FALSE
  )
  mtext(opt$right_title, side = 3, line = 0, cex = 1.1, font = 2)

  if (show_bootstraps && !is.null(tree2_rotated$node.label)) {
    nl <- tree2_rotated$node.label
    nl[nl == ""] <- NA
    nodelabels(nl, frame = "none", cex = boot_cex, font = boot_font, col = boot_col, adj = c(-0.1, -0.4))
  }
  
  draw_family_brackets(tree2_rotated, side = "right")

  # --- SHARED MAIN TITLE ---
  mtext(opt$main_title, outer = TRUE, side = 3, cex = 1.4, font = 2, line = 1)
}

## -------------------------
## SAVE PNG AND PDF
## -------------------------
png_out <- paste0(output_base, ".png")
pdf_out <- paste0(output_base, ".pdf")

message("\nSaving figures...")
png(png_out, width = png_width, height = png_height, res = png_res)
draw_mirror()
invisible(dev.off())
message(" -> Saved PNG: ", png_out)

pdf(pdf_out, width = pdf_width, height = pdf_height)
draw_mirror()
invisible(dev.off())
message(" -> Saved PDF: ", pdf_out)
message("\nAll done!")
