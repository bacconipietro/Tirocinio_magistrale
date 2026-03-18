#!/usr/bin/env Rscript

# ============================================================
# Root-to-Tip Distance Analysis
# Usage: Rscript root_to_tip.R <path_to_tree.nwk> [output_basename]
# Output: - <output_basename>_distances.tsv
#         - <output_basename>_plot.pdf
# If output_basename is omitted, the tree filename is used.
# ============================================================

# --- Check and load required packages -----------------------
required_packages <- c("ape")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Required package not found:", pkg,
               "\nInstall it with: install.packages('", pkg, "')", sep = ""))
  }
}
library(ape)

# --- Parse command-line arguments ---------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("No input file provided.\nUsage: Rscript root_to_tip.R <path_to_tree.nwk> [output_basename]")
}

tree_file <- args[1]

if (!file.exists(tree_file)) {
  stop(paste("File not found:", tree_file))
}

# --- Derive output file names --------------------------------
# If a second argument is provided, use it as the output basename.
# Otherwise, fall back to the tree filename.
if (length(args) >= 2) {
  out_basename <- args[2]
  # Use the tree file's directory unless the user provided a path with folders
  if (dirname(out_basename) == ".") {
    out_dir <- dirname(tree_file)
  } else {
    out_dir      <- dirname(out_basename)
    out_basename <- basename(out_basename)
  }
} else {
  out_basename <- tools::file_path_sans_ext(basename(tree_file))
  out_dir      <- dirname(tree_file)
}

tsv_out <- file.path(out_dir, paste0(out_basename, "_distances.tsv"))
pdf_out <- file.path(out_dir, paste0(out_basename, "_plot.pdf"))

# --- Load tree ----------------------------------------------
cat("Reading tree from:", tree_file, "\n")
tree <- read.tree(tree_file)

if (is.null(tree)) {
  stop("Could not parse the tree file. Make sure it is a valid Newick (.nwk) format.")
}

cat("Tree loaded:", Ntip(tree), "tips,", Nnode(tree), "internal nodes\n")

# --- Calculate root-to-tip distances ------------------------
# node.depth.edgelength returns distances for ALL nodes (tips first, then internal).
# We keep only the first Ntip(tree) entries, which correspond to the tips.
all_distances <- node.depth.edgelength(tree)
tip_distances <- all_distances[1:Ntip(tree)]
names(tip_distances) <- tree$tip.label

# --- Save TSV -----------------------------------------------
df <- data.frame(
  tip_label = names(tip_distances),
  root_to_tip_distance = tip_distances,
  row.names = NULL
)
df <- df[order(df$root_to_tip_distance, decreasing = TRUE), ]

write.table(df, file = tsv_out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Distances saved to:", tsv_out, "\n")

# --- Generate plot ------------------------------------------
# Tips are ordered by distance for a cleaner visual
df$tip_label <- factor(df$tip_label, levels = df$tip_label)

pdf(pdf_out, width = 10, height = max(6, Ntip(tree) * 0.25))

# -- Margin: extra left space for long tip labels
par(mar = c(5, max(8, max(nchar(tree$tip.label)) * 0.5), 4, 2))

barplot(
  df$root_to_tip_distance,
  names.arg  = df$tip_label,
  horiz      = TRUE,
  las        = 1,
  col        = "#4C9BE8",
  border     = "#2C6FAC",
  xlab       = "Root-to-Tip Distance",
  main       = paste("Root-to-Tip Distances\n", basename(tree_file)),
  cex.names  = max(0.5, 1 - Ntip(tree) / 200),  # shrink labels for large trees
  cex.axis   = 0.85
)

# Add a mean line
abline(v = mean(df$root_to_tip_distance), col = "#E84C4C", lty = 2, lwd = 1.5)
legend("bottomright",
       legend = paste("Mean =", round(mean(df$root_to_tip_distance), 5)),
       col    = "#E84C4C",
       lty    = 2,
       lwd    = 1.5,
       bty    = "n",
       cex    = 0.85)

dev.off()
cat("Plot saved to:", pdf_out, "\n")

# --- Summary stats ------------------------------------------
cat("\n--- Summary ---\n")
cat("Min distance  :", round(min(tip_distances), 6), "\n")
cat("Max distance  :", round(max(tip_distances), 6), "\n")
cat("Mean distance :", round(mean(tip_distances), 6), "\n")
cat("Median        :", round(median(tip_distances), 6), "\n")
cat("Done.\n")
