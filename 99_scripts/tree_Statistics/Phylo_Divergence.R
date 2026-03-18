#!/usr/bin/env Rscript

library(ape)

# 1. YOUR IQ-TREE PARSER (To safely read your trees with bootstraps)
read_iqtree_nwk <- function(file_path) {
  raw_txt <- paste(readLines(file_path, warn = FALSE), collapse = "")
  mod_txt <- gsub("\\):([0-9\\.eE-]+)\\[([^]]+)\\]", ")\\2:\\1", raw_txt)
  return(ape::read.tree(text = mod_txt))
}

# 2. LOAD YOUR TREES (Put your actual file paths here)
tree1_path <- "mantids_nu_genaf_ML_MP_fam_rt.nwk"
tree2_path <- "mantids_nu_genaf_ML_W3_MP_fam_rt.nwk"

tree1 <- read_iqtree_nwk(tree1_path)
tree2 <- read_iqtree_nwk(tree2_path)

# 3. RUN THE TEXT COMPARISON
message("\n===========================================")
message("       PHYLOGENETIC TREE COMPARISON        ")
message("===========================================\n")

# This calculates the differences and saves the stats
comparison <- comparePhylo(tree1, tree2, force.rooted = TRUE)

# This prints the text statistics directly to your terminal
print(comparison)

# 4. SAVE THE PLOT TO PDF
# We use a large PDF so your 189 tips are readable!
output_pdf <- "Branch_Divergence_Map.pdf"
pdf(output_pdf, width = 20, height = 15)

# We run it again with plot = TRUE. 
# We add cex = 0.4 so the tip labels aren't giant blurry blobs
comparePhylo(tree1, tree2, plot = TRUE, force.rooted = TRUE, cex = 0.4)

invisible(dev.off())
message("\n===========================================")
message("Success! Red branches show topology differences.")
message("Saved map as: ", output_pdf)
