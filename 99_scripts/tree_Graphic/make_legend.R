#!/usr/bin/env Rscript

# ---------------------------------------------------------
# STANDALONE LEGEND GENERATOR FOR PHYLOGENETIC TREES
# ---------------------------------------------------------

output_filename <- "Family_Color_Legend.png"

# Define the exact dictionary of your families and colors
family_palette <- c(
  "Acanthopidae"     = "#B22222",
  "Amelidae"         = "#0078D7",
  "Amorphoscelidae"  = "#556B2F",
  "Coptopterygidae"  = "#CD5C5C",
  "Deroplatyidae"    = "#0000CD",
  "Empusidae"        = "#8A2BE2",
  "Eremiaphilidae"   = "#008B8B",
  "Galinthiadidae"   = "#00008B",
  "Gonypetidae"      = "#808000",
  "Haaniidae"        = "#D2691E",
  "Hymenopodidae"    = "#C71585",
  "Leptomantellidae" = "#2E8B57",
  "Mantidae"         = "#4B0082",
  "Metallyticidae"   = "#8B0000",
  "Miomantidae"      = "#4169E1",
  "Nanomantidae"     = "#228B22",
  "Rivetinidae"      = "#008080",
  "Thespidae"        = "#B8860B",
  "Toxoderidae"      = "#4682B4",
  "Outgroups"        = "#708090",
)

# Set up the PNG canvas (adjust width/height if you need it wider or taller)
png(output_filename, width = 1500, height = 2000, res = 300)

# Create a completely blank canvas with no margins
par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")

# Draw the legend right in the center
legend(x = "center",
       legend = names(family_palette),
       fill = family_palette,       # Creates the colored boxes
       border = "black",            # Puts a crisp black border around each color box
       cex = 1.2,                   # Text size
       bty = "n",                   # Removes the big bounding box around the whole legend
       title = "Families",          # Title of the legend
       title.font = 2)              # Makes the title bold

invisible(dev.off())
message("Legend saved successfully as: ", output_filename)
