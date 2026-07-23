## 

#Set the working directory to wherever you saved the data files.
#Edit this path to match YOUR machine!
setwd("C:/Users/Utente/OneDrive - Alma Mater Studiorum Università di Bologna/Desktop/wings")

# Load the libraries we will use throughout the tutorial
library(ape)        # core phylogenetics functions in R
library(phytools)   # extra phylogenetic tools (simulation, plotting, ASR, ...)
library(geiger)     # model fitting and utilities (name.check, fitDiscrete, ...)

# Here we will test this hypothesis using a
# comprehensive timetree of Mantodea and a binary trait coding wing
# presence/absence:
#       1  = wings macropterous 
#       2  = wings brachypterous
#       3  = wings apterous
# Read the phylogeny (Newick format)
tree <- read.tree("GBM2_clean_noOutgroup.tre")

# Read the trait data. row.names = 1 tells R that the first column contains
# species names that will be used as row identifiers.
data <- read.csv("wings_recoded.csv", row.names = 1)
# Tree and data table are almost always slightly mismatched (e.g. species in
# one but not the other). name.check() spots discrepancies and returns a list:
#   $tree_not_data  -- tips on the tree with no trait info
#   $data_not_tree  -- trait rows without a corresponding tip
chk <- name.check(tree, data)
# Prune the tree to keep only species that also appear in the data ...
tree.pruned <- drop.tip(tree, chk$tree_not_data)

# ... and prune the data to keep only species that are also on the tree.
# drop = FALSE keeps the result as a data.frame even if a single column remains.
data.pruned <- data[!(rownames(data) %in% chk$data_not_tree), , drop = FALSE]

# Create a NAMED factor vector of trait states. The names MUST match tip labels
# -- many downstream functions silently rely on this.
wings <- setNames(as.factor(data.pruned[, "wings"]), rownames(data.pruned))


# A discrete model is just a transition rate matrix Q telling us how likely
# the trait is to move between states per unit of branch length.
# We will compare four biologically meaningful hypotheses:
#
#   ER        -- Equal Rates: gains and losses happen at the same rate
#   ARD       -- All Rates Different: gain and loss rates are independent
#   Loss-only -- only 1 -> 2 transitions allowed (lose wings, never regain)
#   Gain-only -- only 2 -> 1 transitions allowed (gain wings, never lose)
#
# "All models are wrong, but some are useful." -- George Box


# --- Model 1: Equal Rates (ER) -----------------------------------------------
ER <- fitDiscrete(tree.pruned, data.pruned, model = "ER")
plot(ER)

# --- Model 2: All Rates Different (ARD) --------------------------------------
ARD <- fitDiscrete(tree.pruned, data.pruned, model = "ARD")
plot(ARD)

# --- Model 3: Loss-only (custom Q matrix) ------------------------------------
# This model corresponds to STRICT Dollo's law: it implies the trait could only ever be lost. 
model.loss <- matrix(0, 3, 3)
model.loss[1, 2] <- 1   # M -> B
model.loss[2, 3] <- 1   # B -> A  (same "1" = same rate; use 2 here instead if you want separate rates per step)
rownames(model.loss) <- colnames(model.loss) <- levels(wings)
loss_only <- fitDiscrete(tree.pruned, data.pruned, model = model.loss)
plot(loss_only)

# --- Model 4: Gain-only (custom Q matrix) ------------------------------------
# Biologically this is a bit unusual here ... wings can be gained but never lost.
model.gain <- matrix(0, 3, 3)
model.gain[3, 2] <- 1   # A -> B
model.gain[2, 1] <- 1   # B -> M
rownames(model.gain) <- colnames(model.gain) <- levels(wings)
gain_only <- fitDiscrete(tree.pruned, data.pruned, model = model.gain)
plot(gain_only)


# 1.3  Compare the four models
lnL <- setNames(
  c(ER$opt$lnL, ARD$opt$lnL, loss_only$opt$lnL, gain_only$opt$lnL),
  c("Equal Rates", "All Rates Different", "Loss-only", "Gain-only")
)
lnL
# Equal Rates All Rates Different           Loss-only           Gain-only 
#          -97.43808           -80.42821           -86.02903           -99.46997 
aicc <- setNames(c(ER$opt$aicc, ARD$opt$aicc, loss_only$opt$aicc, gain_only$opt$aicc),c("Equal Rates", "All Rates Different", "Loss-only", "Gain-only"))
aicw(aicc)
#                         fit      delta            w
# Equal Rates         196.8977 23.5771690 4.507049e-06
# All Rates Different 173.3205  0.0000000 5.937580e-01
# Loss-only           174.0796  0.7590708 4.062369e-01
# Gain-only           200.9615 27.6409513 5.908172e-07
fitARD <- ace(wings, tree.pruned, model = "ARD", type = "discrete")

# Recompute branch lengths so they grow only mildly with depth. This is purely
# COSMETIC -- the inference itself is unchanged; we only do it so the fan plot
# below is readable.
no_br_tree <- compute.brlen(tree.pruned, power = 0.4)
# Plot the tree in a fan layout and overlay node / tip pies showing state probabilities.
plot.phylo(no_br_tree,
           type           = "fan",
           lwd            = 1,
           edge.width     = 1,
           label.offset   = 2,
           show.tip.label = FALSE,
           no.margin      = TRUE)

# Pies at internal nodes: reconstructed probabilities for each state
nodelabels(node = 1:tree.pruned$Nnode + Ntip(tree.pruned),
pie  = fitARD$lik.anc,
cex  = 0.2)
# Pies at tips: observed states (weight 1.0 on the observed state)
tiplabels(pie = to.matrix(wings[tree.pruned$tip.label], levels(wings)),
cex = 0.2)

state_levels <- levels(wings)
state_cols   <- rainbow(length(state_levels))
 
# Add legend to the plot
legend("topleft",
        legend = c("Mscropterous", "Brachypterous","Apterous"),
        pch    = 21,
        pt.bg  = state_cols,
        pt.cex = 1.5,
        bty    = "n",
        title  = "Wings")
fitARD
