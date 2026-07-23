###PCM-wings Mantodea - female pipeline

library(ape)        # core phylogenetics functions in R
library(phytools)   # extra phylogenetic tools (simulation, plotting, ASR, ...)
install.packages("geiger")
library(geiger)

setwd("/home/STUDENTI/pietro.bacconi/Tirocinio_magistrale/04_PCM-wings/00_female")
tree <- read.tree("your_tree.tre")
data <- read.csv("your_dataset.csv", row.names = 1)
# Create a NAMED factor vector of trait states. The names MUST match tip labels
# -- many downstream functions silently rely on this.
dimorphism <- setNames(as.factor(data[, "dimorphism"]), rownames(data))

# A discrete model is just a transition rate matrix Q telling us how likely
# the trait is to move between states per unit of branch length.
# We will compare four biologically meaningful hypotheses:
#
#   ER        -- Equal Rates: gains and losses happen at the same rate
#   ARD       -- All Rates Different: gain and loss rates are independent
#   Loss-only -- only 2 -> 1 transitions allowed (lose dimorphism, never regain)
#   Gain-only -- only 1 -> 2 transitions allowed (gain dimorphism, never lose)

# --- Model 1: Equal Rates (ER) -----------------------------------------------
ER <- fitDiscrete(tree, data, model = "ER")
plot(ER)
dev.off()

# --- Model 2: All Rates Different (ARD) --------------------------------------
ARD <- fitDiscrete(tree, data, model = "ARD")
plot(ARD)
dev.off()
# --- Model 3: Loss-only (custom Q matrix) ------------------------------------
# This model corresponds to STRICT Dollo's law: it implies the trait could only ever be lost. 
model.loss <- matrix(c(0,1,0,0),2,2)
rownames(model.loss) <- colnames(model.loss) <- levels(dimorphism)
loss_only <- fitDiscrete(tree, data, model = model.loss)
plot(loss_only)
dev.off()
# --- Model 4: Gain-only (custom Q matrix) ------------------------------------
# Biologically this is a bit unusual here ... dimorphism can be gained but never lost.
model.gain <- matrix(c(0,1,0,0),2,2)
rownames(model.gain) <- colnames(model.gain) <- levels(dimorphism)
gain_only <- fitDiscrete(tree, data, model = model.gain)
plot(gain_only)
dev.off()
# Compare the four models
lnL <- setNames(
  c(ER$opt$lnL, ARD$opt$lnL, loss_only$opt$lnL, gain_only$opt$lnL),
  c("Equal Rates", "All Rates Different", "Loss-only", "Gain-only")
)

lnL

aicc <- setNames(c(ER$opt$aicc, ARD$opt$aicc, loss_only$opt$aicc, gain_only$opt$aicc),c("Equal Rates", "All Rates Different", "Loss-only", "Gain-only"))
aicw(aicc)


## Plot ancestral trait based on best fitted model
fitARD <- ace(dimorphism, tree, model = "ARD", type = "discrete")
fitER <- ace(dimorphism, tree, model = "ER", type = "discrete")
# Recompute branch lengths so they grow only mildly with depth. This is purely
# COSMETIC -- the inference itself is unchanged; we only do it so the fan plot
# below is readable.
no_br_tree <- compute.brlen(tree, power = 0.4)
# Plot the tree in a fan layout and overlay node / tip pies showing state probabilities.
plot.phylo(no_br_tree,
           type           = "fan",
           lwd            = 1,
           edge.width     = 1,
           label.offset   = 0.02,                            
           show.tip.label = TRUE,
           cex=0.25,                        
           no.margin      = TRUE)                         

# Pies at internal nodes: reconstructed probabilities for each state
nodelabels(node = 1:tree$Nnode + Ntip(tree),
pie  = fitARD$lik.anc,
cex  = 0.2)
# Pies at tips: observed states (weight 1.0 on the observed state)
tiplabels(pie = to.matrix(dimorphism[tree$tip.label], levels(dimorphism)),
cex = 0.2)

state_levels <- levels(dimorphism)
state_cols   <- rainbow(length(state_levels))

# Add legend to the plot
legend("topleft",
        legend = c("Macropterous", "Brachypterous","Apterous"),
        pch    = 21,
        pt.bg  = state_cols,
        pt.cex = 1.5,
        bty    = "n",
        title  = "Wings")
dev.off()
# If you want to plot "loss_only" model introduce these lines:

# Fit the loss-only model with fitMk instead of ace
fitLoss <- fitMk(tree, dimorphism, model = model.loss)
fitLoss
# Get marginal ancestral state probabilities (equivalent to ace()$lik.anc)
ancLoss <- ancr(fitLoss)
ancLoss$ace   # matrix: rows = internal nodes, columns = states (1,2,3)
# Pies at internal nodes: reconstructed probabilities for each state
nodelabels(node = 1:tree$Nnode + Ntip(tree),
           pie  = ancLoss$ace,
           cex  = 0.2)










