sensiPhy 27-07-2026 Paterno et al. 2018

install.packages("sensiPhy")
library(sensiPhy)

setwd("/home/STUDENTI/pietro.bacconi/Tirocinio_magistrale/04_PCM-wings/05_sensiPhy/00_sd")
library(sensiPhy)
library(ape)

# 1. Load your data -----------------------------------------------------
tree <- read.tree("aa_iln2_noHPD.tre")
wings_df <- read.csv("binary_pres-abs_wing_female.csv", row.names = 1)

# 2. Build a named factor vector (sensiPhy discrete functions want this) --
wings <- as.factor(as.character(wings_df$wings))
names(wings) <- rownames(wings_df)

# Sanity check: tip labels must match names(wings) exactly
setdiff(names(wings), tree$tip.label)   # should be character(0)
setdiff(tree$tip.label, names(wings))   # should be character(0)

# 3. IMPORTANT: use the SAME model/transform as your main analysis -------
# sensiPhy's discrete functions call geiger::fitDiscrete under the hood.
# Whatever model you already fit (ER, SYM, ARD, or "meristic" if you 
# treated the 3 states as an ordered stepwise reduction of wings) and
# whatever tree transform (none, lambda, kappa, delta, EB) you used 
# should go here too, so the sensitivity results are comparable to your
# original estimate.

model     <- "ARD"     # <- replace with what you actually used
transform <- "none"    # <- replace with what you actually used

# 4. Influential species (leave-one-out) ----------------------------------
influ_wings <- influ_discrete(data = wings, phy = tree, model = model, transform = transform, cutoff = 2, n.cores = 2, track = TRUE)
saveRDS(influ_wings, file = "influ_wings.rds")



summary(influ_wings)
sensi_plot(influ_wings)

# 5. Sample size sensitivity ----------------------------------------------
samp_wings <- samp_discrete(data = wings, phy = tree,
                             model = model, transform = transform,
                             breaks = seq(.1, .3, .1),   # remove 10/20/30% of tips
                             n.sim = 100, n.cores = 2, track = TRUE)

saveRDS(samp_wings, file = "samp_wings.rds")
summary(samp_wings)
sensi_plot(samp_wings)
