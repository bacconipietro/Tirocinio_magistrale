## ============================================================ 
Computing Rate of speciation/extinction Mantodea using .... methods 28-07-2026

Both BiSSE and MuSSE are birth-death models of trait-dependent diversification. The logic is the same for both; MuSSE is just BiSSE generalized to >2 states.

Input:
A time-calibrated phylogeny (your timetree)
A trait vector, one state per tip (NA allowed for missing data)

## BiSSE (binary state speciation and extinction),
## MuSSE (a multistate extension of BiSSE)
## QuaSSE (quantitative state speciation and extinction).

Output:
A fitted set of rates (via find.mle, maximum likelihood) or posterior distributions (via mcmc)
You then compare a model where rates depend on the trait vs. a null model where they don't (likelihood ratio test / AIC) to ask: does wing morphology predict diversification?

Important caveat from the paper itself: these methods are correlative only — a significant result means wing state is associated with rate differences, not necessarily causal (could be a correlated trait driving it instead).
## =============================================================

## =============================================================
## Diversitree pipeline: wing trait vs speciation/extinction in mantids
## Two datasets:
##   A) 3-state wing trait (macropterous / brachypterous / apterous) -> MuSSE
##   B) 2-state wing trait (winged / wingless)                       -> BiSSE
## =============================================================

install.packages("diversitree")  # run once
library(diversitree)
library(ape)
library(geiger)
setwd("/home/STUDENTI/pietro.bacconi/Tirocinio_magistrale/04_PCM-WINGS/04_Diversitree/MuSSE")

## -------------------------------------------------------------
## PART A: 3-state trait -> MuSSE
## -----------------------------------------------------
## --- Load data ---

tree.A <- read.tree("aa_iln2_noHPD.tre")           
## use this instead if newick format

traits.A <- read.csv("female_wings_recoded.csv", row.names = 1, stringsAsFactors = FALSE)

# traits.A should have one column, e.g. "wing", with values 1,2,3
# 1 = macropterous, 2 = brachypterous, 3 = apterous
states.A <- traits.A$wing
names(states.A) <- rownames(traits.A)

## Sanity checks
name.check(tree.A, states.A)          # (from geiger, optional) or manually:
setdiff(tree.A$tip.label, names(states.A))
setdiff(names(states.A), tree.A$tip.label)

## --- Build likelihood functions ---
# Null model: no state-dependent diversification (all lambda equal, all mu equal)
lik.musse.full <- make.musse(tree.A, states.A, k = 3)
argnames(lik.musse.full)

# Constrain to a null model where lambda1=lambda2=lambda3 and mu1=mu2=mu3
lik.musse.null <- constrain(lik.musse.full,lambda2 ~ lambda1, lambda3 ~ lambda1, mu2 ~ mu1, mu3 ~ mu1)

## --- Sensible starting parameters ---
p.musse <- starting.point.musse(tree.A, k = 3)

## --- ML fits ---
fit.musse.null <- find.mle(lik.musse.null, p.musse[argnames(lik.musse.null)])
fit.musse.full <- find.mle(lik.musse.full, p.musse)

## --- Compare models: does wing state affect diversification? ---
anova(fit.musse.full, null = fit.musse.null)

## --- (Optional) Bayesian version ---
prior.musse <- make.prior.exponential(2 / (log(Ntip(tree.A)) / max(branching.times(tree.A))))
samples.musse <- mcmc(lik.musse.full, fit.musse.full$par, nsteps = 10000, w = 1, prior = prior.musse)

library(ggplot2)  # or use base hist()
d1 <- density(samples.musse$lambda1)
d2 <- density(samples.musse$lambda2)
d3 <- density(samples.musse$lambda3)

xr <- range(d1$x, d2$x, d3$x)
yr <- range(d1$y, d2$y, d3$y)

plot(d1, col = "blue", main = "Speciation rate by wing state",
     xlim = xr, ylim = yr, xlab = "Speciation rate (lambda)")
lines(d2, col = "darkgreen")
lines(d3, col = "red")
legend("topright", legend = c("macropterous","brachypterous","apterous"), col = c("blue","darkgreen","red"), lty = 1)

