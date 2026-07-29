## -------------------------------------------------------------
## PART B: 2-state trait -> BiSSE
## -------------------------------------------------------------
setwd("/home/STUDENTI/pietro.bacconi/Tirocinio_magistrale/04_PCM-wings/BiSSE")
tree.B  <- read.tree("GBM2_clean.tre")         # adjust if it's the same tree/subset
traits.B <- read.csv("binary_pres-abs_wing_female.csv", row.names = 1, stringsAsFactors = FALSE)
# traits.B$wing should be coded 0 (wingless) / 1 (winged)
states.B <- traits.B$wings
names(states.B) <- rownames(traits.B)
name.check(tree.B, states.B)

## --- Build likelihood functions ---
lik.bisse.full <- make.bisse(tree.B, states.B)
argnames(lik.bisse.full)   # lambda0 lambda1 mu0 mu1 q01 q10

# Null: state-independent diversification
lik.bisse.null <- constrain(lik.bisse.full, lambda1 ~ lambda0, mu1 ~ mu0)

## --- Starting params & ML fits ---
p.bisse <- starting.point.bisse(tree.B)
fit.bisse.null <- find.mle(lik.bisse.null, p.bisse[argnames(lik.bisse.null)])
fit.bisse.full <- find.mle(lik.bisse.full, p.bisse)

## --- Likelihood ratio test ---
anova(fit.bisse.full, null = fit.bisse.null)

## --- (Optional) Bayesian version ---
prior.bisse <- make.prior.exponential(2 / (log(Ntip(tree.B)) / max(branching.times(tree.B))))
samples.bisse <- mcmc(lik.bisse.full, fit.bisse.full$par, nsteps = 10000, w = 1, prior = prior.bisse)

d0 <- density(samples.bisse$lambda0)
d1 <- density(samples.bisse$lambda1)

xr <- range(d0$x, d1$x)
yr <- range(d0$y, d1$y)

plot(d0, col = "blue", main = "Speciation rate by wing state", xlim = xr, ylim = yr, xlab = "Speciation rate (lambda)")
lines(d1, col = "red")
legend("topright", legend = c("dimorphism (0)","no dimorphism (1)"), col = c("blue","red"), lty = 1)




d0 <- density(samples.bisse$lambda0)
d1 <- density(samples.bisse$lambda1)

xr <- range(d0$x, d1$x)
yr <- range(d0$y, d1$y)

plot(d0, col = "blue", main = "Speciation rate by wing state", xlim = xr, ylim = yr, xlab = "Speciation rate (lambda)")
lines(d1, col = "red")
legend("topright", legend = c("wingless (0)","winged (1)"), col = c("blue","red"), lty = 1)

## -------------------------------------------------------------
## Notes / things to watch for (per FitzJohn 2012)
## -------------------------------------------------------------
# - If your tree is incomplete (not all extant mantid species sampled),
#   consider the "unresolved clade" or "skeleton tree" corrections:
#   make.bisse(tree, states, sampling.f = c(f0, f1))  # proportion sampled per state
#   make.musse(tree, states, k = 3, sampling.f = c(f1, f2, f3))
#
# - For MuSSE with 3 states you already have 12 free parameters in the full model;
#   if convergence is poor, consider constraining extinction (mu) to be equal
#   across states first, and only let lambda vary — mirrors the paper's approach
#   of "free lambda main effects, single mu/q intercepts."
#
# - Always check MCMC trace/mixing (e.g. via coda::traceplot) before trusting
#   posterior intervals, and drop an appropriate burn-in (paper used first 500/10000).
