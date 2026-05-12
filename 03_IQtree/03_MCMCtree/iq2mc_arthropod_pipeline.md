# IQ2MC Timetree Pipeline for Arthropod Phylogenies

A reproducible workflow for inferring divergence times on deep arthropod phylogenies (≤ 350 Mya) using **IQ-TREE 3** and **MCMCTree** via the **IQ2MC** integration. Works for both mitochondrial protein-coding data and nuclear BUSCO data.

**Time convention:** 1 time unit = 100 Mya throughout. So a node at 200 Mya has age `2.0`, and the root (≤ 350 Mya) has age `≤ 3.5`.

---

## 0. Setup

### 0.1 Project structure

```bash
mkdir -p ~/timetree_project/{01_data,02_iqtree,03_calibrations,04_iq2mc,05_mcmctree}
cd ~/timetree_project
```

Place your alignment in `01_data/`. Pick **one** of the two scenarios below depending on your data.

### 0.2 Dependencies

```bash
# IQ-TREE 3.0.1+ (provides --dating mcmctree)
iqtree3 --version

# PAML / MCMCTree 4.10.7+
mcmctree -h 2>&1 | head -5

# R for diagnostics
Rscript -e 'library(ape); cat("ape:", as.character(packageVersion("ape")), "\n")'
```

---

## 1. Choose your scenario

### Scenario A — Mitochondrial protein-coding data (insects)

You have a concatenated alignment of 13 mitochondrial protein-coding genes, **already translated to amino acids**, with a partition file separating each gene. Working at the AA level avoids saturation issues at 3rd codon positions for deep insect divergences.

```bash
# Expected files
01_data/mt_aa.phy           # PHYLIP alignment (or .fasta)
01_data/mt_partitions.nex   # NEXUS partition file, one block per gene
```

Example partition file (`mt_partitions.nex`):

```
#nexus
begin sets;
  charset COX1 = 1-512;
  charset COX2 = 513-740;
  charset COX3 = 741-998;
  charset CYTB = 999-1378;
  charset ND1  = 1379-1696;
  charset ND2  = 1697-2032;
  charset ND3  = 2033-2148;
  charset ND4  = 2149-2592;
  charset ND4L = 2593-2684;
  charset ND5  = 2685-3258;
  charset ND6  = 3259-3422;
  charset ATP6 = 3423-3650;
  charset ATP8 = 3651-3704;
end;
```

### Scenario B — BUSCO amino acid data

You have hundreds of single-copy ortholog alignments (one per BUSCO gene), concatenated into a supermatrix with a partition file.

```bash
01_data/busco_aa.phy           # concatenated supermatrix
01_data/busco_partitions.nex   # one charset per BUSCO gene
```

For BUSCO data with thousands of genes, partitioning every gene gets unwieldy. Either:
- **Merge** with ModelFinder (`-m MFP+MERGE`) — practical for ~50–500 genes
- **Pre-cluster** by rate or compositional similarity, then partition the clusters
- **Concatenate without partitioning** and use a profile mixture model (LG+G4+C60)

The pipeline below handles all three approaches.

---

## 2. Step 1 — Best-fit model + ML tree (IQ-TREE 3)

### 2.1 Mitochondrial AA data

```bash
cd 02_iqtree

# Edge-linked partitions (-p): one set of branch lengths, per-partition rate scaler.
# MFP+MERGE: search for the best partitioning scheme by merging similar partitions.
# Use mtInv4 / mtART / LG family models — relevant for invertebrate mt proteins.
# -B 1000: ultrafast bootstrap for branch support.
# -T AUTO: auto-detect threads.

iqtree3 \
  -s ../01_data/mt_aa.phy \
  -p ../01_data/mt_partitions.nex \
  -m MFP+MERGE \
  -mset mtInv,mtART,mtMet,mtZOA,LG \
  -B 1000 \
  -T AUTO \
  --prefix mt_ml
```

**Output to keep:**
- `mt_ml.treefile` — the ML tree with branch lengths in subs/site
- `mt_ml.best_scheme.nex` — the merged partition scheme actually used
- `mt_ml.iqtree` — the full report (BIC scores, model fit)

### 2.2 BUSCO AA data — three options

```bash
cd 02_iqtree

# Option B1: Edge-linked partitions with model merging (recommended for ~50-500 genes)
iqtree3 \
  -s ../01_data/busco_aa.phy \
  -p ../01_data/busco_partitions.nex \
  -m MFP+MERGE \
  -B 1000 \
  -T AUTO \
  --prefix busco_ml_partitioned

# Option B2: Concatenated with profile mixture (recommended for very deep divergences)
iqtree3 \
  -s ../01_data/busco_aa.phy \
  -m LG+G4+C60 \
  -B 1000 \
  -T AUTO \
  --prefix busco_ml_mixture

# Option B3: Edge-unlinked (-Q): each partition gets its own branch lengths
# Use only when you have abundant data per partition and expect heterogeneous evolution
iqtree3 \
  -s ../01_data/busco_aa.phy \
  -Q ../01_data/busco_partitions.nex \
  -m MFP \
  -B 1000 \
  -T AUTO \
  --prefix busco_ml_unlinked
```

**Pick one** for downstream analysis. Compare BIC scores in the `.iqtree` files; lower is better.

---

## 3. Step 2 — Diagnostics on the ML tree

Before adding calibrations, extract metrics that will inform every prior in the control file.

```bash
cd ../03_calibrations

# Replace TREE_PATH with your actual file:
#   ../02_iqtree/mt_ml.treefile         (Scenario A)
#   ../02_iqtree/busco_ml_partitioned.treefile  (Scenario B1)

TREE=../02_iqtree/mt_ml.treefile

Rscript -e '
library(ape)
tr <- read.tree(commandArgs(TRUE)[1])

# ---- Tree geometry ----
nd <- node.depth.edgelength(tr)
tip_d <- nd[1:Ntip(tr)]

cat("=== Tree geometry (subs/site) ===\n")
cat("Max root-to-tip:    ", round(max(tip_d), 4), "\n")
cat("Mean root-to-tip:   ", round(mean(tip_d), 4), "\n")
cat("Median root-to-tip: ", round(median(tip_d), 4), "\n")
cat("SD root-to-tip:     ", round(sd(tip_d), 4), "\n")
cat("Total tree length:  ", round(sum(tr$edge.length), 4), "\n")
cat("Number of tips:     ", Ntip(tr), "\n")
cat("Number of branches: ", nrow(tr$edge), "\n\n")

# ---- Branch length distribution ----
cat("=== Branch length distribution ===\n")
print(round(summary(tr$edge.length), 5))
cat("Min nonzero branch: ", min(tr$edge.length[tr$edge.length > 0]), "\n")
cat("N near-zero (<1e-6):", sum(tr$edge.length < 1e-6), "\n")
cat("N very long (>1.0): ", sum(tr$edge.length > 1.0), "\n\n")

# ---- Clock-likeness ----
cv_tips <- sd(tip_d) / mean(tip_d)
cat("=== Clock-likeness diagnostic ===\n")
cat("CV of root-to-tip:  ", round(cv_tips, 3), "\n")
cat("Interpretation:\n")
cat("  < 0.05  -> nearly clock-like; clock=1 (strict) is OK\n")
cat("  0.05-0.30 -> moderate variation; use clock=2 (independent rates)\n")
cat("  > 0.30  -> strong variation; clock=2 with diffuse sigma2 prior\n\n")

# ---- Suggested rgene_gamma (uses YOUR assumed root age) ----
# 100 Mya = 1 time unit, so root_age_units = root_Mya / 100
ROOT_MYA <- 350      # <-- EDIT THIS to your expected root age in Mya
root_units <- ROOT_MYA / 100

mean_rate <- mean(tip_d) / root_units
cat("=== Rate prior calibration ===\n")
cat("Assumed root age:   ", ROOT_MYA, "Mya =", root_units, "time units\n")
cat("Implied mean rate:  ", round(mean_rate, 4), "subs/site per 100 My\n")

# Pick alpha=2 (diffuse), solve for beta
alpha <- 2
beta  <- alpha / mean_rate
cat("Suggested rgene_gamma = ", alpha, round(beta, 1), "1\n")
cat("  (alpha=2, beta=", round(beta,1), ", a=1 for variable rates among partitions)\n\n")

# ---- Suggested sigma2_gamma from CV ----
sigma2_mean <- log(1 + cv_tips^2)
cat("=== Clock-relaxation prior ===\n")
cat("Implied sigma^2 mean: ", round(sigma2_mean, 4), "\n")
if (sigma2_mean < 0.05) {
  cat("Suggested sigma2_gamma = 1 20 1   (mean=0.05)\n")
} else if (sigma2_mean < 0.15) {
  cat("Suggested sigma2_gamma = 1 10 1   (mean=0.10)\n")
} else {
  cat("Suggested sigma2_gamma = 2 10 1   (mean=0.20)\n")
}
' "$TREE"
```

**Save the output** — these are the values you will plug into the control file in step 6.

---

## 4. Step 3 — Add fossil calibrations to the tree

MCMCTree reads calibrations as labels on internal nodes in the Newick tree.

### 4.1 Calibration syntax (1 unit = 100 Mya)

| Calibration type | Syntax | Example (200–250 Mya) |
|---|---|---|
| Soft minimum + maximum | `B(min, max)` | `'B(2.0, 2.5)'` |
| Soft minimum only | `L(min)` | `'L(2.0)'` |
| Soft maximum only | `U(max)` | `'U(2.5)'` |
| Gamma density | `G(α, β)` | `'G(186, 82.7)'` (mean ≈ 2.25) |
| Skew normal | `SN(loc, scale, shape)` | `'SN(2.2, 0.15, 4)'` |

The default tail probability for soft bounds is 2.5%. To make a bound nearly hard: `B(2.0, 2.5, 1e-300, 1e-300)`.

### 4.2 Common insect/arthropod calibrations (in 100-My units)

These are illustrative — **always check primary literature** for your specific clade. References below.

| Node | Calibration | Justification |
|---|---|---|
| Root (Insecta crown, if used) | `B(3.0, 4.5)` | Devonian Rhynie chert insects (~410 Ma); Cambrian arthropod max |
| Holometabola crown | `B(3.05, 3.50)` | Carboniferous fossils (Misnaepidaisidae) |
| Diptera crown | `B(2.40, 2.90)` | Triassic dipterans |
| Hymenoptera crown | `B(2.30, 2.85)` | Triassic xyelid sawflies |
| Lepidoptera crown | `B(2.00, 2.85)` | Earliest Triassic Lepidoptera-like wings |
| Coleoptera crown | `B(2.85, 3.15)` | Permian Tshekardocoleidae |

### 4.3 Editing the tree

Open your `.treefile` in **FigTree** (`File > Open`), find each internal node, click it, and add the calibration string in the **Node Label** field. Then `File > Export Trees > Newick`. Or do it programmatically:

```bash
# Open the ML tree and a text editor side by side
cp ../02_iqtree/mt_ml.treefile mt_calibrated.tre
nano mt_calibrated.tre
```

Manually annotate. Final tree should look something like:

```
(((Drosophila, Anopheles)'B(0.8, 2.4)', (Apis, Nasonia)'B(2.3, 2.85)'), (Tribolium, Coleomegilla)'B(2.85, 3.15)')'B(3.0, 4.5)';
```

The first line of the file should have the **number of taxa and 1 (for one tree)**, in PAML format:

```bash
# Add the PAML header
ntips=$(grep -o ',' mt_calibrated.tre | wc -l)
ntips=$((ntips + 1))
echo "$ntips 1" > mt_calibrated_paml.tre
cat mt_calibrated.tre >> mt_calibrated_paml.tre
```

### 4.4 Verify the calibrated tree is readable

```bash
Rscript -e '
library(ape)
tr <- read.tree("mt_calibrated.tre")
cat("Tips:", Ntip(tr), "\n")
cat("Internal nodes:", Nnode(tr), "\n")
cat("Node labels (calibrations):\n")
print(tr$node.label[tr$node.label != ""])
'
```

---

## 5. Step 4 — Compute gradient + Hessian (IQ-TREE 3 with `--dating mcmctree`)

This is the IQ2MC step. IQ-TREE re-evaluates the model at the calibrated tree topology (with branch lengths re-optimized) and writes the gradient vector and Hessian matrix to `in.BV`, plus a starter `mcmctree.ctl`.

```bash
cd ../04_iq2mc

# Scenario A: mitochondrial AA, partitioned
iqtree3 \
  -s ../01_data/mt_aa.phy \
  -p ../02_iqtree/mt_ml.best_scheme.nex \
  -te ../03_calibrations/mt_calibrated.tre \
  --dating mcmctree \
  -T AUTO \
  --prefix mt_iq2mc

# Scenario B2: BUSCO AA, concatenated mixture model
iqtree3 \
  -s ../01_data/busco_aa.phy \
  -m LG+G4+C60 \
  -te ../03_calibrations/busco_calibrated.tre \
  --dating mcmctree \
  -T AUTO \
  --prefix busco_iq2mc
```

**Output:**
- `mt_iq2mc.in.BV` — the gradient + Hessian file MCMCTree needs (one block per partition)
- `mt_iq2mc.mcmctree.ctl` — a draft control file (we will heavily edit this)
- `mt_iq2mc.treefile` — the tree with re-optimized branch lengths

Move the BV file and the calibrated tree into the MCMCTree directory:

```bash
cp mt_iq2mc.in.BV ../05_mcmctree/in.BV
cp ../03_calibrations/mt_calibrated_paml.tre ../05_mcmctree/calibrated.tre
```

---

## 6. Step 5 — Build the MCMCTree control file

Now combine everything: the diagnostic numbers from step 3, the calibrated tree from step 4, and the BV file from step 5.

```bash
cd ../05_mcmctree

cat > mcmctree.ctl << 'EOF'
          seed = -1                  * -1 = clock-time random seed; use positive int for reproducibility
       seqfile = /dev/null           * not used when usedata=2 (likelihood from in.BV)
      treefile = calibrated.tre      * Newick tree with calibration labels
      mcmcfile = mcmc.txt            * posterior samples (read into Tracer/R)
       outfile = mcmc.out            * summary of posterior

         ndata = 13                  * NUMBER OF PARTITIONS in your in.BV
                                     * Scenario A (13 mt genes): 13
                                     * Scenario B concatenated:   1
                                     * Scenario B partitioned:    N partitions

       usedata = 2 in.BV             * 2 = approximate likelihood from in.BV (FAST)
                                     * Use 0 first for prior-only check (see step 7)

         clock = 2                   * 2 = independent log-normal rates per branch
                                     * Use clock=1 only if CV root-to-tip < 0.05
                                     * Use clock=3 (autocorrelated) if you expect rate inheritance

       RootAge = 'B(3.0, 3.5)'       * SOFT BOUND on root in 100-My units
                                     * Required if no fossil calibration on root node
                                     * 3.5 = 350 Mya (your stated max depth)

         model = 0                   * irrelevant when usedata=2 (model already baked into in.BV)
         alpha = 0                   * idem
         ncatG = 5                   * idem
     cleandata = 0                   * keep gaps as missing data

       BDparas = 1 1 0               * birth, death, sampling for node-age prior
                                     * 1 1 0 is the standard default for species trees
                                     * If prior on times is U-shaped, try 2 2 0.1

   rgene_gamma = 2 6 1               * <-- FROM STEP 3: alpha=2, beta=mean_rate^-1*alpha, a=1
                                     * Mean rate = 2/6 = 0.33 subs/site/100My
                                     * EDIT beta based on your diagnostic output

  sigma2_gamma = 1 10 1              * <-- FROM STEP 3: prior on log-rate variance
                                     * 1 10 1 -> mean sigma^2 = 0.1 (moderate clock relaxation)

      finetune = 1: .1 .1 .1 .1 .1 .1 .1
                                     * 1: = auto-tune step sizes during burnin
                                     * Order: times, mu, sigma2, rates, mixing, paras, FossilErr

         print = 1                   * 1 = standard output; 2 = also print per-branch rates (large file)
        burnin = 50000               * discard first 50k iterations
      sampfreq = 50                  * sample every 50th iteration after burnin
       nsample = 20000               * total post-burnin samples
                                     * Total iterations = burnin + sampfreq*nsample = 1,050,000
EOF
```

### 6.1 Quick reference: how each metric maps to a control variable

| Diagnostic | Control variable | How it's set |
|---|---|---|
| Tree depth in Mya | `RootAge` and time unit | Use 100 My units; RootAge in those units |
| Mean root-to-tip / root age | `rgene_gamma` β | β = α / mean_rate |
| CV of root-to-tip | `clock` and `sigma2_gamma` | CV > 0.05 → clock=2; mean σ² ≈ ln(1+CV²) |
| Number of partitions | `ndata` | Count blocks in in.BV (= partitions in step 5) |
| (Prior-only check) | `BDparas` | Adjust if prior on node ages is biologically silly |

---

## 7. Step 6 — Prior-only check (mandatory before the real run)

Run MCMCTree with `usedata = 0` first. This samples from the joint prior on times *as actually constructed by MCMCTree* from your calibrations + BDparas + RootAge. The constructed prior is often subtly different from what you'd intuit from the calibrations alone.

```bash
# Make a copy with usedata=0
sed 's/usedata = 2 in.BV/usedata = 0/' mcmctree.ctl > mcmctree_prior.ctl

# Run a short chain
mcmctree mcmctree_prior.ctl > prior_run.log 2>&1

# Inspect the prior on key node ages
mv mcmc.txt mcmc_prior.txt

Rscript -e '
mp <- read.table("mcmc_prior.txt", header=TRUE)
# Time columns are named t_n<NodeNumber>; the first one (t_n<Ntip+1>) is the root
time_cols <- grep("^t_n", colnames(mp), value=TRUE)
cat("Prior summary on calibrated nodes (in 100-My units):\n\n")
for (col in time_cols[1:min(10,length(time_cols))]) {
  q <- quantile(mp[[col]], c(0.025, 0.5, 0.975))
  cat(sprintf("%-10s  median=%.3f   95%% CI: (%.3f, %.3f)   = %.0f-%.0f Mya\n",
              col, q[2], q[1], q[3], q[1]*100, q[3]*100))
}
'
```

**Sanity questions to ask:**
- Is the root prior peaked anywhere reasonable (not 50 Mya, not 800 Mya)?
- Are calibrated nodes centered close to where you set the calibration?
- Are uncalibrated nodes between their parent and child (i.e. ancestor-descendant order respected)?

If the prior is wildly off, don't run the data — fix calibrations or BDparas first.

---

## 8. Step 7 — Run the real chain (twice)

```bash
# Run 1
mkdir run1 && cd run1
cp ../mcmctree.ctl ../in.BV ../calibrated.tre .
mcmctree mcmctree.ctl > run1.log 2>&1 &

# Run 2 (different random seed)
cd ..
mkdir run2 && cd run2
cp ../mcmctree.ctl ../in.BV ../calibrated.tre .
sed -i 's/seed = -1/seed = -2/' mcmctree.ctl
mcmctree mcmctree.ctl > run2.log 2>&1 &

cd ..
wait
```

### 8.1 Monitor acceptance proportions during burnin

In each `run*.log`, you'll see lines like:

```
-50%  0.31 0.28 0.34 0.29 0.27  ...  -lnL=-95291.4
```

The first numbers are acceptance proportions for: **times, μ, σ², rates, mixing, parameters, FossilErr**.

- All in **0.20–0.40** → great
- Below **0.15** or above **0.70** → bad, stop and re-tune `finetune`

With `finetune = 1: ...` the program auto-tunes during burnin, which usually fixes things.

---

## 9. Step 8 — Convergence and posterior diagnostics

```bash
Rscript -e '
m1 <- read.table("run1/mcmc.txt", header=TRUE)
m2 <- read.table("run2/mcmc.txt", header=TRUE)

# Discard first 10% as extra burnin
m1 <- m1[round(0.1*nrow(m1)):nrow(m1), ]
m2 <- m2[round(0.1*nrow(m2)):nrow(m2), ]

# ESS for time variables (rough, no autocorrelation correction here)
# For real ESS use coda::effectiveSize
library(coda)
mc1 <- as.mcmc(m1[, grep("^t_n", colnames(m1))])
mc2 <- as.mcmc(m2[, grep("^t_n", colnames(m2))])

ess1 <- effectiveSize(mc1)
ess2 <- effectiveSize(mc2)

cat("=== ESS summary ===\n")
cat("Run 1: min ESS =", round(min(ess1)), " mean ESS =", round(mean(ess1)), "\n")
cat("Run 2: min ESS =", round(min(ess2)), " mean ESS =", round(mean(ess2)), "\n")
cat("Target: all ESS > 200\n\n")

# Run-to-run congruence: compare posterior medians
med1 <- apply(m1[, grep("^t_n", colnames(m1))], 2, median)
med2 <- apply(m2[, grep("^t_n", colnames(m2))], 2, median)

cat("=== Between-run congruence ===\n")
cat("Max abs diff in posterior medians:", round(max(abs(med1-med2)), 4), "(time units)\n")
cat("Max relative diff:", round(max(abs(med1-med2)/pmax(med1,med2)), 3), "\n")
cat("Should be < 0.02 (= 2 My on 100 My scale) for well-converged runs\n\n")

# Posterior on rate-variation parameter
if ("sigma2_1" %in% colnames(m1)) {
  cat("=== sigma^2 posterior (clock relaxation) ===\n")
  cat("Run 1 mean:", round(mean(m1$sigma2_1), 3), "\n")
  cat("Run 2 mean:", round(mean(m2$sigma2_1), 3), "\n")
  cat("If much higher than prior mean (0.1), increase prior and re-run\n")
}
'
```

### 9.1 Visual inspection (Tracer or R)

```bash
# Open in Tracer for trace plots, ESS, marginal densities
tracer run1/mcmc.txt run2/mcmc.txt &

# Or in R:
Rscript -e '
library(coda)
m1 <- read.table("run1/mcmc.txt", header=TRUE)
mc <- as.mcmc(m1[, grep("^t_n", colnames(m1))[1:6]])  # first 6 nodes
pdf("trace_plots.pdf", width=10, height=8)
plot(mc)
dev.off()
cat("Trace plots written to trace_plots.pdf\n")
'
```

---

## 10. Step 9 — Sensitivity analysis (don't skip)

The IQ2MC paper shows that with sparse calibrations, model choice can shift dates by tens of Myr. Always re-run with at least one alternative model to bound your uncertainty:

```bash
# Re-do step 5 with a simpler model
cd ~/timetree_project/04_iq2mc
iqtree3 \
  -s ../01_data/mt_aa.phy \
  -m mtInv+G4 \
  -te ../03_calibrations/mt_calibrated.tre \
  --dating mcmctree \
  --prefix mt_iq2mc_simple

# Run MCMCTree with the simpler in.BV in a separate folder, compare medians
```

If the simple-model and complex-model medians differ by < 5% across all calibrated nodes, your dates are robust to model choice. If they differ by > 20% on any node, report both and discuss the discrepancy.

---

## 11. Step 10 — Final output

The posterior tree (median ages with 95% HPD) is in `run1/FigTree.tre` (MCMCTree writes this automatically). View it in FigTree:

```bash
figtree run1/FigTree.tre &
```

Show node bars (`Node Bars > 95%HPD`) to display credible intervals. Convert ages back to Mya by multiplying by 100.

---

## Appendix A — Common errors and fixes

| Symptom | Cause | Fix |
|---|---|---|
| `Cannot read tree file` | Missing `Ntips Ntrees` header line | Add `42 1` (or your tip count) as first line |
| Acceptance prop = 0.00 for one parameter | finetune step size too large | Decrease that finetune value 5×, re-run |
| ESS < 50 after long run | Poor mixing or weak data | Increase nsample × 5, or simplify model |
| Posterior σ² >> prior mean | Underestimated clock relaxation | Increase prior mean (e.g. `1 5 1`) and re-run |
| Run 1 and Run 2 give different medians | Not converged | Run longer, check trace plots for trends |
| Node ages stuck at calibration bound | Bound is too tight | Widen calibration or use Cauchy `L(min, p, c)` |

## Appendix B — Time conversion cheatsheet (100-My units)

| Geological boundary | Mya | Units |
|---|---|---|
| Permo-Triassic | 252 | 2.52 |
| Carboniferous–Permian | 299 | 2.99 |
| Devonian–Carboniferous | 359 | 3.59 |
| Silurian–Devonian | 419 | 4.19 |

## Appendix C — References

- **IQ2MC paper:** Demotte et al. 2025. A new framework to infer phylogenetic time trees using IQ-TREE 3 and MCMCTree. (preprint)
- **MCMCTree:** Yang Z. 2007. PAML 4: Phylogenetic analysis by maximum likelihood. *MBE* 24:1586–1591.
- **Approximate likelihood:** dos Reis M, Yang Z. 2011. *MBE* 28:2161–2172.
- **Soft bounds:** Yang Z, Rannala B. 2006. *MBE* 23:212–226.
- **Cauchy bounds:** Inoue J, Donoghue PCJ, Yang Z. 2010. *Syst Biol* 59:74–89.
- **Insect calibrations:** Misof B et al. 2014. Phylogenomics resolves the timing and pattern of insect evolution. *Science* 346:763–767.
