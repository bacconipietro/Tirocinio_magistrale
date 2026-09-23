setwd("/home/STUDENTI/pietro.bacconi/Tirocinio_magistrale/04_PCM-wings/03_combined")
library(phytools)
Loading required package: ape
Loading required package: maps
librWarning messages:
1: package ‘ape’ was built under R version 4.4.3 
2: package ‘maps’ was built under R version 4.4.3 
library(geiger)
tree <- read.tree("aa_iln2_noHPD.tre")
fem <- read.csv("female_wings_recoded_bin.csv", row.names = 1)  # where winged (macropterous+brachypterous) = 1 and wingless(apterous)  = 2
male <- read.csv("male_wings_recoded_bin.csv", row.names = 1)   # same binary pattern
x <- setNames(as.character(fem$wings), rownames(fem))
y <- setNames(as.character(male$wings), rownames(male))
length(x)
[1] 189
length(y)
[1] 189
fit_dep <- fitPagel(tree, x, y)
print(fit_dep)
