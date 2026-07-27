##Installing packages for BioGeoBEARS pipeline


#####Aggiorniamo la versione di R del terminale 
conda create -n biogeo 
conda activate biogeo
conda install -c conda-forge r-base=4.4
conda install -c conda-forge r-xml r-rnexml r-phylobase

screen -S biogeobears_DEC_7

R

#####Installiamo i pacchetti necessari
install.packages(c(
  "ape", "phangorn", "optimx", "GenSA", "FD", "snow", "numDeriv",
  "Rcpp", "remotes", "rexpokit", "cladoRcpp", "plotrix", "gdata",
  "phylobase", "minqa", "fdrtool", "statmod", "spam", "MultinomialCI",
  "gtools", "scatterplot3d",
  "phytools", "expm", "devtools", "SparseM", "httr", "stringr"
))

install.packages(GenSA)


#####Installiamo BioGeoBEARS
library(devtools)
devtools::install_github(repo="nmatzke/BioGeoBEARS", INSTALL_opts="--byte-compile")

# Open libraries
library(devtools)
library(phytools)
library(BioGeoBEARS)
setwd("/home/STUDENTI/pietro.bacconi/Tirocinio_magistrale/05_BioGeoBEARS")

## load the tree
tree<-read.tree("GBM2_clean_noOutgroup_speciesNames.tre")
## read data from file
mantid.data<-getranges_from_LagrangePHYLIP(lgdata_fn="tmp_7areas.txt")


## Plot Tip Geo Distribution
tmp<-mantid.data@df

## re-name the columns of the matrix to correspond
## to our geographical areas
colnames(tmp)<-c("R.Nearctic", "R.Neotropical", "R.Palearctic", "R.Afrotropical", "R.Madagascan", "R.Indo-Malaysian", "R.Asutralasian")
## set the colors we’ll use for plotting
colors<-setNames(replicate(ncol(tmp),setNames(c("white","darkgreen"),0:1),simplify=FALSE),colnames(tmp))



## Image settings
pdf("mantid_geodistribution.pdf", width=14, height=24)
## 1st value = bottom margin, 2nd value = left margin, 3rd value = top margin, 4th value = right margin
par(mar=c(1,1,2,4))
## graph presence/absence using plotTree.datamatrix
object<-plotTree.datamatrix(tree,tmp,fsize=0.5,yexp=1.1,header=TRUE,xexp=1.25,colors=colors,sep=0.05)
## add a legend
legend("topleft",c("species absent","species present"),pch=22,pt.bg=c("white","darkgray"),pt.cex=1.3,cex=0.8,bty="n")
## save pdf or png, couple command wiht pdf()/png()  
dev.off()


## set the maximum number of areas that a single species is allowed to occupy.
max_range_size<-5
bgb_run<-define_BioGeoBEARS_run(num_cores_to_use=12,max_range_size=max_range_size,trfn="./GBM2_clean_noOutgroup_speciesNames.tre", return_condlikes_table=TRUE)
## update definition of list element geogfn
bgb_run$geogfn<-"./tmp_7areas.txt"
## check if object is complete and there are no issues
check_BioGeoBEARS_run(bgb_run)
## Now we can optimize our DEC model using maximum likelihood.
DEC.fit<-bears_optim_run(bgb_run)

## --- Plot Ancestral Geographic State
pdf("DEC_ancestral_states.pdf", width=12, height=40)
par(mar=c(4.1,0.1,3.1,1.5), cex=0.8)
plot_BioGeoBEARS_results(DEC.fit,
                          analysis_titletxt="DEC model",
                          plotlegend=FALSE,
                          tipcex=0.6,
                          statecex=0.4,
                          splitcex=0.65,
                          pie_tip_statecex=0.6,
                          root.edge=FALSE,        # <- remove unrelevant root data
                          plot_max_age=85)        # <- maxe range for x axis in Ma
## --- reconstruction legend input ---
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn = DEC.fit$inputs$geogfn)
areas <- getareas_from_tipranges_object(tipranges)
max_range_size <- DEC.fit$inputs$max_range_size
possible_ranges_list_txt <- areas_list_to_states_list_new(areas, maxareas=max_range_size,split_ABC=FALSE,include_null_range=DEC.fit$inputs$include_null_range)
colors_matrix <- get_colors_for_numareas(length(areas))
states_list_0based_index <- cladoRcpp::rcpp_areas_list_to_states_list(areas, maxareas=max_range_size, include_null_range=DEC.fit$inputs$include_null_range)
colors_list_for_states <- mix_colors_for_states(colors_matrix, states_list_0based_index,plot_null_range=DEC.fit$inputs$include_null_range)
## --- prune the legend for only extant states ---
relprobs_matrix <- DEC.fit$ML_marginal_prob_each_state_at_branch_top_AT_node
MLstates_all <- get_ML_states_from_relprobs(relprobs_matrix, possible_ranges_list_txt, returnwhat="states", if_ties="takefirst")
used_states <- unique(MLstates_all)
keep_idx <- possible_ranges_list_txt %in% used_states
legend_txt <- possible_ranges_list_txt[keep_idx]
legend_colors <- colors_list_for_states[keep_idx]
## --- draw pruned legend where you wish ---
colors_legend(legend_txt, legend_colors,
              location="left",
              make_blank_plot_first=FALSE,
              legend_ncol=4,          
              legend_cex=0.7)
dev.off()
