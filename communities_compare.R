#### ---------------- LOAD SAVED OBJECTS -------------------------------------------------
rm(list = ls())
setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
source("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/graph_utility_functions.R")

#load("PPI2.RData")
load("ppi2_igraph.rda")
load("ppi2_communities.rda")
load("ppi2_clique.rda")
load("num_edges_ppi2.rda")

ppi2_infomap_copied = ppi2_infomap

##### basic built-in compare functions

# Possible values: ‘vi’ is the variation of information (VI) metric of Meila (2003), 
# ‘nmi’ is the normalized mutual information measure proposed by Danon et al. (2005), 
# ‘split.join’ is the split-join distance of can Dongen (2000), 
# ‘rand’ is the Rand index of Rand (1971),
# ‘adjusted.rand’ is the adjusted Rand index by Hubert and Arabie (1985).
compare(ppi2_infomap, ppi2_infomap_copied, method = "adjusted.rand") # 1
