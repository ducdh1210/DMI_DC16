#### ---------------- LOAD SAVED OBJECTS -------------------------------------------------
rm(list = ls()); gc()
setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
source("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/graph_utility_functions.R")

#load("PPI2.RData")
load("ppi2_igraph.rda")
load("ppi2_communities.rda")
load("ppi2_clique.rda")

#### ---------------- BASIC INSPECTION -------------------------------------------------

## some test codes to find out commonly used commands
str(communities(ppi2_infomap))
# List of 714
# $ 1  : chr [1:378] "G5" "G21" "G22" "G35" ...
# $ 2  : chr [1:337] "G39" "G49" "G50" "G63" ...
# $ 3  : chr [1:235] "G83" "G410" "G732" "G733" ...
# $ 4  : chr [1:435] "G14" "G71" "G144" "G152" ...

ppi2_community_infoMap = communities(ppi2_infomap)

## get all subgraphs of those modules

# get vector whose names = module index, values = number of vertices in that module
ppi2_infomap_community_size = unlist(lapply(ppi2_community_infoMap, function(module) length(module)))

# get density vector of that module
ppi2_infomap_community_density = getModuleDensity(igraphObject = pp2_igraph, 
                                                  communityObject = ppi2_infomap)
ppi2_infomap_community_density [1:10]
# [1] 0.31785328 0.15947789 0.28110566 0.03691933 0.14796522 0.33583960 0.18402985
# [8] 0.27596945 0.04361809 0.07751010

### suppose we want to use clustering coefficent as the main metrix to 
### select significant modules

ppi2_infomap_community_clusterCoeff = getModuleClusterCoef(igraphObject = pp2_igraph,
                                                           communityObject = ppi2_infomap)
ppi2_infomap_community_clusterCoeff[1:10]
# [1] 0.6654397 0.5544604 0.8145371 0.2400649 0.7024261 0.9070722 0.7817665 0.9580171
# [9] 0.2020234 0.3341065

# Note: if a module only has 2 elements, it means its clustering coefficient = NaN
# such cases should be handled carefully

# find all cases whose clustering coefficients = NaN
sum(is.na(ppi2_infomap_community_clusterCoeff)) # 91 These cases should be excluded
# in downstream analysis

max(ppi2_infomap_community_clusterCoeff[!is.na(ppi2_infomap_community_clusterCoeff)]) # 1
# there should be multiple cases like that, check all of whose cases
length(which(ppi2_infomap_community_clusterCoeff == 1)) # 11
which(ppi2_infomap_community_clusterCoeff == 1)
#  [1] 375 457 591 615 625 627 634 636 643 645 646
# Note: this also is the "name" of the modules, since the modules are named by its index

# get how many vectices contained in each of those modules
ppi2_infomap_community_clusterCoeff[which(ppi2_infomap_community_clusterCoeff == 1)]
# [1] 1 1 1 1 1 1 1 1 1 1 1 --> thus all of those modules only has 1 nodes --> bad cases

#### ---------------- MAKE USE OF UTILITY CODES ON CLIQUE ----------------------

# ppi2_infomap_subgraphs = getAllSubgraphs(igraphObject = pp2_igraph, 
#                                          communityObject = ppi2_infomap)
# ppi2_stat_infomap_subgraph = getSubgraphStat(ppi2_infomap_subgraphs)
# ppi2_stat_infomap_subgraph = getSubgraphStatNoNA(ppi2_stat_infomap_subgraph)
# View(ppi2_infomap_subgraph_stat)

# work with cliq
str(cliq) # list of 58424
ppi2_cliq_membership = getCliqMemebership(igraphObject = pp2_igraph, cliq = ppi2_cliq_unordered)
length(unique(ppi2_cliq_membership)) # [1] 10910 --> thus there are 10910 cliq

# Very important: names of the cliq also imply the relative size order of the clique
# with respect to other cliques, thus cliq 1 is the biggest clique, clique 10910
# is the smallest clique

ppi2_cliq_membership[1:10]
#  G0  G1  G2  G3  G4  G5  G6  G7  G9 G10 
# 173 174 175 176 177 178 179 180 181 182

table(ppi2_cliq_membership)[1:10]
# ppi2_cliq_membership
# 1  2  3  4  5  6  7  8  9 10 
# 73 55 46 45 44 41 32 31 28 26 

length(table(ppi2_cliq_membership)) # [1] 10910

table(table(ppi2_cliq_membership))
# 1     5     6     7     8     9    10    11    12    13    14    15    16    17    19    21    23    25    26    28    31    32 
# 10738    63    25    18    16     6     7     8     3     2     3     1     2     1     1     1     2     2     2     1     1     1 
# 41    44    45    46    55    73 
# 1     1     1     1     1     1 

length(which(table(ppi2_cliq_membership) == 1)) # [1] 10738
# for experimental purpose, let's look at a the biggest cliq, whose size = 73.
# since it is the biggest cliq, it is also the first

# get clique whose size is more than than 1
ppi2_cliqBiggerThanOne = unname(which(table(ppi2_cliq_membership) >1))

ppi2_cliq_ordered = reorderCliqBySize(igraphObject = pp2_igraph, cliq = ppi2_cliq_unordered)
ppi2_cliq_subgraphs = getAllSubgraphs(igraphObject = pp2_igraph, 
                                        clique = ppi2_cliq_ordered[ppi2_cliqBiggerThanOne])

ppi2_stat_cliq_subgraphs = getSubgraphStat(ppi2_cliq_subgraphs)

ppi2_cliq_infomap_subgraphs[[1]]

ppi2_cliq_size1 = unname(which(table(ppi2_cliq_membership)  == 1))
ppi2_cliq_subgraphs_size1 = induced.subgraph(graph = pp2_igraph,
                                             vids = ppi2_cliq_size1)










