#### ---------------- LOAD SAVED OBJECTS ---------------------------------------
rm(list = ls())
setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
source("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/graph_utility_functions.R")

#load("PPI2.RData")
load("ppi2_igraph.rda")
load("ppi2_communities.rda")
load("ppi2_clique.rda")
load("num_edges_ppi2.rda")
load("ppi2_compare_summary.rda")

ppi2_infomap_copied = ppi2_infomap

#### ---------------- BASIC COMPARE -------------------------------------------

# Possible values: ‘vi’ is the variation of information (VI) metric of Meila (2003), 
# ‘nmi’ is the normalized mutual information measure proposed by Danon et al. (2005), 
# ‘split.join’ is the split-join distance of can Dongen (2000), 
# ‘rand’ is the Rand index of Rand (1971),
# ‘adjusted.rand’ is the adjusted Rand index by Hubert and Arabie (1985).
compare(ppi2_infomap, ppi2_infomap_copied, method = "adjusted.rand") # 1
compare(ppi2_infomap, ppi2_infomap_copied, method = "nmi") # [1] 1
compare(ppi2_infomap, ppi2_infomap_copied, method = "rand") # [1] 1
compare(ppi2_infomap, ppi2_infomap_copied, method = "vi") # [1] 0
compare(ppi2_infomap, ppi2_infomap_copied, method = "split.join") # [1] 0

# create a modified object for testing purpose

ppi2_infomap_changed = ppi2_infomap
attributes(ppi2_infomap)
# $names
# [1] "membership" "codelength" "names"      "vcount"     "algorithm"  "modularity" --> mean we can do things like ppi2_infomap$codelength
# $class
# [1] "communities" --> means we can do things like communities(ppi2_infomap)

ppi2_infomap
# IGRAPH clustering infomap, groups: 714, mod: 0.44
# + groups:
#   $`1`
# [1] "G5"     "G21"    "G22"    "G35"    "G36"    "G131"   "G147"   "G209"   "G337"   "G346"   "G417"   "G437"   "G489"   "G672"  
# [15] "G682"   "G709"   "G734"   "G735"   "G736"   "G737" .....

#### -----------PLOT NUM OF MODULES VS MODULE SIZE ------------------------------

ppi2_module_sizes = list(cc1 = as.numeric(table(membership(ppi2_infomap))), 
                        cc2 = as.numeric(table(membership(ppi2_louvain))), 
                        cc3 = as.numeric(table(membership(ppi2_walktrap))), 
                        cc4 = as.numeric(table(membership(ppi2_fastgreedy))),
                        cc5 = as.numeric(table(membership(ppi2_labelPropagation)))
                        ) 
par(xpd=T, mar=par()$mar+c(0,3,0,0))
boxplot(ppi2_module_sizes,  horizontal = TRUE, ylim = c(0,50), xaxt = 'n', yaxt = 'n')
axis(1, at = seq(from = 0, to = 50, by = 2), labels = as.character(seq(from = 0, to = 50, by = 2)))
y_labels = c(paste("infomap",as.character(max(membership(ppi2_infomap))), sep = "-"),
             paste("louvain",as.character(max(membership(ppi2_louvain))), sep = "-"),
             paste("walktrap",as.character(max(membership(ppi2_walktrap))), sep = "-"),
             paste("fastgreedy",as.character(max(membership(ppi2_fastgreedy))), sep = "-"),
             paste("labelprop",as.character(max(membership(ppi2_labelPropagation))), sep = "-") 
            )
axis(2, at = 1:length(y_labels), las = 1, labels = y_labels)

par(xpd=F)
# boxplot(as.numeric(table(membership(ppi2_infomap))), horizontal = TRUE, ylab = as.character(max(membership(ppi2_infomap))) )
boxplot(as.numeric(table(membership(ppi2_louvain))), horizontal = TRUE, 
        xaxt = 'n' )
axis(1, at = seq(from = 0, to = 1000, by = 30), labels = as.character(seq(from = 0, to = 1000, by = 30)))
axis(2, at = 1, labels = as.character(max(membership(ppi2_louvain))) )

summary(getModuleVertexSize(igraphObject = pp2_igraph, communityObject = ppi2_infomap))
summary(getModuleVertexSize(igraphObject = pp2_igraph, communityObject = ppi2_louvain))
summary(getModuleVertexSize(igraphObject = pp2_igraph, communityObject = ppi2_walktrap))
summary(getModuleVertexSize(igraphObject = pp2_igraph, communityObject = ppi2_fastgreedy))
summary(getModuleVertexSize(igraphObject = pp2_igraph, communityObject = ppi2_labelPropagation))

ppi2_stat_alg_matrix

#### ---------------- COMPARE RESULT PPI2  -------------------------------------

ppi2_stat_alg_matrix

#### infomap
ppi2_subgraphs_infomap = getAllSubgraphs(igraphObject = pp2_igraph, communityObject = ppi2_infomap)
ppi2_stat_infomap = getSubgraphStat(list_of_subgraphs = ppi2_subgraphs_infomap, minVertexCount = 3)
ppi2_infomap_size = getModuleVertexSize(igraphObject = pp2_igraph, communityObject = ppi2_infomap)
summary(getModuleVertexSize(igraphObject = pp2_igraph, communityObject = ppi2_infomap))

# we can either use ppi2_stat_infomap or ppi2_subgraphs_infomap to compute the weighted module
#ppi2_stat_module_infomap = getCommunityWeightedStat(stat_dataframe = ppi2_stat_infomap, communityObject = ppi2_infomap)
ppi2_stat_module_infomap = getCommunityWeightedStat(list_of_subgraphs = ppi2_subgraphs_infomap, 
                                                    communityObject = ppi2_infomap)

#### fast-greedy
ppi2_subgraphs_fastgreedy = getAllSubgraphs(igraphObject = pp2_igraph, communityObject = ppi2_fastgreedy)
ppi2_stat_module_fastgreedy = getCommunityWeightedStat(list_of_subgraphs = ppi2_subgraphs_fastgreedy, 
                                                       communityObject = ppi2_fastgreedy)

#### walk_trap
ppi2_subgraphs_walktrap = getAllSubgraphs(igraphObject = pp2_igraph, communityObject = ppi2_walktrap)
ppi2_stat_module_walktrap = getCommunityWeightedStat(list_of_subgraphs = ppi2_subgraphs_walktrap,
                                                     communityObject = ppi2_walktrap)

### louvain
ppi2_subgraphs_louvain = getAllSubgraphs(igraphObject = pp2_igraph, communityObject = ppi2_louvain)
ppi2_stat_module_louvain = getCommunityWeightedStat(list_of_subgraphs = ppi2_subgraphs_louvain,
                                                    communityObject = ppi2_louvain)

### multilevel --> similar result as louvain
# ppi2_subgraphs_multilevel = getAllSubgraphs(igraphObject = pp2_igraph, communityObject = ppi2_multilevel)
# ppi2_stat_module_multilevel = getCommunityWeightedStat(list_of_subgraphs = ppi2_subgraphs_multilevel)

### label propagation
ppi2_subgraphs_labelPropagation = getAllSubgraphs(igraphObject = pp2_igraph, communityObject = ppi2_labelPropagation)
ppi2_stat_module_labelPropagation = getCommunityWeightedStat(list_of_subgraphs = ppi2_subgraphs_labelPropagation,
                                                    communityObject = ppi2_labelPropagation)

### combine result into a single matrix
ppi2_stat_alg_matrix = rbind(ppi2_stat_module_infomap, ppi2_stat_module_fastgreedy, 
                             ppi2_stat_module_walktrap, ppi2_stat_module_louvain,
                             ppi2_stat_module_labelPropagation)
rownames(ppi2_stat_alg_matrix) = c("infomap", "fastgreedy", "walktrap", 
                                   "louvain", "label_propagation")
colnames(ppi2_stat_alg_matrix) = c("density", "clusCoef", "avg_simmiliarity", "modularity")

# save result
save(ppi2_stat_alg_matrix, file = "ppi2_compare_summary.rda")

#### -------------- CREATE COMPARE MATRICES ------------------------------------

ppi2_stat_alg_matrix
#                     density   clusCoef   avg_simmiliarity modularity
# infomap           0.25319805 0.3908631       0.13978319  0.4370835
# fastgreedy        0.02784615 0.4827973       0.01298743  0.4127827
# walktrap          0.13368894 0.4310397       0.06726763  0.3913991
# louvain           0.03702952 0.5355259       0.01839634  0.4961280
# label_propagation 0.03067184 0.3716713       0.01497208  0.1716730

mclust::adjustedRandIndex(ppi2_infomap$membership, ppi2_fastgreedy$membership) # [1] 0.03747002
mclust::adjustedRandIndex(ppi2_infomap$membership, ppi2_walktrap$membership) # [1] 0.05862351
mclust::adjustedRandIndex(ppi2_infomap$membership, ppi2_louvain$membership) # [1] 0.1197415
mclust::adjustedRandIndex(ppi2_fastgreedy$membership, ppi2_walktrap$membership) # [1] 0.1142185
mclust::adjustedRandIndex(ppi2_fastgreedy$membership, ppi2_louvain$membership) # [1] 0.1405325
mclust::adjustedRandIndex(ppi2_walktrap$membership, ppi2_louvain$membership) # [1] 0.1075633

#### --------------- RETRIEVE NUMBER OF FRONTIER EDGES -------------------------

subgraph1 = ppi2_subgraphs_infomap[[1]]
vcount(subgraph1) # 378
ecount(subgraph1) # 22468

V(subgraph1)$name [1:10]
#[1] "G5"   "G21"  "G22"  "G35"  "G36"  "G131" "G147" "G209" "G337" "G346"

E(subgraph1)
# + 22648/22648 edges (vertex names):
# [1] G5--G35   G5--G209  G5--G346  G5--G709  G5--G1964 G5--G2177 G5--G2181 G5--G2275 G5--G2422

intra_edges = E(subgraph1) [ from(V(subgraph1)[1]) ]
# intra_edges
# + 89/22648 edges (vertex names):
#   [1] G5--G35   G5--G209  G5--G346  G5--G709  G5--G1964 G5--G2177 G5--G2181 G5--G2275 G5--G2422
idx = match(V(subgraph1)[1]$name, V(pp2_igraph)$name)
all_edges = E(pp2_igraph) [ from(idx) ]

num_frontier_edges = length(all_edges) - length(intra_edges)

list_of_subgraphs = ppi2_subgraphs_infomap
length(list_of_subgraphs) # 714
min_vertex_count = 94

size = getSubgraphSize(list_of_subgraphs)

selected_subg_indices = c()
for (i in 1:length(list_of_subgraphs)){
  subg = list_of_subgraphs[[i]]
  if (vcount(subg) >= min_vertex_count){
    selected_subg_indices = append(selected_subg_indices, i)
  }
}

list_of_subgraphs[selected_subg_indices];

sum_frontier_edges = c()
# loop through each community
for (i in 1:length(list_of_subgraphs)){
  subg = list_of_subgraphs[[i]]
  # loop through each vertex in the community
  total_num_frontier_edges = 0
  for (j in 1:vcount(subg)) { 
    intra_edges = E(subg) [ from(V(subg)[j]) ]
    idx = match(V(subg)[1]$name, V(subg)$name)
    all_edges = E(pp2_igraph) [ from(idx) ]
    num_frontier_edges = length(all_edges) - length(intra_edges)
    total_num_frontier_edges = total_num_frontier_edges + num_frontier_edges
  }
  sum_frontier_edges = append(sum_frontier_edges,total_num_frontier_edges)
}

#### ---------------- SHORT TUTORIAL -------------------------------------------

g <- make_graph("Zachary")
sg <- cluster_spinglass(g)
le <- cluster_leading_eigen(g)
compare(sg, le, method="rand")
compare(membership(sg), membership(le))




