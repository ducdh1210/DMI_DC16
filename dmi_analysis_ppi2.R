# References:
# http://www.r-bloggers.com/summary-of-community-detection-algorithms-in-igraph-0-6/
# http://stackoverflow.com/questions/9471906/what-are-the-differences-between-community-detection-algorithms-in-igraph
# https://www.synapse.org/#!Synapse:syn6156761/wiki/400652

#### ---------------- load library -------------------------------------------------
library('gplots'); library('ggplot2'); library('WGCNA'); library('igraph')
library('clue'); library('ProNet')
#### ---------------- load saved image --------------------------------------------
setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
source("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/graph_utility_functions.R")
rm(list = ls());  gc();

#load("PPI2.RData")
load("ppi2_igraph.rda")
load("ppi2_communities.rda")
load("ppi2_clique_objects.rda")


#save(ppi2_fastgreedy, ppi2_walktrap, ppi2_labelPropagation, ppi2_infomap, ppi2_louvain, ppi2_multilevel, file = "ppi2_communities.rda")
# #### ---------------- import file -------------------------------------------------
setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
rm(list = ls())
download_path = "/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge"
downloaded_folder <- "DREAM11_subchallenge1_initial_networks"
file_names = paste(download_path,"/",downloaded_folder,"/",list.files(downloaded_folder)
                   ,sep = "")
ppi_2 = read.table(file = file_names[2], header = TRUE, sep = "\t")


# #### ---------------- PPI_2 -------------------------------------------------
# to free up memory, remove all variables that are not ppi_2
rm(list = ls()[-which(ls() == "ppi_2" )])
temp_pp2 = ppi_2
colnames(temp_pp2) = c("V1","V2","weight"); head(temp_pp2)
temp_pp2$V1 = paste("G",temp_pp2$V1,sep = "")
temp_pp2$V2 = paste("G",temp_pp2$V2,sep = "")
pp2_igraph = graph.data.frame(d = temp_pp2, directed = F);
#pp2_igraph = graph.data.frame(d = temp_pp2, directed = T);
rm(temp_pp2)

#### ---------------- WCGNA -------------------------------------------

pp2_edgeList = get.edgelist(pp2_igraph)

# get adjacency matrix of type dgCMatrix
pp2_adjMatrix = get.adjacency(pp2_igraph, attr = "weight")
pp2_adjMatrix@Dim # [1] 12420 12420

# convert dgCMatrix to normal matrix object
pp2_adjM = as.matrix(pp2_adjMatrix)


# Cluster gene expression profiles; the flashClust function from
# the authors of WGCNA is another options for larger datasets.
# For input, we use the reciprocal of the adjacency matrix; hierarchical
# clustering works by comparing the _distance_ between objects instead of the
# _similarity_.
ppi2_gene_tree = hclust(as.dist(1 - pp2_adjM), method="average")


gc()
# we will use the cuttreeDynamicTree method to break apart the hc dendrogram
# into separate modules
ppi2_wgcna_modules = cutreeDynamicTree(dendro=ppi2_gene_tree, minModuleSize=3,
                                       deepSplit=TRUE)
gc()
# assign a color to each module for easier visualization and referencing
module_colors = labels2colors(ppi2_wgcna_modules);

nlevels(as.factor(module_colors)) # [1] 339 ==> thus 338 modules 

V(pp2_igraph)$color = col2hex(module_colors)

file_name = "/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/subchallenge1_ppi2.graphml"
write.graph(graph = pp2_igraph, file = file_name, format='graphml') # able to display in Cytoscape

#### ---------------- MCODE ---------------------------------------------------------------
library(ProNet)

mcode_result = mcode(pp2_igraph,vwp=0.05,haircut=TRUE,fluff=FALSE,fdt=0.8,loops=FALSE)
rm(result)

#### ---------------- FAST GREEDY ---------------------------------------------------------

library(igraph)
# This function tries to find dense subgraph, also called communities in graphs 
# via directly optimizing a modularity score.
ppi2_fastgreedy = fastgreedy.community(pp2_igraph, weights = E(pp2_igraph)$weight)
table(membership(ppi2_fastgreedy)) # 125 modules

attributes(ppi2_fastgreedy); communities(ppi2_fastgreedy)

ppi2_fastgreedy_den = as.dendrogram(ppi2_fastgreedy)

# check if fastgreedy.community works EXACTLY the same as cluster_fast_greedy
ppi2_fastgreedy2 = cluster_fast_greedy(graph = pp2_igraph, weights = E(pp2_igraph)$weight)
sum(ppi2_fastgreedy$membership - ppi2_fastgreedy2$membership) # [1] 0 --> confirmed that 
# this alg works exactly the same


# Possible values: ‘vi’ is the variation of information (VI) metric of Meila (2003), 
# ‘nmi’ is the normalized mutual information measure proposed by Danon et al. (2005), 
# ‘split.join’ is the split-join distance of can Dongen (2000), 
# ‘rand’ is the Rand index of Rand (1971),
# ‘adjusted.rand’ is the adjusted Rand index by Hubert and Arabie (1985).
compare(ppi2_fastgreedy, ppi2_fastgreedy2, method = "adjusted.rand") # 1
compare(ppi2_fastgreedy, ppi2_fastgreedy2, method = "nmi") # 1
compare(ppi2_fastgreedy, ppi2_fastgreedy2, method = "vi") # 0
compare(ppi2_fastgreedy, ppi2_fastgreedy2, method = "split.join") # 0

class(ppi2_fastgreedy_den)
plot_dendrogram(ppi2_fastgreedy_den)

plot(rev(ppi2_fastgreedy$modularity), xlab = 'Number of clusters', ylab = 'Modularity value')
which.max(rev(ppi2_fastgreedy$modularity)) # 83

module_matrix_fastgreedy = getModuleMatrix(igraph_object = pp2_igraph, 
                                           community_membership = membership(ppi2_fastgreedy))

### Test out different evaluation metrics

modularity(x = ppi2_fastgreedy) # [1] 0.4127827
modularity(x = ppi2_fastgreedy, membership = ppi2_fastgreedy2$membership, weights = V(pp2_igraph)$weight) # [1] 0.4127827
mod.matrix(graph = pp2_igraph,  membership = ppi2_fastgreedy2$membership, weights = V(pp2_igraph)$weight)
mod.matrix(ppi2_fastgreedy)
ppi2_fastgreedy_module_density = getModuleDensity(igraphObject = pp2_igraph, communityObject = ppi2_fastgreedy)
ppi2_fastgreedy_module_clusterCoeff = getModuleClusterCoef(igraphObject = pp2_igraph, communityObject = ppi2_fastgreedy)
ppi2_fastgreedy_simmiliarty = getModuleSimmilarity(igraphObject = pp2_igraph, communityObject = ppi2_fastgreedy, method = "jaccard")

#### ---------------- WALK TRAP ----------------------------------------------------------
library(igraph)
# This function tries to find densely connected subgraphs, also called communities
# in a graph via random walks. The idea is that short random walks tend to stay 
# in the same community.
# http://arxiv.org/abs/physics/0512106
ppi2_walktrap = walktrap.community(pp2_igraph, weights = E(pp2_igraph)$weight)
attributes(ppi2_walktrap)
membership(ppi2_walktrap) 
which.max(rev(ppi2_walktrap$modularity)) # 2364
table(membership(ppi2_walktrap)) # 2364 modules

#### ---------------- LEADING EIGENVECTOR -------------------------------------------
library(igraph)
# This function tries to find densely connected subgraphs in a graph by calculating the 
# leading non-negative eigenvector of the modularity matrix of the graph.
# https://arxiv.org/abs/physics/0605087
ppi2_leadingEigenvector = leading.eigenvector.community(pp2_igraph, weights = E(pp2_igraph)$weight)
module_matrix_leadingEigenvector = getModuleMatrix(igraph_object = pp2_igraph, 
                                                   community_membership = membership(ppi2_leadingEigenvector))

#### ---------------- LABEL PROPODAGATION -------------------------------------------
library(igraph)
# This is a fast, nearly linear time algorithm for detecting community structure in networks. 
# In works by labeling the vertices with unique labels and then updating the labels 
# by majority voting in the neighborhood of the vertex.
# https://arxiv.org/pdf/0709.2938.pdf
ppi2_labelPropagation = label.propagation.community(pp2_igraph, weights = E(pp2_igraph)$weight)
attributes(ppi2_labelPropagation)
membership(ppi2_labelPropagation) 
table(membership(ppi2_labelPropagation)) # 92 modules

#### ---------------- INFOMAP ---------------------------------------------------------
library(igraph)
# Find community structure that minimizes the expected description length
# of a random walker trajectory
# http://arxiv.org/abs/0707.0609
ppi2_infomap = infomap.community(pp2_igraph, e.weights = E(pp2_igraph)$weight)
membership(ppi2_infomap) 
table(membership(ppi2_infomap)) # 711 modules

#### ---------------- MULTILEVEL ---------------------------------------------------------
library(igraph)
# This function implements the multi-level modularity optimization algorithm for 
# finding community structure, see references below. It is based on the modularity
# measure and a hierarchial approach.
# https://arxiv.org/abs/0803.0476
ppi2_multilevel = multilevel.community(pp2_igraph, weights = E(pp2_igraph)$weight)
membership(ppi2_multilevel) 
table(membership(ppi2_multilevel)) # 61 modules

#### ---------------- OPTIMAL  ---------------------------------------------------------
library(igraph)
# This function implements the multi-level modularity optimization algorithm for 
# finding community structure, see references below. It is based on the modularity
# measure and a hierarchial approach.
# https://arxiv.org/abs/0803.0476
ppi2_optimal = optimal.community(pp2_igraph, weights = E(pp2_igraph)$weight) # overflow ==> CRASHED!
ppi2_optimal = optimal.community(pp2_igraph, weights = NULL) # also CRASHED!

#### ---------------- LUOVAIN  ---------------------------------------------------------
library(igraph)
# This function implements the multi-level modularity optimization algorithm for finding 
# community structure, see VD Blondel, J-L Guillaume, R Lambiotte and E Lefebvre: 
# Fast unfolding of community hierarchies in large networks, http://arxiv.org/abs/arXiv:0803.0476 for the details.

# It is based on the modularity measure and a hierarchial approach. 
# Initially, each vertex is assigned to a community on its own. 
# In every step, vertices are re-assigned to communities in a local, greedy way: 
# each vertex is moved to the community with which it achieves the highest contribution to modularity. 
# When no vertices can be reassigned, each community is considered a vertex on its own, 
# and the process starts again with the merged communities. 
# The process stops when there is only a single vertex left or 
# when the modularity cannot be increased any more in a step.

# https://arxiv.org/abs/0803.0476
ppi2_louvain = cluster_louvain(pp2_igraph, weights = E(pp2_igraph)$weight) 
class(ppi2_louvain)# [1] "communities"
length(unique(ppi2_louvain$membership)) # return 61 modules
table(ppi2_louvain$membership)

plot(rev(ppi2_louvain$modularity), xlab = 'Number of clusters', ylab = 'Modularity value')
which.max(rev(ppi2_louvain$modularity)) # 83

test = as.dendrogram(ppi2_louvain)

#### ---------------- MAXIMAL CLIQUE -------------------------------------------


pp2_igraph_contracted <- contract(graph = pp2_igraph, membership_clique, vertex.attr.comb=toString)

# original igraph object
vcount(pp2_igraph); ecount(pp2_igraph)
# [1] 12420
# [1] 397308

# cliquedified igraph object
vcount(pp2_igraph_contracted); ecount(pp2_igraph_contracted)
# [1] 10910
# [1] 397308

# look more into: http://blog.revolutionanalytics.com/2015/08/contracting-and-simplifying-a-network-graph.html
pp2_igraph_contracted_simplified = simplify(pp2_igraph_contracted); 
vcount(pp2_igraph_contracted_simplified); ecount(pp2_igraph_contracted_simplified)
# [1] 10910
# [1] 194841
is_simple(pp2_igraph_contracted_simplified) # TRUE


#### ---------------- MAXIMAL CLIQUE (2)----------------------------------------

# obtain list of igraph.vs objects
ppi2_clique_list = getMaximalClique(igraphObject = pp2_igraph,
                                    minCliqueSize = 5, 
                                    degreeThreshold = 10)
# NOTE that in this list, different cliques can share similar vertices 

# obtain vector of clique membership
ppi2_clique_membership = getCliqMemebership(igraphObject = pp2_igraph,
                                            cliq = ppi2_clique_list,
                                            isCliqOrdered = F)

save(ppi2_clique_list, ppi2_cliq_membership, file = "ppi2_clique_objects.rda")
load("ppi2_clique_objects.rda")

length(ppi2_clique_list) # 58424
length(ppi2_clique_list[[1]])
clique_sizes = sapply(ppi2_clique_list, function(clique) length(clique))
summary(clique_sizes)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.0     8.0    18.0    19.8    28.0    73.0 

# some counting statistic
max(ppi2_clique_membership) # 10910 --> number of predicted modules
table(table(ppi2_clique_membership))
# 1     5     6     7     8     9    10    11    12    13    14    15    16    17    19    21 
# 10738    63    25    18    16     6     7     8     3     2     3     1     2     1     1     1 
# 23    25    26    28    31    32    41    44    45    46    55    73 
# 2     2     2     1     1     1     1     1     1     1     1     1 

# meaning: 10738 cliques only has 1 elements --> trivial

modularity(x = pp2_igraph, 
           membership = ppi2_clique_membership, 
           weights = E(pp2_igraph)$weight) # [1] 0.05431305

ppi2_clique_subgraph_list = getAllSubgraphs(igraphObject = pp2_igraph,
                                            clique = ppi2_clique_list, 
                                            cliqueMembership = ppi2_clique_membership)

getCommunityWeightedStat(list_of_subgraphs = ppi2_clique_subgraph_list)

pp2_contracted = contract(graph = pp2_igraph, ppi2_clique_membership, vertex.attr.comb=toString)

pp2_clique_simplified1 = simplify(graph = pp2_contracted, 
                                  remove.multiple = TRUE, 
                                  remove.loops = FALSE)

pp2_clique_simplified1_louvain = cluster_louvain(pp2_igraph_simplified1, 
                                                 weights = E(pp2_igraph_simplified1)$weight) 

max(membership(pp2_clique_simplified1_louvain)) # 61

table(table(membership(pp2_clique_simplified1_louvain)))

modularity(x = pp2_clique_simplified1, membership = pp2_clique_simplified1_louvain)

pp2_clique_simplified2 = simplify(graph = pp2_contracted, 
                                  remove.multiple = TRUE, 
                                  remove.loops = TRUE)

pp2_clique_simplified2_louvain = cluster_louvain(pp2_clique_simplified2, 
                                                 weights = E(pp2_clique_simplified2)$weight) 

max(membership(pp2_clique_simplified2_louvain)) # 63
table(table(membership(pp2_clique_simplified2_louvain)))



##### ------- How to deal with size 1 clique -----------------------------------

# Strategy 1: 
#   1. combine 10738 nodes from those cliques into a graph 
#   2. apply community detection algorithms on that graph
#   3. integrate the newly found modules with the clique
# 
# Strategy 2:
#   1. For each nodes belonging to one of 10738 cliques, find which other (bigger)
#      cliques it can be moved into. Question: How to decide which cliques should 
#      it moves into?

# Strategy 3: Combine Louvain result with clique
# General questions: How to
# 1. Merge smaller modules
# 2. Divide up big modules
# 3. Order of operation for (1) and (2): sequential or simultanous?


# get indices of clique whose size = 1
ppi2_cliq_size_one = names(table(ppi2_clique_membership)[which(table(ppi2_clique_membership) == 1)])
ppi2_alone_vertices = names(ppi2_clique_membership[which(ppi2_clique_membership %in% ppi2_cliq_size_one)])
ppi2_combined_alone_vertices_graph = induced.subgraph(graph = pp2_igraph, vids = ppi2_alone_vertices)
graph.density(ppi2_combined_alone_vertices_graph) # 0.001914204, it's lower than 
# pp2_graph's density, which is 0.005151682
transitivity(ppi2_combined_alone_vertices_graph, type = "global") # 0.4160488, actually higher than 
# pp2_graph's clustering coefficient, which is 0.3625158

test = infomap.community(graph = ppi2_combined_alone_vertices_graph, 
                         e.weights = E(ppi2_combined_alone_vertices_graph)$weight)

modularity(test) # 0.5332497
modularity(x = ppi2_combined_alone_vertices_graph, membership = membership(test)) # 0.3974415
max(membership(test))
test2 = getAllSubgraphs(igraphObject = ppi2_combined_alone_vertices_graph,communityObject = test) 
getCommunityWeightedStat(list_of_subgraphs = test2)

















