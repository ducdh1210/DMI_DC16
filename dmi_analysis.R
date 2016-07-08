# References:
# http://www.r-bloggers.com/summary-of-community-detection-algorithms-in-igraph-0-6/
# http://stackoverflow.com/questions/9471906/what-are-the-differences-between-community-detection-algorithms-in-igraph


#### ---------------- load library -------------------------------------------------
library('gplots')
library('ggplot2')
library('WGCNA')
library('igraph')
library('clue')

#### ---------------- load saved image --------------------------------------------
setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
rm(list = ls()); load("PPI2.RData")

# #### ---------------- import file -------------------------------------------------
# setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
# rm(list = ls())
# download_path = "/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge"
# downloaded_folder <- "DREAM11_subchallenge1_initial_networks"
# file_names = paste(download_path,"/",downloaded_folder,"/",list.files(downloaded_folder)
#                   ,sep = "")
# ppi_2 = read.table(file = file_names[2], header = TRUE, sep = "\t")
# 
# #### ---------------- PPI_2 -------------------------------------------------
# # to free up memory, remove all variables that are not ppi_2
# rm(list = ls()[-which(ls() == "ppi_2" )])
# temp_pp2 = ppi_2
# colnames(temp_pp2) = c("V1","V2","weight"); head(temp_pp2)
# temp_pp2$V1 = paste("G",temp_pp2$V1,sep = "")
# temp_pp2$V2 = paste("G",temp_pp2$V2,sep = "")
# pp2_igraph = graph.data.frame(d = temp_pp2, directed = F);
# rm(temp_pp2)

  
# temp_pp2[] = lapply(temp_pp2, as.character)
# length(base::unique(df_ppi2$X0)) # [1] 9650
# length(base::unique(df_ppi2$X1)) # [1] 12093   
# length(intersect(base::unique(df_ppi2$X0), base::unique(df_ppi2$X1))) # 9323
# length(union(base::unique(df_ppi2$X0), base::unique(df_ppi2$X1))) # 12420

pp2_edgeList = get.edgelist(pp2_igraph)

# get adjacency matrix of type dgCMatrix
pp2_adjMatrix = get.adjacency(pp2_igraph, attr = "weight")
pp2_adjMatrix@Dim # [1] 12420 12420

# convert dgCMatrix to normal matrix object
pp2_adjM = as.matrix(pp2_adjMatrix)

# What is maximum clique size
clique.number(pp2_igraph) # too slow to run, have to stop before completion

#### ---------------- WCGNA -------------------------------------------

# Cluster gene expression profiles; the flashClust function from
# the authors of WGCNA is another options for larger datasets.
# For input, we use the reciprocal of the adjacency matrix; hierarchical
# clustering works by comparing the _distance_ between objects instead of the
# _similarity_.
gene_tree = hclust(as.dist(1 - pp2_adjM), method="average")
gc()
# we will use the cuttreeDynamicTree method to break apart the hc dendrogram
# into separate modules
module_labels = cutreeDynamicTree(dendro=gene_tree, minModuleSize=15,
                                   deepSplit=TRUE)
gc()
# assign a color to each module for easier visualization and referencing
module_colors = labels2colors(module_labels);

nlevels(as.factor(module_colors)) # [1] 339 ==> thus 338 modules 

V(pp2_igraph)$color = col2hex(module_colors)

file_name = "/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/subchallenge1_ppi2.graphml"
write.graph(graph = pp2_igraph, file = file_name, format='graphml') # able to display in Cytoscape

#### ---------------- MCODE ---------------------------------------------------------------
library(ProNet)

result = mcode(pp2_igraph,vwp=0.05,haircut=TRUE,fluff=FALSE,fdt=0.8,loops=FALSE)
mcode_result = result; 
rm(result)

#### ---------------- FAST GREEDY ---------------------------------------------------------

library(igraph)
# This function tries to find dense subgraph, also called communities in graphs 
# via directly optimizing a modularity score.
ppi2_fastgreedy = fastgreedy.community(pp2_igraph, weights = E(pp2_igraph)$weight)
attributes(ppi2_fastgreedy)
modularity(ppi2_fastgreedy) # [1] 0.4127827
communities(ppi2_fastgreedy)

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
table(membership(ppi2_fastgreedy)) # 125 modules

module_matrix_fastgreedy = getModuleMatrix(igraph_object = pp2_igraph, 
                                   community_membership = membership(ppi2_fastgreedy))

ppi2_fastgreedy_module_density = getModuleDensity(igraphObject = pp2_igraph, communityObject = ppi2_fastgreedy)

getModuleDensity = function(igraphObject = NULL, communityObject = NULL){
  module_density = sapply(sort(unique(membership(communityObject))), function(assigned_cluster) {
    subg <-induced.subgraph(igraphObject, which(membership(communityObject)==assigned_cluster)) #membership id differs for each cluster
    #ecount(subg)/ecount(g)
    graph.density(subg)
  })
  return(module_density)
}

ppi2_fastgreedy_module_density = sapply(sort(unique(membership(ppi2_fastgreedy))), function(assigned_cluster) {
  subg <-induced.subgraph(pp2_igraph, which(membership(ppi2_fastgreedy)==assigned_cluster)) #membership id differs for each cluster
  #ecount(subg)/ecount(g)
  graph.density(subg)
})

getModuleClusterCoef = function(igraphObject = NULL, communityObject = NULL){
  module_clusterCoeff = sapply(sort(unique(membership(communityObject))), function(assigned_cluster) {
    subg <-induced.subgraph(igraphObject, which(membership(communityObject)==assigned_cluster)) #membership id differs for each cluster
    #ecount(subg)/ecount(g)
    transitivity(subg, type = "global")
  })
  return(module_clusterCoeff)
}

mm_module_clusCoef = sapply(sort(unique(membership(mm_zachary))), function(assigned_cluster) {
  subg <-induced.subgraph(g_zachary, which(membership(mm)==assigned_cluster)) #membership id differs for each cluster
  #ecount(subg)/ecount(g)
  transitivity(subg, type = "global")
})

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



#### ---------------- UTILITY FUNCTIONS -------------------------------------------

getModuleMatrix = function(igraph_object, community_membership){
  # get total number of vertices
  size = length(V(igraph_object)$name)
  # initialize the square module matrix to be returned
  module_matrix = matrix(FALSE, nrow = size, ncol = size)
  diag(module_matrix) = TRUE
  rownames(module_matrix) = V(igraph_object)$name
  colnames(module_matrix) = rownames(module_matrix)
  # foreach vertice (row), get the module it belongs to, and search for other vertices
  # belonging to the same module, then set the cell = TRUE in such cases
  for (row_index in rownames(module_matrix)){
    assigned_module = community_membership[row_index]
    vertices_in_same_module =  names(which(community_membership == assigned_module))
    module_matrix[row_index, vertices_in_same_module] = TRUE
  }
  return(module_matrix)
}


