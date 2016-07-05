# References:
# http://www.r-bloggers.com/summary-of-community-detection-algorithms-in-igraph-0-6/
# http://stackoverflow.com/questions/9471906/what-are-the-differences-between-community-detection-algorithms-in-igraph

setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")

library('gplots')
library('ggplot2')
library('knitr')
library('limma')
library('reshape2')
library('RColorBrewer')
library('WGCNA')
library(igraph)
library(clue)

rm(list = ls())
load(".RData")
# class(ppi_1); head(ppi_1)
# 
# working_ppi1 = ppi_1
# colnames(working_ppi1) = c("V1","V2","weight"); head(working_ppi1)
# 
# # create igraph object 
# pp1_graphObj = graph.data.frame(working_ppi1, directed = F);
# class(pp1_graphObj)
# 
# E(pp1_graphObj) #2232404/2232404 edges 
# V(pp1_graphObj) #17397/17397 vertices
# head(E(pp1_graphObj)$weight)
# 
# # plot the graph
# plot.igraph(pp1_graphObj) # takes long time
# groups <- membership(cluster_louvain(pp1_graphObj))
# communities <- communities(cluster_louvain(pp1_graphObj))
# 
# # create file of type graphml to be exported to Cytoscape
# file_name = "/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/subchallenge1_ppi1.graphml"
# write.graph(graph = pp1_graphObj, file = file_name, format='graphml') # could not export
# 
# library(RCy3)
# 
# pp1_graphNEL = igraph.to.graphNEL(pp1_graphObj)
# cw <- CytoscapeWindow ('vignette', graph=pp1_graphObj, overwrite=TRUE)
# displayGraph (cw)

#### ---------------- PPI_2 -------------------------------------------------
# to free up memory, remove all variables that are not ppi_2
rm(list = ls()[-which(ls() == "ppi_2" )])
temp_pp2 = ppi_2
colnames(temp_pp2) = c("V1","V2","weight"); head(temp_pp2)
temp_pp2$V1 = paste("G",temp_pp2$V1,sep = "")
temp_pp2$V2 = paste("G",temp_pp2$V2,sep = "")

pp2_igraph = graph.data.frame(d = temp_pp2, directed = F);

  
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

#### ---------------- WCGNA -------------------------------------------

# Cluster gene expression profiles; the flashClust function from
# the authors of WGCNA is another options for larger datasets.
# For input, we use the reciprocal of the adjacency matrix; hierarchical
# clustering works by comparing the _distance_ between objects instead of the
# _similarity_.
gene_tree <- hclust(as.dist(1 - pp2_adjM), method="average")
gc()
# we will use the cuttreeDynamicTree method to break apart the hc dendrogram
# into separate modules
module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=15,
                                   deepSplit=TRUE)
gc()
# assign a color to each module for easier visualization and referencing
module_colors <- labels2colors(module_labels);

nlevels(as.factor(module_colors)) # [1] 339 ==> thus 338 modules 

V(pp2_igraph)$color = col2hex(module_colors)

file_name = "/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/subchallenge1_ppi2.graphml"
write.graph(graph = pp2_igraph, file = file_name, format='graphml') # able to display in Cytoscape

#### ---------------- MCODE ---------------------------------------------------------------
library(ProNet)

result <- mcode(pp2_igraph,vwp=0.05,haircut=TRUE,fluff=FALSE,fdt=0.8,loops=FALSE)
mcode_result = result; 
rm(result)

#### ---------------- FAST GREEDY ---------------------------------------------------------
library(igraph)
# This function tries to find dense subgraph, also called communities in graphs 
# via directly optimizing a modularity score.
ppi2_fastgreedy <- fastgreedy.community(pp2_igraph, weights = E(pp2_igraph)$weight)
attributes(ppi2_fastgreedy)
modularity(ppi2_fastgreedy) # [1] 0.4127827
communities(ppi2_fastgreedy)

table(membership(ppi2_fastgreedy)) # 125 modules

module_matrix_fastgreedy = getModuleMatrix(igraph_object = pp2_igraph, 
                                   community_membership = membership(ppi2_fastgreedy))

#### ---------------- WALK TRAP ----------------------------------------------------------
library(igraph)
# This function tries to find densely connected subgraphs, also called communities
# in a graph via random walks. The idea is that short random walks tend to stay 
# in the same community.
# http://arxiv.org/abs/physics/0512106
ppi2_walktrap <- walktrap.community(pp2_igraph, weights = E(pp2_igraph)$weight)
attributes(ppi2_walktrap)
membership(ppi2_walktrap) 
table(membership(ppi2_walktrap)) # 2364 modules

#### ---------------- LEADING EIGENVECTOR -------------------------------------------
library(igraph)
# This function tries to find densely connected subgraphs in a graph by calculating the 
# leading non-negative eigenvector of the modularity matrix of the graph.
# https://arxiv.org/abs/physics/0605087
ppi2_leadingEigenvector <- leading.eigenvector.community(pp2_igraph, weights = E(pp2_igraph)$weight)
module_matrix_leadingEigenvector = getModuleMatrix(igraph_object = pp2_igraph, 
                                           community_membership = membership(ppi2_leadingEigenvector))

#### ---------------- LABEL PROPODAGATION -------------------------------------------
library(igraph)
# This is a fast, nearly linear time algorithm for detecting community structure in networks. 
# In works by labeling the vertices with unique labels and then updating the labels 
# by majority voting in the neighborhood of the vertex.
# https://arxiv.org/pdf/0709.2938.pdf
ppi2_labelPropagation <- label.propagation.community(pp2_igraph, weights = E(pp2_igraph)$weight)
attributes(ppi2_labelPropagation)
membership(ppi2_labelPropagation) 
table(membership(ppi2_labelPropagation)) # 92 modules

#### ---------------- INFOMAP ---------------------------------------------------------
library(igraph)
# Find community structure that minimizes the expected description length
# of a random walker trajectory
# http://arxiv.org/abs/0707.0609
ppi2_infomap <- infomap.community(pp2_igraph, e.weights = E(pp2_igraph)$weight)
membership(ppi2_infomap) 
table(membership(ppi2_infomap)) # 711 modules

#### ---------------- MULTILEVEL ---------------------------------------------------------
library(igraph)
# This function implements the multi-level modularity optimization algorithm for 
# finding community structure, see references below. It is based on the modularity
# measure and a hierarchial approach.
# https://arxiv.org/abs/0803.0476
ppi2_multilevel <- multilevel.community(pp2_igraph, weights = E(pp2_igraph)$weight)
membership(ppi2_infomap) 
table(membership(ppi2_infomap)) # 711 modules

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


