setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")

library('gplots')
library('ggplot2')
library('knitr')
library('limma')
library('reshape2')
library('RColorBrewer')
library('WGCNA')
library(igraph)

load(".RData")
class(ppi_1); head(ppi_1)

working_ppi1 = ppi_1
colnames(working_ppi1) = c("V1","V2","weight"); head(working_ppi1)

# create igraph object 
pp1_graphObj = graph.data.frame(working_ppi1, directed = F);
class(pp1_graphObj)

E(pp1_graphObj) #2232404/2232404 edges 
V(pp1_graphObj) #17397/17397 vertices
head(E(pp1_graphObj)$weight)

# plot the graph
plot.igraph(pp1_graphObj) # takes long time
groups <- membership(cluster_louvain(pp1_graphObj))
communities <- communities(cluster_louvain(pp1_graphObj))

# create file of type graphml to be exported to Cytoscape
file_name = "/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/subchallenge1_ppi1.graphml"
write.graph(graph = pp1_graphObj, file = file_name, format='graphml') # could not export

library(RCy3)

pp1_graphNEL = igraph.to.graphNEL(pp1_graphObj)
cw <- CytoscapeWindow ('vignette', graph=pp1_graphObj, overwrite=TRUE)
displayGraph (cw)

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

# http://stackoverflow.com/questions/14849835/how-to-calculate-adjacency-matrices-in-r

# require(reshape2)
# m = as.matrix(dcast(temp_pp2, V1 ~ V2, value.var = "weight", fill = 0 ))
# require(intergraph)


#### -----------------------------------------------------------------
class(signal_directed)
working_signal_directed = signal_directed
colnames(working_signal_directed) = c("V1","V2","weight"); head(working_signal_directed)
signal_directed_graphObj = graph.data.frame(working_signal_directed, directed = T)

E(signal_directed_graphObj) #21825
V(signal_directed_graphObj) #5254
head(E(signal_directed_graphObj)$weight)

plot.igraph(signal_directed_graphObj) # takes long time, runable!

file_name = "/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/subchallenge1_signalDirected.graphml"
write.graph(graph = signal_directed_graphObj, file = file_name, format='graphml') 


signalDirected_graphNEL = igraph.to.graphNEL(signal_directed_graphObj)
cw <- CytoscapeWindow ('vignette', graph=pp1_graphObj, overwrite=TRUE)
displayGraph (cw)
