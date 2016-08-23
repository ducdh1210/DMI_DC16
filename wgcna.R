setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
require(igraph); require(WGCNA)
allowWGCNAThreads()

## -------------Dataset 1 ----------------------------------------------------------------
load("ppi1.Rda")

ppi1_edgeList = get.edgelist(ppi1_igraph)

# get adjacency matrix of type dgCMatrix
ppi1_adjMatrix = get.adjacency(ppi1_igraph, attr = "weight")
ppi1_adjMatrix@Dim # [1] 12420 12420

# convert dgCMatrix to normal matrix object
ppi1_adjM = as.matrix(ppi1_adjMatrix)

gc()
# Cluster gene expression profiles; the flashClust function from
# the authors of WGCNA is another options for larger datasets.
# For input, we use the reciprocal of the adjacency matrix; hierarchical
# clustering works by comparing the _distance_ between objects instead of the
# _similarity_.
ppi1_gene_tree = hclust(as.dist(1 - ppi1_adjM), method="average")
gc()
# we will use the cuttreeDynamicTree method to break apart the hc dendrogram
# into separate modules
ppi1_wgcna_modules = cutreeDynamicTree(dendro=ppi1_gene_tree, minModuleSize=3,
                                       deepSplit=TRUE)
names(ppi1_wgcna_modules) = V(ppi1_igraph)$name

save(ppi1_wgcna_modules, file = "wgcna_modules.rda")

## -------------Dataset 2 ----------------------------------------------------------------

load("ppi2_igraph.rda")
load("ppi2_communities.rda")

adjMatrix = get.adjacency(pp2_igraph, attr = "weight")
adjM = as.matrix(adjMatrix); gc()
gene_tree = hclust(as.dist(1 - adjM), method="average"); gc()
wgcna_modules = cutreeDynamicTree(dendro=gene_tree, 
                                  minModuleSize=3,
                                  deepSplit=TRUE)

names(wgcna_modules) = V(pp2_igraph)$name
table(table(wgcna_modules))

modularity(x = pp2_igraph, membership = ppi2_louvain$membership, weights = E(pp2_igraph)$weight)

mclust::adjustedRandIndex(wgcna_modules, ppi2_louvain$membership) # 0.01
mclust::adjustedRandIndex(wgcna_modules, ppi2_walktrap$membership) # 0.006277767
mclust::adjustedRandIndex(wgcna_modules, ppi2_fastgreedy$membership) # 0.003222752
mclust::adjustedRandIndex(wgcna_modules, ppi2_infomap$membership) # [1] 0.1066305

gc()

## -------------Dataset 3 ----------------------------------------------------------------
rm(list = ls()); gc()
load("signaling3.Rda")
adjMatrix = get.adjacency(signaling3_igraph, attr = "weight")
adjM = as.matrix(adjMatrix)
gene_tree = hclust(as.dist(1 - adjM), method="average")
wgcna_modules = cutreeDynamicTree(dendro=gene_tree, minModuleSize=3,
                                  deepSplit=TRUE)

## -------------Dataset 4 ----------------------------------------------------------------
rm(list = ls()); gc()
load("coexpr4.Rda")
adjMatrix = get.adjacency(coexpr4_igraph, attr = "weight")
adjM = as.matrix(adjMatrix)
gene_tree = hclust(as.dist(1 - adjM), method="average")
wgcna_modules = cutreeDynamicTree(dendro=gene_tree, minModuleSize=3,
                                  deepSplit=TRUE)
table(table(wgcna_modules))

## -------------Dataset 5 ----------------------------------------------------------------
rm(list = ls()); gc()
load("cancer5.Rda")
adjMatrix = get.adjacency(cancer5_igraph, attr = "weight")
adjM = as.matrix(adjMatrix)
gene_tree = hclust(as.dist(1 - adjM), method="average"); gc()
wgcna_modules = cutreeDynamicTree(dendro=gene_tree, minModuleSize=3,
                                  deepSplit=TRUE)
table(table(wgcna_modules))

## -------------Dataset 6 ----------------------------------------------------------------

load("homology6.Rda")
# get adjacency matrix of type dgCMatrix
homo6_adjMatrix = get.adjacency(homology6_igraph, attr = "weight")
homo6_adjMatrix@Dim 
gc()

# convert dgCMatrix to normal matrix object
homo6_adjM = as.matrix(homo6_adjMatrix)
homo6_gene_tree = hclust(as.dist(1 - homo6_adjM), method="average")
gc()

homo6_wgcna_modules = cutreeDynamicTree(dendro=homo6_gene_tree, minModuleSize=3,
                                       deepSplit=TRUE)
names(homo6_wgcna_modules) = V(homology6_igraph)$name

save(homo6_wgcna_modules, file = "wgcna_modules6.rda")





