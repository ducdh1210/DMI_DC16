#### ---------------- BASIC SETUP --------------------------------------------

setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
rm(list = ls())

download_path = "/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge"
downloaded_folder <- "DREAM11_subchallenge1_initial_networks"
file_names = paste(download_path,"/",downloaded_folder,"/",list.files(downloaded_folder)
                   ,sep = "")

ppi1 = read.table(file = file_names[1], header = TRUE, sep = "\t")

#### ---------------- LOAD LIBRARIES -------------------------------------------------

library('gplots'); library('ggplot2'); library('WGCNA'); library('igraph'); library('clue'); library(ProNet)

#### ---------------- CREATE IGRAPH OBJECT -------------------------------------------------

colnames(ppi1) = c("V1","V2","weight"); head(ppi1)
ppi1$V1 = paste("P",ppi1$V1,sep = "")
ppi1$V2 = paste("P",ppi1$V2,sep = "")
ppi1_igraph = graph.data.frame(d = ppi1, directed = F);
save(ppi1,ppi1_igraph, file = "ppi1.Rda")

#### ---------------- FAST GREEDY ---------------------------------------------------------
ptm = proc.time()
ppi1_fastGreedy = fastgreedy.community(graph = ppi1_igraph, 
                                          weights = E(ppi1_igraph)$weight); proc.time() - ptm
# user  system elapsed 
# 375.596   0.176 378.403 

#### ---------------- INFO MAP ---------------------------------------------------------
ptm = proc.time()
ppi1_infoMap = infomap.community(graph = ppi1_igraph, 
                                    e.weights = E(ppi1_igraph)$weight); proc.time() - ptm
# user  system elapsed 
# 406.864   0.320 408.083
# save(ppi1_fastGreedy,ppi1_infoMap, file = "ppi1_communities.Rda")

#### ---------------- WALK TRAP ---------------------------------------------------------
ptm = proc.time()
ppi1_walkTrap = walktrap.community(graph = ppi1_igraph, 
                                 weights = E(ppi1_igraph)$weight); proc.time() - ptm
# user  system elapsed 
# 325.424   0.392 326.185 

#### ---------------- LABEL PROPAGRATION ---------------------------------------------------------
ptm = proc.time()
ppi1_labelPropagation = label.propagation.community(graph = ppi1_igraph, 
                                                    weights = E(ppi1_igraph)$weight); proc.time() - ptm
# user  system elapsed 
# 3.080   0.012   3.087

#### ---------------- LOUIVAIN ---------------------------------------------------------
ptm = proc.time()
ppi1_louvain = cluster_louvain(ppi1_igraph, 
                               weights = E(ppi1_igraph)$weight); proc.time() - ptm
# user  system elapsed 
# 8.152   0.052   8.665 

save(ppi1_fastGreedy, ppi1_infoMap, ppi1_walkTrap, 
     ppi1_labelPropagation, ppi1_louvain, 
     file = "ppi1_communities.Rda")

#### ---------------- SPIN GLASS ---------------------------------------------------------
ptm = proc.time()
ppi1_spinglass = spinglass.community(graph = ppi1_igraph, 
                                        weights = E(ppi1_igraph)$weight); proc.time() - ptm
