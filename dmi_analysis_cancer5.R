#### ---------------- LOAD SAVED OBJECTS --------------------------------------------
setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")

load("cancer5.Rda")
load("cancer5_communities.Rda")

#### ---------------- BASIC SETUP --------------------------------------------

setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
rm(list = ls())

download_path = "/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge"
downloaded_folder <- "DREAM11_subchallenge1_initial_networks"
file_names = paste(download_path,"/",downloaded_folder,"/",list.files(downloaded_folder)
                   ,sep = "")

cancer5 = read.table(file = file_names[5], header = TRUE, sep = "\t")

#### ---------------- LOAD LIBRARIES -------------------------------------------------

library('gplots'); library('ggplot2'); library('WGCNA'); library('igraph'); library('clue'); library(ProNet)

#### ---------------- CREATE IGRAPH OBJECT -------------------------------------------------

colnames(cancer5) = c("V1","V2","weight"); head(cancer5)
cancer5$V1 = paste("G",cancer5$V1,sep = "")
cancer5$V2 = paste("G",cancer5$V2,sep = "")
cancer5_igraph = graph.data.frame(d = cancer5, directed = F);
vcount(cancer5_igraph) # 14679
ecount(cancer5_igraph) # 999999
save(cancer5,cancer5_igraph, file = "cancer5.Rda")

#### ---------------- FAST GREEDY ---------------------------------------------------------
ptm = proc.time()
cancer5_fastGreedy = fastgreedy.community(graph = cancer5_igraph, 
                                            weights = E(cancer5_igraph)$weight); 
cancer5_fastGreedy_time = proc.time() - ptm

#### ---------------- INFO MAP ---------------------------------------------------------
ptm = proc.time()
cancer5_infoMap = infomap.community(graph = cancer5_igraph, 
                                      e.weights = E(cancer5_igraph)$weight); 
cancer5_infoMap_time = proc.time() - ptm

save(cancer5_fastGreedy,cancer5_infoMap, file = "cancer5_communities.Rda")

#### ---------------- WALK TRAP ---------------------------------------------------------
ptm = proc.time()
cancer5_walkTrap = walktrap.community(graph = cancer5_igraph, 
                                        weights = E(cancer5_igraph)$weight); 
cancer5_walkTrap_time = proc.time() - ptm

#### ---------------- LABEL PROPAGRATION ---------------------------------------------------------
ptm = proc.time()
cancer5_labelPropagation = label.propagation.community(graph = cancer5_igraph, 
                                                         weights = E(cancer5_igraph)$weight); 
cancer5_labelPropagation_time = proc.time() - ptm

#### ---------------- LOUIVAIN ---------------------------------------------------------
ptm = proc.time()
cancer5_louvain = cluster_louvain(cancer5_igraph, 
                                    weights = E(cancer5_igraph)$weight); 
cancer5_louvain_time = proc.time() - ptm

save(cancer5_fastGreedy, cancer5_infoMap, cancer5_walkTrap, 
     cancer5_labelPropagation, cancer5_louvain, 
     file = "cancer5_communities.Rda")

table(cancer5_fastGreedy$membership)
table(table(cancer5_fastGreedy$membership))
