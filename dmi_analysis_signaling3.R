#### ---------------- BASIC SETUP --------------------------------------------

setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
rm(list = ls())

download_path = "/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge"
downloaded_folder <- "DREAM11_subchallenge1_initial_networks"
file_names = paste(download_path,"/",downloaded_folder,"/",list.files(downloaded_folder)
                   ,sep = "")

signaling3 = read.table(file = file_names[3], header = TRUE, sep = "\t")

#### ---------------- LOAD LIBRARIES -------------------------------------------------

library('gplots'); library('ggplot2'); library('WGCNA'); library('igraph'); library('clue'); library(ProNet)

#### ---------------- CREATE IGRAPH OBJECT -------------------------------------------------

colnames(signaling3) = c("V1","V2","weight"); head(signaling3)
signaling3$V1 = paste("G",signaling3$V1,sep = "")
signaling3$V2 = paste("G",signaling3$V2,sep = "")
signaling3_igraph = graph.data.frame(d = signaling3, directed = T);
vcount(signaling3_igraph) # 5254
ecount(signaling3_igraph) # 21825
save(signaling3,signaling3_igraph, file = "signaling3.Rda")

#### ---------------- INFO MAP ---------------------------------------------------------
ptm = proc.time()
signaling3_infoMap = infomap.community(graph = signaling3_igraph, 
                                      e.weights = E(signaling3_igraph)$weight); 
signaling3_infoMap_time = proc.time() - ptm
# user  system elapsed 
# 32.212   0.008  33.710 

#### ---------------- WALK TRAP ---------------------------------------------------------
ptm = proc.time()
signaling3_walkTrap = walktrap.community(graph = signaling3_igraph, 
                                        weights = E(signaling3_igraph)$weight); 
signaling3_walkTrap_time = proc.time() - ptm
# user  system elapsed 
# 3.192   0.032   3.914

save(signaling3_infoMap, signaling3_walkTrap, 
     file = "signaling3_communities.Rda")


