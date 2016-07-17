#### ---------------- BASIC SETUP --------------------------------------------

setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
rm(list = ls())

download_path = "/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge"
downloaded_folder <- "DREAM11_subchallenge1_initial_networks"
file_names = paste(download_path,"/",downloaded_folder,"/",list.files(downloaded_folder)
                   ,sep = "")

coexpr4 = read.table(file = file_names[4], header = TRUE, sep = "\t")

#### ---------------- LOAD LIBRARIES -------------------------------------------------

library('gplots'); library('ggplot2'); library('WGCNA'); library('igraph'); library('clue'); library(ProNet)

#### ---------------- CREATE IGRAPH OBJECT -------------------------------------------------

colnames(coexpr4) = c("V1","V2","weight"); head(coexpr4)
coexpr4$V1 = paste("G",coexpr4$V1,sep = "")
coexpr4$V2 = paste("G",coexpr4$V2,sep = "")
coexpr4_igraph = graph.data.frame(d = coexpr4, directed = F);
save(coexpr4,coexpr4_igraph, file = "coexpr4.Rda")

#### ---------------- FAST GREEDY ---------------------------------------------------------
ptm = proc.time()
coexpr4_fastGreedy = fastgreedy.community(graph = coexpr4_igraph, 
                                          weights = E(coexpr4_igraph)$weight); proc.time() - ptm
# user  system elapsed 
# 201.056   0.024 201.642

#### ---------------- INFO MAP ---------------------------------------------------------
ptm = proc.time()
coexpr4_infoMap = infomap.community(graph = coexpr4_igraph, 
                                    e.weights = E(coexpr4_igraph)$weight); proc.time() - ptm
# user  system elapsed 
# 85.424   0.072  85.852 
save(coexpr4_fastGreedy,coexpr4_infoMap, file = "coexpr4_communities.Rda")

#### ---------------- SPIN GLASS (NOT RUNABLE!)---------------------------------------------------------
ptm = proc.time()
coexpr4_spinglass = spinglass.community(graph = coexpr4_igraph, 
                                      weights = E(coexpr4_igraph)$weight); proc.time() - ptm



