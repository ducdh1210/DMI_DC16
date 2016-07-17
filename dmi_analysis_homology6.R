#### ---------------- BASIC SETUP --------------------------------------------

setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
rm(list = ls())

download_path = "/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge"
downloaded_folder <- "DREAM11_subchallenge1_initial_networks"
file_names = paste(download_path,"/",downloaded_folder,"/",list.files(downloaded_folder)
                   ,sep = "")

homology6 = read.table(file = file_names[6], header = TRUE, sep = "\t")

#### ---------------- LOAD LIBRARIES -------------------------------------------------

library('gplots'); library('ggplot2'); library('WGCNA'); library('igraph'); library('clue'); library(ProNet)

#### ---------------- CREATE IGRAPH OBJECT -------------------------------------------------

colnames(homology6) = c("V1","V2","weight"); head(homology6)
homology6$V1 = paste("G",homology6$V1,sep = "")
homology6$V2 = paste("G",homology6$V2,sep = "")
homology6_igraph = graph.data.frame(d = homology6, directed = F);
vcount(homology6_igraph) # 10475
ecount(homology6_igraph) # 8506428
homology6_igraph = simplify(homology6_igraph)
vcount(homology6_igraph) # 10475
ecount(homology6_igraph) # 4223606
save(homology6,homology6_igraph, file = "homology6.Rda")

#### ---------------- FAST GREEDY ---------------------------------------------------------
ptm = proc.time()
homology6_fastGreedy = fastgreedy.community(graph = homology6_igraph, 
                                       weights = E(homology6_igraph)$weight); 
homology6_fastGreedy_time = proc.time() - ptm
# user  system elapsed 
# 712.224   0.376 714.187

#### ---------------- INFO MAP ---------------------------------------------------------
ptm = proc.time()
homology6_infoMap = infomap.community(graph = homology6_igraph, 
                                 e.weights = E(homology6_igraph)$weight); 
homology6_infoMap_time = proc.time() - ptm
# user  system elapsed 
# 506.560   0.336 507.330 

#save(homology6_fastGreedy,homology6_infoMap, file = "homology6_communities.Rda")

#### ---------------- WALK TRAP ---------------------------------------------------------
ptm = proc.time()
homology6_walkTrap = walktrap.community(graph = homology6_igraph, 
                                   weights = E(homology6_igraph)$weight); 
homology6_walkTrap_time = proc.time() - ptm
# user  system elapsed 
# 283.508   0.796 284.300


#### ---------------- LABEL PROPAGRATION ---------------------------------------------------------
ptm = proc.time()
homology6_labelPropagation = label.propagation.community(graph = homology6_igraph, 
                                                    weights = E(homology6_igraph)$weight); 
homology6_labelPropagation_time = proc.time() - ptm
# user  system elapsed 
# 7.108   0.064   7.164 

#### ---------------- LOUIVAIN ---------------------------------------------------------
ptm = proc.time()
homology6_louvain = cluster_louvain(homology6_igraph, 
                               weights = E(homology6_igraph)$weight); 
homology6_louvain_time = proc.time() - ptm
# user  system elapsed 
# 7.984   0.100   8.075  

save(homology6_fastGreedy, homology6_infoMap, homology6_walkTrap, 
     homology6_labelPropagation, homology6_louvain, 
     file = "homology6_communities.Rda")


#### ---------------- SPIN GLASS ---------------------------------------------------------
ptm = proc.time()
homology6_spinglass = spinglass.community(graph = homology6_igraph, 
                                     weights = E(homology6_igraph)$weight); proc.time() - ptm
