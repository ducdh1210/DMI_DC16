#### ------------- Load libaries ----------------
library(igraph)
library(stinepack)
library(Hmisc)
library(pracma)

#### ------------- Change directory ----------------
setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/WG-Cluster-Applied/WG-Cluster-Sources")


#### ------------- Load input ---------------
# file_name = "2_ppi_anonym.txt"
file_name = "3_signal_anonym_directed.txt"
input.db = read.table(file = file_name, header = FALSE, sep = "\t")
colnames(input.db) = c("NodeA", "NodeB", "Weight")
input.db$NodeA = paste("G",input.db$NodeA,sep = "")
input.db$NodeB = paste("G",input.db$NodeB,sep = "")

node_names = unique(c(input.db$NodeA,input.db$NodeB))
node.prob.db = data.frame(Node = node_names, 
                          Probability = 0.5, 
                          Standard.error = 0)


#### --------- function uploaded ---------------------------
print("Function uploading and output-folder creation.")
source("./functions.r")
# set the seed to make the output reproducible
set.seed(1)

#### --------- graphical summary of the input data---------------------------
jpeg("./summary_plot.jpg", width = 780, height = 480)
par(mfrow=c(2,2))
hist(input.db[,3], main="Edge weight", xlab="Weight", col="lightgray")
hist(as.numeric(node.prob.db[,2]), main="Node probability", xlab="Probability", col="lightgray")
hist(as.numeric(node.prob.db[,3]), main="Standard error on node probability", xlab="Standard error", col="lightgray")
errbar(1:dim(node.prob.db)[1], as.numeric(node.prob.db[,2]), 
       as.numeric(node.prob.db[,2]) + as.numeric(node.prob.db[,3]),
       as.numeric(node.prob.db[,2]) - as.numeric(node.prob.db[,3]),
       xlab="Node ID", ylab="Probability")
title("Error bars plot of node probability")
dev.off()

#### --------- create graph ---------------------------
print("Graph construction.")
nodeA <- as.vector(input.db[,1]) # length = 397308
nodeB <- as.vector(input.db[,2]) # 
w.g <- as.numeric(as.vector(input.db[,3]))
actors <- data.frame(name=unique(c(nodeA, nodeB))) #12420  rows
relations <- data.frame(from=nodeA,
                        to=nodeB,
                        w=w.g)
directedness = TRUE;
g <- graph.data.frame(relations, directed=directedness, vertices=actors)
# assign weight to edges
E(g)$weight <- relations$w
# assign names to node labels
V(g)$label <- as.vector(actors$name)
# save the graph
write.graph(g, file="Network_to_analyze.graphml", format="graphml")

#### --------- Graph clustering by edge weight ---------------------------
# Determine number of clusters:
# we chose the minimum number of clusters that zeroes the within gourps sum of squares
mydata <- as.matrix(E(g)$weight)
# This is the number of clusters
print("Estimation of the optimal number of clusters.")
source("functions.r")
### correspond to section 2.2.1 - Detection of subgraphs
ptm = proc.time()
# Important: cluster.max is selected to be 64 by experiment to avoid
# error message "more cluster centers than distinct data points"
# a systematic way to determine the number of cluster.max would be needed 
# In brief, for now value for cluster.max should be choosen case by case
nr.clusters <- find_N_cluster(my.data = mydata,
                              cluster.max = 64,
                              iter.max = 1000); proc.time() - ptm 
# takes too long to run, have to stop in the middle
paste("Number of subgraphs (nr.cluster) is: ", nr.clusters, sep="")
#########################
# Now, use k-means to find clusters of edges w.r.t. edges weight
print("K-means clustering of the graph into sub-graphs of similar edge weight.")
weight.cluster <- kmeans(E(g)$weight, nr.clusters, algorithm="Lloyd", iter.max=10000)
class(weight.cluster)

# matrix whose first column is the vectors of edges weights
# and the second column is the cluster to which an edge belongs
edges.weights.cl <- cbind(E(g)$weight,weight.cluster$cluster)


# Extraction of sub-networks by clustering edge weights
# This is the list of sub-graphs

# Extraction of sub-network w.r.t. clusters of edge weights
# n is the number of cluster
# weight.clusters is the matrix whose first column is the vectors of edges weights
# and the second column is the cluster to which an edge belongs
list.of.subgraph <- extract.subn(graph = g, 
                                 weight.clusters = edges.weights.cl, 
                                 n = nr.clusters)
class(list.of.subgraph) # list
length(list.of.subgraph) 
analysis.summary <- analyze.subgraph(list.of.subgraph)

write(c("Subgraph", "Is.connected", "Nr.connected.components", "Nr.nodes", "Nr.edges", "Mean.weight"), 
      file="./analysis.summary.txt", ncolumns=6, sep="\t")
write(t(as.matrix(analysis.summary)), file="./analysis.summary.txt", 
      ncolumns=6, sep="\t", append=T)

##############################################################
# ANALYSIS OF CONNECTED COMPONENTS: ENTROPY AND SIGNIFICANCE #
##############################################################

# extract connected components from each sub-graph
print("Detection of connected components in each sub-graph.")

# this function generates the graphs of the connected components of a sub-graph
# the output is a list of lists of graphs (i.e. for each sub-graph in the input 
# list we have a list of connected components.) saved in .graphml format and a list
# the entropies for each conencted components
connected_components <- decompose_subgraph(list.of.subgraph, node.prob.db) # runable, yay!
sum(sapply(X = connected_components, function(subgraph) length(subgraph))) # 491


# this function saves as GRAPHML format the connected components of each sub-graph 
save_cc_component(list.of.subgraph)

# Entropy of connected components
print("Calculation of connected components's entropy.")

C_entropies <- community_entropy(connected_components)

# generation of random connected components
print("Generation of random connected components.")
system.time(source("./random_cc_components.r")) # runable, about 10 mins

# generate SIF formats sub.graphs and conencted components
print("Data rearrangments.")
system.time(source("./SIFgenerator.r")) # quick, run in 13 seconds

# convolution analysis
system.time(source("./addition.r")) # quick, run in 36 seconds

print("Convolution od mean edge weight and connected component entropy.")
system.time(source("./convolve.r")) # run very quickly



