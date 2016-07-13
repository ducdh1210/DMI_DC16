#############################################################################
#                                                                           #
# WG-Cluster (Weighted Graph CLUSTERing)                                    #
# Copyright (C) 2015, 2016 Paola Lecca, Angela Re                           #
#                                                                           #
# This program is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation, version 2 of the licence.                   #
#                                                                           #
# This program is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with this program.  If not, see <http://www.gnu.org/licenses/>      #
#                                                                           #
# Authors' contacts:                                                        #
#                                                                           #
# Paola Lecca,                                                              #
# Dept. of Mathematics, University of Trento,                               #
# via Sommarive 14, 38123 Trento, Italy                                     # 
# E-mail: paola.lecca@unitn.it                                              #
#                                                                           #
# Angela Re                                                                 # 
# Centre for Integrative Biology, University of Trento                      #
# via Sommarive 9, 38123 Trento, Italy                                      #
# E-mail: angela.re@unitn.it                                                #
#                                                                           #
#############################################################################



library(igraph)
library(stinepack)

#######################################################################################
# Prepare the sub-directories to stores the graphs
mainDir <- "./"
subDir1 <- "SubGraphs_GraphML"
subDir2 <- "ConnectedComponents_GraphML"
subDir3 <- "ConnectedComponents_Entropy"

subDir <- c(subDir1, subDir2, subDir3)



for (i in 1:length(subDir))
{
  if (file.exists(paste(mainDir, subDir[i], "/", sep = "/", collapse = "/"))) {
    cat("subDir exists in mainDir and is a directory")
  } 
  else if (file.exists(paste(mainDir, subDir[i], sep = "/", collapse = "/"))) {
    cat("subDir exists in mainDir but is a file")
    # you will probably want to handle this separately
  } 
  else {
    cat("subDir does not exist in mainDir - creating")
    dir.create(file.path(mainDir, subDir[i]))
  }
}
#######################################################################################


# function to to round to the next order of magnitude in R
log10_ceiling <- function(x) {
  (10^(ceiling(log10(x))))/10
}

# if you want to visualize the output of log10_ceiling function
# log10_ceiling (abs(stinemanSlopes(1:150,wss,scale=FALSE)))

#min(log10_ceiling (abs(stinemanSlopes(1:150,wss,scale=FALSE))))
#mean(log10_ceiling (abs(stinemanSlopes(1:150,wss,scale=FALSE)))

m.clusters <- 250 # defau
find_N_cluster <- function(my.data) 
{		
  wss <- (nrow(my.data)-1)*sum(apply(my.data,2,var))
  # the seed of random generator is fixed to initialize the k-means
  # always with the same centroids
  set.seed(1)
  for (i in 2:m.clusters) {
    wss[i] <- sum(kmeans(my.data, centers=i, algorithm="Lloyd", iter.max=10000)$withinss)
  }
  slopes.magnitude <- log10_ceiling (abs(stinemanSlopes(1:m.clusters,wss,scale=FALSE)))
  tr <- slopes.magnitude[1]*0.0001
  min.no <-  match(tr, slopes.magnitude)
  # find the first index of min.no.zero in the wss vector
  min.no
}
					 
					 
# Extraction of sub-network w.r.t. clusters of edge weights
# n is the number of cluster
# weight.clusters is the matrix whose first column is the vectors of edges weights
# and the second column is the cluster to which an edge belongs
extract.subn <- function(graph, weight.clusters, n)
{
  mat <- list()
  length(mat) <- n
  #
  glist <- list()
  length(glist) <- n
  #
  for (i in 1:n)
  {
    mat[[i]] <- which(weight.clusters[,2]!=i)
    glist[[i]] <- delete.edges(graph, E(g)[   as.vector(mat[[i]]) ]    )
    #
    # save each subgraph in graphml format
    #g.fileout <- paste("./", subDir[1], "/", "subgraph_", i, ".graphml", sep="")
    #write.graph(glist[[i]], file=g.fileout, format="graphml")
  }
  glist
}
						
						
analyze.subgraph <- function(list_of_graph)
{                
  len <- length(list_of_graph)
  is.connected <- array("", len )
  n.connected.components <- array(0, len )
  n.nodes <- array(0, len )
  n.edges <- array(0, len )
  clss <- list()
  length(clss) <- len
  #
  mean.weight <- array(0, len )
  #
  new_list_of_graph <- list()
  length(new_list_of_graph) <- len
  #
  for (i in 1:len)
  {
    # remove nodes with null degree to avoid the computeation of single-node
    # connected components
    new_list_of_graph[[i]]<- delete.vertices(list_of_graph[[i]], which(degree(list_of_graph[[i]]) == 0))
    #
    is.connected[i] <- is.connected(new_list_of_graph[[i]], mode="weak")
    clss[[i]] <- clusters( new_list_of_graph[[i]], mode="weak")
    n.connected.components [i] <- no.clusters(new_list_of_graph[[i]], mode="weak") 
    #
    n.edges[i] <- ecount(new_list_of_graph[[i]])
    n.nodes[i] <- vcount(new_list_of_graph[[i]])
    #
    # color the clusters
    V(new_list_of_graph[[i]])$color=clss[[i]]$membership
    #
    # save each subgraph in graphml format
    g.fileout <- paste("./", subDir[1], "/", "subgraph_", i, ".graphml", sep="")
    write.graph(new_list_of_graph[[i]], file=g.fileout, format="graphml")
    #
    # average of edge weights
    mean.weight[i] <- mean(E(new_list_of_graph[[i]])$weight)
    
  }
  subg <- 1:length(is.connected)
  res <- data.frame(subg, is.connected, n.connected.components, n.nodes, n.edges, mean.weight )
  res
}
						
						



# this function computes the weighted degree of the ndoes of connected components
# it takes as an input the list of connected component of subgraphs
#
cc_node_stregth <- function(g.list)
{  
  d.w <- list()
  length(d.w) <- length(g.list)
  for (i in 1:length(g.list))
  {
   d.w[[i]] <- graph.strength(g.list[[i]])
   }
  d.w
}

# function to combine lists of lists of matrixes
attach.list <- function(list1, list2)
{
  res <- list()
  length(res) <- length(list1)
  for (i in 1:length(list1))
  {
    for (j in 1:length(list1[[i]]))
    {
      res[[i]][[j]] <- cbind(list1[[i]][[j]], list2[[i]][[j]])
    }
  }
  res
}



# this function generates the graphs of the connected components of a sub-graph
# the output is a list of lists of graphs (i.e. for each sub-graph in the input 
# list we have a list of connected components.) saved in .graphml format and a list
# the entropies for each conencted components
decompose_subgraph <- function(list_of_subgraph, node.prob)
{                
  
  len <- length(list_of_subgraph)
  #
  new_list_of_subgraph <- list()
  length(new_list_of_subgraph ) <- len
  #
  prob.db <- list()
  length(prob.db) <- len
  #
  n.connected.components <- array(0, len )
  #
  graphs <- list()
  length(graphs) <- len
  #
  strengths <- list()
  length(strengths) <- len
  #
  for (i in 1:len)
  {
    # remove nodes with null degree to avoid the computeation of single-node
    # connected components
    new_list_of_subgraph[[i]]<- delete.vertices(list_of_subgraph[[i]], which(degree(list_of_subgraph[[i]]) == 0))
    #
    n.connected.components[i] <- no.clusters(new_list_of_subgraph[[i]], mode="weak")                            
    # list of list of graph
    graphs[[i]] <- decompose.graph(new_list_of_subgraph[[i]], mode = "weak", max.comps=n.connected.components[i])
    
    strengths[[i]] <- cc_node_stregth(graphs[[i]])
    
    prob.db[[i]] <- assign.probability.to.node(graphs[[i]], node.prob) 
    
  }
  # connected components with the information of node probability and node strenght
  res<- attach.list(prob.db, strengths)
  res
}

# Assignment of probabilities to nodes.
# The inputs of this function are the list of connectened components and the dataset
# of node names and node probabilities.
# The output is a list of dataframes. Each dataframe containes the names and the probabilities
# of the nodes of the corresponding to the connected components.
assign.probability.to.node <- function(list_of_cc_components, node.probability.db)
{
  len <- length(list_of_cc_components)
  cc_component_nodes_weights <- list()
  length(cc_component_nodes_weights) <- len
  #
  #
  # i ranges over the number of connected components to be analyzed
  for (i in 1:len)
  {
    # names of the nodes in the i-th connected component
    cc_component_nodes_names <- as.vector(V(list_of_cc_components[[i]])$label)
    lenss <- length(cc_component_nodes_names)
    #
    # initialize an element of the list as an empty dataframe
    cc_component_nodes_weights[[i]] <- data.frame(node.name=character(lenss), 
                                                  node.probability=numeric(lenss), 
                                                  node.prob.err=numeric(lenss),
                                                  stringsAsFactors=FALSE)
    #
    # j ranges over the number of nodes in a connected component
    for (j in 1:lenss)
    {
      #
      if(cc_component_nodes_names[j] %in% as.vector(node.probability.db[,1]))
      {
        pos <- which(as.vector(node.probability.db[,1]) == cc_component_nodes_names[j])
        cc_component_nodes_weights[[i]]$node.probability[j] <- node.probability.db[pos,2]
        cc_component_nodes_weights[[i]]$node.prob.err[j] <- node.probability.db[pos,3]
        cc_component_nodes_weights[[i]]$node.name[j] <- cc_component_nodes_names[j]
      }
      #
      if((cc_component_nodes_names[j] %in% as.vector(node.probability.db[,1]))==FALSE)
      { 
        print(i);
        print(j);
        print(cc_component_nodes_names[j]);
        stop(paste("Error: a node in a subgraph does not belong to the whole network.", i, j, sep=" "))
      }   
    }
  }
  #cc_component_nodes_weights
  cc_component_nodes_weights
}
							  
							  
# this function saves as GRAPHML format the connected components of each sub-graph 
save_cc_component <- function(list_of_subgraph)
{                
  
  len <- length(list_of_subgraph)
  #
  new_list_of_subgraph <- list()
  length(new_list_of_subgraph ) <- len
  #
  n.connected.components <- array(0, len )
  #
  graphs <- list()
  length(graphs) <- len
  #
  for (i in 1:len)
  {
    # remove nodes with null degree to avoid the computeation of single-node
    # connected components
    new_list_of_subgraph[[i]]<- delete.vertices(list_of_subgraph[[i]], which(degree(list_of_subgraph[[i]]) == 0))
    #
    n.connected.components[i] <- no.clusters(new_list_of_subgraph[[i]], mode="weak")                            
    # list of list of graph
    graphs[[i]] <- decompose.graph(new_list_of_subgraph[[i]], mode = "weak", max.comps=n.connected.components[i])
    
    for(k in 1:length(graphs[[i]]))
    { 
      file.out <- paste("./", subDir[2], "/", "Subgraph_", i, "_component_", k, ".graphml", sep="")
      write.graph(graphs[[i]][[k]], file=file.out, format="graphm")
    }
    
  }
}
	

# compute error on node probability
cc_std_err <- function(prob, prob.err, node_strength)
{
  res <- (prob.err/abs(node_strength))*(abs(1 + log(prob)))
  res 
}

# function to calculate entropy of a connected components
c.entropy <- function(prob, prob.err, node_strength)
                        {
                          cc_entropy <- - sum( ((prob*log(prob))/ abs(node_strength) ), na.rm = T) 
                          cc_entropy_err <- sqrt(sum( (cc_std_err(prob, prob.err, node_strength))^2))
                          res <- paste(cc_entropy, "+/-", cc_entropy_err, sep=" ")
                          res
                        }


double.unlist <- function(community.entropy.res)
{ 
  mat <- list()
  length(mat) <- length(community.entropy.res)
  #
  for (i in 1:length(community.entropy.res))
  {
    mat[[i]] <- do.call(rbind, community.entropy.res[[i]])
  }
  res <- do.call(rbind, mat)
  res
}

community_entropy <- function(connected.components)
                             {
                              CE <- list()
                              length(CE) <- length(connected.components)
                              #
                              row <- list()
                              length(row) <- length(connected.components)
                              #
                              for (i in 1:length(connected.components))
                                  {
                                   for (j in 1:length(connected.components[[i]]))
                                       {
                                         CE[[i]][[j]] <- c.entropy(as.numeric(connected_components[[i]][[j]][,2]), # node probability
                                                                   as.numeric(connected_components[[i]][[j]][,3]), # node probability error
                                                                   as.numeric(connected_components[[i]][[j]][,4]))  # node strength
                                         row[[i]][[j]] <- c(i, j, CE[[i]][[j]])
                                         
                                       }
                                   }
                              res <- double.unlist(row)
                              write(c("Subgraph", "Connected.Component", "Entropy"), file="./ConnectedComponents_Entropy/Entropy_Connected_components.txt", ncol=3)
                              write(t(res), file="./ConnectedComponents_Entropy/Entropy_Connected_components.txt", ncol=3, append=TRUE)
                              res 
                              
                              }
                      





  
