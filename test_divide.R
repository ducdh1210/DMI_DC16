# --------------------------- divide big graph into small graphs -----------------
rm(list = ls());  gc();
setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
source("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/graph_utility_functions.R")
source("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/divide_subgraph.R")

# load all graphs from 6 datasets
load("ppi1.Rda")
# load("ppi2_igraph.rda")
# load("signaling3.Rda")
# load("coexpr4.Rda")
# load("cancer5.Rda")
load("homology6.Rda")
# load libraries
library('igraph'); library('rlist'); library('reshape')
#library('WGCNA'); library('clue'); library('ProNet')

# set this to increase the recursion depth
options(expressions = 500000) # default value is 5000

# graph to work with
graph = ppi1_igraph

# ---------- find all connected components -----------------------------------------------

global_good_list = list()
initial_list_of_big_graphs = list()

# find connected_components
connected_components = components(graph)
inital_list_of_graphs = decompose.graph(graph)

# for (i in sort(unique(connected_components$membership))){
#   list_of_graphs[[i]] = induced.subgraph(graph,
#                                                 which(connected_components$membership==i))
# }

# doing 2 things
# 1: get inital list of big graphs
# 2: add small list into the result class
for (i in 1:length(inital_list_of_graphs)){
  if (vcount(inital_list_of_graphs[[i]]) > 100){
    initial_list_of_big_graphs <<- list.append(initial_list_of_big_graphs, inital_list_of_graphs[[i]])
  }else if(vcount(inital_list_of_graphs[[i]]) >= 3){
    global_good_list <<- list.append(global_good_list, inital_list_of_graphs[[i]])
  }
}

# we do this because we have only 1 big graph 
main_graph = initial_list_of_big_graphs[[1]]

#table(table(connected_components$membership))
# 2     3     4 12325 
# 33     7     2     1 

# # groups of at least 3 connected components 
# selected_connected_component_groups = unname(which(table(connected_components$membership) >= 3))
# # select the node names in from those above groups 
# main_nodes = names(connected_components$membership[
#   which(connected_components$membership %in% selected_connected_component_groups)])
# # create the subgraph out the the above nodes 
# main_graph = induced_subgraph(graph, vids = main_nodes)
# vcount(main_graph) 
# ecount(main_graph) 
# hist(E(main_graph)$weight)
# plot(density(E(main_graph)$weight))

# apply louvain on the subgraphs 
main_graph_communities = cluster_louvain(graph = main_graph,
                                         weights = E(main_graph)$weight )
# main_graph_communities = cluster_infomap(graph = main_graph,
#                                          e.weights = E(main_graph)$weight )
# modularity(main_graph_communities) # 0.4959873

table(main_graph_communities$membership)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18 
# 892  923 2914 5022 1962    7 4617  131  405  353    6  136    2    4    5    4    3    2 

table(table(main_graph_communities$membership))
# 2    3    4    5    6    7  131  136  353  403  899  925 1964 2921 4614 5009 
# 2    2    3    1    1    1    1    1    1    1    1    1    1    1    1    1
# which clusters have more than 100 memberships
# which(table(main_graph_communities$membership) > 100)

# create list of graphs 
louvain_list_of_graphs = getAllSubgraphsFromMembership(igraphObject = main_graph, 
                                               membership = main_graph_communities$membership)


# only for testing purpose
# temp_list_of_graphs_is_connected = unlist(lapply(temp_list_of_graphs, function(subgraph) is.connected(subgraph)))
# temp_list_of_graphs_size = unlist(lapply(temp_list_of_graphs, function(subgraph) vcount(subgraph)))
# names(temp_list_of_graphs_is_connected) = as.character(temp_list_of_graphs_size)

#------------ function to iteratively decompose list_of_graph ---------------------------

require(rlist)
list_of_connected_big_graphs = list()

decomposeGraphsInList = function(graphList){
  for (i in 1:length(graphList)){
    #subgraph = graphList[[i]];
    if (is.connected(graphList[[i]])){
      if (vcount(graphList[[i]]) > 100){
        list_of_connected_big_graphs <<- list.append(list_of_connected_big_graphs, graphList[[i]])
      }else if (vcount(graphList[[i]]) >= 3){
        global_good_list <<- list.append(global_good_list, graphList[[i]])
      }
    }else{ # not connected
      temp_list = decompose.graph(graphList[[i]])
      for (j in 1:length(temp_list)){
        if (vcount(temp_list[[j]]) > 100){ 
          list_of_connected_big_graphs <<- list.append(list_of_connected_big_graphs, temp_list[[j]])
        }else if (vcount(temp_list[[j]]) >= 3){ 
          global_good_list <<- list.append(global_good_list, temp_list[[j]])
        }
      }
    }
  }
}

decomposeGraphsInList(louvain_list_of_graphs) # worked!

# ---------------- recursive bisection on temp_list_of_graphs ----------------------------
list_of_connected_big_graphs
size = sapply(list_of_connected_big_graphs, function(e) vcount(e)); size

# if the big-graph has odd node, remove 1 node whose edge from it has smallest weight


removeOneNode = function(igraphObject){
  node_to_removed = names(sort(degree(igraphObject), decreasing = F)[1])
  nodes_to_retain = names(V(igraphObject))[-which(names(V(igraphObject)) == node_to_removed)]
  graph = induced_subgraph(graph = igraphObject, vids = nodes_to_retain)
  return(graph)
}

# https://www.r-bloggers.com/graph-bisection-in-r/
## Approximate bisection
# returns a bisection of the graph that minimizes the cost using Kernighan/Lin Algorithm
# http://www.eecs.berkeley.edu/~demmel/cs267/lecture18/lecture18.html#link_4.2
# partition<-approximateBisection(weightMatrix)
# weightMatrix is symmetric matrix of size 2Nx2N made of non-negative values.
# partition is a list of two vectors of N indices.
# C.Ladroue

approximateBisection<-function(weightMatrix,mode="matrix",minimumGain=1e-5){
  #   minimumGain<-1e-5 # minimum value for gain, setting it to 0 might lead to infinite loop due to numerical inaccuracy
  
  N<-dim(weightMatrix)[1] # number of elements
  m<-N/2
  
  # start off with a random partition
  A<-sample(1:N,N/2,replace=FALSE)
  B<-(1:N)[-A]
  
  maxGain<-Inf
  while(maxGain>minimumGain){
    DA<-rowSums(weightMatrix[A,B])-rowSums(weightMatrix[A,A])+diag(weightMatrix[A,A])
    DB<-rowSums(weightMatrix[B,A])-rowSums(weightMatrix[B,B])+diag(weightMatrix[B,B])
    unmarkedA<-1:m
    unmarkedB<-1:m
    markedA<-rep(0,m)
    markedB<-rep(0,m)
    gains<-rep(0,m)
    for(k in 1:m){
      # find best pair from remainder
      # with 2 loops, slow but easy on memory
      if(mode=='2loops'){
        bestGain<--Inf
        besti<-0
        bestj<-0
        for(i in unmarkedA)
          for(j in unmarkedB){
            gain<-DA[i]+DB[j]-2*weightMatrix[A[i],B[j]]
            if(gain>bestGain) {bestGain<-gain; besti<-i;bestj<-j}
          }
        #           mark the best pair
        unmarkedA<-unmarkedA[-which(unmarkedA==besti)]
        unmarkedB<-unmarkedB[-which(unmarkedB==bestj)]
        markedA[k]<-besti
        markedB[k]<-bestj
      }              
      # with one matrix, much faster but builds a matrix as large as weightMatrix
      if(mode=='matrix'){
        dimension<-m+1-k
        fasterGain<-matrix(DA[unmarkedA],nrow=dimension,ncol=dimension,byrow=FALSE)+
          matrix(DB[unmarkedB],nrow=dimension,ncol=dimension,byrow=TRUE)-
          2*weightMatrix[A[unmarkedA],B[unmarkedB]]
        # mark the best pair
        best<-arrayInd(which.max(fasterGain),.dim=c(dimension,dimension))
        besti<-unmarkedA[best[1]]
        bestj<-unmarkedB[best[2]]
        bestGain<-fasterGain[best]
        markedA[k]<-unmarkedA[best[1]]
        markedB[k]<-unmarkedB[best[2]]
        unmarkedA<-unmarkedA[-best[1]]
        unmarkedB<-unmarkedB[-best[2]]
      }
      # record gain
      gains[k]<-bestGain
      
      # update D for unmarked indices 
      DA[unmarkedA]<-DA[unmarkedA]+2*weightMatrix[A[unmarkedA],A[besti]]-2*weightMatrix[A[unmarkedA],B[bestj]]
      DB[unmarkedB]<-DB[unmarkedB]+2*weightMatrix[B[unmarkedB],B[bestj]]-2*weightMatrix[B[unmarkedB],A[besti]]
    }
    gains<-cumsum(gains)
    bestPartition<-which.max(gains)
    maxGain<-gains[bestPartition]
    if(maxGain>minimumGain){ 
      # swap best pairs
      A1<-c(A[-markedA[1:bestPartition]],B[markedB[1:bestPartition]])
      B1<-c(B[-markedB[1:bestPartition]],A[markedA[1:bestPartition]])
      A<-A1
      B<-B1
    }
  }
  list(A,B)
}

# -------------------------------------------------------------------------------

# bisection procedure:
# first check if the graph has odd size; if yes, then remove 1 node
#     then perform bisection
# 
#     check if the two bisected graphs has size < 100
#         if yes, then add to the good graph
#         if no, then call recursive bisection

recursiveBisection = function(igraphObject){
  # base case 
  if(vcount(igraphObject) <= 100 & vcount(igraphObject) > 2){
    global_good_list <<- append(global_good_list, igraphObject)
    return()
  }
  
  # if the graph size is odd, remove a node to make it even 
  if (vcount(igraphObject) %% 2 == 1){
    igraphObject = removeOneNode(igraphObject)
  }
  
  # build weighted matrix from graph
  adjMatrix = get.adjacency(igraphObject, attr = "weight")
  # convert dgCMatrix to normal matrix object
  adjM = as.matrix(adjMatrix); rm(adjMatrix); gc();
  
  # call approximateBisection
  bisection_result_asssignment = approximateBisection(adjM);
  graph1 = induced_subgraph(graph = igraphObject, vids = bisection_result_asssignment[[1]]);
  graph2 =  induced_subgraph(graph = igraphObject, vids = bisection_result_asssignment[[2]]);
  bisected_graph_list = list(graph1, graph2)

  # loop through each graph in the list bisection_result_list
  for (i in 1:2){
    # if the graph is connected
    bisected_graph = bisected_graph_list[[i]]
    if (is.connected(bisected_graph)){
      if (vcount(bisected_graph) >100) {
        # if the graph is big size (and connected from above condiction) --> bisect it
        recursiveBisection(bisected_graph)
      }else if (vcount(bisected_graph) >3 ){
        # if the graph is good size (and connected from above condiction) --> add to good list
        global_good_list <<- list.append(global_good_list, bisected_graph)
      }
    }else{ # if the graph is not connected, first decompose it into connected component
      temp_list = decompose.graph(bisected_graph)
      # loop through each connected component
      for (j in 1:length(temp_list)){
        if (vcount(temp_list[[j]]) > 100){ # if the size is big, then bisecting it
          recursiveBisection(temp_list[[j]])
        }else if (vcount(temp_list[[j]]) >= 3){ # if the size is good, then
          global_good_list <<- list.append(global_good_list, temp_list[[j]])
        }
      }
    }
  }
}


ptm = proc.time()
# recursiveBisection(list_of_connected_big_graphs[[4]])
for (i in 1:length(list_of_connected_big_graphs)){
  recursiveBisection(list_of_connected_big_graphs[[i]])
}
ptm = proc.time() - ptm

global_good_list
save(global_good_list, file = "result/ppi1_result.rda")

# --------------------------------------------------------------------------------------

length(list_of_graphs)
list_size = unlist(lapply(temp_list_of_graphs, function(e) vcount(e)))

# list_size
# [1]  899  925 2921 5009 1964    7 4614  131  403  353    6  136    2    4    4    3    5    4    3    2
# create lists to retain and manipulate results
global_good_list = list()
list_of_small_graph = list()
list_of_big_graphs = list()

# doing 2 things
# 1: get inital list of big graphs
# 2: add small list into the result class
for (i in 1:length(list_of_graphs)){
  if (vcount(list_of_graphs[[i]]) > 100){
    list_of_big_graphs <<- list.append(list_of_big_graphs, list_of_graphs[[i]])
  }else if(vcount(list_of_graphs[[i]]) > 2){
    global_good_list <<- list.append(global_good_list, list_of_graphs[[i]])
  }
}

subgraph_stat = getSubgraphStat(list_of_big_graphs)

test_graph = list_of_big_graphs[[5]]

test_adjMatrix = get.adjacency(test_graph, attr = "weight")
test_adjMatrix@Dim # [1] 12420 12420
# convert dgCMatrix to normal matrix object
test_adjM = as.matrix(test_adjMatrix); rm(test_adjMatrix); gc()

ptm = proc.time()
bisection_result = approximateBisection(ppi1_adjM)
proc.time() - ptm

# https://www.r-bloggers.com/graph-bisection-in-r/
## Approximate bisection
# returns a bisection of the graph that minimizes the cost using Kernighan/Lin Algorithm
# http://www.eecs.berkeley.edu/~demmel/cs267/lecture18/lecture18.html#link_4.2
# partition<-approximateBisection(weightMatrix)
# weightMatrix is symmetric matrix of size 2Nx2N made of non-negative values.
# partition is a list of two vectors of N indices.
# C.Ladroue

approximateBisection<-function(weightMatrix,mode="matrix",minimumGain=1e-5){
  #   minimumGain<-1e-5 # minimum value for gain, setting it to 0 might lead to infinite loop due to numerical inaccuracy
  
  N<-dim(weightMatrix)[1] # number of elements
  m<-N/2
  
  # start off with a random partition
  A<-sample(1:N,N/2,replace=FALSE)
  B<-(1:N)[-A]
  
  maxGain<-Inf
  while(maxGain>minimumGain){
    DA<-rowSums(weightMatrix[A,B])-rowSums(weightMatrix[A,A])+diag(weightMatrix[A,A])
    DB<-rowSums(weightMatrix[B,A])-rowSums(weightMatrix[B,B])+diag(weightMatrix[B,B])
    unmarkedA<-1:m
    unmarkedB<-1:m
    markedA<-rep(0,m)
    markedB<-rep(0,m)
    gains<-rep(0,m)
    for(k in 1:m){
      # find best pair from remainder
      # with 2 loops, slow but easy on memory
      if(mode=='2loops'){
        bestGain<--Inf
        besti<-0
        bestj<-0
        for(i in unmarkedA)
          for(j in unmarkedB){
            gain<-DA[i]+DB[j]-2*weightMatrix[A[i],B[j]]
            if(gain>bestGain) {bestGain<-gain; besti<-i;bestj<-j}
          }
        #           mark the best pair
        unmarkedA<-unmarkedA[-which(unmarkedA==besti)]
        unmarkedB<-unmarkedB[-which(unmarkedB==bestj)]
        markedA[k]<-besti
        markedB[k]<-bestj
      }              
      # with one matrix, much faster but builds a matrix as large as weightMatrix
      if(mode=='matrix'){
        dimension<-m+1-k
        fasterGain<-matrix(DA[unmarkedA],nrow=dimension,ncol=dimension,byrow=FALSE)+
          matrix(DB[unmarkedB],nrow=dimension,ncol=dimension,byrow=TRUE)-
          2*weightMatrix[A[unmarkedA],B[unmarkedB]]
        # mark the best pair
        best<-arrayInd(which.max(fasterGain),.dim=c(dimension,dimension))
        besti<-unmarkedA[best[1]]
        bestj<-unmarkedB[best[2]]
        bestGain<-fasterGain[best]
        markedA[k]<-unmarkedA[best[1]]
        markedB[k]<-unmarkedB[best[2]]
        unmarkedA<-unmarkedA[-best[1]]
        unmarkedB<-unmarkedB[-best[2]]
      }
      # record gain
      gains[k]<-bestGain
      
      # update D for unmarked indices 
      DA[unmarkedA]<-DA[unmarkedA]+2*weightMatrix[A[unmarkedA],A[besti]]-2*weightMatrix[A[unmarkedA],B[bestj]]
      DB[unmarkedB]<-DB[unmarkedB]+2*weightMatrix[B[unmarkedB],B[bestj]]-2*weightMatrix[B[unmarkedB],A[besti]]
    }
    gains<-cumsum(gains)
    bestPartition<-which.max(gains)
    maxGain<-gains[bestPartition]
    if(maxGain>minimumGain){ 
      # swap best pairs
      A1<-c(A[-markedA[1:bestPartition]],B[markedB[1:bestPartition]])
      B1<-c(B[-markedB[1:bestPartition]],A[markedA[1:bestPartition]])
      A<-A1
      B<-B1
    }
  }
  list(A,B)
}

# getSubgraphStat(list_of_subgraphs = list_of_big_graphs)
# # 
# # m9 = list_of_big_graphs[[9]]
# # m9_louvain = cluster_louvain(graph = m9, weights = E(m9)$weight)
# # m9.1 = induced_subgraph(m9, vids = which(m9_louvain$membership == 1))
# # m9.1_louvain = cluster_louvain(m9.1, weights = E(m9.1)$weight)
# # divide_subgraph(list_of_big_graphs)
# # global_good_list
# 
# divide_subgraph(list_of_big_graphs = list_of_big_graphs)
# 
# 
# # there will be some graphs whose size < 3, thus remove those graphs from the list
# good_list_size = sapply(global_good_list, function(module) vcount(module))
# index_to_be_removed = which(good_list_size < 3)
# global_good_list = global_good_list[-index_to_be_removed]
# 
# # use reshape::melt to create membership vector
# global_good_list_node_names = lapply(global_good_list, function(module) names(V(module)))
# module_assignment = melt(global_good_list_node_names)
# membership = module_assignment$L1
# names(membership) = module_assignment$value
# 
# # since we have removed some graphs whose size < 3 
# # (implying we removed some nodes), thus we need to create a new induced graph 
# # from the retained nodes in order to test the modularity
# induced_graph_from_main = induced_subgraph(graph = main_graph, 
#                                            vids = names(membership))
# name_order = V(induced_graph_from_main)$name
# membership = membership[name_order]
# 
# # modularity(induced_graph_from_main, membership = membership, weights = E(induced_graph_from_main)$weights) # 0.1396317
# # save(membership, file = "result/signal3_infomap.rda")



cut.dend <- function(graph, comm, maxsize) {
  vc <- vcount(graph)
  steps <- 0
  sizes <- c(rep(1, vc), rep(0, vc-1))
  msize <- 1
  while (msize <= maxsize) {
    act <- comm$merges[steps+1,]
    sizes[steps+vc+1] <- sizes[act[1]+1] + sizes[act[2]+1]
    if (sizes[steps+vc+1] > msize) { msize <- sizes[steps+vc+1] }
    steps <- steps + 1
  }
  steps-1
}

comm <- fastgreedy.community(test_graph, weights = E(test_graph)$weight)
steps <- cut.dend(test_graph, comm, 10)

## this is the right number of steps
test = cutat(comm, steps)

table(table(test))
















# --------------------------------------------------------------------------------------------------



rm(list = ls());  gc();
setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")
source("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/graph_utility_functions.R")
source("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/divide_subgraph.R")

# load libraries
library('igraph'); library('rlist'); library('reshape')
library('WGCNA'); library('clue'); library('ProNet'); 
# load all graphs from 6 datasets
load("ppi1.Rda"); load("result/temp_wcgna_modules_ppi1.RData")
# load("ppi2_igraph.rda")
# load("signaling3.Rda")
# load("coexpr4.Rda")
# load("cancer5.Rda")
# load("homology6.Rda")

graph = ppi1_igraph 

# remove small subgraphs out of the big graph
connected_components = components(graph)

#table(table(connected_components$membership))
# 2     3     4 12325 
# 33     7     2     1 

# groups of at least 3 connected components 
selected_connected_component_groups = unname(which(table(connected_components$membership) >= 3))

# select the node names in from those above groups 
main_nodes = names(connected_components$membership[
  which(connected_components$membership %in% selected_connected_component_groups)])

# create the subgraph out the the above nodes 
main_graph = induced_subgraph(graph, vids = main_nodes)

# get modules using wcgna
wcgna_modules = get_wcgna_modules(igraphObject = main_graph)
#save(wcgna_modules, file = "result/temp_wcgna_modules_ppi1.RData")

table(table(wcgna_modules))

# create list of graphs 
list_of_graphs = getAllSubgraphsFromMembership2(igraphObject = main_graph, 
                                                membership = wcgna_modules)

# create lists to retain and manipulate results
global_good_list = list()
list_of_small_graph = list()
list_of_big_graphs = list()

gc()

# doing 2 things
# 1: get inital list of big graphs
# 2: add small list into the result class
for (i in 1:length(list_of_graphs)){
  if (vcount(list_of_graphs[[i]]) > 100){
    list_of_big_graphs <<- list.append(list_of_big_graphs, list_of_graphs[[i]])
  }else if(vcount(list_of_graphs[[i]]) > 2){
    global_good_list <<- list.append(global_good_list, list_of_graphs[[i]])
  }
}

divide_subgraph(list_of_big_graphs)
divide_subgraph2(list_of_big_graphs)
