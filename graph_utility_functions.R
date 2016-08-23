#### ---------------- LOAD LIBRARIES-------------------------------------------
require(igraph)

#### ---------------- GET LIST OF SUBGRAPHS -------------------------------------------

# OUTPUT: a list of subgraphs either reprsenting modules from community objects 
#         or cliques from the clique list
getAllSubgraphs = function(igraphObject = NULL, 
                           communityObject = NULL, 
                           clique = NULL, 
                           cliqueSorted = FALSE, 
                           cliqueMembership = NULL,
                           minCliqueSize = 5){
  list_of_all_subgraphs  = list();
  if ((is.null(communityObject) & is.null(clique)) || 
      (!is.null(communityObject) & !is.null(clique))){
    print("either communityObject or clique must be null, and only 1 of them")
    break;
  }else{
    # if communityObject is provided
    if(!is.null(communityObject)){
      for (i in sort(unique(membership(communityObject)))){
        list_of_all_subgraphs[[i]] = induced.subgraph(igraphObject,
                                                      which(membership(communityObject)==i))
      }
    }else{
      # if cliq is provided, first check if it is sorted by size; if not, then 
      #  call reorderCliqBySize to return the ordered clique list
      if (!cliqueSorted){
        clique = reorderCliqBySize(igraphObject = igraphObject, cliq = clique)
      }
      clique = clique[unname(which(table(cliqueMembership) >= minCliqueSize))]
      for (i in 1:length(clique)){
        list_of_all_subgraphs[[i]] = induced.subgraph(igraphObject,
                                                      names(clique[[i]]))
      }
    }
  }

  return(list_of_all_subgraphs)
}

getAllSubgraphsFromMembership = function(igraphObject = NULL, membership = NULL){
  list_of_all_subgraphs = list()
  for (i in sort(unique(membership))){
    list_of_all_subgraphs[[i]] = induced.subgraph(igraphObject,
                                                  which(membership==i))
  }
  return(list_of_all_subgraphs)
}


# get subgraph sizes given a list of subgraphs
# to run this on clique: step 1: call getMaximalClique, step 2: call getAllSubgraphs
# step 3: call this function
getSubgraphSize = function(list_of_subgraphs){
  size = sapply(list_of_subgraphs, function(subgraph) vcount(subgraph))
  return(size)
}


#### ---------------- CLIQUE FUNCTIONS -------------------------------------------

# INPUT:  (1) igraph object 
#         (2) minCliqueSize: minimum clique size to returned
#         (3) maxCliqueSize: maximum clique size to returned 
#         (4) verticesDegree: vectors whose element are the node degree of each 
#             vertex in the original graph. If not provided, the function will compute it
#         (5) degreeThreshold: for computational purpose, it is better to subset 
#             the original graph while calling max_clique(). degreeThreshold is used 
#             to select only vertices whose degree is higher than a set threshold
# OUTPUT: list of cliques whose mininum size specified by minCliqueSize
getMaximalClique = function(igraphObject = NULL, 
                            minCliqueSize = 5, 
                            maxCliqueSize = NULL, 
                            verticesDegree = NULL,
                            degreeThreshold = 10)
{
  # if the vector verticesDegree is not provided, then calling getVerticesDegree()
  if(is.null(verticesDegree)){
    verticesDegree = degree(igraphObject)
  }
  # only get vertices whose degrees are higher than the provided threshold
  idx = which(verticesDegree > degreeThreshold)
  s = V(igraphObject)[idx]
  # obtain cliques with provided parameters
  cliq = max_cliques(igraphObject, 
                     min = minCliqueSize, 
                     max = maxCliqueSize, 
                     subset=s) 
  return(cliq)
}

# INPUT: (1) igraph object (2) a list of cliq (e.g return from max.cliq(..))
# OUTPUT: a list of cliqs ordered in decreasing order based on the  number of 
#         vertices contained in each cliq (first cliq has the biggest size, last
#         cliq has the smallest size)
reorderCliqBySize = function(igraphObject,cliq){
  cliq_sizes = sapply(cliq, function(vertices) length(vertices))
  # name the vector cliq_size by the cliq index order 
  names(cliq_sizes) = 1:length(cliq_sizes)
  # reorder the names based on the decreasing order of the size
  indices_ordered_by_size = as.numeric(names(sort(cliq_sizes,decreasing = TRUE)))
  # reorder the clique by new order
  cliq_ordered_by_size = cliq[indices_ordered_by_size]
  return(cliq_ordered_by_size)
}

# INPUT: (1) igraph object (2) a list of cliq (e.g return from max.cliq(..)) 
#        (3) whether the clique list object is ordered by size  
# OUTPUT: a vector whose each names is a vertex in the graph and each corresponding
#         values are the cliq the vertex belong to. 
# NOTE: in the cliq list, a vertex can belong to multiple cliques. However, 
#       in the returned vector, each vertex only belongs to the its maximal clique 
getCliqMemebership = function(igraphObject = NULL, 
                              cliq = NULL, 
                              isCliqOrdered = FALSE){
  if (!isCliqOrdered){
    cliq = reorderCliqBySize(igraphObject = igraphObject, cliq = cliq)
  }
  
  membership_clique = numeric(length = length(V(igraphObject)))
  names(membership_clique) = V(igraphObject)$name
  
  c = 1
  
  for (i in 1:length(cliq)) {
    # a clique
    cCliq = cliq[[i]]
    # vertex in clique
    idx = V(igraphObject)$name[cCliq]
    # if none of the members of this clique are already in a clique
    if (sum(membership_clique[idx]) == 0 ) {
      # the members of this clique are in group c
      membership_clique[idx] = c
      c = c + 1
    }
  }
  
  #anything not in a clique is in it's own group
  max_membership = max(membership_clique)+1
  for (i in 1:length(V(igraphObject))){
    if (membership_clique[i] == 0){
      membership_clique[i] = max_membership
      max_membership = max_membership + 1
    }
  }
  return(membership_clique)
}


# get total number of edges in a community which are not intra-community edges
# this code takes really long time to run
getNumFrontierEdges = function(igraph_object, list_of_subgraphs, vertices_degree){
  names(vertices_degree) = names(V(igraph_object))
  sapply(list_of_subgraphs)
  
  subgraph = list_of_subgraphs[[1]]
  for (i in 1:vcount(subgraph)) { 
    intra_edges = E(subgraph) [ from(V(subgraph)[i]) ]
    idx = match(V(subgraph)[i]$name, V(igraph_object)$name)
    all_edges = E(igraph_object) [ idx ]
    num_inter_edges = length(all_edges) - length(intra_edges)
  }
  
}

#### ---------------- NETWORK STATISTIC-RELATED FUNCTIONS -------------------------------------------

# To run this function, call getAllSubgraphs(..) to obtain list_of_subgraphs
getSubgraphStat = function(list_of_subgraphs, minVertexCount = 3){
  module_vertex_count = sapply(list_of_subgraphs, function(subg) vcount(subg) )
  module_edge_count = sapply(list_of_subgraphs, function(subg) ecount(subg))
  module_densities = sapply(list_of_subgraphs, function(subg) graph.density(subg))
  module_clusCoeff = sapply(list_of_subgraphs, function(subg) transitivity(subg, type = "global"))
  # module_simmilarity = sapply(list_of_subgraphs, function(subg) mean(similarity(subg, method = "jaccard"))/2)
  
  # stat = data.frame(vertex_count = module_vertex_count, edge_count = module_edge_count,
  #                   density = module_densities, clusCoef = module_clusCoeff,
  #                   avg_simmiliarity =  module_simmilarity)   
  
  stat = data.frame(vertex_count = module_vertex_count, edge_count = module_edge_count,
                    density = module_densities, clusCoef = module_clusCoeff)  
  
  rownames(stat) = paste("M",1:nrow(stat),sep = "")
  
  # if(!allowMissingValues){
  #   stat = stat[which(!is.na(stat$clusCoef)),]
  # }
  stat = stat[which(stat$vertex_count >= minVertexCount),]
  return(stat)
}

getSubgraphStatNoMissingValues = function(igraphObject = NULL, stat = NULL){
  if (!is.null(igraphObject)){
    stat = getSubgraphStat(igraphObject)
  }
  stat = stat[which(!is.na(stat$clusCoef)),]
  return(stat)
}

# INPUT: either stat_dataframe OR list_of_subgraphs
# OUTPUT: vector of weighted module statistic
getCommunityWeightedStat = function(stat_dataframe = NULL, list_of_subgraphs = NULL, communityObject = NULL){
  if ((is.null(stat_dataframe) & is.null(list_of_subgraphs)) || 
      (!is.null(stat_dataframe) & !is.null(list_of_subgraphs))){
    print("either stat_dataframe or clique list_of_subgraphs be null, and only 1 of them")
    break;
  }
  # if provided list_of_subgraphs
  if (!is.null(list_of_subgraphs)){
    stat_dataframe = getSubgraphStat(list_of_subgraphs = list_of_subgraphs, 
                                     minVertexCount = 3)
  }
  # else if provided stat_dataframe
  if (!is.null(stat_dataframe)){
    vertex_weighted = stat_dataframe$vertex_count/sum(stat_dataframe$vertex_count)
    weighted_density = sum(vertex_weighted*stat_dataframe$density)
    weighted_clusCoef = sum(vertex_weighted*stat_dataframe$clusCoef) 
    weighted_aveSimmiliarity = sum(vertex_weighted*stat_dataframe$avg_simmiliarity) 
  }
  # if provided communityObject
  if (!is.null(communityObject)){
    weighted_stat = c(weighted_density, weighted_clusCoef, 
                      weighted_aveSimmiliarity, modularity(communityObject))
  }else{
    weighted_stat = c(weighted_density, weighted_clusCoef,
                     weighted_aveSimmiliarity, NA)
  }
  
  return(weighted_stat)
}

getSubgraphStatSummary = function(stat_dataframe){
  stat_dataframe = stat_dataframe[which(!is.na(stat_dataframe$clusCoef)),]
  
  summary = rbind(summary(stat_dataframe$vertex_count), summary(stat_dataframe$edge_count),
                  summary(stat_dataframe$density), summary(stat_dataframe$clusCoef),
                  summary(stat_dataframe$avg_simmiliarity))
  rownames(summary) = c("vertex_count", "edge_count", "density",
                        "clusCoef", "avg_simmiliarity")
  return(summary)
}

getModuleMatrix = function(igraph_object, community_membership){
  # get total number of vertices
  size = length(V(igraph_object)$name)
  # initialize the square module matrix to be returned
  module_matrix = matrix(FALSE, nrow = size, ncol = size)
  diag(module_matrix) = TRUE
  rownames(module_matrix) = V(igraph_object)$name
  colnames(module_matrix) = rownames(module_matrix)
  # foreach vertice (row), get the module it belongs to, and search for other vertices
  # belonging to the same module, then set the cell = TRUE in such cases
  for (row_index in rownames(module_matrix)){
    assigned_module = community_membership[row_index]
    vertices_in_same_module =  names(which(community_membership == assigned_module))
    module_matrix[row_index, vertices_in_same_module] = TRUE
  }
  return(module_matrix)
}

#### ---------------- OLD FUNCTIONS -------------------------------------------
getModuleDensity = function(igraphObject = NULL, communityObject = NULL){
  module_density = sapply(sort(unique(membership(communityObject))), function(assigned_cluster) {
    subg <-induced.subgraph(igraphObject, which(membership(communityObject)==assigned_cluster)) #membership id differs for each cluster
    #ecount(subg)/ecount(g)
    graph.density(subg)
  })
  return(module_density)
}

getModuleClusterCoef = function(igraphObject = NULL, communityObject = NULL){
  module_clusterCoef = sapply(sort(unique(membership(communityObject))), function(assigned_cluster) {
    subg <-induced.subgraph(igraphObject, which(membership(communityObject)==assigned_cluster)) #membership id differs for each cluster
    transitivity(subg, type = "global")
  })
  return(module_clusterCoef)
}

getModuleSimmilarity= function(igraphObject = NULL, communityObject = NULL, method = "jaccard"){
  module_similarity = sapply(sort(unique(membership(communityObject))), function(assigned_cluster) {
    subg <-induced.subgraph(igraphObject, which(membership(communityObject)==assigned_cluster)) #membership id differs for each cluster
    similarity(subg, method = method)
  })
  return(module_similarity)
}

getModuleVertexSize= function(igraphObject = NULL, communityObject = NULL){
  module_similarity = sapply(sort(unique(membership(communityObject))), function(assigned_cluster) {
    subg <-induced.subgraph(igraphObject, which(membership(communityObject)==assigned_cluster)) #membership id differs for each cluster
    vcount(subg)
  })
  return(module_similarity)
}






























