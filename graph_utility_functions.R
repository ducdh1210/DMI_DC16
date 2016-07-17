#### ---------------- UTILITY FUNCTIONS -------------------------------------------

require(igraph)
require(WGCNA)

#### ---------------- GET LIST OF SUBGRAPHS -------------------------------------------

# note that if use cliq, the cliq elements must be sorted by size
getAllSubgraphs = function(igraphObject, communityObject = NULL, clique = NULL, cliqueSorted = TRUE){
  list_of_all_subgraphs  = list();
  if ((is.null(communityObject) & is.null(clique)) || 
      (!is.null(communityObject) & !is.null(clique))){
    print("either communityObject or clique must be null, and only 1 of them")
    break;
  }else{
    if(!is.null(communityObject)){
      for (i in sort(unique(membership(communityObject)))){
        list_of_all_subgraphs[[i]] = induced.subgraph(igraphObject,
                                                      which(membership(communityObject)==i))
      }
    }else{
      for (i in 1:length(clique)){
        list_of_all_subgraphs[[i]] = induced.subgraph(igraphObject,
                                                      names(clique[[i]]))
      }
    }
  }
  # for (i in sort(unique(membership(communityObject)))){
  #     list_of_all_subgraphs[[i]] = induced.subgraph(igraphObject,
  #                                                   which(membership(communityObject)==i))
  # }
  return(list_of_all_subgraphs)
}


#### ---------------- CLIQUE FUNCTIONS -------------------------------------------


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

getCliqMemebership = function(igraphObject,cliq, isCliqOrdered = FALSE){
  
  if (!isCliqOrdered){
    cliq_sizes = sapply(cliq, function(vertices) length(vertices))
    # name the vector cliq_size by the cliq index order 
    names(cliq_sizes) = 1:length(cliq_sizes)
    # reorder the names based on the decreasing order of the size
    indices_ordered_by_size = as.numeric(names(sort(cliq_sizes,decreasing = TRUE)))
    # reorder the clique by new order
    cliq = cliq[indices_ordered_by_size]
  }
  
  membership_clique = numeric(length = length(V(igraphObject)))
  names(membership_clique) = V(igraphObject)$name
  
  c= 1
  
  for (i in 1:length(cliq)) {
    # a clique
    cCliq = cliq_ordered_by_size[[i]]
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

# get number of edges for each vertex in the graph, 
# INPUT: igraph object
# OUPUT: vector showing number of edges for each vertex
# note: taking long time to run

getVerticesDegree = function(igraphObject){
  for (i in 1:vcount(igraphObject)) { 
    edges = E(igraphObject) [ from(V(igraphObject)[i]) ]
    num_edges = append(num_edges, length(edges))
  }
  return(num_edges)
}

# OUTPUT: list of cliques whose mininum size specified by minCliqueSize
# Note: verticesDegree can be obtained by calling getVerticesDegree()
getMaxClique = function(igraphObject = NULL, 
                        minCliqueSize = 5, 
                        maxCliqueSize = NULL, 
                        verticesDegree = NULL,
                        degreeThreshold = 10)
{
  # if the vector verticesDegree is not provided, then calling getVerticesDegree()
  if(is.null(verticesDegree)){
    verticesDegree = getVerticesDegree(igraphObject)
  }
  # only get vertices whose degrees are higher than the provided threshold
  idx = which(verticesDegree > degreeThreshold)
  s = V(igraphObject)[idx]
  # obtain cliques with provided parameters
  cliq = max_cliques(pp2_igraph, 
                     min = minCliqueSize, 
                     max = maxCliqueSize, 
                     subset=s) 
  return(cliq)
}

#### ---------------- STAT-RELATED FUNCTIONS -------------------------------------------

# To run this function, call getAllSubgraphs(..) to obtain list_of_subgraphs
getSubgraphStat = function(list_of_subgraphs){
  module_vertex_count = sapply(list_of_subgraphs, function(subg) vcount(subg) )
  module_edge_count = sapply(list_of_subgraphs, function(subg) ecount(subg))
  module_densities = sapply(list_of_subgraphs, function(subg) graph.density(subg))
  module_clusCoeff = sapply(list_of_subgraphs, function(subg) transitivity(subg, type = "global"))
  module_simmilarity = sapply(list_of_subgraphs, function(subg) mean(similarity(subg, method = "jaccard"))/2)
  
  stat = data.frame(vertex_count = module_vertex_count, edge_count = module_edge_count,
                    density = module_densities, clusCoef = module_clusCoeff,
                    avg_simmiliarity =  module_simmilarity)   
  rownames(stat) = paste("M",1:nrow(stat),sep = "")
  return(stat)
}

getSubgraphStatSummary = function(stat, allowMissingValues = TRUE){
  # missing values coming from the cases where 
  if (!allowMissingValue){
    stat = stat[which(!is.na(stat$clusCoef)),]
  }
  summary = rbind(summary(stat$vertex_count), summary(stat$edge_count),
                  summary(stat$density), summary(stat$clusCoef),
                  summary(stat$avg_simmiliarity))
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






























