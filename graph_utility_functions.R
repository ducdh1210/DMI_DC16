#### ---------------- UTILITY FUNCTIONS -------------------------------------------

require(igraph)
require(WGCNA)

# 
getAllSubgraphs = function(igraphObject, communityObject){
  list_of_all_subgraphs  = list();
  for (i in sort(unique(membership(communityObject)))){
      list_of_all_subgraphs[[i]] = induced.subgraph(igraphObject,
                                                    which(membership(communityObject)==i))
  }
  return(list_of_all_subgraphs)
}

getSubgraphStat = function(list_of_subgraphs){
  module_vertex_count = sapply(list_of_subgraphs, function(subg) ecount(subg) )
  module_edge_count = sapply(list_of_subgraphs, function(subg) vcount(subg))
  module_densities = sapply(list_of_subgraphs, function(subg) graph.density(subg))
  module_clusCoeff = sapply(list_of_subgraphs, function(subg) transitivity(subg, type = "global"))
  module_simmilarity = sapply(list_of_subgraphs, function(subg) mean(similarity(subg, method = "jaccard"))/2)
  
  stat = data.frame(vertex_count = module_vertex_count, edge_count = module_edge_count,
                    density = module_densities, clusCoef = module_clusCoeff,
                    avg_simmiliarity =  module_simmilarity)                   
  return(stat)
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

