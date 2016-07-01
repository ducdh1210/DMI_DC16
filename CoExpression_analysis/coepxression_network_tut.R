# http://khughitt.github.io/2016-iscb-dc-rsg-workshop-presentation/#77
# https://github.com/iscb-dc-rsg/2016-summer-workshop/blob/master/3B-Hughitt-RNASeq-Coex-Network-Analysis/tutorial/README.md
setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge/CoExpression_analysis")

########## some unrelated simulation ##########################################
gene_per_cluster = 15
num_times = 10

nvals = gene_per_cluster * num_times

# highly co-expressed cluster; low-high expression

cluster1 = matrix(rep(1:num_times, gene_per_cluster) + rnorm(nvals, sd = 0.25),
                  nrow = gene_per_cluster, byrow = TRUE)
cluster2 = matrix(rep(num_times:1, gene_per_cluster) + rnorm(nvals, sd = 0.75),
                  nrow = gene_per_cluster, byrow = TRUE)
cluster3 = matrix(sample(1:2, nvals, replace = TRUE) + rnorm(nvals, sd = 1),
                  byrow = TRUE)


########## tutorial starting here ##########################################

library('gplots')
library('ggplot2')
library('knitr')
library('limma')
library('reshape2')
library('RColorBrewer')
library('WGCNA')

# Make sure results are reproducible
set.seed(1)

samples = read.csv('sample_metadata.csv')
kable(samples)

# load RNA-Seq raw count
raw_counts = read.csv('count_table.csv', row.names = 1)
head(raw_counts)

dim(raw_counts) # [1] 20956    19

# For gene annotations, we can use the Bioconductor Homo.sapiens OrganismDb package. 
# This meta-package combines several Human-specific annotation packages, providing gene- and transcript-level details.

library('Homo.sapiens')

# OrganismDb packages can be queried in a manner similar to querying a database. 
# You have to specify one or more gene identifiers ('keys'), 
# along with the type of the identifier ('key type'), 
# and one or more fields that you are interested in querying.

keytypes(Homo.sapiens)
columns(Homo.sapiens)

# To query the package, you use the select() function, e.g.:
genids = head(keys(Homo.sapiens, keytype='ENSEMBL'), 2)
select(Homo.sapiens, keytype='ENSEMBL', keys=genids, 
       columns=c('ALIAS', 'TXCHROM', 'TXSTART', 'TXEND'))

####  Sample check

# First, it is always wise to check the quality of your samples before continuing with an analysis
# like this. There are many such checks that one can (and should) perform, starting at the 
# level of read quality (e.g. FastQC).

# Here, we will just do a quick check using a sample-correlation heatmap.

# add a colorbar along the heatmap with sample condition
num_conditions <- nlevels(samples$condition)
pal <- colorRampPalette(brewer.pal(num_conditions, "Set1"))(num_conditions)
cond_colors <- pal[as.integer(samples$condition)]

raw_counts[] = lapply(raw_counts,as.numeric)

heatmap.2(cor(raw_counts), RowSideColors=cond_colors,
          trace='none', main='Sample correlations (raw)')


### Low count filtering

# Now that we are satisfied with the quality of our samples, 
# we will want to filter our data to keep only the genes 
# which will be informative during differential expression analysis 
# and network construction.

# Remove all rows with less than n counts across all samples, where n=#samples
low_count_mask <- rowSums(raw_counts) < ncol(raw_counts)

sprintf("Removing %d low-count genes (%d remaining).", sum(low_count_mask), 
        sum(!low_count_mask))

dim(raw_counts)
raw_counts = raw_counts[which(!low_count_mask),]

### Log2 transformation

# Most of the methods developed for co-expression network analysis and network inference were written for use with microarray data, including WGCNA!
  
#  Attempting to apply a method such as this to discrete-count RNA-Seq data will not work out well.

# There are a number of methods for working around this, in effect, making RNA-Seq data "look" more like microarray data, but the simplest thing is just to log the data. This will transform our discrete, over-dispersed counts to a more Poisson-like continuous distribution.


log_counts <- log2(raw_counts + 1)

x = melt(as.matrix(log_counts))

colnames(x) = c('gene_id', 'sample', 'value')
ggplot(x, aes(x=value, color=sample)) + geom_density()

heatmap.2(cor(log_counts), RowSideColors = cond_colors,
          trace = 'none', main = 'Sample correlations (log2-transformed')


# As you can see, after the low-count filtering and log-transformation, the samples within each condition are starting to behave better.

#### Remove non differentially-expressed genes

# Next, we will perform a series of differential expression contrasts, and use the results to further filter out genes for which there is not a significant amount of variance.

# If you were just interested in performing differential expression analysis, this may not be the most appropriate approach. In this case, you may want to consider additional steps such as quantile normalization and/or mean-variance adjustment with voom.

# first, let's remove any genes with _zero_ variance since these are not
# going to help us, and may cause problems with some of the models
log_counts <- log_counts[apply(log_counts, 1, var) > 0,]

# create design matrix for differential expression analysis;
# if you wanted to account for batch here, you could simply include a batch
# term in the linear model at this step, e.g.:
# mod <- model.matrix(~0+samples$condition+samples$batch)

mod <- model.matrix(~0+samples$condition)

# make model terms easier to work with
colnames(mod) <- levels(samples$condition)

fit <- lmFit(log_counts, design=mod)

# generate a list of all possible pairwise contrasts
condition_pairs <- t(combn(levels(samples$condition), 2))                                                                                                                               

comparisons <- list()                                                                                                                                          
for (i in 1:nrow(condition_pairs)) {                                                                                                                                     
  comparisons[[i]] <- as.character(condition_pairs[i,])                                                                                                      
}    

# vector to store de genes
sig_genes <- c()

# iterate over the contrasts, and perform a differential expression test for
# each pair
for (conds in comparisons) {
  # generate string contrast formula, "infLM24 - infLM4"
  contrast_formula <- paste(conds, collapse=' - ')
  
  contrast_mat <- makeContrasts(contrasts=contrast_formula, levels=mod)
  contrast_fit <- contrasts.fit(fit, contrast_mat)
  eb <- eBayes(contrast_fit)
  
  # Grab highly ranked genes; this is a pretty stringent p-value cutoff, but
  # it serves to limit the total number of genes we will use for this
  # tutorial
  sig_genes <- union(sig_genes, 
                     rownames(topTable(eb, number=Inf, p.value=0.005)))
}

# Filter out genes which were not differentially expressed for any contrast
log_counts <- log_counts[rownames(log_counts) %in% sig_genes,]

### Co-expression network construction

#'
#' Similarity measure which combines elements from Pearson correlation and
#' Euclidean distance.
#' 
cordist <- function(dat) {
  cor_matrix  <- cor(t(dat))
  
  dist_matrix <- as.matrix(dist(dat, diag=TRUE, upper=TRUE))
  dist_matrix <- log1p(dist_matrix)
  dist_matrix <- 1 - (dist_matrix / max(dist_matrix))
  
  sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix)/ 2)
}

sim_matrix <- cordist(log_counts)

# Let's see what our similarity matrix looks like at this point. 
# Because the heatmap.2 function (which includes a biclustering step) can be
# pretty slow, we will use a sub-sample of our data --
# for visualization purposes this is fine.

heatmap_indices <- sample(nrow(sim_matrix), 500)

heatmap.2(t(sim_matrix[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          labRow=NA, labCol=NA, 
          trace='none', dendrogram='row',
          xlab='Gene', ylab='Gene',
          main='Similarity matrix',
          density.info='none', revC=TRUE)

### Construct adjacency matrix

adj_matrix <- adjacency.fromSimilarity(sim_matrix, power=12, type='signed')

# Delete similarity matrix to free up memory
rm(sim_matrix)
gc()

# Convert to matrix
gene_ids <- rownames(adj_matrix)

adj_matrix <- matrix(adj_matrix, nrow=nrow(adj_matrix))
rownames(adj_matrix) <- gene_ids
colnames(adj_matrix) <- gene_ids

# same plot as before
heatmap.2(t(adj_matrix[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          labRow=NA, labCol=NA, 
          trace='none', dendrogram='row',
          xlab='Gene', ylab='Gene',
          main='Adjacency matrix',
          density.info='none', revC=TRUE)

### Co-expression module detection

# Cluster gene expression profiles; the flashClust function from
# the authors of WGCNA is another options for larger datasets.
# For input, we use the reciprocal of the adjacency matrix; hierarchical
# clustering works by comparing the _distance_ between objects instead of the
# _similarity_.
gene_tree <- hclust(as.dist(1 - adj_matrix), method="average")

# we will use the cuttreeDynamicTree method to break apart the hc dendrogram
# into separate modules
module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=15,
                                   deepSplit=TRUE)

# assign a color to each module for easier visualization and referencing
module_colors <- labels2colors(module_labels)

### Exporting network

export_network_to_graphml <- function (adj_mat, filename=NULL, weighted=TRUE,
                                       threshold=0.5, max_edge_ratio=3,
                                       nodeAttr=NULL, nodeAttrDataFrame=NULL,
                                       edgeAttributes=NULL, verbose=FALSE) {
  library('igraph')
  
  # Determine filename to use
  if (is.null(filename)) {
    filename='network.graphml'
  }
  
  # TODO 2015/04/09
  # Add option to rescale correlations for each module before applying
  # threshold (this is simpler than the previous approach of trying to
  # determine a different threshold for each module)
  #
  # Still, modules with very low correlations should be given somewhat
  # less priority than those with very high correlations.
  
  #module_colors <- unique(nodeAttrDataFrame$color)
  #module_genes <- which(nodeAttrDataFrame$color == color)
  #module_adjmat <- adj_mat[module_genes,]
  #num_genes <- length(module_genes)
  
  # Adjust threshold if needed to limit remaining edges
  max_edges <- max_edge_ratio * nrow(adj_mat)
  
  edge_to_total_ratio <- max_edges / length(adj_mat)
  edge_limit_cutoff <- as.numeric(quantile(abs(adj_mat), 1 - edge_to_total_ratio))
  
  # Also choose a minimum threshold to make sure that at least some edges
  # are left
  min_threshold <- as.numeric(quantile(abs(adj_mat), 0.9999))
  
  threshold <- min(min_threshold, max(threshold, edge_limit_cutoff))
  
  # Remove edges with weights lower than the cutoff
  adj_mat[abs(adj_mat) < threshold] <- 0
  
  # Drop any genes with no edges (TODO: Make optional)
  orphaned <- (colSums(adj_mat) == 0)
  adj_mat <- adj_mat[!orphaned, !orphaned]
  
  # Also remove annotation entries
  if (!is.null(nodeAttr)) {
    nodeAttr <- nodeAttr[!orphaned]
  }
  if (!is.null(nodeAttrDataFrame)) {
    nodeAttrDataFrame <- nodeAttrDataFrame[!orphaned,]
  }
  
  # Keep track of non-positive edges and rescale to range 0,1
  is_zero     <- adj_mat == 0
  is_negative <- adj_mat < 0
  
  adj_mat <- (abs(adj_mat) - threshold) / (max(adj_mat) - threshold)
  adj_mat[is_zero] <- 0
  adj_mat[is_negative] <- -adj_mat[is_negative]
  
  if (verbose) {
    message(sprintf("Outputting matrix with %d nodes and %d edges", 
                    nrow(adj_mat), sum(adj_mat > 0)))
  }
  
  # Create a new graph and add vertices
  # Weighted graph
  if (weighted) {
    g <- graph.adjacency(adj_mat, mode='undirected', weighted=TRUE, diag=FALSE)
  } else {
    adj_mat[adj_mat != 0] <- 1
    g <- graph.adjacency(adj_mat, mode='undirected', diag=FALSE)
  }
  
  # Add single node annotation from vector
  if (!is.null(nodeAttr)) {
    g <- set.vertex.attribute(g, "attr", value=nodeAttr)
  }
  
  # Add node one or more node annotations from a data frame
  if (!is.null(nodeAttrDataFrame)) {
    for (colname in colnames(nodeAttrDataFrame)) {
      g <- set.vertex.attribute(g, colname, value=nodeAttrDataFrame[,colname])
    }
  }
  
  edge_correlation_negative <- c()
  
  # neg_correlations[edge_list]
  edge_list <- get.edgelist(g)
  
  for (i in 1:nrow(edge_list)) {
    from <- edge_list[i, 1]    
    to   <- edge_list[i, 2]    
  }
  
  # Save graph to a file
  write.graph(g, filename, format='graphml')
  
  # return igraph
  return(g)
}

# use OrganismDb to retrieve gene annotations
gene_info <- select(Homo.sapiens, keytype='ENSEMBL', keys=rownames(log_counts),
                    columns=c('TXCHROM', 'TXSTRAND', 'GENENAME'))

colnames(gene_info) <- c('gene_id', 'description', 'chr', 'strand')

# for now, just grab the description for the first transcript
gene_info <- gene_info[!duplicated(gene_info$gene_id),]

gene_info <- cbind(gene_info, module=module_colors)

# Include RGB versions of module colors for better assignment in Cytoscape
gene_info$color_rgb <- col2hex(gene_info$module)

# first, it's a good idea to check the distribution of edges weights in our
# correlation matrix. This will help us choose a reasonable cutoff for
# exporting the network.
g <- export_network_to_graphml(adj_matrix, filename='~/network.graphml',
                               threshold=0.4, nodeAttrDataFrame=gene_info)

class(g) # [1] "igraph" ==> checked!





















































