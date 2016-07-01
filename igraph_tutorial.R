setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge")

######################################################################################
# http://www.r-bloggers.com/network-visualization-in-r-with-the-igraph-package/

require(igraph)
bsk<-read.table("http://www.dimiter.eu/Data_files/edgesdata3.txt", sep='\t', dec=',', header=T)#specify the path, separator(tab, comma, ...), decimal point symbol, etc.


# Transform the table into the required graph format:
bsk.network<-graph.data.frame(bsk, directed=F,) #the 'directed' attribute specifies whether the edges are directed
# or equivelent irrespective of the position (1st vs 2nd column). For directed graphs use 'directed=T'

class(bsk.network)

# Inspect the data:

V(bsk.network) #prints the list of vertices (people)
E(bsk.network) #prints the list of edges (relationships)
degree(bsk.network) #print the number of edges per vertex (relationships per people)

# First try. We can plot the graph right away but the results will usually be unsatisfactory:
plot(bsk.network)

#Subset the data. If we want to exclude people who are in the network only tangentially (participate in one or two relationships only)
# we can exclude the by subsetting the graph on the basis of the 'degree':

bad.vs<-V(bsk.network)[degree(bsk.network)<3] #identify those vertices part of less than three edges
bsk.network<-delete.vertices(bsk.network, bad.vs) #exclude them from the graph

# Plot the data.Some details about the graph can be specified in advance.
# For example we can separate some vertices (people) by color:

V(bsk.network)$color<-ifelse(V(bsk.network)$name=='CA', 'blue', 'red') #useful for highlighting certain people. Works by matching the name attribute of the vertex to the one specified in the 'ifelse' expression

# We can also color the connecting edges differently depending on the 'grade': 

E(bsk.network)$color<-ifelse(E(bsk.network)$grade==9, "red", "grey")

# or depending on the different specialization ('spec'):

E(bsk.network)$color<-ifelse(E(bsk.network)$spec=='X', "red", ifelse(E(bsk.network)$spec=='Y', "blue", "grey"))

# Note: the example uses nested ifelse expressions which is in general a bad idea but does the job in this case
# Additional attributes like size can be further specified in an analogous manner, either in advance or when the plot function is called:

V(bsk.network)$size<-degree(bsk.network)/10#here the size of the vertices is specified by the degree of the vertex, so that people supervising more have get proportionally bigger dots. Getting the right scale gets some playing around with the parameters of the scale function (from the 'base' package)

# Note that if the same attribute is specified beforehand and inside the function, the former will be overridden.
# And finally the plot itself:
par(mai=c(0,0,1,0)) 			#this specifies the size of the margins. the default settings leave too much free space on all sides (if no axes are printed)
plot(bsk.network,				#the graph to be plotted
     layout=layout.fruchterman.reingold,	# the layout method. see the igraph documentation for details
     main='Organizational network example',	#specifies the title
     vertex.label.dist=0.5,			#puts the name labels slightly off the dots
     vertex.frame.color='blue', 		#the color of the border of the dots 
     vertex.label.color='black',		#the color of the name labels
     vertex.label.font=2,			#the font of the name labels
     vertex.label=V(bsk.network)$name,		#specifies the lables of the vertices. in this case the 'name' attribute is used
     vertex.label.cex=1			#specifies the size of the font of the labels. can also be made to vary
)

# Save and export the plot. The plot can be copied as a metafile to the clipboard, or it can be saved as a pdf or png (and other formats).
# For example, we can save it as a png:
png(filename="org_network.png", height=800, width=600) #call the png writer
#run the plot
dev.off() #dont forget to close the device
#And that's the end for now.

#########################################################################################################
## http://www.shizukalab.com/toolkits/sna/weighted-edgelists

el=read.csv(file.choose()) # read  the 'el.with.weights.csv' file
el[,1]=as.character(el[,1]) #Because the vertex IDs in this dataset are numbers, we make sure igraph knows these should be treated as characters. Otherwise, it'll create problems (see page on data import)
el[,2]=as.character(el[,2])
el=as.matrix(el) #igraph needs the edgelist to be in matrix format
g=graph.edgelist(el[,1:2]) #We first greate a network from the first two columns, which has the list of vertices
E(g)$weight=as.numeric(el[,3]) #We then add the edge weights to this network by assigning an edge attribute called 'weight'.
plot.igraph(g,edge.label = E(g)$weight)

el=read.csv(file.choose()) # read the 'el.with.weights.csv' file 
g2=graph.data.frame(el,directed = F)
plot.igraph(g2, edge.label = E(g2)$weight)

#########################################################################################################
library(igraph)
load(".RData")
class(ppi_1); head(ppi_1)

working_ppi1 = ppi_1
colnames(working_ppi1) = c("V1","V2","weight"); head(working_ppi1)
pp1_graphObj = graph.data.frame(working_ppi1, directed = F)

E(pp1_graphObj) #2232404/2232404 edges 
V(pp1_graphObj) #17397/17397 vertices
head(E(pp1_graphObj)$weight)

plot.igraph(pp1_graphObj) # takes long time
groups <- membership(cluster_louvain(pp1_graphObj))
communities <- communities(cluster_louvain(pp1_graphObj))



class(signal_directed)
working_signal_directed = signal_directed
colnames(working_signal_directed) = c("V1","V2","weight"); head(working_signal_directed)
signal_directed_graphObj = graph.data.frame(working_signal_directed, directed = T)

E(signal_directed_graphObj) #21825
V(signal_directed_graphObj) #5254
head(E(signal_directed_graphObj)$weight)

plot.igraph(signal_directed_graphObj) # takes long time, runable!



