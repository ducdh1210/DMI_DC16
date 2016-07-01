library(ProNet)
nodes<-data.frame(c("1855","1856","1857"))
network<-construction(input=nodes,db="Biogrid",species="human",ID.type="Entrez Gene",hierarchy=1)

net1<-extraction(network, mode="sample", sample.number=20)
net2<-extraction(network, mode="exact", nodes=1:20)
net3<-assemble(net1, net2, mode="union")

# Functional modules can be achieved by clustering
cluster(network, method="MCODE", plot=TRUE, layout="fruchterman.reingold")

# --------------------- network construction and operation ---------------------------------------
iavPath <-file.path(system.file("example",package="ProNet"),"iav.txt")
iav <- read.table(iavPath, header=TRUE, sep="\t")
dim(iav) # 179   5

g1 <- construction(iav[,c("Gene_name_1","Gene_name_2")],local.net=TRUE)
sp <- unique(cbind(c(as.vector(iav[,"Gene_name_1"]),as.vector(iav[,"Gene_name_2"])),
                   c(as.vector(iav[,"Adscription_1"]),as.vector(iav[,"Adscription_2"]))))
dim(sp) # 114 2

summary(g1)

## Second, the non-local network between host proteins that interact with nodes having at least IAV protein
## partners was constructed based on the integrated Biogrid database.
hostPath <-file.path(system.file("example",package="ProNet"),"host.txt")
host <- read.table(hostPath, header=TRUE, sep="\t")
g2 <- construction(input=as.data.frame(unique(host[,"Protein.name"])),
                     hierarchy=1,db="HPRD",species="human",ID.type="Gene symbol")
summary(g2)

## The expanded IAV-host network was then abtained by integrating these two networks.
hprd <- construction(db="HPRD",ID.type= c("Gene symbol"))
id <- match(unique(c(V(g1)$name,V(g2)$name)),V(hprd)$name)
gtemp <- induced.subgraph(hprd, id[!is.na(id)])
g3 <- assemble(g1,gtemp,mode="union")
summary(g3)

# --------------------- network construction and operation ---------------------------------------

color <- rep(1,vcount(g3))
color[V(g3)$species=="DHP of IAV"] <- "red"
color[V(g3)$species=="IAV protein"] <- "black"
color[is.na(V(g3)$species)] <- "green"
visualization(g3,node.size=3,node.fill.color=color,node.label="",edge.color="gray")
legend("topleft",col=c("black","red","green"),
          legend=c("virus","human_direct","human_indirect"),pch=19)

# Or subcellular localization based visualization:

V(g3)$expression<-rexp(vcount(g3),1)
location(g3,species=c("human"),vertex.size=3,vertex.label.cex=0.5,
           vertex.color="expression",xlim=c(-1,1),ylim=c(-1,1))

# --------------------- topological analysis ---------------------------------------

topology(g3,simple.parameters=TRUE)
tp <- topology(g2,degree.distribution=TRUE)

head(as.data.frame(tp))

# Other topological parameters like clustering coefficient, betweeness, shortest path, eigenvector centrality,
# connectivity and closeness can be obtained similarly by changing the default setting of the parameters to be
# TRUE.

tp <- topology(g2,shortest.paths=TRUE)
head(as.data.frame(tp))


# --------------------- Topological comparison of networks ---------------------------------------
net.comparing(g3,hprd,topology.parameters=TRUE)


# ---------------------Network Clustering ---------------------------------------

# Several graph based network clustering algorithms were integrated into the package, 
# such as the FN (A Clauset et al., 2004), linkcomm (Kalinka et al., 2011), 
# MCL (van Dongen SM, 2000) and MCODE (Bader GD et al., 2003) methods. 

# There are 7 clusters found by the FN method, and the number of nodes in each cluster is also shown.

result <- cluster(g3, method="FN")
clusters <- rep(1, vcount(g3))
for(i in 1:vcount(g3)){
  clusters[i] <- result[[i]]
}

clusters <- as.factor(clusters)
table(clusters)

# MCODE method can be perfomed using the individual mcode module. 11 clusters were found, with the
# largest containing 77 elements. Scores of each cluster were also shown.

result <- mcode(g3,vwp=0.05,haircut=TRUE,fluff=FALSE,fdt=0.8,loops=FALSE)

summary(result$COMPLEX)

result$score

cluster1<-induced.subgraph(g3,result$COMPLEX[[1]])
summary(cluster1)

visualization(cluster1,node.size=4,node.label=V(cluster1)$name,node.label.color="blue")







