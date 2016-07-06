library(igraph)
library(ggplot2)
library(grid)
g=read.graph(file='~/Documents/MyLocal/NakedMole/FManalysi/Network/PlotingNet/Net_female.gml',format='gml')
g=read.graph(file='~/Documents/MyLocal/NakedMole/FManalysi/Network/PlotingNet/Net_male.gml',format='gml')
g=as.undirected(g)
g=simplify(g)
V(g)$degree=degree(g)
V(g)$between=betweenness(graph=g)
fc <- edge.betweenness.community(g)
memb <- community.to.membership(g,fc$merges,which.max(fc$modularity))
V(g)$memb=memb$membership
V(g)['GNRH1']$memb
NMR.eusociality.uniqReads.rpkm$Name[match(c(V(g)[V(g)$memb==V(g)['GNRH1']$memb]$name),NMR.eusociality.uniqReads.rpkm$Symbol)]
V(g)[V(g)$memb==V(g)['GNRH1']$memb]
V(g)[V(g)$memb==V(g)['GNRH1']$memb]$degree
ModuleA=memb$membership %in% c(V(g)['GNRH1']$memb)
mark.groups=V(g)[ModuleA]
E(g)[to('GNRH1')]
lay <- layout.fruchterman.reingold(g)
E(g)$color <- "darkgrey"
col=rep('yellow',length(V(g)))
col[V(g)$Color==1]='#F8B39B'
col[V(g)$Color==2]='lightblue'
col[V(g)$Color==3]='yellow'
col[is.na(col)] <- "yellow"
NS=rep('circle',length(V(g)))
Lable_Color='black'
par(bg = 'white')
par(ps=26)

Size=1+7*(V(g)$degree)/max(V(g)$degree)
Size=1+7*(V(g)$between)/max(V(g)$between)
Vertex.label.cex=0.1#(V(Sg)$degree)^1.9/max(V(Sg)$degree)^1.9+0.1
plot(g, layout=lay, vertex.size=Size,asp=FALSE,
     vertex.label.cex=Vertex.label.cex,
     vertex.label.color=Lable_Color,
     vertex.label.font=1,
     vertex.color=col,
     vertex.frame.color=col,
     vertex.shape=NS,
     edge.color=E(g)$color,
     edge.curved=F,
     edge.color='black',
     edge.width=0.9,
     mark.groups=mark.groups,
     mark.col=rainbow(length(mark.groups), alpha=0.05))
V(g)[order(V(g)$degree,decreasing=T)]$degree
V(g)[order(V(g)$between,decreasing=T)]$between
V(g)['NPY']$degree
V(g)['NPY']$between
V(g)['CRH']$degree
V(g)['CRH']$between

