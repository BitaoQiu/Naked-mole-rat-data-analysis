g=read.graph(file='~/Documents/MyLocal/NakedMole/FManalysi/Network/PlotingNet/Net_female.gml',format='gml')
g=read.graph(file='~/Documents/MyLocal/NakedMole/FManalysi/Network/PlotingNet/Net_male.gml',format='gml')
g=as.undirected(g)
g=simplify(g)
fc <- edge.betweenness.community(g)
memb <- community.to.membership(g,fc$merges,which.max(fc$modularity))
V(g)$memb=memb$membership
V(g)[V(g)$memb==V(g)['GNRH1']$memb]
g2=induced.subgraph(g,V(g)[V(g)$memb==V(g)['GNRH1']$memb])
lay <- layout.fruchterman.reingold(g2)
V(g2)$degree=degree(g2)
V(g2)$between=betweenness(graph=g2)
col=rep('yellow',length(V(g2)))
col[V(g2)$Color==1]='#F8B39B'
col[V(g2)$Color==2]='lightblue'
col[V(g2)$Color==3]='yellow'
Size=1+7*(V(g2)$degree)/max(V(g2)$degree)
Size=1+7*(V(g2)$between)/max(V(g2)$between)

plot(g2, layout=lay,
     vertex.size=Size,
     asp=FALSE,
     vertex.label.font=1,
     vertex.label.cex=0.5,
     vertex.color=col,
     vertex.frame.color=col,
     edge.curved=F,
     edge.color='black',
     edge.width=0.9)

