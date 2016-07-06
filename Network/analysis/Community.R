#Detect community based on modularity
fc <- edge.betweenness.community (g)
memb <- community.to.membership(g,fc$merges,which.max(fc$modularity))
lay <- layout.auto(g)
jg3 <- graph.empty(n=vcount(g))
#colbar <- rainbow(5)
col[V(g)$Color==1]='palevioletred3'
col[V(g)$Color==2]='lightblue'
col[is.na(col)] <- "yellow"
Size=5*V(g)$degree/max(V(g)$degree)
Vertex.label.cex=0.5*V(g)$degree/max(V(g)$degree)
par(bg = 'white')
plot(g, layout=lay, vertex.size=Size,asp=FALSE,
     vertex.label.cex=Vertex.label.cex,
     vertex.label.color='black',vertex.label.family='Arials',
     vertex.color=col,vertex.frame.color=col,
     mark.groups=V(g)[(memb$membership==3)],
     vertex.label.degree=pi/2,
     edge.color='black',
     edge.width=0.5)

Memb0=V(g)$name[memb$membership==0]
write.table(Memb0,file='Sub/Female_0',row.names=F,col.names=F,quote=F)
