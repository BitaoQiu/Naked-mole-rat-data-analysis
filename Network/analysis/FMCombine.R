library(igraph)
BF=read.csv('/home/enjoycloud/Documents/MyLocal/NakedMole/FManalysi/BF.D',sep='\t',header=T)
BM=read.csv('/home/enjoycloud/Documents/MyLocal/NakedMole/FManalysi/BM.txt',sep='\t',header=T)
#BF=BM
BF=toupper(BF[BF$X16.S!='-',])
BM=toupper(BM[BM$X16.S!='-',])
Mus_selc=Mus_action[(Mus_action$A1 %in% c(BF,BM)) & (Mus_action$A2 %in% c(BF,BM)),]
Mus_selc=Mus_selc[Mus_selc$A1!=Mus_selc$A2 & Mus_selc$combined_score>=400,]
Mus_selc2=Mus_action2[(Mus_action2$A1 %in% c(BF,BM)) & (Mus_action2$A2 %in% c(BF,BM)),]
Mus_selc2=Mus_selc2[Mus_selc2$mode=='binding' |
                      (Mus_selc2$mode!='binding' & Mus_selc2$a_is_acting==1),]
Mus_selc2=Mus_selc2[Mus_selc2$A1!=Mus_selc2$A2,]
colnames(Mus_selc)[1:5]=colnames(Mus_selc2[,c(1:4,6)])
Mus_selc2=Mus_selc2[rownames(Mus_selc2) %in% rownames(unique(Mus_selc2[,c(1,2)])),]
Mus_selc=Mus_selc[rownames(Mus_selc) %in% rownames(unique(Mus_selc[,c(1,2)])),]
useM=!( (paste(Mus_selc$A1,Mus_selc$A2) %in% paste(Mus_selc2$A1, Mus_selc2$A2)) | 
          (paste(Mus_selc$A2,Mus_selc$A1) %in% paste(Mus_selc2$A1, Mus_selc2$A2)))
Mus_selc3=rbind(Mus_selc[useM,c(1:5)],Mus_selc2[,c(1:4,6)])
#g3_M=graph.data.frame(Mus_selc3,directed=T)
g=graph.data.frame(Mus_selc3,directed=F)
V(g)$degree=degree(g)
clusters(g)
V(g)[order(V(g)$degree)]$degree
V(g)[order(V(g)$degree)]
length(E(g))
V(g)$color[V(g)$name %in% BF]='red'
V(g)$color[V(g)$name %in% BM]='blue'
Over <- read.table("~/Dropbox/Writing_NMR/POSTER/Brain/NETWORK/Over.txt", header=T, quote="\"")
Over$X16.Symbol=toupper(Over$X16.Symbol)
V(g)$color[V(g)$name %in% Over$X16.Symbol]='yellow'
V(g)['GNRH1']$color='green'
V(g)['LHB']$color='green'
g=simplify(g)
fc <- edge.betweenness.community (g)
memb <- community.to.membership(g,fc$merges,which.max(fc$modularity))
lay <- layout.fruchterman.reingold(g)
#jg3 <- graph.empty(n=vcount(g))
#colbar <- rainbow(5)
#col[V(g)$Color==1]='palevioletred3'
#col[V(g)$Color==2]='lightblue'
#col[is.na(col)] <- "yellow"
markgroup=list(V(g)[memb$membership==0],V(g)[memb$membership==1],V(g)[memb$membership==2],
               V(g)[memb$membership==3],V(g)[memb$membership==4])
Size=3
Vertex.label.cex=0.5*V(g)$degree/max(V(g)$degree)
par(bg = 'white')
plot(g, layout=lay, vertex.size=Size,asp=FALSE,
     vertex.label.cex=Vertex.label.cex,
     vertex.label.color='black',vertex.label.family='Arials',
     vertex.color=V(g)$color,
     edge.arrow.size=0,
     vertex.frame.color=V(g)$color,
     mark.groups=markgroup,
     vertex.label.degree=pi/2,
     edge.color='black',
     edge.width=0.5)

