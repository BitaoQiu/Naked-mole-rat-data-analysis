library(igraph)
library(ggplot2)
library(grid)
g=read.graph(file='PlotingNet/Net_female.gml',format='gml')
g=read.graph(file='PlotingNet/Net_male.gml',format='gml')
#lay=as.matrix(read.table("Female_net.lay"))
#lay=as.matrix(read.table("Male_net.lay"))
g=simplify(g,remove.multiple=F)
fc <- edge.betweenness.community (g)
memb <- community.to.membership(g,fc$merges,which.max(fc$modularity))
lay <- layout.auto(g)#fruchterman.reingold(g,niter=2000)
E(g)$color <- "darkgrey"
E(g)[ mode!="0" ]$color <- "purple"
E(g)[ action=="activation" ]$color <- "red"
E(g)[ action=="inhibition" ]$color <- "blue"
E(g)[ mode=="binding" ]$color <- "black"
#E(g)[159]$color='purple'
col=rep('yellow',length(V(g)))
col[V(g)$Color==1]='#F8B39B'
col[V(g)$Color==2]='lightblue'
col[is.na(col)] <- "yellow"
Sg=simplify(g,remove.multiple=T,remove.loops=T)
V(Sg)$degree=degree(Sg)
AT=rep(0,length(E(g)))
AT[E(g)$color!='darkgrey' & E(g)$color!='black']=2
AT[E(g)$color=='purple']=0
LTY=rep(5,length(E(g)))
LTY[E(g)$color!='darkgrey' ]=1
NS=rep('circle',length(V(g)))
Lable_Color='black'
#Lable_Color[V(g)$name %in% Over$X16.Symbol]='black'
#col[V(g)$name %in% Over$X16.Symbol]='green'
#NS[V(g)$name=='GNRH1']='csquare'
par(bg = 'white')
ModuleA=memb$membership %in% c(11)
ModuleB=memb$membership %in% c(1)#0,11)
mark.groups=list(V(g)[ModuleA],V(g)[ModuleB])
#Module=memb$membership %in% c(0,3)
#mark.groups=list(V(g)[Module])
par(ps=26)
Size=2+5*(V(Sg)$degree)/max(V(Sg)$degree)
Vertex.label.cex=(V(Sg)$degree)^1.9/max(V(Sg)$degree)^1.9+0.1
#Vertex.label.cex[V(g)$name=='NPY']=0.5
#Vertex.label.cex[V(g)$name=='GNRH1']=1
#pdf('Male_Net.pdf',width=10,height=8)
#pdf('Female_Net.pdf',width=10,height=8)
plot(g, layout=lay, vertex.size=Size,asp=FALSE,
     vertex.label.cex=Vertex.label.cex,
     vertex.label.color=Lable_Color,
     vertex.label.font=1,
     vertex.color=col,
     vertex.frame.color=col,
     vertex.shape=NS,
     edge.arrow.size=0.5,
     edge.arrow.width=1.2,
     edge.arrow.mode=AT,
     edge.color=E(g)$color,
     edge.lty=LTY,
     edge.curved=F,
     mark.groups=mark.groups,
     mark.col=rainbow(length(mark.groups), alpha=0.05),
     edge.color='black',
     edge.width=0.9)
#dev.off()
#write.table(V(g)$name[Module],file='Male_Module_Other',quote=F,row.names=F,col.names=F)
Name='SST'
#lay[V(g)$name=='NPY']=c(lay[V(g)$name=='NPY',1]+2,lay[V(g)$name=='NPY',2]-2)
#lay[V(g)$name=='GNRH1']=c(lay[V(g)$name=='GNRH1',1]+2,lay[V(g)$name=='GNRH1',2]+2)
#lay[V(g)$name==Name]=c(lay[V(g)$name==Name,1]-2,lay[V(g)$name==Name,2]-4)


legend('bottomright',c('Up regulated in Breeders',"Down Regulated in Breeders",'GnRH',
                       'Activation',"Inhibition","Expression",
                       'Binding','Functional Interaction'),
       col=c("#F8B39B","lightblue",'yellow','red','blue','purple','black','darkgrey'),
       lty=c(NA,NA,NA,1,1,1,1,5),
       lwd=c(NA,NA,NA,rep(2,5)),
       pch=c(rep(19,3),rep(NA,5)),
       cex=0.5,pt.cex=1)
TitleName='Male'
title(main=paste('Functional Interaction Network:\nBrain DEGs in ',TitleName,sep=''),cex=0.9,font=2)
#dev.off()
#write.table(lay,"Female_net.lay",quote=F)
#write.table(lay,"Male_net.lay",quote=F)
#Memb0=V(g)$name[memb$membership==0]
#write.table(Memb0,file='Sub/Female_0',row.names=F,col.names=F,quote=F)