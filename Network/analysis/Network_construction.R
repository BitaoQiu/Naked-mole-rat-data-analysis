library(igraph)
## Load the mouse gene interaction scores from strings database.
action=read.csv('Link_score.txt',header=T,sep=' ')
#Load annotation between NMR gene ID and string gene ID.
Mus.PS=read.csv('anno.str',header=T,sep='\t')
Mus.PS$Gene=toupper(Mus.PS$V3)
Mus.PS=Mus.PS[!is.na(Mus.PS$Gene),]
probe=match(action$protein1,Mus.PS$X7.Subject_id)
action$A1=Mus.PS$Gene[probe]
probe=match(action$protein2,Mus.PS$X7.Subject_id)
action$A2=Mus.PS$Gene[probe]
Mus_action=action[,c(11,12,3:10)]
##################################################################################
#Load list of candidate genes
BF=read.csv('../SE_2v2/DEG/BF_ensemble.GeneID',sep='\t',header=F)
BM=read.csv('../SE_2v2/DEG/BM_ensemble.GeneID',sep='\t',header=F)
BF=BM
BF$V1=toupper(BF[BF$V1!='-',])
Mus_selc=Mus_action[(Mus_action$A1 %in% BF$V1) & (Mus_action$A2 %in% BF$V1),]
Mus_selc=Mus_selc[Mus_selc$A1!=Mus_selc$A2 & Mus_selc$combined_score>=400,]
colnames(Mus_selc)[1:5]=colnames(Mus_selc2[,c(1:4,6)])
Mus_selc2=Mus_selc2[rownames(Mus_selc2) %in% rownames(unique(Mus_selc2[,c(1,2)])),]
Mus_selc=Mus_selc[rownames(Mus_selc) %in% rownames(unique(Mus_selc[,c(1,2)])),]
useM=!( (paste(Mus_selc$A1,Mus_selc$A2) %in% paste(Mus_selc2$A1, Mus_selc2$A2)) | 
          (paste(Mus_selc$A2,Mus_selc$A1) %in% paste(Mus_selc2$A1, Mus_selc2$A2)))
############## Combine#########
Mus_selc3=rbind(Mus_selc[useM,c(1:5)],Mus_selc2[,c(1:4,6)])
#g3_M=graph.data.frame(Mus_selc3,directed=T)
g=graph.data.frame(Mus_selc3,directed=T)
V(g)$degree=degree(g)
clusters(g)
V(g)[order(V(g)$degree)]$degree
V(g)[order(V(g)$degree)]
length(E(g))
E(g)$color <- "grey"
E(g)[ mode %in% names(table(Mus_selc2$mode))]$color <- "purple"
E(g)[ action=="activation" ]$color <- "red"
E(g)[ action=="inhibition" ]$color <- "blue"
E(g)[ mode=="ptmod" ]$color <- "purple"
E(g)[ mode=="binding" ]$color <- "black"
BF_info <- read.delim(header=F,"~/Documents/MyLocal/NakedMole/FManalysi/SE_2v2/BF_ensemble")[,c(11,16)]
BM_info <- read.delim(header=F,"~/Documents/MyLocal/NakedMole/FManalysi/SE_2v2/BM_ensemble")[,c(11,16)]
Shared=BM_info[paste(BM_info$V11,BM_info$V16) %in% paste(BF_info$V11,BF_info$V16)  & BM_info$V16!='-'
        ,2]
#############################################
V(g)$Color=BF_info$V11[match(V(g)$name,toupper(BF_info$V16))]
V(g)[V(g)$name %in% Shared]$Color=3
V(g)[is.na(V(g)$Color)]$Color=3
write.graph(g,'PlotingNet/Net_female.gml',format='gml')
V(g)$Color=BM_info$V11[match(V(g)$name,toupper(BM_info$V16))]
write.graph(g,'PlotingNet/Net_male.gml',format='gml')
