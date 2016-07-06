#Bootstrap to simulation the degree distribution of FI network
#and candidate genes (CRH and NPY)
library(igraph)
BF <- read.table("~/Documents/MyLocal/NakedMole/FManalysi/SE_2v2/DEG/BM_ensemble.GeneID", quote="\"")
B_back <- read.table("~/Documents/MyLocal/NakedMole/FManalysi/SE_2v2/Background/BM_brain.back.gene", quote="\"")
Trans=c()
NPY_list2=c()
Edge_list2=c()
for(i in seq(1,10,1)){Sample_Set=toupper(c(as.character(B_back[sample(row.names(B_back),size=nrow(BF)),1])))
                     Sample_Net=Mus_action[(Mus_action$A1 %in% Sample_Set) & (Mus_action$A2 %in% Sample_Set),]
                     Sample_Net=Sample_Net[Sample_Net$A1!=Sample_Net$A2 & Sample_Net$combined_score>=400,]
                     Sample_g=graph.data.frame(Sample_Net,directed=F)
                     Sample_g=simplify(Sample_g)
                     V(Sample_g)$degree=degree(Sample_g)
                     Trans=c(Trans,transitivity(Sample_g))}
                    
                     #   if(!'NPY' %in% Sample_Net$A1) X1=0 else X1=V(Sample_g)['NPY']$degree;
                     #X2=length(E(Sample_g))
                     #  NPY_list2=c(NPY_list2,X1)}
                     #Edge_list2=c(Edge_list2,X2)}
CRH_list2=c()
for(i in seq(1,1000,1)){Sample_Set=toupper(c(as.character(B_back[sample(row.names(B_back),size=nrow(BF)-1),1]),"CRH"))
                        Sample_Net=Mus_action[(Mus_action$A1 %in% Sample_Set) & (Mus_action$A2 %in% Sample_Set),]
                        Sample_Net=Sample_Net[Sample_Net$A1!=Sample_Net$A2 & Sample_Net$combined_score>=400,]
                        Sample_g=graph.data.frame(Sample_Net,directed=F)
                        Sample_g=simplify(Sample_g)
                        V(Sample_g)$degree=degree(Sample_g)
                        if(!'CRH' %in% Sample_Net$A1) X1=0 else X1=V(Sample_g)['CRH']$degree;
                        X2=length(E(Sample_g))
                        CRH_list2=c(CRH_list2,X1)
                        Edge_list2=c(Edge_list2,X2)}
Edge_list=c()
Node_list=c()
anno=read.delim("~/Documents/MyLocal/NakedMole/Common/anno", 
                header=F)
BF=read.delim("../BF.D",header=T)
BF_back=toupper(anno$V2[match(B_back$V1,anno$V1)])
sample_size=nrow(BF)-2
for(i in seq(1,10000,1)){Sample_Set=toupper(c(sample(BF_back,size=sample_size)));
                      Sample_Net=Mus_action[(Mus_action$A1 %in% Sample_Set) & (Mus_action$A2 %in% Sample_Set),];
                      Sample_Net=Sample_Net[Sample_Net$A1!=Sample_Net$A2 & Sample_Net$combined_score>=400,];
                      Sample_g=graph.data.frame(Sample_Net,directed=F);
                      Sample_g=simplify(Sample_g);
                      V(Sample_g)$degree=degree(Sample_g);
                      X1=length(E(Sample_g));
                      X2=length(V(Sample_g));
                      Edge_list=c(Edge_list,X1);
                      Node_list=c(Node_list,X2)}
Edge2_list=c()
Node2_list=c()
B_back <- read.table("~/Documents/MyLocal/NakedMole/FManalysi/BM_brain-vs-WM_brain.compare4.refine.finale.BackGround_list", quote="\"")
BM=read.delim("../BM.txt",header=T)
BM_back=toupper(anno$V2[match(B_back$V1,anno$V1)])
sample_size=nrow(BM)-2
for(i in seq(1,10000,1)){Sample_Set=toupper(c(sample(BM_back,size=sample_size)));
                     Sample_Net=Mus_action[(Mus_action$A1 %in% Sample_Set) & (Mus_action$A2 %in% Sample_Set),];
                     Sample_Net=Sample_Net[Sample_Net$A1!=Sample_Net$A2 & Sample_Net$combined_score>=400,];
                     Sample_g=graph.data.frame(Sample_Net,directed=F);
                     Sample_g=simplify(Sample_g);
                     V(Sample_g)$degree=degree(Sample_g);
                     X1=length(E(Sample_g));
                     X2=length(V(Sample_g));
                     Edge2_list=c(Edge2_list,X1);
                     Node2_list=c(Node2_list,X2)}                  
save(Edge2_list,Edge_list,Node2_list,Node_list,file='Simulation.Rdata')
