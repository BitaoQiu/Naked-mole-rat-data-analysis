#Bootstrap simulation to test the egde of FI network.
library(igraph)
anno=read.delim("/home/enjoycloud/Documents/MyLocal/NakedMole/Gonad//SE2V2/BACK/Testes_Back.gene", header=F)
anno=read.delim("/home/enjoycloud/Documents/MyLocal/NakedMole/Common/anno", header=F)
GeneSet1='/home/enjoycloud/Documents/MyLocal/NakedMole/FManalysi/SE_2v2/DEG/BM_ensemble.GeneID'
GeneSet2='/home/enjoycloud/Documents/MyLocal/NakedMole/Meth2/SAMPLE/MEDIPS/TEST/DMG/GO_geneID/MB_DMG_all_ID.txt'
BM=read.csv(GeneSet1,sep='\t',header=F)
BF=read.csv(GeneSet2,sep='\t',header=F)
BF=toupper(BF[BF!='-'])
BM=toupper(BM[BM!='-'])
Round=10
Edge_list=rep(0,Round)
Tran_list=rep(0,Round)
Node_list=rep(0,Round)
sample_size=length(BM)
Background=anno$V1
for(i in seq(1,Round,1)){Sample_Set=toupper(sample(Background,size=sample_size));
                     Sample_Net=Mus_action[(Mus_action$A1 %in% c(BF,Sample_Set)) &
                                             (Mus_action$A2 %in% c(BF,Sample_Set)),];
                     Sample_Net=Sample_Net[Sample_Net$A1!=Sample_Net$A2 & 
                                             Sample_Net$combined_score>=400,];
                     Sample_Net=Sample_Net[order(Sample_Net[,c('A1','A2')]),]
                     Sample_Net=(unique(Sample_Net[,c('A1','A2')]))
                     Sample_g=graph.data.frame(Sample_Net,directed=F);
                     Sample_g=simplify(Sample_g);
                     Edge_list[i]=length(E(Sample_g));
                     Node_list[i]=length(V(Sample_g));
                     Tran_list[i]=transitivity(Sample_g)
}

#save(Edge2_list,Edge_list,Node2_list,Node_list,file='Simulation.Rdata')

