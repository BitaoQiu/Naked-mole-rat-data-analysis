mouse1=read.table('19-tissues-expr/testes-1.expr',header=T)
mouse2=read.table('19-tissues-expr/testes-2.expr',header=T)
mouse1$FPKM2=mouse2$FPKM[match(mouse1$gene_id,mouse2$gene_id)]
mouse=mouse1[,c('gene_id','FPKM','FPKM2')]
mouse[is.na(mouse)]=0
anno=read.delim('mart_export.txt',header = T)
names(anno)[2]='Ref'
anno=anno[!anno$Ref=='',]
table(mouse$gene_id %in% anno$Ref)
naked_mole=RPKM[,grep('BM_testis',names(RPKM))]
naked_mole$name=RPKM$GeneID
naked_mole$ID=anno_2$V2[match(naked_mole$name,anno_2$V1)]
mouse$ID=anno$Associated.Gene.Name[match(mouse$gene_id,anno$Ref)]
mouse$ID = toupper(mouse$ID)
naked_mole$ID = toupper(naked_mole$ID)
table(mouse$ID %in% naked_mole$ID)
table(apply(naked_mole[,c(2,3)],MARGIN = 1,FUN = min)>1)
naked_mole_F=naked_mole[apply(naked_mole[,c(2,3)],MARGIN = 1,FUN = min)>1
                        ,c(5,2,3)]
mouse_F=mouse[,c(4,2,3)]
mouse_F$mean=rowMeans(mouse_F[,c(2,3)])
naked_mole_F$mean=rowMeans(naked_mole_F[,c(2,3)])
mouse_F2=mouse_F[mouse_F$ID %in% naked_mole_F$ID,]
naked_mole_F2=naked_mole_F[naked_mole_F$ID %in%
                             mouse_F$ID,]
naked_mole_F2$ID=factor(naked_mole_F2$ID)
mouse_F2$ID=factor(mouse_F2$ID)
naked_mole_F2_sum=aggregate(mean~ID,data = naked_mole_F2,FUN =
                            function(x) c(Testis=sum(x)))
mouse_F2_sum=aggregate(mean~ID,data = mouse_F2,FUN =
                              function(x) c(Testis=sum(x)))

mouse_F2_sum$mean_n=naked_mole_F2_sum$mean[match(mouse_F2_sum$ID,naked_mole_F2_sum$ID)]

boxplot(mouse_F2_sum[,c(2,3)],outline = F)

human_altas$hgnc_symbol=factor(human_altas$hgnc_symbol)
human_altas_sum=aggregate(testes~hgnc_symbol,data = human_altas,FUN =
            function(x) c(Testis=sum(x)))
human_altas_sum=human_altas_sum[human_altas_sum$testes>0.5,]
mouse_F2_sum=mouse_F2_sum[mouse_F2_sum$ID %in% 
               human_altas_sum$hgnc_symbol,]
human_altas_sum=human_altas_sum[human_altas_sum$hgnc_symbol
                                %in% mouse_F2_sum$ID,]
mouse_F2_sum$mean_h=human_altas_sum$testes[match(
  mouse_F2_sum$ID,human_altas_sum$hgnc_symbol)]

Ratbodymap_Gene_FPKM_v2=Ratbodymap_Gene_FPKM_v2[,c(1,grep('Tst_M_021',names(Ratbodymap_Gene_FPKM_v2)))]
Ratbodymap_Gene_FPKM_v2$median=apply(Ratbodymap_Gene_FPKM_v2[,-1],
                                     1,FUN=function(x) 
                                       return(median(x)))
Ratbodymap_sum=aggregate(median~GeneSymbol,data = Ratbodymap_Gene_FPKM_v2,FUN =
                            function(x) c(Testis=sum(x)))
Ratbodymap_sum$GeneSymbol=toupper(Ratbodymap_sum$GeneSymbol)
Ratbodymap_sum=Ratbodymap_sum[Ratbodymap_sum$GeneSymbol
                                %in% mouse_F2_sum$ID,]
mouse_F2_sum=mouse_F2_sum[mouse_F2_sum$ID %in% 
                            Ratbodymap_sum$GeneSymbol,]
mouse_F2_sum$mean_R=Ratbodymap_sum$median[match(
  mouse_F2_sum$ID,Ratbodymap_sum$GeneSymbol)]

Chimpanzee_Testis=Chimpanzee_Ensembl57_TopHat_UniqueReads[,UseID]
Chimpanzee_Testis$ID=Chimpanzee_anno$Associated.Gene.Name[match(Chimpanzee_Testis$GeneID,Chimpanzee_anno$Ensembl.Gene.ID)]
Chimpanzee_sum=aggregate(Chimpanzee_Testis_Male~ ID,data = Chimpanzee_Testis,FUN =
                           function(x) c(Testis=sum(x)))
Chimpanzee_sum=Chimpanzee_sum[Chimpanzee_sum$ID
                              %in% mouse_F2_sum$ID,]
mouse_F2_sum=mouse_F2_sum[mouse_F2_sum$ID %in% 
                            Chimpanzee_sum$ID,]
mouse_F2_sum$mean_C=Chimpanzee_sum$Chimpanzee_Testis_Male[match(
  mouse_F2_sum$ID,Chimpanzee_sum$ID)]

library(preprocessCore)
mouse_F_normalize=normalize.quantiles(as.matrix(
  mouse_F2_sum[,-1]))
rownames(mouse_F_normalize)=mouse_F2_sum$ID
boxplot(mouse_F_normalize,outline = F)
mouse_sperm=mouse_F_normalize[rownames(mouse_F_normalize) %in% sperm2$V1,]

boxplot(mouse_sperm,outline = F,
        names=c('mouse','NMR','Human','Rat','Chimpan'))
