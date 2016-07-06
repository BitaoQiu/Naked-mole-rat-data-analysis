library(DESeq2)
FF=c(3:4,9:10)
MM=c(12,13,15,16)
FF2=c(3,9,4,10)
MM2=c(12,15,13,16)
FF3=c(3,10,4,9)
MM3=c(12,16,15,13)
data_name=paste(rep(c('female','male'),3),rep(c(1,2,3),each = 2),sep ='_')
lib=read.table('/Volumes/BACK1/NMR/Volumes/BACK1/NakedMole/Common/lib.toNorm',header=T)
Data_list=list(FF,MM,FF2,MM2,FF3,MM3)
Female=lib[,FF]
plot(hclust(dist(t(log(Female+1)))))
head(Female)
for (i in seq(1,6,1)){
  Female=lib[,Data_list[[i]]];
  rownames(Female)=lib$GeneID;
  colData=data.frame(condition=c(rep('breeder',2),rep('nbreeder',2)),Type=rep('Female',4));
  rownames(colData)=names(Female);
  dds <- DESeqDataSetFromMatrix(countData = Female,
                                colData = colData,
                                design = ~ condition);
  colData(dds)$condition <- factor(colData(dds)$condition,
                                   levels=c("breeder","nbreeder"));
  dds <- DESeq(dds);
  res <- results(dds);
  res <- res[order(res$padj),];
  Out_name=data_name[i];
  write.table(as.data.frame(res),file=paste(Out_name,'D3',sep='.'),quote=F,sep='\t')}                      #head(res)
