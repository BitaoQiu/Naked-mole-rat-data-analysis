for i in */;
do
cd $i;
LenStr=${#i}-1
/share/raid12/qiubitao/bin/bin/python3.3 ../Get_Compare.py edgeR/*.OtherNorm.DESeq.txt Cuffdiff/gene_exp.diff DEseq/*.SelfNorm.DESeq.txt *_*.D2 >${i:0:$LenStr}.compare4_3
cd ..
done
