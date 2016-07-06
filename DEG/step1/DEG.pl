#!/usr/bin/perl -w
#Zhang Pei: zhangpei@genomics.org.cn
use strict;
use Getopt::Long;
GetOptions(
		'ReadNum=s'=>\our$ReadNum,
		'BamList=s'=>\our$BamList,
		'Lib=s'=>\our$Lib,
		'Genome=s'=>\our$Genome,
		'Gtf=s'=>\our$Gtf,
		'LibraryType=s'=>\our$LibraryType,
#		'Analysis=s'=>\our$Analysis,
#		'Method=s'=>\our$Method,
#		'DESeq=s'=>\our$DESeq,
		'LibrarySize=s'=>\our$LibrarySize
#		'Dispersion=s'=>\our$Dispersion
		);
my$help=<<"HELP";
	--ReadNum: <str> Read number count file
	--BamList: <str> File that contains bam file path
	--Lib: <str> file control DEG calculation
	--Genome: <str> Genome file
	--Gtf: <str> Gtf file
	--LibraryType: <str> library type, only for Cuffdiff, default fr-unstranded, others are fr-firststrand and fr-secondstrand
#	--Analysis: <str> analysis type: Pairwise analysis: Pairwise; ANOVA-like analysis: ANOVA; Time-series: TimeSeries; Count data transformations: VSD (Just stabilizing variance)
#	--Method: <str> Methods used for DEG detection. Pairwise analysis: DESeq, Cuffdiff, edgeR; ANOVA-like analysis: DESeq, edgeR; 
#	--DESeq: <str> Wald or LRT; default: Wald for pairwise or one-way anova like analysis, LRT for multi-way anova like analysis
#	--LibrarySize: <str> library size factors
HELP
if(!$ReadNum || !$BamList || !$Lib || !$Genome || !$Gtf){die "$help\n"}
if(!$LibraryType){$LibraryType='fr-unstranded'}
my$Dispersion;
our%LSize;
our%LFactor;
our%LFactor2;
if($LibrarySize){
	my$TTSize;
	my$TTnum;
	open IN,"$LibrarySize";
	<IN>;
	while(<IN>){
		chomp;
		my@A=split(/\t/);
		$LSize{$A[0]}=$A[3];
		$LFactor{$A[0]}=$A[2];
		$TTSize+=$A[1];
		$TTnum++;
	}
	close IN;
	my$avg=$TTSize/$TTnum;
	foreach my$key(keys %LSize){
		$LFactor2{$key}=$LSize{$key}/$avg;
	}
}
our$LSize="";
our$LFactor="";
our$LFactor2="";
my%num;
open IN,"$ReadNum";
my$id=<IN>;
chomp$id;
my@id=split(/\t/,$id);
while(<IN>){
	chomp;
	my@A=split(/\t/);
	for(my$i=1;$i<=$#A;$i++){
		$num{$A[0]}{$id[$i-1]}=$A[$i];
	}
}
close IN;
my%bam;
open IN,"$BamList";
while(<IN>){
	chomp;
	my@A=split(/\s+/);
	$bam{$A[0]}=$A[1];
}
close IN;
open IN,"$Lib";
my%tag;
my%compare;
while(<IN>){
	chomp;
	my@A=split(/\s+/);
	if($#A ==1){
		push @{$tag{$A[0]}},$A[1];
	}else{
		$compare{$A[0]}++;
	}
}
close IN;
foreach my$key(keys %compare){
	mkdir "$key";
	chdir "$key";
	open OUT,">$key.num";
	open OUT1,">$key.design";
	my$tag1;
	my$tag2;
	my$num1=0;
	my$num2=0;
	if($key=~/^(.+?)-vs-(.+?)$/){$tag1=$1;$tag2=$2}
	my$LSize="";
	my$LFactor="";
	my$LFactor2="";
	print OUT "GeneID";
	print OUT1 "Sample\n";
	foreach my$ele(@{$tag{$tag1}}){
		print OUT "\t$ele";
		print OUT1 "$ele\t$tag1\n";
		$num1++;
		if($LibrarySize){
			$LSize=$LSize."$LSize{$ele}".",";
			$LFactor=$LFactor."$ele"."="."$LFactor{$ele}".",";
			$LFactor2=$LFactor2."$ele"."="."$LFactor2{$ele}".",";
		}
	}
	foreach my$ele(@{$tag{$tag2}}){
		print OUT "\t$ele";
		print OUT1 "$ele\t$tag2\n";
		$num2++;
		if($LibrarySize){
			$LSize=$LSize."$LSize{$ele}".",";
			$LFactor=$LFactor."$ele"."="."$LFactor{$ele}".",";
			$LFactor2=$LFactor2."$ele"."="."$LFactor2{$ele}".",";
		}
	}
	if($LibrarySize){
		chop($LSize);
		chop($LFactor);
		chop($LFactor2);
	}
	print OUT "\n";
	foreach my$gene(keys %num){
		print OUT "$gene";
		foreach my$ele(@{$tag{$tag1}}){
			print OUT "\t$num{$gene}{$ele}";
		}
		foreach my$ele(@{$tag{$tag2}}){
			print OUT "\t$num{$gene}{$ele}";
		}
		print OUT "\n";
	}
	close OUT;
	close OUT1;
	if($LibrarySize){
		&edgeR($key,$num1,$num2,$tag1,$tag2,$LFactor);
		&DEseq($key,$tag1,$tag2,$LFactor2);
		&baySeq($key,$num1,$num2,$tag1,$tag2,$LSize);
	}else{
		&edgeR($key,$num1,$num2,$tag1,$tag2);
		&DEseq($key,$tag1,$tag2);
		&baySeq($key,$num1,$num2,$tag1,$tag2);
	}
	mkdir "Cuffdiff";
	chdir "Cuffdiff";
	open OUT,">Cuffdiff.sh";
	print OUT "export LD_LIBRARY_PATH=/ifs5/PC_PA_UN/ANIMAL/USER/GROUP2/zhangpei/bin/lib/boost_1_57_0/lib:\$LD_LIBRARY_PATH\n";
	print OUT "export LD_LIBRARY_PATH=/ifs5/PC_PA_UN/ANIMAL/USER/GROUP2/zhangpei/bin/lib/zlib-1.2.8/lib:\$LD_LIBRARY_PATH\n";
	print OUT "export LD_LIBRARY_PATH=/ifs5/PC_PA_UN/ANIMAL/USER/GROUP2/zhangpei/bin/samtools-1.1/bin/lib:\$LD_LIBRARY_PATH\n";
	print OUT "export PATH=/ifs5/PC_PA_UN/ANIMAL/USER/GROUP2/zhangpei/bin/cufflinks-2.2.1.Linux_x86_64:\$PATH\n";
	my$bam1="";
	my$bam2="";
	foreach my$ele(@{$tag{$tag1}}){
		$bam1.="$bam{$ele},";
	}
	foreach my$ele(@{$tag{$tag2}}){
		$bam2.="$bam{$ele},";
	}
	chop$bam1;
	chop$bam2;
	print OUT "/ifs5/PC_PA_UN/ANIMAL/USER/GROUP2/zhangpei/bin/cufflinks-2.2.1.Linux_x86_64/cuffdiff -b $Genome -u -N -p 8 --library-type $LibraryType -L $tag1,$tag2 $Gtf $bam1 $bam2";
	close OUT;
	chdir "../";
	open OUT, ">run.sh";
	print OUT "cd edgeR; sh edgeR.sh; cd ../\n";
	print OUT "cd DEseq; sh DEseq.sh; cd ../\n";
	print OUT "cd baySeq; sh baySeq.sh; cd ../\n";
	print OUT "cd Cuffdiff; sh Cuffdiff.sh; cd ../\n";
	chdir "../";
}
sub edgeR{
	my$compare=shift;
	my$num1=shift;
	my$num2=shift;
	my$tag1=shift;
	my$tag2=shift;
	mkdir "edgeR";
	chdir "edgeR";
	open OUT, ">edgeR.R";
	my$value;
	for(my$i=1;$i<=$num1;$i++){
		$value.="\"$tag1\",";
	}
	for(my$i=1;$i<=$num2;$i++){
		$value.="\"$tag2\",";
	}
	chop$value;
	my$edger;
	if($LibrarySize){
		 $edger=<<"EDGER";
library( edgeR )
contsTab<-read.delim("../$compare.num",row.names=1)
cound<-c($value)
libsizes<-c($_[-1])
dgl<-DGEList(counts = contsTab, group = cound,norm.factors=libsizes)
cpm.d<-cpm(dgl)
dgl<-dgl[rowSums(cpm.d > 1) >=3,]
dgl <- estimateTagwiseDisp( dgl )
edgerResTBVSW <- exactTest( dgl, pair=c("$tag1","$tag2"))
edgerResTBVSW\$table\$Padj <- p.adjust( edgerResTBVSW\$table\$PValue, method="BH")
write.table(edgerResTBVSW\$table,"$tag1-vs-$tag2.OtherNorm.DESeq.txt",sep="\\t")
EDGER
	}else{
	$edger=<<"EDGER";
library( edgeR )
contsTab<-read.delim("../$compare.num",row.names=1)
cound<-c($value)
dgl<-DGEList(counts = contsTab, group = cound)
cpm.d<-cpm(dgl)
dgl<-dgl[rowSums(cpm.d > 1) >=3,]
dgl<-calcNormFactors(dgl, method="upperquartile")
dgl <- estimateTagwiseDisp( dgl )
edgerResTBVSW <- exactTest( dgl, pair=c("$tag1","$tag2"))
edgerResTBVSW\$table\$Padj <- p.adjust( edgerResTBVSW\$table\$PValue, method="BH")
write.table(edgerResTBVSW\$table,"$tag1-vs-$tag2.OtherNorm.DESeq.txt",sep="\\t")
EDGER
	}
	print OUT "$edger\n";
	close OUT;
	open OUT, ">edgeR.sh";
	print OUT "/ifs5/PC_PA_UN/ANIMAL/USER/GROUP2/zhangpei/bin/R-3.2.3/bin/R -f edgeR.R";
	close OUT;
	chdir "../";
}
sub DEseq{
	my$compare=shift;
	my$tag1=shift;
	my$tag2=shift;
	mkdir "DEseq";
	chdir "DEseq";
	open OUT,">DEseq.R";
	my$deseq;
	$deseq=<<"DESEQ";
library( DESeq2 )
contsTab<-read.delim("../$compare.num",row.names=1)
desi<-read.table("../$compare.design")
dds<-DESeqDataSetFromMatrix(contsTab,desi,design =~ Sample)
dds<-DESeq(dds)
res<-results(dds)
write.table(res,"$tag1-vs-$tag2.SelfNorm.DESeq.txt",sep="\\t")
DESEQ
	print OUT "$deseq\n";
	close OUT;
	open OUT, ">DEseq.sh";
	print OUT "/ifs5/PC_PA_UN/ANIMAL/USER/GROUP2/zhangpei/bin/R-3.2.3/bin/R -f DEseq.R";
	close OUT;
	chdir "../";
}
sub baySeq{
	my$compare=shift;
	my$num1=shift;
	my$num2=shift;
	my$tag1=shift;
	my$tag2=shift;
	mkdir 'baySeq';
	chdir "baySeq";
	open OUT, ">baySeq.R";
	my$value;
	my$NDE;
	my$DE;
	for(my$i=1;$i<=$num1;$i++){
		$value.="\"$tag1\",";
		$NDE.="1,";
		$DE.="1,";
	}
	for(my$i=1;$i<=$num2;$i++){
		$value.="\"$tag2\",";
		$NDE.="1,";
		$DE.="2,";
	}
	chop$value;
	chop$NDE;
	chop$DE;
	my$bayseq;
	if($LibrarySize){
	$bayseq=<<"BAYSEQ";
library(baySeq)
library(snow)
data<-read.delim("../$compare.num",row.names=1)
replicates<-c($value)
lib<-c($_[-1])
groups<-list(NDE=c($NDE),DE=c($DE))
data<-as.matrix(data)
CD<-new("countData", data=data, replicates=replicates, groups=groups,libsizes = as.integer(lib))
cDPair <- getPriors.NB(CD, cl = NULL)
cDPair <- getLikelihoods.NB(cDPair, nullData = TRUE, cl = NULL)
cD.TPs <- getTPs(cDPair, group=2, TPs = 1:100)
num<-nrow(cDPair\@data)
bayseq_de = topCounts(cDPair, group=2, number=num)
write.table(bayseq_de, file="$tag1-vs-$tag2.SelfNorm.baySeq.txt",sep="\\t")
BAYSEQ
	}else{
	$bayseq=<<"BAYSEQ";
library(baySeq)
library(snow)
data<-read.delim("../$compare.num",row.names=1)
replicates<-c($value)
groups<-list(NDE=c($NDE),DE=c($DE))
data<-as.matrix(data)
CD<-new("countData", data=data, replicates=replicates, groups=groups)
libsizes(CD) <- getLibsizes(CD, estimationType = "quantile")
cDPair <- getPriors.NB(CD, cl = NULL)
cDPair <- getLikelihoods.NB(cDPair, nullData = TRUE, cl = NULL)
cD.TPs <- getTPs(cDPair, group=2, TPs = 1:100)
num<-nrow(cDPair\@data)
bayseq_de = topCounts(cDPair, group=2, number=num)
write.table(bayseq_de, file="$tag1-vs-$tag2.SelfNorm.baySeq.txt",sep="\\t")
BAYSEQ
	}
	open OUT, ">baySeq.R";
	print OUT "$bayseq\n";
	close OUT;
	open OUT,">baySeq.sh";
	print OUT "/ifs5/PC_PA_UN/ANIMAL/USER/GROUP2/zhangpei/bin/R-3.2.3/bin/R -f baySeq.R";
	close OUT;
	chdir "../";
}
