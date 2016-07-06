import sys
EdgR=open(sys.argv[1]).readlines()
Cufflink=open(sys.argv[2]).readlines()
DESeq=open(sys.argv[3]).readlines()
DESeq2=open(sys.argv[4]).readlines()
def Store_info(File,Name,Pval,Database,Method):
	for i in File[1:]:
		line=i.split('\t')
		GeneID=line[Name].replace('"','')
		Pvalue=line[Pval].replace('"','')
		if GeneID not in Database:
			Database[GeneID]={}
		Database[GeneID][Method]=Pvalue
	return Database
DEG_data={}
File_list=(EdgR,Cufflink,DESeq,DESeq2)
Name_P_list=[[0,3],[0,11],[1,7],[0,5]]
Method_list=('EdgR','Cufflink','DESeq','DESeq2')


for i in range(len(File_list)):
	DEG_data=Store_info(File_list[i],Name_P_list[i][0],Name_P_list[i][1],DEG_data,Method_list[i])
print("ID",'EdgR','Cufflink','DESeq','DESeq2',sep='\t')
for i in DEG_data:
	print(i,end='')
	for method in Method_list:
		print('\t',end='')
		if method not in DEG_data[i]:
			print('NA',end='')
		else:
			print(DEG_data[i][method],end='')
	print('')

	
