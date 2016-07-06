import sys
#RPKM=open(sys.argv[1],'r').readlines()
RPKM=open(r'NMR.eusociality.uniqReads.rpkm.ann','r').readlines()
File=open(sys.argv[1],'r')
CompareA=sys.argv[1].split('-vs-')[0]
CompareB=sys.argv[1].split('-vs-')[1].split('.')[0]
for i in RPKM[0:30]:
    if i.startswith('#') and 'col' in i:
        Position,Name=i.split(': ')[0],i.split(':')[1].strip('\n').strip(' ')
       #print(Name)
        if CompareA==Name:
            A1=int(Position.split('-')[0].strip('# col'))
            A2=int(Position.split('-')[1].strip('col'))
        if CompareB==Name:
            B1=int(Position.split('-')[0].strip('# col'))
            B2=int(Position.split('-')[1].strip('col'))
Gene_Dict={}
for i in File:
    line=i.split('\t')
    if line[0] not in Gene_Dict:
        Gene_Dict[line[0]]={}
    Gene_Dict[line[0]]['P_val']=line[1:]
    M=0
    for i in line[1:]:
        if i=='NA' or i=='NA\n':
            pass
        elif float(i)<0.05:
            M+=1
    if M>=2:
        Gene_Dict[line[0]]['Sig']='1'
    else:
        Gene_Dict[line[0]]['Sig']='0'


for i in RPKM:
    if not i.startswith('#'):
        line=i.split('\t')
        if line[0] in Gene_Dict:
            C1=line[A1-1:A2]
            C2=line[B1-1:B2]
            for i in range(len(C1)):
                C1[i]=float(C1[i])
            for i in range(len(C2)):
                C2[i]=float(C2[i])
            C1A=sum(C1)/len(C1)
            C2A=sum(C2)/len(C2)
            if min(C1A,C2A)<1:
                Fold=max(C1A,C2A)
            else:
                Fold=max(C1A,C2A)/min(C1A,C2A)
            Gene_Dict[line[0]]['C1']=C1
            Gene_Dict[line[0]]['C2']=C2
            Gene_Dict[line[0]]['C1A']=C1A
            Gene_Dict[line[0]]['C2A']=C2A
            Gene_Dict[line[0]]['Fold']=Fold
            Gene_Dict[line[0]]['Symbol']=line[30].strip('\n')
            Gene_Dict[line[0]]['Descript']=line[31].strip('\n')
            if C2A>C1A:
                Gene_Dict[line[0]]['Direct']='Up'
            else:
                Gene_Dict[line[0]]['Direct']='Down'
Col=0
Col=Col+1
print(str(Col)+'.'+'GeneID',end='\t')
for i in range(len(C1)):
    Col=Col+1
    print(str(Col)+'.'+CompareA+'_'+str(i+1),end='\t')
for i in range(len(C2)):
    Col=Col+1
    print(str(Col)+'.'+CompareB+'_'+str(i+1),end='\t')
print(str(Col+1)+'.'+CompareA+'_avg',str(Col+2)+'.'+CompareB+'_avg',
      str(Col+3)+'.'+'Fold_C',str(Col+4)+'.'+'Direction',
      str(Col+5)+'.''edgR',str(Col+6)+'.'+'Cufflink',
      str(Col+7)+'.'+'DESeq2',str(Col+8)+'.'+'Sig',
      str(Col+9)+'.'+'Symbol',str(Col+10)+'.'+'Description',sep='\t')

for i in Gene_Dict:
    if 'C1' in Gene_Dict[i]:
        print(i,end='\t')
        for j in Gene_Dict[i]['C1']:
            print(round(j,2),end='\t')
        for j in Gene_Dict[i]['C2']:
            print(round(j,2),end='\t')
        print(round(Gene_Dict[i]['C1A'],2),round(Gene_Dict[i]['C2A'],2),
              round(Gene_Dict[i]['Fold'],4),Gene_Dict[i]['Direct'],sep='\t',end='\t')
        for j in Gene_Dict[i]['P_val']:
            if j.strip('\n')=='NA':
                print(j.strip('\n'),end='\t')
            else:
                print(round(float(j.strip('\n')),4),end='\t')
        print(Gene_Dict[i]['Sig'],Gene_Dict[i]['Symbol'],
              Gene_Dict[i]['Descript'],sep='\t')


                

                
    
 

