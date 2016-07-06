#Change Y and N1+N2 to modify the criteria.
a=open(r'../Heter_glaber.v1.7.correct.fa','r')
Out=open(r'Ref.Pos.V2','w')
N1=0 #Relative position on pseudo chromosome.
N2=500 # Length of scaffold (N1+N2: End of scaffold on pseudo chromosome.) Minimal Length: N*500.
ID='ID'
X=0 #Number of tested scaffold
Chromo='PseudoC_'
Y=0 #Chromosome Name
while True:
    line=a.readline()
    Chromo_ID=Chromo+str(Y)
    if line=='' :#or Y>1: #Stop when over 24 pseudochromosomes are generated.
        print(ID,Chromo_ID,str(N1),str(N1+N2),str(N2),sep='\t',file=Out)
        break
    if '>' in line:
        if X!=0:
            print(ID,Chromo_ID,str(N1),str(N1+N2),str(N2),sep='\t',file=Out)
            N1=N1+N2
        ID=line.strip('>').rstrip('\n')
        N2=500
        if N1+N2>= 1e7: #2664766285/20: #split into 20 chromosome (New Pseudochrom, when the start over is too long)
            Y+=1
            N1=0
        X+=1
    else:
        N2+=len(line.rstrip('\n'))
Out.close()
