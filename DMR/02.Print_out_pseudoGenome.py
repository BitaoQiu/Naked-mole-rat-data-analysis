Ref=open(r'Ref.Pos.V2','r').readlines()
Genome=open(r'../Heter_glaber.v1.7.correct.fa','r')
def break_line(seq,n,file):
    Nline=len(seq)//n
    for i in range(Nline):
        print(seq[i*n:(i+1)*n],file=file)
    if len(seq)%n!=0:
        print(seq[Nline*n:],file=file)
    return
def Print_seq(Chrome):
    Out_file=open('BS_genome/'+Chrome,'w')
    print('>'+Chrome,file=Out_file)
    break_line(Chrome_seq[Chrome],80,Out_file)
    Out_file.close()
    return 
Scaf_chrome={}
Chrome_seq={}
X=0#Tested Scaffolds
#Build scaffold name reference
for i in Ref:
    Scaf_N=i.split('\t')[0]
    Scaf_chrome[Scaf_N]=i.split('\t')[1]
while True:
    line=Genome.readline()
    if line=='':
        Print_seq(Chromo)
        break
    if '>' in line:
        if X!=0:
            Chrome_seq[Chromo]+=500*'N'
        Name=line.lstrip('>').rstrip('\n')
        if Name not in Scaf_chrome.keys():
            Print_seq(Chromo)
            break
        else:
            if Scaf_chrome[Name] not in Chrome_seq.keys():
                if X!=0:
                    Print_seq(Chromo)
                Chrome_seq={}
                Chromo=Scaf_chrome[Name]
                Chrome_seq[Chromo]=''
            Chromo=Scaf_chrome[Name]
        X=1
    else:
        Chrome_seq[Chromo]+=line.rstrip('\n')

