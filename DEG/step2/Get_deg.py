#Get DEG
import sys
File=open(sys.argv[1],'r').readlines()
m=0
for i in File[0].split('\t'):
    if '.Sig' in i:
        Pos=m
        break
    else:
        m+=1
print(File[0],end='')
for i in File[1:]:
	if int(i.split('\t')[Pos])==1:
		print(i,end='')

    
        
