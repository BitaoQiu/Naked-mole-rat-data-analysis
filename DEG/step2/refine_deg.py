import sys
File=open(sys.argv[1],'r').readlines()
Header=File[0]
PA1=1
for pos in range(1,len(Header.split('\t'))):
	if 'Direction' in Header.split('\t')[pos]:
		Dp=pos
	if Header.split('\t')[pos].split('.')[1][0:-2]==Header.split('\t')[1].split('.')[1][0:-2]:
		PA2=pos
		PB1=pos+1
	elif Header.split('\t')[pos].split('.')[1][0:-2]==Header.split('\t')[PB1].split('.')[1][0:-2]:
		PB2=pos
Out_S=open(sys.argv[2],'w')
Out_filter=open(sys.argv[3],'w')
def mean(List):
	Sum=0
	for i in List:
		Sum+=float(i)
	return Sum/len(List)
def Compare(A1,A2): #Compare the median#
	nA1=round(len(A1)/2)-1
	if len(A1)%2==0:
		nA2=nA1+1
	else:
		nA2=nA1
	nB1=round(len(A2)/2)-1
	if len(A2)%2==0:
		nB2=nB1+1
	else:
		nB2=nB1
	
	if (sorted(A1)[nA1]+sorted(A1)[nA2])/2 < (sorted(A2)[nB1]+sorted(A2)[nB2])/2 and min(A1)<min(A2):
		return True
	else:
		return False
print(Header,end='',file=Out_filter)
print(Header,end='',file=Out_S)
for line in File[1:]:
	i=line.split('\t')
	Dir=i[Dp]
	A1=[]
	for j in range(PA2-1,(PA2+1)):# Change to 2v2
		A1.append(float(i[j]))
	A2=[]
	for j in range(PB2-1,(PB2+1)): #Change to 2v2
		A2.append(float(i[j]))
	if mean(A1)<mean(A2):#Dir=='Up':
		if Compare(A1,A2):
			print(line,end='',file=Out_S)
		else:
			print(line,end='',file=Out_filter)
	if mean(A1)>mean(A2):#Dir=='Down':
		if Compare(A2,A1):
			print(line,end='',file=Out_S)
		else:
			print(line,end='',file=Out_filter)
Out_filter.close()
Out_S.close()

