import sys
matrix = sys.argv[1]
AL = 119
HT = 201
print 'ID\tAL-HT_af\tAL_af\tHT_af\tEAS_af'
f=open(matrix,"r").readlines()
for line in f:
	line = line.strip().split()
	if line[0] == "ID":
		continue
	AL_ID = sum([int(line[i]) for i in range(1,AL+1)])/float(AL)
	HT_ID = sum([int(line[i]) for i in range(AL+1,AL+HT+1)]) / float(HT)
	AT_ID = sum([int(line[i]) for i in range(1,AL+HT+1)])/float(AL+HT)
	EAS_ID = sum([int(line[i]) for i in range(AL+HT+1,len(line))]) / float(len(line)-AL-HT-1)
	print line[0]+'\t'+str(AT_ID)+'\t'+str(AL_ID)+'\t'+str(HT_ID)+'\t'+str(EAS_ID)
