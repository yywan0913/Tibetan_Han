import sys
gff = sys.argv[1]
bed = sys.argv[2]
windows = 100000 if sys.argv[3]==None else int(sys.argv[3])
d = {}
with open(gff,"r") as f:
	for line in f:
		line = line.strip().split('\t')
		if '#' in line[0] or 'Contig' in line[0] or 'Un' in line[0] or 'un' in line[0]:
			continue
		if 'gene' != line[2]:
			continue
		chrom = line[1]
		start = int(line[2])
		mode = start/windows
                sv_type = line[4]
                #print(sv_type)
                #print(str(mode)+"="+str(line[1])+"/"+str(windows))
		if chrom not in d:
                        d[chrom]={}
                else:
                        if mode not in d[chrom]:
                                d[chrom][mode]=1
                        else:
                                d[chrom][mode]+=1
		#d[chrom].setdefault(mode,[]).append(length)
                #print(line)



f= open(bed,"r").readlines()
for line in f:
	line = line.strip().split('\t')
	chrom = line[0]
	mode = int(line[1])/windows
        #print(str(mode)+"="+str(line[1])+"/"+str(windows))
	if mode in d[chrom]:
                count=d[chrom][mode]
	else:
		count = 0
	print line[0]+'\t'+line[1]+'\t'+str(int(line[2])-1)+'\t'+str(count)
		
