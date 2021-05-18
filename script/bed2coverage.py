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
		chrom = line[0]
		start = int(line[1])
		end = int(line[2])
		mode = start/windows
                #print(str(mode)+"="+str(line[1])+"/"+str(windows))
		length = end-start+1
		d.setdefault(chrom,{})
		d[chrom].setdefault(mode,[]).append(length)



f= open(bed,"r").readlines()
for line in f:
	line = line.strip().split('\t')
	chrom = line[0]
	mode = int(line[1])/windows
        #print(str(mode)+"="+str(line[1])+"/"+str(windows))
	if mode in d[chrom]:
		coverage = float(sum(d[chrom][mode]))/(int(line[2])-int(line[1]))*100
		coverage = 100 if coverage>100 else coverage
	else:
		coverage = 0
	print line[0]+'\t'+line[1]+'\t'+str(int(line[2])-1)+'\t'+str(coverage)
		
