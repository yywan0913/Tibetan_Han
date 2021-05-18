import sys

def getsamplevcfinfo(File):
	d ={}
	d.setdefault('samples',{})
	n = 0
	f=open(File,"r")
	for line in f.readlines():
		with open(line.strip()) as sf:
			n +=1
			for sfline in sf:
				if '##' in sfline:
					if n==1:d.setdefault('header',[]).append(sfline.strip())
					continue
				if '#' in sfline:
					sample = sfline.strip().split()[9]
					samplegroup = sfline.strip().split()[9:]
					d['samples'].setdefault(sample,{})
					d['samples_'+sample]=samplegroup
					d.setdefault('samplegroupone',[]).append(sample)
					for k in samplegroup : d.setdefault('samplegroup',[]).append(k)
					continue
				line0=sfline.strip().split('\t')
				ID = line0[2]
				d['samples'][sample][ID] = [k.split(':')[0] for k in line0[9:]]  #line0[-1].split(':')[0]
		sf.close()
	f.close()
	return d


def readmerge(File):
	dmerge = {}
	with open(File) as f:
		for line in f:
			line0=line.strip().split('\t')  ## svID
			dmerge.setdefault(line0[0],{})#.append(line0[1:9])
			sample = line0[-2]
			dmerge[line0[0]].setdefault(sample,[]).append(line0[-8:])
	f.close()
	return dmerge

def get_median(data):
	data = [int(i) for i in data]
	data.sort()
	half = len(data) // 2
	return (data[half] + data[~half]) / 2

def genois(g1,g2):
	if g1=="0|0" or g1 == "0/0":
		return g2
	elif g1=="1|0" or g1 == "0|1" or g1 =="0/1":
		if g2 == "0|0" or g2 == "0/0":
			return g1
		else:
			return g2
	else:
		return g1

def cbindlist(List):
	out = List[0]
	for i in range(1,len(List)):
		for k in range(0,len(List[0])):
			out[k] = genois(out[k],List[i][k])
	return out
	

def getvcf(dmerge,dsamplevcfinfo):
	samples = dsamplevcfinfo['samplegroup']
	groupsampleone = dsamplevcfinfo['samplegroupone']
	if 'header' in dsamplevcfinfo:
		for k in dsamplevcfinfo['header']:
			print k
	print '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(samples)
	for i in sorted(dmerge): ## mergeID
		start = []
		end = []
		gt = []
		svlength = []
		SUPP_VEC=''
		svtype = ''
		N = 0
		for j in groupsampleone:
			if j in dmerge[i]:
				N += 1
				gtindex = []
				for index in dmerge[i][j]:
					chr1 = index[0]
					chr2 = index[2]
					svtype = index[4]
					start.append(index[1])
					end.append(index[3])
					svlength.append(index[5])
					gtindex.append(dsamplevcfinfo['samples'][j][index[-1]])
				for gtindexout in cbindlist(gtindex): gt.append(gtindexout) #
				SUPP_VEC+='1'
			else:
				for jsamples in dsamplevcfinfo['samples_'+j]: gt.append('0/0')
				SUPP_VEC+='0'
		if N != len(groupsampleone):
			continue
		if svtype == "TRA":svlength='0'
		if svtype == "INS":END = str(get_median(start)+1)
		AVGLEN = str(round(sum([int(k) for k in svlength])/len(svlength),2))
		INFO = "SUPP=%s;SUPP_VEC=%s;AVGLEN=%s;SVTYPE=%s;CHR2=%s;END=%s;CIPOS=0,0;CIEND=0,0;STRANDS=++" %(str(len(start)),SUPP_VEC,AVGLEN,svtype,chr2,str(get_median(end)))
		print '\t'.join([chr1,str(get_median(start)),i,'N',svtype,'.','pass',INFO,'GT'])+'\t'+'\t'.join(gt)

def main():
	mergevcf = sys.argv[1]
	vcflist = sys.argv[2]	
	dmerge = readmerge(mergevcf)
	dsamplevcfinfo=getsamplevcfinfo(vcflist)
	getvcf(dmerge,dsamplevcfinfo)

if __name__ == "__main__":
	main()
			
		
