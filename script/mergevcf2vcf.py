#!/sfs-grand-med-research/home/xiaoyh/software/miniconda2/bin/python2
import os
import re
import sys
import time
import argparse
from glob import glob

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass

def getsamplevcfinfo(File,prefix):
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
					sample = os.path.basename(sfline.strip().split()[9])
					#if 'ngmlr/' in sample and not prefix:sample=sample.replace('ngmlr/',prefix+'.')
					d['samples'].setdefault(os.path.basename(sample),{})
					continue
				line0=sfline.strip().split('\t')
				ID = line0[2]
				if 'RNAMES=' in line0[7]:
					ReadsName = line0[7].split('RNAMES=')[1].split(';')[0]
				else:
					ReadsName = ''
				#gt = line0[9].split(':')[0]
				gt =  line0[9]
				d['samples'][sample][ID] = [ReadsName,gt]
				#d['samples'][sample][ID] = line0[-1].split(':')[0]
		sf.close()
	f.close()
	return d

def getsamples(File,change):
	d={}
	d.setdefault('samples',{})
	d.setdefault('S',[])
	f=open(File,"r")
	for line in f.readlines():
		line0 = line.strip().split()
		if change:
			sample = line0[0]+'.'+line0[1]+'.ngmlr.hg19.merged.bam'
		else :
			sample = os.path.basename(line0[1])
		d['samples'][sample] = line0[2]
		d.setdefault('S',[]).append(sample)
	f.close()
	return d

def readmerge(File,prefix):
	dmerge = {}
	with open(File) as f:
		for line in f:
			line0=line.strip().split('\t')  ## svID
			dmerge.setdefault(line0[0],{})#.append(line0[1:9])
			#if 'ngmlr/' in line0[7]: 
			#	line0[7]=line0[7].replace('ngmlr/',prefix+'.')
			sample = os.path.basename(line0[7])
			dmerge[line0[0]][sample] = line0[1:]
	f.close()
	return dmerge

def get_median(data):
	data = [int(i) for i in data]
	data.sort()
	half = len(data) // 2
	return (data[half] + data[~half]) / 2

def getvcf(dmerge,dsamples,dsamplevcfinfo,outfile):
	#samples = sorted(dsamples)
	if not dsamples:
		samples = dsamplevcfinfo['samples'].keys()
		newsamples = samples
	else:
		samples = dsamples['S']
		newsamples = [dsamples['samples'][k] for k in samples]
	out1 = ''
	#out2 = ''
	for k in dsamplevcfinfo['header']:
		out1 += k+'\n'
		#out2 += k+'\n'
	out1 += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(newsamples)+'\n'
	#out2 += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(newsamples)+'\t'+'\t'.join(newsamples)+'\n'
	for i in sorted(dmerge): ## mergeID
		start = []
		end = []
		gt = []
		svlength = []
		SUPP_VEC=''
		svtype = ''
		RNAMES = []
		for j in samples:
			if j in dmerge[i]:
				chr1 = dmerge[i][j][0]
				chr2 = dmerge[i][j][2]
				svtype = dmerge[i][j][4]
				start.append(dmerge[i][j][1])
				end.append(dmerge[i][j][3])
				svlength.append(dmerge[i][j][5])
				ID = dmerge[i][j][-1]
				#gtsample = dsamplevcfinfo['samples'][j][dmerge[i][j][-1]]
				gtsample = dsamplevcfinfo['samples'][j][ID][1]
				gt.append(gtsample)
				RNAMES.append(dsamplevcfinfo['samples'][j][ID][0])
				if '1' in gtsample.split(':')[0] or './.' in gtsample.split(':')[0] :
					SUPP_VEC+='1'
				else:
					SUPP_VEC+='0'
			else:
				gt.append('0/0:NaN:0')
				RNAMES.append('')
				SUPP_VEC+='0'
		if svtype == "TRA":svlength='0'
		if svtype == "INS":END = str(get_median(start)+1)
		if svlength==[]:
			continue
		if svtype == "INS" or svtype == "TRA":
			AVGLEN = str(int(sum([int(k) for k in svlength])/float(len(svlength))))
		else:
			AVGLEN = str(int(get_median(end) - get_median(start)))
		INFO = "SUPP=%s;SUPP_VEC=%s;SVLEN=%s;SVTYPE=%s;CHR2=%s;END=%s;CIPOS=NA,NA;CIEND=NA,NA;STRANDS=++" %(str(len(start)),SUPP_VEC,AVGLEN,svtype,chr2,str(get_median(end)))
		out1 += '\t'.join([chr1,str(get_median(start)),i,'N',svtype,'.','pass',INFO,'GT:DR:DV'])+'\t'+'\t'.join(gt)+'\n'
		#out2 += '\t'.join([chr1,str(get_median(start)),i,'N',svtype,'.','pass',INFO,'GT'])+'\t'+'\t'.join(gt)+'\t'+'\t'.join(RNAMES) +'\n'
		
	print >>open(outfile,"w"),out1.strip()
	#print >>open(outfile+'.tmp',"w"),out2.strip()

def main(args):
	mergevcf = args[0].merge
	vcflist = args[0].vcflist
	sampleinfo = args[0].sampleinfo
	change = args[0].changesample
	outfile = args[0].outfile
	if change == False or change == "False":change = False
	if sampleinfo != None:
		prefix = os.popen("awk 'NR==1{print $1}' %s" %(sampleinfo)).readlines()[0].strip()
		dsamples = getsamples(sampleinfo,change)
	else:
		prefix = False
		dsamples = False
	dmerge = readmerge(mergevcf,prefix)
	dsamplevcfinfo=getsamplevcfinfo(vcflist,prefix)
	getvcf(dmerge,dsamples,dsamplevcfinfo,outfile)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		formatter_class=HelpFormatter,
                description='''
	description:
	run:
		mergevcf2vcf.py --merge allmerge --vcflist allvcflist --sampleinfo sampleinfo -o out.vcf
		mergevcf2vcf.py --merge allmerge --vcflist allvcflist --sampleinfo sampleinfo --changesample False -o out.vcf
	''')

	parser.add_argument('-m', '--merge', metavar='FILE', required=True, help='set the input vcfmerge file.')
	parser.add_argument('-v', '--vcflist', metavar='FILE', required=True, help='set the input vcfmergelist file')
	parser.add_argument('-s', '--sampleinfo', metavar='FILE', help='set the sample information file. not must')
	parser.add_argument('-c', '--changesample', metavar='True/False', default=True,
		help='T:change(prefix+sampleinfo[1]+"ngmlr.hg19.merged.bam"-->sampleinfo[2]); F:sampleinfo[1]--sampleinfo[2]')
	parser.add_argument('-o', '--outfile', metavar='FILE', help='set the outfile.')
	args = parser.parse_known_args()
	main(args)
			
	
