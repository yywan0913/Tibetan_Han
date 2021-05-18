import sys
vcf = sys.argv[1]

with open(vcf) as f:
	for line in f:
		if '##' in line:
			continue
		line0 = line.strip().split('\t')
		if '#' in line[0]:
			print line0[2]+'\t'+'\t'.join(line0[9:])
			continue
		value = ['1' if '1' in i.split(':')[0] else '0' for i in line0[9:]]
		print line0[2]+'_'+line0[4]+'\t'+'\t'.join(value)
f.close()
