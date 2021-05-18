library(data.table)
library(reshape2)
args = commandArgs(TRUE)
outpdf = args[1]
region.num = as.numeric(args[2])
if(is.null(region.num)) region.num = 3

chromfile = "input/chrom.len"
svfile = "input/sv.fst.fordraw.xls"
snpfile = "input/snp.fst.fordraw.filter.xls"
indelfile = "input/indel.fst.fordraw.xls"
#outpdf = "union.manhattan.pdf"
genegfffile = "input/hg19.bed.gtf.gz"
genegff = gettextf("gzip -dc %s",genegfffile)
posfile = "pos.txt"

options(stringsAsFactors=F)
chrom.data = read.table(chromfile,header=F,sep="\t")
data1 = fread(svfile,header=T,sep="\t")
data1 = as.data.frame(data1)
data2 = fread(snpfile,header=T,sep="\t")
data2 = as.data.frame(data2)
data3 = fread(indelfile,header=T,sep="\t")
data3 = as.data.frame(data3)
gene.data = fread(genegff,header=F,sep="\t")
gene.data = as.data.frame(gene.data)
gene.data[,1] = gsub('chr','',gene.data[,1])

pos.data = read.table(posfile,header=F,sep="\t",fill=F)

var.color =c("#E64B35","#00A087","#3C5488")
pch = 15:17


if(grepl('.png$',outpdf)) png(outpdf,width=24,height=30, units = "in",res=400)
if(grepl('.pdf$',outpdf)) pdf(outpdf,width=24,height=30)

if(region.num==3){
	layout(mat=matrix(c(1,1,1,2,3,4,5,7,9,6,8,10),ncol=3,byrow=T),height=c(10,10,2,8))
}
if(region.num==2){
	layout(mat=matrix(c(1,1,2,3,4,6,5,7),ncol=2,byrow=T),height=c(11,11,2,8))
}

# up 
source('source/source.up.manhattan.r')
par(mar=c(3,5,3,4))
print(11)
up.manhattan(chrom.data,sv.data=data1,snp.data=data2,indel.data=data3,col=var.color,pos.data=pos.data)

# median
source('source/source.median.manhattan.r')
par(mar=c(3,5,3,4))
for(i in 1:region.num){
	pos.datai = pos.data[i,]
	print('22')
	median.manhanttan(pos.datai,gene.data=gene.data,sv.data=data1,snp.data=data2,indel.data=data3,color=var.color,pch=pch)
}
# bottom 
source('source/source.bottom.r')
AL="input/AL"
HT="input/HT"
Tibetan= scan(AL,what="")
Han = scan(HT,what="")
gene.col = c("#009CFFFF","#FF5500FF","#D20000FF","#9E819BFF","#763400FF","#FFC600FF","#54FFAAFF",'yellow')
set.seed(123)
for(i in 1:region.num){
	posi.file = pos.data[i,7]
	data = read.table(posi.file,sep="\t",check.names=F,header=T)
	sampleAL = intersect(Tibetan,colnames(data))
	sampleHT = intersect(Han,colnames(data))
	sample20AL = sample(sampleAL,30)
	sample20HT = sample(sampleHT,30)
	matchsampleAL = match(sample20AL,colnames(data))
	matchsampleHT = match(sample20HT,colnames(data))
	data.AL = data[,matchsampleAL]
	data.HT = data[,matchsampleHT]
	data.info = data[,1:3]
	PlotRegion(data.info,data.AL,data.HT,gene.data,gene.col)
}
dev.off()
