readvcf <- function(file){
	library(data.table)
	df <- fread(file)
	df <- as.data.frame(df)
	df <- df[,-c(3:9)]
	colnames(df)[1]="CHROM"
	colnames(df)=gsub('(.*).ngmlr.*','\\1',colnames(df))
	df[,-(1:2)] = apply(df[,-(1:2)],1:2, gsub ,pattern=':.*',replacement='')
	dfdata <- df[,-(1:2)]
	dfdata <- ifelse(dfdata=="1/1"|dfdata=="0/1",1,0)
	data = data.frame(df[,1:2],dfdata,check.names=F)
	return(data)
}
library(data.table)
args=commandArgs(TRUE)
vcfFile = args[1]
samplelist=args[2]
#outdir = args[3]
prefix=args[3]
outtre = paste0(prefix,".tre")
#vcf = readvcf(vcfFile)
vcf = fread(vcfFile,header=T)
vcf = as.data.frame(vcf)
data = vcf[,-1]
MAF = apply(data,1,sum)
data = data[MAF>1,]
hc = hclust(dist(t(data),method = "euclidean"), method = "complete") ## complete  ave  ## euclidean  binary
merge = hc$merge
labels = hc$labels
n = nrow(merge)
list = list()
for(i in 1:n){
	if(merge[i,1]<0&merge[i,2]<0)list[[i]]=paste("(",labels[-1*merge[i,1]],",",labels[-1*merge[i,2]],")",sep="");
	if(merge[i,1]<0&merge[i,2]>0)list[[i]]=paste("(",labels[-1*merge[i,1]],",",list[[merge[i,2]]],")",sep="");
	if(merge[i,1]>0&merge[i,2]<0)list[[i]]=paste("(",list[[merge[i,1]]],",",labels[-1*merge[i,2]],")",sep="");
	if(merge[i,1]>0&merge[i,2]>0)list[[i]]=paste("(",list[[merge[i,1]]],",",list[[merge[i,2]]],")",sep="")
}
tre = paste(list[[n]],";",sep="")
write.table(tre,outtre,col.names=F,row.names=F,quote=F)

options(stringsAsFactors=F)
library(ggtree)
pdf(paste0(outtre,'.pdf'),width=10,height=10)
tree<-read.tree(outtre)
samplegroup = read.table(samplelist,sep="\t",header=F)
samplegroup = samplegroup[match(tree$tip.label,samplegroup[,1]),]
group = samplegroup[,2]
cls = split(samplegroup[,1],group)
tree = groupOTU(tree,cls)
print(tree$group)
ggtree(tree,layout="fan",ladderize=F,branch.length="none",aes(color=group))+geom_tiplab(aes(angle=angle),size=2)+theme(legend.position='right')
dev.off()
