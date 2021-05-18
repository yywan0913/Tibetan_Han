args = commandArgs(TRUE)
k2order = args[1] # k2.structure.order.xls
FSTmatrix = args[2] # FST0.1.AL--HT.weir.fst.matrix.xls
outpdf = args[3] # structure.fst0.1.diff.complete.pdf
options(stringsAsFactors=F)
sampleorder = read.table(k2order,sep="\t",header=F)
sample = sampleorder[,1]
data = read.table(FSTmatrix,header=T,sep="\t",check.names=F)
data = data[,match(sample,colnames(data))]
ALcounts = apply(data[,1:119],1,sum)/119 + 1e-10
HTcounts = apply(data[,120:ncol(data)],1,sum)/201 + 1e-10
rate = ALcounts/HTcounts

data = data[rate>1.5|rate<2/3,]
Rate = rate[rate>1.5|rate<2/3]
data = data[order(Rate),]
col = "#FC4E07"
#color = c("#00AFBB", "#E7B800")

library(pheatmap)
pdf(outpdf,width=10)
pheatmap(mat=as.matrix(data),show_rownames=F,show_colnames=F,cluster_cols=F,color=c('white',col),legend=F,treeheight_row=60)
dev.off()

