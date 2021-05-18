library(corrplot)
library(psych)
library(RColorBrewer)
args = commandArgs(TRUE)
input = args[1]
outpdf = args[2]
df = read.table(input,header=T,check.names=F)
df = df[c(13,8,7,6,5,12,11,10,9,4,14,15,16,17,18,19,20,21,22,23)]
colnames(df) = c('Gene Density','SD(genome)','LTR(genome)','LINE(genome)','SINE(genome)','SD(SV)','LTR(SV)','LINE(SV)','SINE(SV)','ALL SV Density','DEL Density','DUP Density','INS Density','INV Density','TRA Density','Enhancer','Promoter','HiC(boundary)','LAMIN(boundary)','LAMIN(domains)')
cor = cor(df)
cor_p = corr.p(cor,n=nrow(df),alpha=.05)
#cor[lower.tri(cor_p$p)] = cor_p$p[lower.tri(cor_p$p)]
#write.table(cor,"TE-SV.corAndP.spearman.xls",col.names=T,row.names=T,quote=F,sep="\t")
res1 = cor.mtest(df,conf.level=.95)
pdf(outpdf,width=15,height=15)
#corrplot.mixed(cor,p.mat=res1$p,lower = "number",upper = "circle",
#               tl.pos = "lt",lower.col = "black",insig='label_sig',sig.level=c(.001,.01,.05),pch.cex=0.6)
diag(cor)=1
diag(res1$p)=1
col1 = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
col = colorRampPalette(c(rgb(0,0,0.8),rgb(0.8,0,0)))(100)

corrplot(cor, p.mat=res1$p,method = "circle", type = "upper", tl.pos = "lt",insig='label_sig',sig.level=c(.001,.01,.05),pch.cex=0.6,diag=TRUE,col=col,pch.col="white")
corrplot(cor, add = TRUE, type = "lower", method = "number", diag = FALSE, tl.pos = "n", cl.pos = "n",col=col)
dev.off()

