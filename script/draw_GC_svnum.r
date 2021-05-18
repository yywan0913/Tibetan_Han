args = commandArgs(TRUE)
bedgc = args[1]
ALHT.DEL = args[2]
ALHT.INS = args[3]
EAS.DEL = args[4]
EAS.INS = args[5]
outpdf = args[6]

options(stringsAsFactors=F)
ref=read.table("windows.bed.gc.txt",header=F,sep="\t")
ref[,2] = ref[,2]-1
ref[,1] = paste0('chr',ref[,1])
colnames(ref) = c('chr','start','end','gc')
ref = ref[,-3]
AL.DEL = read.table("DEL.splide.counts.xls",header=T,sep="\t")
AL.DEL = AL.DEL[,-3]
mer = merge(ref,AL.DEL,by=c("chr",'start'),all=T)
colnames(mer)[ncol(mer)] = 'AL-DEL'

AL.INS = read.table("INS.splide.counts.xls",header=T,sep="\t")
AL.INS = AL.INS[,-3]
mer = merge(mer,AL.INS,by=c("chr",'start'),all=T)
colnames(mer)[ncol(mer)] = 'AL-INS'

EAS.DEL = read.table("EAS/DEL.splide.counts.xls",header=T,sep="\t")
EAS.DEL = EAS.DEL[,-3]
mer = merge(mer,EAS.DEL,by=c("chr",'start'),all=T)
colnames(mer)[ncol(mer)] = 'EAS-DEL'

EAS.INS = read.table("EAS/INS.splide.counts.xls",header=T,sep="\t")
EAS.INS = EAS.INS[,-3]
mer = merge(mer,EAS.INS,by=c("chr",'start'),all=T)
colnames(mer)[ncol(mer)] = 'EAS-INS'

mer[is.na(mer)] = 0

out = mer[,-c(1:2)]
class = cut(out$gc,breaks=seq(0,1,by=0.05))
out$gc = as.character(class)
library(reshape2)
mel = melt(out,id="gc")
mel[is.na(mel)]="(0,0.05]"
dca=dcast(mel,gc~variable,sum)   ## col=5

nr = nrow(dca)
rect.width = 0.4
pdf('DEL--INS.splide.gc.distribution.bar.pdf')
plot(1,type="n",xlab="",ylab="",xlim=c(0,nrow(dca)+1),ylim=c(0,max(apply(dca[,2:3],1,sum))),bty="l",xaxt="n")
rect(1:nr-rect.width,0,1:nr,dca[,2],col="darkred",border=NA)
rect(1:nr-rect.width,dca[,2],1:nr,dca[,2]+dca[,3],col="darkgreen",border=NA)

rect(1:nr,0,1:nr+rect.width,dca[,4],col=rgb(50/255,50/255,50/255),border=NA)
rect(1:nr,dca[,4],1:nr+rect.width,dca[,4]+dca[,5],col="grey",border=NA)

axis(1,1:nr,dca[,1])
legend("topleft",legend=c('AL-HT_DEL','AL-HT_INS','EAS-DEL','EAS-INS'),col=c('darkred','darkgreen',rgb(50/255,50/255,50/255),'grey'),bty="n",pch=15,cex=1.5)
dev.off()

pdf('DEL--INS.splide.gc.distribution.heatmap.pdf')
par(mar=c(4,6,4,4),las=1,mgp=c(3,0.5,0))
plot(1,type="n",xlab="",ylab="",xlim=c(0,4),ylim=c(0,nr),bty="n",xaxt="n",yaxt="n")
data = dca[,-1]
data = apply(data,2,function(x)(x-mean(x))/sd(x))
DEL.col = colorRampPalette(c('yellow','red'))(nr*2)
del = c(data[,1],data[,3])
del.col = DEL.col[as.numeric(as.factor(del))]
INS.col = colorRampPalette(c('blue','green'))(nr*2)
ins = c(data[,2],data[,4])
ins.col = INS.col[as.numeric(as.factor(ins))]
rect(0,0:(nr-1),1,1:nr,col=del.col[1:nr],border=NA)
rect(1,0:(nr-1),2,1:nr,col=del.col[(nr+1):(nr*2)],border=NA)
rect(2,0:(nr-1),3,1:nr,col=ins.col[1:nr],border=NA)
rect(3,0:(nr-1),4,1:nr,col=ins.col[(nr+1):(nr*2)],border=NA)
#legend(par('usr')[2],par('usr')[4],legend=c())
axis(3,seq(0.5,4,by=1),c('AL-HT_DEL','EAS-DEL','AL-HT_INS','EAS-INS'),tick=F)
axis(2,1:nr-0.5,dca[,1],tick=F)
dev.off()
