args =commandArgs(TRUE)
file = args[1]
outpdf = args[2]
#file = "data.input.nogap.xls"
library(circlize)
data = read.table(file,sep="\t",header=T)
main.data = data[,5:8]
#main.data = apply(main.data,2,function(x)(x-mean(x))/sd(x))

## hg19
plotcytoBand = function(chrom,bm=1000000,start,d=0.15,D=0.5,cex=1.6,par=T,y1,y2){
        #file = system.file(package = "circlize","extdata", "cytoBand.txt")
	file = "cytoBand.txt"
        options(stringsAsFactors=F)
        cyto = read.table(file,header=F,sep="\t")
        colnames(cyto)<-c("chr","start","end","name","type")
        #cyto$start = cyto$start / bm
        #cyto$end = cyto$end / bm
        cyto$Color = "white"  # white
        cyto[cyto$type=="gneg",]$Color<-rgb(255,255,255, maxColorValue=255)
        cyto[cyto$type=="gpos25",]$Color<-rgb(200,200,200, maxColorValue=255)
        cyto[cyto$type=="gpos50",]$Color<-rgb(200,200,200, maxColorValue=255)
        cyto[cyto$type=="gpos75",]$Color<-rgb(130,130,130, maxColorValue=255)
        cyto[cyto$type=="gpos100",]$Color<-rgb(200,200,200, maxColorValue=255)
        cyto[cyto$type=="acen",]$Color<-"yellow" # red centromere 
        cyto[cyto$type=="stalk",]$Color<-rgb(100,127,164, maxColorValue=255) # repeat regions 
        cyto[cyto$type=="gvar",]$Color<-rgb(220,220,220, maxColorValue=255) # indented region
	
	cyto$chr = gsub('chr','',cyto$chr)
	cyto = cyto[cyto$chr==as.character(chrom),,drop=F]
	bei = max(cyto$end)
	cyto$start = cyto$start/bei*bm+start
	cyto$end = cyto$end/bei*bm+start

	acen.w = which(cyto$type=="acen")
	cyto.up = cyto[1:(acen.w[1]-1),,drop=F]
	cyto.acen1 = cyto[acen.w[1],,drop=F]
	cyto.acen2 = cyto[acen.w[2],,drop=F]
	cyto.down = cyto[(acen.w[2]+1):nrow(cyto),,drop=F]
	rect(min(cyto.up$start),y1,max(cyto.up$end),y2,col=NA,border="black",lwd=2)
	rect(cyto.up$start,y1,cyto.up$end,y2,col=cyto.up$Color,border=NA)
	
	rect(min(cyto.down$start),y1,max(cyto.down$end),y2,col=NA,border="black",lwd=2)
	rect(cyto.down$start,y1,cyto.down$end,y2,col=cyto.down$Color,border=NA)

	polygon(x=c(cyto.acen1$start[1],cyto.acen1$end[1],cyto.acen1$start[1],cyto.acen1$start[1]),y=c(y1,(y1+y2)/2,y2,y1),col=2,border=2)
	polygon(x=c(cyto.acen2$end[1],cyto.acen2$start[1],cyto.acen2$end[1],cyto.acen2$end[1]),y=c(y1,(y1+y2)/2,y2,y1),col=2,border=2)

}

rmsk.max = max(data[,5])
SINE.genome =  data[,5]/rmsk.max
LINE.genome = data[,6]/rmsk.max
LTR.genome = data[,7]/rmsk.max
SD.genome = data[,8]/rmsk.max
SD.genome[SD.genome>1]=1

chromlen = table(data[,1])
fenmu = min(max(data[,9]),max(data[,10]),max(data[,11]),max(data[,12]))
sine.sv = data[,9]/fenmu
sine.sv[sine.sv>1]=1
line.sv = data[,10]/fenmu
line.sv[line.sv>1]=1
ltr.sv = data[,11]/fenmu
ltr.sv[ltr.sv>1]=1
sd.sv=data[,12]/fenmu
sd.sv[sd.sv>1]=1
allsv = data[,4]/fenmu
allsv[allsv>1]=1
gene.density = data[,13]/max(data[,13])

svtype.data = data[,14:18]
svtype.data = svtype.data/fenmu
del.data = ifelse(svtype.data[,1]>1,1,svtype.data[,1])
dup.data = ifelse(svtype.data[,2]>1,1,svtype.data[,2])
ins.data = ifelse(svtype.data[,3]>1,1,svtype.data[,3])
inv.data = ifelse(svtype.data[,4]>1,1,svtype.data[,4])
tra.data = ifelse(svtype.data[,5]>1,1,svtype.data[,5])

els.data = data[,19:20]
els.data[,1] = els.data[,1]/max(els.data[,1])
els.data[,2] = els.data[,2]/max(els.data[,2])
hicdamin.data = data[,21:23]
hicdamin.data[,1] = hicdamin.data[,1]/max(hicdamin.data[,1])
hicdamin.data[,2] = hicdamin.data[,2]/max(hicdamin.data[,2])
hicdamin.data[,3] = hicdamin.data[,3]/max(hicdamin.data[,3])
##heatmap.col = c(rgb(31/255,49/255,94/255),rgb(34/255,103/255,179/255),rgb(63/255,148/255,198/255),rgb(148/255,197/255,224/255),rgb(212/255,229/255,245/255),rgb(251/255,250/255,251/255),rgb(252/255,220/255,200/255),rgb(247/255,165/255,130/255),rgb(219/255,95/255,79/255),rgb(180/255,30/255,48/255),rgb(106/255,9/255,36/255))
col.breaks = c(min(main.data),-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,max(main.data))
##f = colorRamp2(breaks = col.breaks, colors = heatmap.col)
##heatmap.data = f(unlist(main.data))
##heatmap.data = matrix(heatmap.data,nc=ncol(main.data))

#chrom_col=c("1"="#a2b9bc","2"="#b2ad7f","3"="#878f99","4"="#6b5b95","5"="#feb236","6"="#d64161","7"="#ff7b25","8"="#d6cbd3","9"="#eca1a6","10"="#bdcebe","11"="#ada397","12"='#66c2a5',"13"='#fc8d62',"14"='#8da0cb',"15"='#e78ac3',"16"='#a6d854',"17"='#ffd92f',"18"='#e5c494',"19"='#b3b3b3',"20"="#20854E","21"="#7876B1","22"="#6F99AD","23"="#FFDC91","24"="#EE4C97")
##TE.col = c('#E64B35', '#4DBBD5', '#00A087', '#3C5488')
#cluster_col = c("#b2ad7f","#878f99",'#e78ac3',"#feb236","#d64161","#ff7b25","#d6cbd3")
#cluster_col = c("#009CFF", "#FF5500", "#763400","#9E819B","#D20000", "#FFC600")
chrom_col=c("1"="#a2b9bc","2"="#b2ad7f","3"="#878f99","4"="#6b5b95","5"="#feb236","6"="#d64161","7"="#ff7b25","8"="#d6cbd3","9"="#eca1a6","10"="#bdcebe","11"="#ada397","12"='#66c2a5',"13"='#fc8d62',"14"='#8da0cb',"15"='#e78ac3',"16"='#a6d854',"17"='#ffd92f',"18"='#e5c494',"19"='#b3b3b3',"20"="#20854E","21"="#7876B1","22"="#6F99AD","23"="#FFDC91","24"="#EE4C97")

pdf(outpdf,width=30,height=20)
layout(mat=matrix(c(1,2,3),nc=1),height=c(1.5,9,2))


## k label (2)
par(mar=c(1,2,5,35)) # 1,8
plot(1,type="n",xlim=c(0.5,nrow(main.data)+10.5),ylim=c(0,2),xlab="",ylab="",axes=F,xaxs="i")
#chr.uniq = unique(data[,1])
start = 0.5
start.index = 1
k = 22
chrom.seq = c()
bgdata = matrix(0,nc=2,nr=11)
for(i in 1:k){
	cluster.num = i
	nx = sum(data[,1]==cluster.num)
	chrom.seq = c(chrom.seq,start+nx)
	#rect(start,0.1,start+nx,1.9,col=chrom_col[i],border=NA)
	if(i%%2==1) {
		bgdata[(i+1)/2,1] = start
		bgdata[(i+1)/2,2] = start+nx
		rect(start,1.1,start+nx,1.9,col='grey80',border=NA)
		text((start+start+nx)/2,1.5,cluster.num,col="black",cex=3)
	}else{
		rect(start,1.1,start+nx,1.9,col='grey30',border=NA)
		text((start+start+nx)/2,1.5,cluster.num,col="white",cex=3)
	}
	
	plotcytoBand(chrom=i,bm=nx,start,d=0.15,D=0.5,cex=1.6,par=T,y1=0.1,y2=0.9)

	start = start+nx
	start.index = start+0.5
}
require("RColorBrewer")
line.col = colorRampPalette(brewer.pal(9, "Paired"))(8)
line.col[1] = 'purple'
line.col[4] = 'Cyan1'
par(mar=c(8,2,0,35),cex.axis=5,las=1,mgp=c(0.5,0.1,0),cex.lab=5,font.lab=2)
y2 = 7+4+2+3+4
plot(1,type="n",xlim=c(0.5,nrow(main.data)+13.5),ylim=c(-0.5,y2),ylab="",axes=F,xlab="",xaxs="i",yaxs="i")
rect(bgdata[,1],-0.5,bgdata[,2],y2,col=rgb(204/255,204/255,204/255,0.6),border=NA)
lines(1:nrow(data),y2-1+gene.density,col = line.col[1])
lines(1:nrow(data),y2-2+SD.genome,col = line.col[2])
lines(1:nrow(data),y2-3+LTR.genome,col = line.col[2])
lines(1:nrow(data),y2-4+LINE.genome,col = line.col[2])
lines(1:nrow(data),y2-5+SINE.genome,col = line.col[2])
lines(1:nrow(data),y2-6+sd.sv,col = line.col[3])
lines(1:nrow(data),y2-7+ltr.sv,col = line.col[3])
lines(1:nrow(data),y2-8+line.sv,col = line.col[3])
lines(1:nrow(data),y2-9+sine.sv,col = line.col[3])
lines(1:nrow(data),y2-10+allsv,col = line.col[4])
lines(1:nrow(data),y2-11+del.data,col = line.col[4]) # DEL
lines(1:nrow(data),y2-12+dup.data,col = line.col[4]) # DUP
lines(1:nrow(data),y2-13+ins.data,col = line.col[4]) # INS
lines(1:nrow(data),y2-14+inv.data,col = line.col[4]) # INV
lines(1:nrow(data),y2-15+tra.data,col = line.col[4]) # TRA
lines(1:nrow(data),y2-16+els.data[,1],col = line.col[5])
lines(1:nrow(data),y2-17+els.data[,2],col = line.col[6])
lines(1:nrow(data),y2-18+hicdamin.data[,1],col = line.col[7])
lines(1:nrow(data),y2-19+hicdamin.data[,2],col = line.col[8])
lines(1:nrow(data),y2-20+hicdamin.data[,3],col = line.col[8])
#chrom.seq = chrom.seq[-length(chrom.seq)]
#abline(v=chrom.seq,col="black",lwd=0.1)
par(xpd=T)
axis(4,at=(y2-1):0+0.5,c('Gene Density','SD(genome)','LTR(genome)','LINE(genome)','SINE(genome)','SD(SV)','LTR(SV)','LINE(SV)','SINE(SV)','ALL SV Density','DEL Density','DUP Density','INS Density','INV Density','TRA Density','Enhancer','Promoter','HiC(boundary)','LAMIN(boundary)','LAMIN(domains)'),tick=F)

par(mar=c(0,1,2.5,35))
plot(1,type="n",axes=F,xlab="",ylab="")
legend('top',legend=c('Gene Density','Repeat Coverage(genome)','Repeat (SV Density)','SV Density','Enhancer','Promoter','HiC-boundary','LAMIN'),col=line.col,bty="n",ncol=4,lty=1,cex=3,lwd=4)
dev.off()
