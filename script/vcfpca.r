args = commandArgs(TRUE)
groupvcf = args[1]
samplelist = args[2]
#genome1000groupFile="/export/home/xiaoyh/Pipeline/SVanalysis/database/genome1000groupSamplelist"
genome1000groupFile = args[3]
outdir = args[4] 
if(is.na(outdir)) outdir ="./"

library(data.table)
Groupvcf = fread(groupvcf,header=T)
Groupvcf = as.data.frame(Groupvcf)
data = Groupvcf[,-1]
MAF = apply(data,1,function(x)min(table(x)/length(x)))
data = data[MAF>1/ncol(data),]
Groupvcfpos = Groupvcf[MAF>1/ncol(data),1]

genome1000group = read.table(genome1000groupFile,sep="\t",header=T)
colnames(genome1000group) = c('V1','V2')
samplegroup = read.table(samplelist,sep="\t",header=F)
Group = rbind(samplegroup,genome1000group)
intersect = intersect(colnames(data),Group$V1)
data = data[,match(intersect,colnames(data))]
Group = Group[match(intersect,Group$V1),]

group_col = rainbow(length(unique(Group[,2])),s=0.6,v=0.6)
col = group_col[as.numeric(as.factor(Group[,2]))]
group = levels(as.factor(Group[,2])) 
pch=c(15,16,17,18,25,13,10,8)
Pch = pch[as.numeric(as.factor(Group[,2]))]
group_pch = pch[1:length(group)]
pcat <- prcomp(data, scale = TRUE)$rotation
rownames(pcat) = colnames(data)
write.table(pcat,paste0(outdir,'/ALL.pca.xls'),sep="\t",col.names=T,row.names=T,quote=F)
 
plot2d=function(x,y,labx=NULL,laby=NULL,col,group,group_col,outpdf,samples,pch,group_pch){
	pdf(outpdf,height=10,width=12)
	par(las=1,mar=c(4,4,4,6))
	min_x = min(x)
	max_x = max(x)
	min_y = min(y)
	max_y = max(y)
	if(is.null(labx)) labx = "PC1"
	if(is.null(laby)) laby = "PC2"
	plot(1,xlim=c(min_x,max_x),ylim=c(min_y,max_y),type="n",xlab=labx,ylab=laby,main="")
	points(x=x,y=y,col=col,pch=pch)
	if(length(samples)<50) text(x,y,samples,cex=0.8,pos=3)
	abline(v=0,col="grey")
	abline(h=0,col="grey")
	if(length(samples)<50) legend(par('usr')[2],sum(par('usr')[3:4])/2,legend=group,col=group_col,pch=group_pch,bty="n",xpd=T) else
		legend(par('usr')[2],par('usr')[4],legend=group,col=group_col,pch=group_pch,bty="n",xpd=T)
	dev.off()
}
deig=eigen(cor(scale(data)))
k = sum(deig$values)
pc1_rate=sprintf("%2.2f",deig$values[1]/k)
pc2_rate=sprintf("%2.2f",deig$values[2]/k)
pc3_rate=sprintf("%2.2f",deig$values[3]/k)
lab1 = paste("PC1",pc1_rate)
lab2 = paste("PC2",pc2_rate)
lab3 = paste("PC3",pc3_rate)
plot2d(x=pcat[,1],y=pcat[,2],col=col,group=group,group_col=group_col,outpdf=paste0(outdir,"/g1000.pca12.pdf"),samples=rownames(pcat),labx=lab1,laby=lab2,group_pch=group_pch,pch=Pch)
plot2d(x=pcat[,1],y=pcat[,3],col=col,group=group,group_col=group_col,outpdf=paste0(outdir,"/g1000.pca13.pdf"),samples=rownames(pcat),labx=lab1,laby=lab3,group_pch=group_pch,pch=Pch)
plot2d(x=pcat[,2],y=pcat[,3],col=col,group=group,group_col=group_col,outpdf=paste0(outdir,"/g1000.pca23.pdf"),samples=rownames(pcat),labx=lab2,laby=lab3,group_pch=group_pch,pch=Pch)


