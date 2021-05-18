
args = commandArgs(TRUE)
input = args[1]
#input = "ALL.pca.info.xls"
outpdf = args[2]
library(data.table)
data = fread(input,sep="\t",header=T)
data = as.data.frame(data)
group_col = c("#009CFF", "#FF5500", "#763400","#9E819B","#D20000", "#FFC600","#54FFAA")
#group_col = apply(col2rgb(group_col)/255,2,function(x)rgb(x[1],x[2],x[3],0.6))
#Group = data$Group
data$Group = as.factor(data$Group)
group.levels = levels(data$Group)
sample.col = group_col[as.numeric(data$Group)]
names(group_col) = group.levels
draw.data = data[,-c(1:2)]

pc.rowdata = fread("usandgenotypes.1000genomeintersect.noXY.matrix.filter",header=T,sep="\t")
pc.rowdata = as.data.frame(pc.rowdata)
pc.rowdata = pc.rowdata[,-1]
MAF = apply(pc.rowdata,1,function(x)min(table(x)/length(x)))
pc.rowdata = pc.rowdata[MAF>1/ncol(pc.rowdata),]
deig=eigen(cor(scale(pc.rowdata)))
k = sum(deig$values)
pc1_rate=sprintf("%2.2f%%",deig$values[1]/k*100)
pc2_rate=sprintf("%2.2f%%",deig$values[2]/k*100)
pc3_rate=sprintf("%2.2f%%",deig$values[3]/k*100)
pc4_rate=sprintf("%2.2f%%",deig$values[4]/k*100)
pc5_rate=sprintf("%2.2f%%",deig$values[5]/k*100)
pc6_rate=sprintf("%2.2f%%",deig$values[6]/k*100)
lab1 = paste("PC1,",pc1_rate)
lab2 = paste("PC2,",pc2_rate)
lab3 = paste("PC3,",pc3_rate)
lab4 = paste("PC4,",pc4_rate)
lab5 = paste("PC5,",pc5_rate)
lab6 = paste("PC6,",pc6_rate)
PClabel = c(lab1,lab2,lab3,lab4,lab5,lab6)

plot1 =function(x,group){
	group.data = split(x,group)
	grouphistinfo = lapply(group.data,hist,plot=F,breaks=10)
	y.max = max(unlist(lapply(grouphistinfo,function(x)max(x$density))))
	plot(1,xlim=c(min(x),max(x)),ylim=c(0,y.max),type="n",xlab="",ylab="",xaxt="n",yaxt="n")
	grid()
	lapply(1:length(grouphistinfo),function(x){
		groupi = names(grouphistinfo)[x]
		grouphistinfoi = grouphistinfo[[x]]
		polygon(x=c(grouphistinfoi$breaks,grouphistinfoi$breaks[length(grouphistinfoi$breaks)]),
			y=c(0,grouphistinfoi$density,0),col = group_col[groupi],border=NA)
	})
}

plot2 = function(x,y,col){
	plot(x,y,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
	grid()
	points(x,y,col=col,pch=16,cex=0.9)
}

main.plot = function(mar=mar){
	for(i in 1:6){
		for(j in i:6){
			par(mar=mar)
		if(i==j){
			plot1(draw.data[,i],group=data$Group)
		}else{
			plot2(x=draw.data[,i],y=draw.data[,j],col=sample.col)
		}
	}
}
}
#outpdf = "g1000.PC1-PC6.pdf"
png(outpdf,width=15,height=15,units="cm",res=400)
##pdf(outpdf,width=15,height=15)
#mat = matrix(c(22,1:6,23,0,7:11,24,0,0,12:15,25,0,0,0,16:18,26,0,0,0,0,19:20,27,0,0,0,0,0,21,0,28,29,30,31,32,33),nc=7)
mat = matrix(c(22,1:6,23,0,7:11,24,0,0,12:15,25,0,0,0,16:18,26,34,34,34,0,19:20,27,34,34,34,0,0,21,0,28,29,30,31,32,33),nc=7)
layout(mat=mat,height = c(2,5,5,5,5,5,5),width=c(5,5,5,5,5,5,2))
# main
mar = c(0.5,0.5,0.5,0.5)
main.plot(mar=mar)
# col pc1 ---pc6  label
for(i in 1:6){
	par(mar=mar)
	plot(1,xlim=c(0,1),ylim=c(0,1),type="n",axes=F,xlab="",ylab="",xaxt="n",yaxt="n")
	rect(0,0,1,1,col="grey80",border="grey80")
	text(0.5,0.5,PClabel[i],cex=1)
}
# row PC1--PC6 label
for(i in 1:6){
        par(mar=mar,srt=90)
        plot(1,xlim=c(0,1),ylim=c(0,1),type="n",axes=F,xlab="",ylab="",xaxt="n",yaxt="n")
        rect(0,0,1,1,col="grey80",border="grey80")
        text(0.5,0.5,PClabel[i],cex=1)
}
par(srt=0)
# legend 
plot(1,xlim=c(0,1),ylim=c(0,1),type="n",axes=F,xlab="",ylab="",xaxt="n",yaxt="n")
legend('center',legend=group.levels,col=group_col,bty="n",pch=16,cex=2.5)

dev.off()
