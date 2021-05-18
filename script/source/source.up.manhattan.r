library(data.table)
library(reshape2)
#output = "Manhattan.VAR.png"

up.manhattan = function(chrom.data,sv.data,snp.data,indel.data,col,pos.data){
	# chrom.data 
	chrom.len = chrom.data[,2]
	chrom.label = chrom.data[,1]
	x.cumsum = cumsum(as.numeric(chrom.len))
	x.axis = x.cumsum-chrom.len/2
	x.label = chrom.label
	data1 = sv.data
	data2 = snp.data
	data3 = indel.data
	# for line 1M one points
	data1.max = max1Mfst(data1)
	data2.max = max1Mfst(data2)
	data3.max = max1Mfst(data3)
	pos.data = pos.data[,1:3]
	for(k in 1:nrow(pos.data)){
		match.chr = match(pos.data[k,1],chrom.label)
		addnum = ifelse(match.chr==1,0,x.cumsum[match.chr-1])
		pos.data[k,2] = pos.data[k,2]+addnum
		pos.data[k,3] = pos.data[k,3]+addnum
	}
	
	# draw
	#col = c("#E64B35","#00A087","#3C5488")
	par(las=1,cex.lab=2.5,cex.main=2.5,cex.axis=1.2,mgp=c(2.5,0.6,0),lab=c(5,6,0))
	cex.p = 0.8
	plot(1,xlab="",ylab="MEAN_FST",main="SV",xaxt="n",ylim=c(-0.2,0.9),type="n",xlim=c(0,max(x.cumsum)))
	rect(pos.data[,2],0.05,pos.data[,3],1,col="#FFFF00B2",border="#FFFF00B2")
	segments(x.cumsum,-0.05,x.cumsum,0.05,col="grey60",lty=1)
	points(data2$x,data2[,3],cex=cex.p,col = data2$color,pch=data2$pch)
	points(data3$x,data3[,3],cex=cex.p,col = data3$color,pch=data3$pch)
	points(data1$x,data1[,3],cex=cex.p+0.1,col = data1$color,pch=data1$pch)
	abline(h=0.05,lty=1,col="grey")
	lines(data1.max$x,data1.max[,2]-0.2,col=col[1])
	lines(data2.max$x,data2.max[,2]-0.15,col=col[2])
	lines(data3.max$x,data3.max[,2]-0.1,col=col[3])
	abline(h=c(-0.1,-0.15,-0.2),lty=1,col="grey60")
	axis(1,x.axis,x.label)
	axis(2,0.05,0.05)
	data1 = data1[data1[,4]!="False",,drop=F]
	text(data1$x,data1[,3],data1[,4],pos=4,cex=2,font=4)
	legend('top',legend=c('SV','SNV','InDel'),col=col,pch=15:17,cex=2,bty="n",ncol=3,lty=1,lwd=2)
}

max1Mfst = function(data){
	data.deal = data[,c('WEIR_AND_COCKERHAM_FST','x')]
	data.deal$x = ceiling(data.deal$x/1000000)
	melt.data = melt(data.deal,id="x")
	melt.data.max = dcast(melt.data,x~variable,max)
	melt.data.max$x = (melt.data.max$x-0.5)*1000000
	melt.data.max[,2] = melt.data.max[,2]/max(melt.data.max[,2])*0.05
	return(melt.data.max)  ## x fst // 2col
}


