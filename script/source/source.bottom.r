
PlotRegion = function(data.info,data.AL,data.HT,gene.data,gene.col){
	HanColor = "#009CFFFF"
	TibetanColor = "#FFC600FF"
	SNPColor = "#FC4E07"
	INSColor = "#9E819BFF"
	chrom = data.info[1,1]
	start = min(data.info[,2])-10
	end = max(data.info[,2])+10
	gene.data = gene.data[gene.data[,1]==chrom,]
	gene.data = gene.data[!(gene.data[,2]>end|gene.data[,3]<start),,drop=F]
	gene.data.gene = gene.data[gene.data[,5]=="gene",,drop=F]
	gene.data.gene = gene.data.gene[order(gene.data.gene[,2]),]
	gene.data.exon = gene.data[gene.data[,5]=="exon",,drop=F]
	#layout(mat=matrix(c(1,2),nc=1),height=c(2,8))
	par(mar=c(0.1,2,4,2))
	par(xpd=T)
	plot(1,xlim=c(start-100,end+100),ylim=c(-1,10),axes=F,type="n",xlab="",ylab="")
	if(nrow(gene.data.gene)>0){
	for(i in 1:nrow(gene.data.gene)){
		genei = gene.data.gene[i,4]
		segments(gene.data.gene[i,2],5,gene.data.gene[i,3],5,lwd=1,col="grey")
		gene.data.exon.i = gene.data.exon[gene.data.exon[,4]==genei,,drop=F]
		rect(gene.data.exon.i[,2],3,gene.data.exon.i[,3],7,col=gene.col[i],border=gene.col[i])
		if(i%%2==1){
			text(sum(max(gene.data.gene[i,2],start),min(gene.data.gene[i,3],end))/2,7,genei,cex=1.2,font=4,pos=3)
		}else{
			text(sum(max(gene.data.gene[i,2],start),min(gene.data.gene[i,3],end))/2,3,genei,cex=1.2,font=4,pos=1)
		}
	}
	}
	# 2
	par(mar=c(4,2,0.1,2),xpd=F)
	if(chrom=="2") start = start-10000 else start=start
	plot(1,xlim=c(start,end),ylim=c(0,ncol(data.AL)+ncol(data.HT)),axes=F,type="n",xlab="",ylab="",yaxs="i")
	d = 0.05
	ins.text.cex = 0.7
	for(pos in 1:ncol(data.HT)){
		i = pos
		rect(start,i-1+d,end,i-d,col=HanColor,border=NA)
		for(j in 1:nrow(data.info)){
			samplevariable = data.HT[j,pos]
			if (samplevariable == 0) next
			alt = data.info[j,3]
			if(grepl('DEL',alt)){
				del.start = as.numeric(strsplit(alt,":|-")[[1]][2])
				del.end = as.numeric(strsplit(alt,"-|_")[[1]][3])
				rect(del.start,i-1+d,del.end,i-d,col="white",border=NA)
				segments(del.start,i-0.5,del.end,i-0.5,col="black")
				segments((del.start+del.end)/2,i-0.4,(del.start+del.end)/2,i-0.6,lwd=2)
			}else if(!grepl('_',alt)){
				snp.start = data.info[j,2]
				rect(snp.start,i-1+d,snp.start+1,i-d,col=SNPColor,border=SNPColor)
			}
		}
		for(j in 1:nrow(data.info)){
			samplevariable = data.HT[j,pos]
			if (samplevariable == 0) next
			alt = data.info[j,3]
			if(grepl('INS',alt)){
				ins.start = as.numeric(strsplit(alt,":|-")[[1]][2])
				ins.len = as.numeric(gsub('.*_','',alt))
				ins.strwidth = strwidth(ins.len,cex=ins.text.cex,family="mono")
				rect(ins.start,i-1+d,ins.start+ins.strwidth,i-d,col=INSColor,border=NA)
				text(ins.start,i-0.5,ins.len,adj=0,cex=ins.text.cex,family="mono",col="white")
			}
		}
	}
	for(pos in 1:ncol(data.HT)){
		i = pos +ncol(data.HT)
		rect(start,i-1+d,end,i-d,col=TibetanColor,border=NA)
		for(j in 1:nrow(data.info)){
			samplevariable = data.AL[j,pos]
			if (samplevariable == 0) next
			alt = data.info[j,3]
			if(grepl('DEL',alt)){
				del.start = as.numeric(strsplit(alt,":|-")[[1]][2])
				del.end = as.numeric(strsplit(alt,"-|_")[[1]][3])
				rect(del.start,i-1+d,del.end,i-d,col="white",border=NA)
				segments(del.start,i-0.5,del.end,i-0.5,col="black")
				segments((del.start+del.end)/2,i-0.4,(del.start+del.end)/2,i-0.6,lwd=2)
			}else if(!grepl('_',alt)){
				snp.start = data.info[j,2]
				rect(snp.start,i-1+d,snp.start+1,i-d,col=SNPColor,border=SNPColor)
			}
		for(j in 1:nrow(data.info)){
			samplevariable = data.AL[j,pos]
			if (samplevariable == 0) next
			alt = data.info[j,3]
			if(grepl('INS',alt)){
				ins.start = as.numeric(strsplit(alt,":|-")[[1]][2])
				ins.len = as.numeric(gsub('.*_','',alt))
				ins.strwidth = strwidth(ins.len,cex=ins.text.cex,family="mono")
				rect(ins.start,i-1+d,ins.start+ins.strwidth,i-d,col=INSColor,border=NA)
				text(ins.start,i-0.5,ins.len,adj=0,cex=ins.text.cex,family="mono",col="white")
			}
		}
		}
	}
}
	
