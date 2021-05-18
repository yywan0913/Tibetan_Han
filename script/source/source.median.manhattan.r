args = commandArgs(TRUE)
library(data.table)
require("RColorBrewer")
#inregion = args[1]
#gene.label = args[2]

options(stringsAsFactors=F)
options(scipen=200)
#pos.data = read.table(inregion,sep="\t",header=F)


median.manhanttan = function(pos.data,gene.data,sv.data,snp.data,indel.data,color,pch){
	gene.label = sv.data[sv.data[,4]!="False",4]
	genelabel = unique(unlist(strsplit(gsub('\\(|\\)','',gene.label),";")))
	sv.data = sv.data
	snp.data = snp.data
	indel.data = indel.data
	colnames(sv.data)[1:3]=colnames(indel.data)[1:3]=colnames(snp.data)[1:3] = c('chrom','Pos','FST')
	regioni.chrom = as.character(pos.data[1])
	regioni.start = as.numeric(pos.data[2])
	regioni.end = as.numeric(pos.data[3])
	yellowi.start = as.numeric(pos.data[5])
	yellowi.end = as.numeric(pos.data[6])
	sv.datai = sv.data[sv.data$chrom==regioni.chrom & sv.data$Pos>=regioni.start &
		sv.data$Pos<=regioni.end&!is.na(sv.data$FST),,drop=F]
	indel.datai = indel.data[indel.data$chrom==regioni.chrom & indel.data$Pos>=regioni.start &
		indel.data$Pos<=regioni.end&!is.na(indel.data$FST),,drop=F]
	snp.datai = snp.data[snp.data$chrom==regioni.chrom & snp.data$Pos>=regioni.start &
		snp.data$Pos <= regioni.end&!is.na(snp.data$FST),,drop=F]
	gene.datai = gene.data[gene.data[,1]==regioni.chrom & gene.data[,2]>=regioni.start &
		gene.data[,3]<=regioni.end,,drop=F]
	gene.datai.gene = gene.datai[gene.datai[,5]=="gene",,drop=F]
	gene.datai.exon = gene.datai[gene.datai[,5]=="exon",,drop=F]
	gene.numi = nrow(gene.datai.gene)
	gene.col = colorRampPalette(brewer.pal(9, "Paired"))(gene.numi)
	main_title = gettextf('%s:%s-%s',regioni.chrom,prettyNum(regioni.start,","),prettyNum(regioni.end,","))
	mcolor = color
	mpch = pch
	y1 = sv.datai$FST; pos1 = sv.datai$Pos
	y3 = indel.datai$FST; pos3 = indel.datai$Pos
	y3 = snp.datai$FST; pos3 = snp.datai$Pos
	ym = max(max(y1),max(y2),max(y3))
	if(ym<0.5) ym = 0.7
	par(las=1,cex.lab=2.5,cex.main=2.5,cex.axis=1.2,mgp=c(2.5,0.6,0),lab=c(5,6,0))
	plot(1,ylim=(c(-0.1,ym+0.1)),xlim=(c(regioni.start,regioni.end)),cex=1,xlab="",ylab="Fst",main=main_title,axes=F,type="n")
	rect(yellowi.start,0,yellowi.end,1,col="#FFFF00B2",border=NA)
	#points(pos3,y3,col=indel.datai$color,pch=mpch[3]) # snp
	#points(pos2,y2,col=snp.datai$color,pch=mpch[2]) # indel
	#points(pos1,y1,col=sv.datai$color,pch=mpch[1],cex=1.5) # sv
	points(pos2,y2,col=color[2],pch=mpch[2],cex=1.5) # snp
	points(pos3,y3,col=color[3],pch=mpch[3],cex=1.5) # indel
	points(pos1,y1,col=color[1],pch=mpch[1],cex=1.5)
	
	if(gene.numi>0)
	{
		for(j in 1:nrow(gene.datai.gene))
		{
			genej = gene.datai.gene[j,4]
			segments(gene.datai.gene[j,2],-0.075,gene.datai.gene[j,3],-0.075,lwd=1,col="grey")
			gene.datai.exon.j = gene.datai.exon[gene.datai.exon[,4]==genej,,drop=F]
			rect(gene.datai.exon.j[,2],-0.1,gene.datai.exon.j[,3],-0.05,col=gene.col[j],border=gene.col[j])
			#if(sum(!is.na(match(genej,genelabel)))>0|j==1|j==nrow(gene.datai.gene))
			if((j==1|j==nrow(gene.datai.gene))&genej!="TTTY14")
			{
				text(sum(gene.datai.gene[j,2:3])/2,-0.1,genej,cex=1,font=4,pos=1,col=gene.col[j])
			}
			if(genej=="CD24"|genej=="EGLN1"|genej=="EPAS1"){
				text(sum(gene.datai.gene[j,2:3])/2,-0.1,genej,cex=1,font=4,pos=1,col=gene.col[j])
			}
		}
	}
	sv.datai = sv.datai[sv.datai[,4]!="False",,drop=F]
	if(nrow(sv.datai)>0){
		text(sv.datai[,2],sv.datai[,3],sv.datai[,4],pos=4,cex=4,font=4)
	}
	axis(1,at=seq(from=regioni.start,to=regioni.end,by=100000))
	axis(2)
	legend("topleft",c("SV","SNV","InDel"),col=mcolor,text.col=mcolor,pch=mpch,cex=2.5)
	
}





