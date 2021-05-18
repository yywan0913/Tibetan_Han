orderK <- function (x){
	mydata<-x
	row.num<-nrow(mydata)
	column.num<-ncol(mydata)
	order.data<-data.frame()
	for (i in 1:row.num){
		max.index<-which(mydata[i,4:column.num]==max(mydata[i,4:column.num]))
		mydata[i,column.num+1]<-max.index+3
	}
	struc.num<-as.factor(mydata[,column.num+1])
	for (j in levels(struc.num)){
		temp<-mydata[mydata[,column.num+1]==as.integer(j),]
		temp<-temp[order(temp[,as.integer(j)]),]
		order.data<-rbind(order.data,temp)
	}
	ordered<-order.data[,-(column.num+1)]
	return(ordered)

}

######################################

library(getopt)
library(scales)
library(RColorBrewer)
h<-function(x){
	cat("Usage: Rscript --vanilla draw_structure_frappe.R -s k_start -e k_end  [-o Prefix_of_output_files (output)] [-h]\n\n")
	cat("Note: Before run this script,please link all K file in current dir,name format must k_2,k_3....\n\n")
	quit(save="no",status=x)
}
opt <- getopt(matrix(c(
    'help',   'h', 0, "logical",
    'prefix', 'o', 1, "character",
    'kstart', 's', 1, "integer",
    'kend',   'e', 1, "integer",
    'samplegroup','p',1,"character"
), ncol=4, byrow=TRUE));

if (! is.null(opt$help)) { h(0) }
file.prefix <- ifelse(is.null(opt$prefix), "output", opt$prefix)
range.start<-ifelse(is.null(opt$kstart),2,opt$kstart)

if (! is.null(opt$kend)) {
    range.end <- opt$kend
} else {
    h(1)
}
if (range.end<=12){
	mycol<-brewer.pal(12,"Paired")
}else{
	mycol<-c(brewer_pal(palette="Paired")(12),brewer_pal(palette="Dark2")(8));
}

palette(mycol)
file.name<-paste(file.prefix,".pdf",sep="")


pdf(file.name,h=8,w=10)
layout(c(1:10));
par(xpd=T,mfcol=c(4,1))
library(RColorBrewer)
palette(mycol)
##################k=2##############
#mylist<-read.table("list",header=F)
if(!is.null(opt$samplegroup)){
	group = read.table(opt$samplegroup,header=F,sep="\t")
	group_col = rainbow(length(unique(group[,2])),s=0.7,v=0.7)
}
for (k in range.start:range.end ){
	fn<-paste("k_",k,sep="")
	if (! file.exists(fn)) { cat(sprintf("sturcture file (%s) was not found!\n", fn)); h(1) }
	mydata<-read.table(fn,header=F)
	ordered.data<-orderK(mydata)
	mylist<-ordered.data[,1]
	ordered.data<-ordered.data[,c(-1,-2,-3)]
	rownames(ordered.data)<-mylist
	ordered.data<-t(ordered.data)
	par(mar=c(5,4,4,5)+0.1)
	bar = barplot(as.matrix(ordered.data[,1:length(mylist)]),main="",col=1:range.end,axisnames=FALSE,cex.axis=0.7,space =c(0,0),ylim=c(0,1),ylab="",border = NA)
	#axis(1,at=seq(0.5,length(mylist)-0.5,1),labels=colnames(ordered.data),col="white",las=2,cex.axis=0.7)
	if(!is.null(opt$samplegroup)) axiscol = group_col[as.numeric(as.factor(group[match(colnames(ordered.data),group[,1]),2]))] else axiscol = "black"
	par(xpd=T,srt=60);text(bar,par('usr')[3],colnames(ordered.data),col=axiscol,adj=1);par(srt=0)
	k.num<-paste("K=",k,sep="")
	legend("top",legend=1:k,pch=15,col=1:range.end,ncol=k,cex=0.8,bty="n",inset=c(0,-0.25))
	mtext(k.num,side=4)

}

dev.off()
q(save="no")
