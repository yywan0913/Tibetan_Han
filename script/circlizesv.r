#!/export/software/bin/Rscript
library('getopt',lib.loc="/home/xiaoyuhui/R/x86_64-pc-linux-gnu-library/3.4")
library(data.table)
spec = matrix(c(
        'help' ,    'h', 0, "logical",
        'groupvcf',  'v', 1, "character",
	#'samplelist', 'l',1,"character",
        #'prefix',    'p', 1, "character",
        'outdir',   'o', 1, "character",
	'windows',  'w', 0,"numeric",
	'splide',   'sp',0,"numeric"
), byrow=TRUE, ncol=4);
opt = getopt(spec);
# define usage function
print_usage <- function(spec=NULL){
        cat(getopt(spec, usage=TRUE));
        cat("
Usage example:
1) Rscript DiffFreqBox.r --groupvcf testmerge.vcf -o /path/outdir/

Options:
--help          NULL            get this help;
--groupvcf      character       the group file [*vcf] [forced];
#--samplelist  character         the group [forced];
#--prefix       character       the group names('A--B') [forced];
--outdir      character       outdir [forced];
--windows     numeric         windows [default:10Mb];
--splide      numeric         splided windows [default:1Mb];
\n")
        q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if(!dir.exists(opt$outdir)) dir.create(opt$outdir)
if(is.null(opt$windows)) opt$windows = 10000000
if(is.null(opt$splide)) opt$splide = 1000000
nsplide = opt$windows/opt$splide
scriptdir = dirname(get_Rscript_filename())
#centromere = paste0(scriptdir,"/hg19.centromere.txt")
#if(!file.exists(centromere)) centromere = "/export/home/xiaoyh/script/windowsdiff/hg19.centromere.txt"

get_centromere = function(centromere,splide){
	df = read.table(centromere,sep="\t",header=F)
	Centromere = list()
	for(i in 1:nrow(df)){
		Centromere[[df[i,1]]] = c(df[i,2]/splide,df[i,3]/splide)
	}
	return(Centromere)
}

readvcf <- function(file){
	df <- fread(file,sep="\t")
	df <- as.data.frame(df)
	s_ = sum(grepl('##',df[,1]))
	if(s_ >=0){
		df <- fread(file,header=T,sep="\t",skip=s_)
		df <- as.data.frame(df)
	}
	svtype = gsub(';.*','',gsub('.*SVTYPE=','',df[,8]))
	svtype[svtype=="BND"]="TRA"
	df <- df[,-c(3:9)]
	colnames(df)[1]="CHROM"
	#colnames(df)=gsub('(.*).ngmlr.*','\\1',colnames(df))
	df[,-(1:2)] = apply(df[,-(1:2)],1:2, gsub ,pattern=':.*',replacement='')
	dfdata <- df[,-(1:2)]
	dfdata <- ifelse(dfdata=="1/1"|dfdata=="0/1"|dfdata=="1",1,0)
	data = data.frame(df[,1:2],svtype=svtype,dfdata,check.names=F)
	return(data)
}

readTRA <- function(file){
	df <- fread(file,sep="\t")
	df <- as.data.frame(df)
	s_ = sum(grepl('##',df[,1]))
	if(s_ >=0){
		df <- fread(file,header=T,sep="\t",skip=s_)
		df <- as.data.frame(df)
	}
	svtype = gsub(';.*','',gsub('.*SVTYPE=','',df[,8]))
	if(is.na(match('TRA',svtype))) return(NULL) else{
		df = df[svtype=="TRA",,drop=F]
		chr2 = gsub(';.*','',gsub('.*CHR2=','',df[,8]))
		chr2pos = gsub(';.*','',gsub('.*;END=',"",df[,8]))
		if(!grepl('chr',df[1,1])) {
			chr1 = paste0('chr',df[,1]) 
			chr2 = paste0('chr',chr2)
		}else {
			chr1 = df[,1]
		}
		data = data.frame(chr=chr1,start=as.numeric(df[,2]),chr2=chr2,chr2pos=as.numeric(chr2pos))
		return(data)
		
	}
}

getboxdata = function(data,nsplide,svtype){
	Box = list()
	chr = unique(data[,1])
	nsample = ncol(data)-3
	#MaxNum = 0
	for(i in chr){   ## chr
		datai = data[data[,1]==i,]
		max = max(datai[,2])-nsplide
		svtypenum = list()
		for(k in svtype) svtypenum[[k]]= 0 
		start = c()
		end = c()
		for(j in 1:max){
			w=which(datai[,2]%in%j:(j+nsplide))
			start = c(start,j*opt$splide)
			end = c(end,(j+nsplide)*opt$splide)
			#if(sum(w)==0) num=c(num,rep(NA,nsample)) else num=c(num,apply(datai[w,-c(1:2)],2,sum))
			if(sum(w)==0) for(k in svtype) svtypenum[[k]]=c(svtypenum[[k]],0) 
			if(sum(w)>0){
				dataij = datai[w,]
				for(k in svtype){
					dataijt = dataij[dataij[,3]==k,-c(1:3),drop=F]
					if(nrow(dataijt)==0) svtypenum[[k]]=c(svtypenum[[k]],0) else
						svtypenum[[k]] = c(svtypenum[[k]],nrow(dataijt))
				}
			}
		}
		for(k in svtype){  ## type
			if(!grepl('chr',i)) chri = paste0('chr',i) else chri=i
			Box[[k]][[i]] = data.frame(chr=chri,start=start,end=end,value1=svtypenum[[k]][-1])
			#Box[[k]][[i]]$value1 = ifelse(Box[[k]][[i]]$value1<=1,0,log(Box[[k]][[i]]$value1,2))
			#MaxNum = max(MaxNum,Box[[k]][[i]]$value1)
		}
	}
	return(Box)
}


getcirclize = function(box,TRAbed,outtiff){
	#col = rainbow(7,s=0.6,v=0.6)
	col = c(rgb(182/255,53/255,188/255),rgb(52/255,70/255,137/255),rgb(189/255,36/255,188/255),rgb(44/255,136/255,62/255),rgb(56/255,35/255,98/255),rgb(204/255,68/255,112/255),rgb(208/255,202/255,106/255))
	svtype = names(box)
	library(circlize)
	tiff(outtiff,height=32,width=32,units="cm",res=600,compression="lzw")
	par(mar = c(4, 4, 4, 4))
	circos.par("track.height" = 0.3)
	circos.par('gap.degree'=4)   ## 
	circos.par(gap.after = 2, start.degree = 90)
	#circos.par('start.degree'=8)  ## jian ge qishi  jiaodu
	#circos.initializeWithIdeogram(plotType = NULL)
	circos.initializeWithIdeogram(chromosome.index = paste0('chr',c(1:22,"X","Y")), 
        	plotType = c("ideogram", "labels"), ideogram.height = 0.03)
	draw_small_circle = function(box,TYPE,svtype,num){
		if(!is.na(match(TYPE,svtype))){
			bed = do.call(rbind,box[[TYPE]])
			circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
				circos.genomicLines(region, value, col=col[num], type = "s", area = TRUE,...)  ## lwd=1 type="l"
			}, track.height=0.1,bg.border=NA,bg.col="grey90")
		}
	}
	draw_small_circle(box,'DEL',svtype,num=1)
	draw_small_circle(box,'INS',svtype,num=2)
	draw_small_circle(box,'DUP',svtype,num=3)
	draw_small_circle(box,'INV',svtype,num=4)
	draw_small_circle(box,'TRA',svtype,num=5)
	if(!is.null(TRAbed)){
		#TRAbed = TRAbed[!duplicated(TRAbed[,1:2]),]
		#TRAbed = TRAbed[!duplicated(TRAbed[,3:4]),]
		TRAbed = TRAbed[TRAbed[,1]!="chrMT"&TRAbed[,3]!="chrMT",]
		chr1chr2 = paste0(TRAbed[,1],"--",TRAbed[,1])
		unichr1chr2 = unique(chr1chr2)
		chr1chr2table = table(chr1chr2)
		TRAcol = rainbow(length(unichr1chr2),s=0.8,v=0.8)
		#TRAcol = rgb(228/255,53/255,39/255)
		chr1chr2alpha = chr1chr2table/max(chr1chr2table)   ## uniq
		tablename = names(chr1chr2table)
		tmpcol = col2rgb(TRAcol)/255
		TRAcolor = c()
		for(j in 1:length(chr1chr2)){
			mj = match(chr1chr2[j],unichr1chr2)
			alpha = chr1chr2alpha[tablename==chr1chr2[j]]
			TRAcolor = c(TRAcolor,rgb(tmpcol[1,mj],tmpcol[2,mj],tmpcol[3,mj],alpha=alpha))
		}
		#TRAcolor = TRAcol[as.numeric(as.factor(chr1chr2))]
		bed1 = data.frame(chr=TRAbed[,1],start=as.numeric(TRAbed[,2]),end=as.numeric(TRAbed[,2])+5000000,value1=1)
		bed2 = data.frame(chr=TRAbed[,3],start=as.numeric(TRAbed[,4]),end=as.numeric(TRAbed[,4])+5000000,value1=1)
		#for(x in 1:nrow(TRAbed)) circos.link(TRAbed[x,1],c(TRAbed[x,2],TRAbed[x,2]+10000),TRAbed[x,3],c(TRAbed[x,4],TRAbed[x,4]+10000),h=0.3,border=NA,col=col[6])
		circos.genomicLink(bed1, bed2, col = TRAcolor, border = NA,lwd=2)
	}
	circos.clear()
	dev.off()
}

options(stringsAsFactors=F)

groupfile = opt$groupvcf
data = readvcf(groupfile)
data[,2] = ceiling(data[,2]/opt$splide)  ## 
svtype = unique(data[,3])
Box = getboxdata(data,nsplide,svtype=svtype)
TRAbed = readTRA(opt$groupvcf)
MaxNum = 1
for(i in svtype){
	svtypeibed = do.call(rbind,Box[[i]])
	MaxNum = max(MaxNum,svtypeibed$value1)
	write.table(svtypeibed,paste0(opt$outdir,"/",i,".splide.counts.xls"),sep="\t",col.names=T,row.names=F,quote=F)
}
write.table(TRAbed,paste0(opt$outdir,"/TRA.bed.xls"),sep="\t",col.names=T,row.names=F,quote=F)

for(i in unique(data[,1])){
	for (k in svtype){
		#Box[[k]][[i]]$value1 = Box[[k]][[i]]$value1/MaxNum
		Box[[k]][[i]]$value1 = ifelse(Box[[k]][[i]]$value1<=1,0,log(Box[[k]][[i]]$value1,2))
		Box[[k]][[i]]$value1 = Box[[k]][[i]]$value1/log(MaxNum,2)
	}
}

outtiff = paste0(opt$outdir,"/splide.circlize.tiff")
getcirclize(Box,TRAbed=TRAbed,outtiff=outtiff)

