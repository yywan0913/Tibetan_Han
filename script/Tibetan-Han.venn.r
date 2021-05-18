library(VennDiagram)
library(data.table)
args = commandArgs(TRUE)
mergefile = args[1]
outpng = args[2]

df = fread(mergefile,sep="\t",header=F)
df = as.data.frame(df)
df = df[df[,7]>=50&df[,7]<=50000,,drop=F]
tab = table(df[,1])
oneID = names(tab[tab==1])

term = unique(df[,8])
for (i in term){
	if (grepl('^Tibetan',i)){
		AL = unique(df[df[,8]==i,1])
	}else if (grepl('^Han',i)){
		HT = unique(df[df[,8]==i,1])
	}
}
x = list(Tibetan=AL,Han=HT)
venn.diagram(x=x, outpng, height = 450, width = 450, resolution =300, imagetype="png", col="transparent", fill=c('cornflowerblue','green'), alpha=0.6, cex=0.45, cat.cex=0.45,label.col="black",fontfamily = "serif", fontface = "bold")
