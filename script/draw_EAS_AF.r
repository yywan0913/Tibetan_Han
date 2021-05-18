args = commandArgs(TRUE)
input = args[1]
output = args[2]
df = read.table(input,sep="\t",header=T)
EAS = df[,5]
AT = df[,2]
AL = df[,3]
HT = df[,4]
c1 = cor.test(EAS,AT);R1 = (c1$estimate)^2;p1=c1$p.value
c2 = cor.test(EAS,AL);R2 = (c2$estimate)^2;p2=c2$p.value
c3 = cor.test(EAS,HT);R3 = (c3$estimate)^2;p3=c3$p.value
p1 = ifelse(p1==0,"<2.2e-16",paste("=",p1))
p2 = ifelse(p2==0,"<2.2e-16",paste("=",p2))
p3 = ifelse(p3==0,"<2.2e-16",paste("=",p3))

#col.raw = c("#00AFBB", "#E7B800", "#FC4E07") 
#col.raw = c(rgb(115/255,82/255,130/255),rgb(253/255,253/255,0),rgb(118/255,181/255,87/255))
##col.raw = c(rgb(115/255,82/255,130/255),'darkred',rgb(118/255,181/255,87/255)) # #DAA520
col.raw = c("#00AFBB", "purple", "#FC4E07")
col = apply(col2rgb(col.raw)/255,2,function(x)rgb(x[1],x[2],x[3],0.8))
pdf(output,width=13,height=15)
mat = matrix(c(2,1,4,0,3,0),nc=2)
layout(mat=mat,width=c(10,2.8),height=c(2,10,3))

# main
par(mar=c(0.1,8,0.5,0.1),cex.lab=2.5,cex.axis=2.2,mgp=c(3.5,1,0),las=1,tck=-0.01,cex.main=2.5,font.lab=2,mgp=c(4,1,0))
plot(EAS,AT,type="n",ylab="Tibetan and/or Han AF",xlab="",ylim=c(0,1),xlim=c(0,1),xaxt="n",main="") # EAS AF
segments(0,0,1,1,lwd=1.2)
## lines(EAS,as.numeric(fitted(z)),col=col[1],lwd=3)
## AL
points(EAS,AL,cex=1.2,col=col[2],pch=18)
#z=lm(AL~EAS)
## lines(EAS,as.numeric(fitted(z)),col=col[2])

## HT
points(EAS,HT,cex=1.2,col=col[3],pch=20)
#z=lm(HT~EAS)
## lines(EAS,as.numeric(fitted(z)),col=col[3])

# ALL
points(EAS,AT,cex=1.2,col=col[1],pch=16)
#z=lm(AT~EAS)

abline(h=1.04)
axis(1)

#  up EAS
par(mar=c(0.1,8,3,0.1),mgp=c(3,0.5,0))
hist(EAS,col="grey60",border="grey60",xaxt="n",xlab="",ylab="",breaks=seq(0,1,by=0.025),main="",axes=F)
legend('center',legend="EAS",bty="n",cex=2.5)
a2=axis(2,tick=F,labels=F)
axis(2,a2,prettyNum(a2,","))

#axis(1)
## right AL-HT
par(mar=c(0.1,0.1,2,1),mgp=c(3,0.5,0))
h1 = hist(AT,breaks=seq(0,1,by=0.05),plot=F)
h2 = hist(AL,breaks=seq(0,1,by=0.05),plot=F)
h3 = hist(HT,breaks=seq(0,1,by=0.05),plot=F)
wid = 0.014

par(las=3)
plot(1,xlab="AL-HT",ylab="",yaxt="n",xaxt="n",col=col[1],type="n",axes=F,xlim=c(0,max(h1$counts,h2$counts,h3$counts)),ylim=c(0,1))

for(i in 1:(length(h1$breaks)-1)){
	rect(0,h1$breaks[i],h1$counts[i],h1$breaks[i]+wid,col=col.raw[1],border=NA)
	rect(0,h1$breaks[i]+wid,h2$counts[i],h1$breaks[i]+wid*2,col=col.raw[2],border=NA)
	rect(0,h1$breaks[i]+wid*2,h3$counts[i],h1$breaks[i]+wid*3,col=col.raw[3],border=NA)
}
legend('right',legend=c('Tibetan&Han','Tibetan','Han'),col=col.raw[1:3],pch=15,bty="n",cex=2.5)
a3=axis(3,tick=F,labels=F)
axis(3,a3,prettyNum(a3,","))
## bottom
par(mar=c(0.1,8,1,0.1),las=1,mgp=c(3,0.5,0))
plot(1,type="n",ylab="",xlab="",ylim=c(0,1),xlim=c(0,1),xaxt="n",axes=F)
legends <- vector('expression',3)
legends[1] <- substitute(expression(paste("EAS--Tibetan&Han: ",R^2 ," = ",rR1,", Pvalue ",p1)),list(rR1=round(R1,4),p1=p1))[2]
legends[2] <- substitute(expression(paste("EAS--Tibetan: ",R^2," = ", rR2,", Pvalue ",p2)),list(rR2=round(R2,4),p2=p2))[2]
legends[3] <- substitute(expression(paste("EAS--Han: ",R^2 ," = ", rR3,", Pvalue ",p3)),list(rR3=round(R3,4),p3=p3))[2]
legend('left',legend=legends,bty="n",col=col.raw,pch=c(16,18,20),cex=2.2)
dev.off()

