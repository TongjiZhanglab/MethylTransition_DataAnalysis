###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol01 <- c("#CE001301","#16557A01","#C7A60901","#87C23201","#00879201","#A14C9401","#15A08C01","#8B7E7501","#1E7CAF01","#EA425F01","#46489A01","#E5003301","#0F231F01","#1187CD01")
cccol05 <- c("#CE001305","#16557A05","#C7A60905","#87C23205","#00879205","#A14C9405","#15A08C05","#8B7E7505","#1E7CAF05","#EA425F05","#46489A05","#E5003305","#0F231F05","#1187CD05")
cccol30 <- c("#CE001330","#16557A30","#C7A60930","#87C23230","#00879230","#A14C9430","#15A08C30","#8B7E7530","#1E7CAF30","#EA425F30","#46489A30","#E5003330","#0F231F30","#1187CD30")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E8003380","#0F231F80","#1187CD80")
options(scipen = 200)
library(beeswarm)

ErrorBar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
HalfErrorBar <- function(x, y, upper, lower=upper, length=0.03,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	arrows(x,y+upper, x, y, angle=90, code=1, length=length, ...)
}
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
u <- read.table("Parameters_u_DropoutDisturbing.txt")
d <- read.table("Parameters_d_DropoutDisturbing.txt")
p <- read.table("Parameters_p_DropoutDisturbing.txt")

###################################################################################################
################################                PLOT               ################################
###################################################################################################
pdf("Fig2B.Parameters_DropoutDisturbing.pdf",width=3.5,height=5)
plot_list <- list(u[,9],u[,8],u[,7],u[,6],u[,5],u[,4],u[,3],u[,2],u[,1])
boxplot(plot_list,boxwex=0.5, col="white",border=cccol[1],outline=F,ylim=c(0,1),names=seq(0.1,0.9,0.1),main="",xlab="Dropout ratio",ylab="Parameter value",las=2);box(lwd=2)
plot_list <- list(d[,9],d[,8],d[,7],d[,6],d[,5],d[,4],d[,3],d[,2],d[,1])
boxplot(plot_list,boxwex=0.5, col="white",border=cccol[2],outline=F,add=T,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot_list <- list(p[,9],p[,8],p[,7],p[,6],p[,5],p[,4],p[,3],p[,2],p[,1])
boxplot(plot_list,boxwex=0.5, col="white",border=cccol[3],outline=F,add=T,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
legend("topleft",c("u","d","p"),col=cccol,pch=15,bty="n")
dev.off()

# pdf("LineParameters_StabilityDisturbing.pdf",width=3.8,height=5)
# plot(1,type="n",xaxt="n",yaxt="n",xlab="Drop-out ratio",ylab="Parameter value",ylim=c(0,1),xlim=c(1,9));box(lwd=2)
# axis(side=1,at=seq(9),labels=seq(0.1,0.9,0.1))
# axis(side=2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
# u <- rev(u)
# n <- apply(u,2,length); mean <- apply(u,2,mean,na.rm=T); sd <- apply(u,2,sd,na.rm=T)
# points(seq(9),mean,lwd=2,col=cccol[1],type="b",pch=".")
# ErrorBar(seq(9),mean, qt(0.99, n-1) * sd / sqrt(n),lwd=2,length=0.03,col=cccol[1])
# d <- rev(d)
# n <- apply(d,2,length); mean <- apply(d,2,mean,na.rm=T); sd <- apply(d,2,sd,na.rm=T)
# points(seq(9),mean,lwd=2,col=cccol[2],type="b",pch=".")
# ErrorBar(seq(9),mean, qt(0.99, n-1) * sd / sqrt(n),lwd=2,length=0.03,col=cccol[2])
# p <- rev(p)
# n <- apply(p,2,length); mean <- apply(p,2,mean,na.rm=T); sd <- apply(p,2,sd,na.rm=T)
# points(seq(9),mean,lwd=2,col=cccol[3],type="b",pch=".")
# ErrorBar(seq(9),mean, qt(0.99, n-1) * sd / sqrt(n),lwd=2,length=0.03,col=cccol[3])
# legend("topleft",c("u","d","p"),col=cccol,pch=15,bty="n")
# dev.off()

pdf("FigR8.ParametersRobustness.Dropout.pdf",width=5,height=4)
par(mar=c(6,4,4,4),lwd=2)
l1 <- rbind(u[,9],d[,9],p[,9])
l2 <- rbind(u[,8],d[,8],p[,8])
l3 <- rbind(u[,7],d[,7],p[,7])
l4 <- rbind(u[,6],d[,6],p[,6])
l5 <- rbind(u[,5],d[,5],p[,5])
l6 <- rbind(u[,4],d[,4],p[,4])
l7 <- rbind(u[,3],d[,3],p[,3])
l8 <- rbind(u[,2],d[,2],p[,2])
l9 <- rbind(u[,1],d[,1],p[,1])
plot_list <- list(as.numeric(cor(l1)),as.numeric(cor(l2)),as.numeric(cor(l3)),as.numeric(cor(l4)),as.numeric(cor(l5)),as.numeric(cor(l6)),as.numeric(cor(l7)),as.numeric(cor(l8)),as.numeric(cor(l9)))
boxplot(plot_list,xaxt="n",col = cccol50[5],boxwex=0.6,border=cccol[5],outline=F,main="Dropout",xlab="",ylab="Parameter correlation",ylim=c(-1,1),las=2);box(lwd=2)
axis(at=seq(length(plot_list)),side=1,labels=seq(0.1,0.9,0.1),las=2)

boxplot(plot_list,xaxt="n",col = cccol50[5],boxwex=0.6,border=cccol[5],outline=F,main="Dropout",xlab="",ylab="Parameter correlation",las=2);box(lwd=2)
axis(at=seq(length(plot_list)),side=1,labels=seq(0.1,0.9,0.1),las=2)

par(mar=c(6,4,4,4),lwd=1)
plot_list <- list(as.numeric(dist(t(l1))),as.numeric(dist(t(l2))),as.numeric(dist(t(l3))),as.numeric(dist(t(l4))),as.numeric(dist(t(l5))),as.numeric(dist(t(l6))),as.numeric(dist(t(l7))),as.numeric(dist(t(l8))),as.numeric(dist(t(l9))))
boxplot(plot_list,xaxt="n",yaxt="n",col = cccol50[5],boxwex=0.6,border=cccol[5],outline=F,main="Dropout",xlab="",ylab="Parameter Distance",ylim=c(0,sqrt(3)),las=2);box(lwd=2)
axis(at=seq(length(plot_list)),side=1,labels=seq(0.1,0.9,0.1),las=2)
axis(at=c(0,sqrt(3)/2,sqrt(3)),side=2,labels=c("0","s3/2","s3"),las=2)
abline(h=seq(0,sqrt(3),sqrt(3)/10),lty=2,col=cccol50[8])

par(mar=c(2,6,4,4),lwd=2)
n <- rbind(length(as.numeric(cor(l1))),length(as.numeric(cor(l2))),length(as.numeric(cor(l3))),length(as.numeric(cor(l4))),length(as.numeric(cor(l5))),length(as.numeric(cor(l6))),length(as.numeric(cor(l7))),length(as.numeric(cor(l8))),length(as.numeric(cor(l9))))
mean <- rbind(mean(as.numeric(cor(l1)),na.rm=T),mean(as.numeric(cor(l2)),na.rm=T),mean(as.numeric(cor(l3)),na.rm=T),mean(as.numeric(cor(l4)),na.rm=T),mean(as.numeric(cor(l5)),na.rm=T),mean(as.numeric(cor(l6)),na.rm=T),mean(as.numeric(cor(l7)),na.rm=T),mean(as.numeric(cor(l8)),na.rm=T),mean(as.numeric(cor(l9)),na.rm=T))
sd <- rbind(sd(as.numeric(cor(l1)),na.rm=T),sd(as.numeric(cor(l2)),na.rm=T),sd(as.numeric(cor(l3)),na.rm=T),sd(as.numeric(cor(l4)),na.rm=T),sd(as.numeric(cor(l5)),na.rm=T),sd(as.numeric(cor(l6)),na.rm=T),sd(as.numeric(cor(l7)),na.rm=T),sd(as.numeric(cor(l8)),na.rm=T),sd(as.numeric(cor(l9)),na.rm=T))
bp <- barplot(mean,width = 1,lwd=2,main="",space=c(0.5,rep(0.3,8)),xlim=c(0,2.4+9+1),ylim=c(0.98,1.003),beside=T,col=cccol50[5],border=cccol[5],xpd=F,las=2,ylab="Pearson correlation coefficient");box(lwd=2)
axis(at=bp,side=1,labels=seq(0.1,0.9,0.1),las=2)
# HalfErrorBar(bp,mean, qt(0.975, n-1) * sd / sqrt(n),col=cccol[7],lwd=2)
HalfErrorBar(bp,mean, 1.96 * sd / 10,col=cccol[5],lwd=2)
# plot_list <- list(as.numeric(cor(l1)),as.numeric(cor(l2)),as.numeric(cor(l3)),as.numeric(cor(l4)),as.numeric(cor(l5)),as.numeric(cor(l6)),as.numeric(cor(l7)),as.numeric(cor(l8)),as.numeric(cor(l9)))
# beeswarm(plot_list,xaxt="n",spacing = 1,col = cccol50[5],pch = 16,method = "swarm",corral = "wrap",add=TRUE)

bp <- barplot(mean,width = 1,lwd=2,main="",space=c(0.5,rep(0.3,8)),xlim=c(0,2.4+9+1),ylim=c(0.7,1.01),beside=T,col=cccol50[5],border=cccol[5],xpd=F,las=2,ylab="Pearson correlation coefficient");box(lwd=2)
axis(at=bp,side=1,labels=seq(0.1,0.9,0.1),las=2)
HalfErrorBar(bp,mean, 1.96 * sd / 10,col=cccol[5],lwd=2)
dev.off()