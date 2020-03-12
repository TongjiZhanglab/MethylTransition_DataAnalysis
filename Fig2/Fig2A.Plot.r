###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)

ErrorBar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
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
pdf("Fig2A.Parameters_DropoutDisturbing.pdf",width=3.5,height=5)
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