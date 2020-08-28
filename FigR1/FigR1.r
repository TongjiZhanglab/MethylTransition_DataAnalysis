###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)
# library(MethylTransition)
# library(beeswarm)

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
Genebody <- read.table("hg_Genebody_parameters.txt",header=T)
Enhancer <- read.table("hg_Enhancer3_parameters.txt",header=T)
###################################################################################################
################################                PLOT               ################################
###################################################################################################

# Genebody
C12parameters <- Genebody[seq(2),]
C24parameters <- Genebody[seq(3,10),]
C48parameters <- Genebody[seq(11,42),]

plot_n <- c()
plot_mean <- c()
plot_sd <- c()

plot_data <- C12parameters
plot_n <- rbind(plot_n,apply(plot_data,2,length))
plot_mean <- rbind(plot_mean,apply(plot_data,2,mean))
plot_sd <- rbind(plot_sd,apply(plot_data,2,sd))
plot_data <- C24parameters
plot_n <- rbind(plot_n,apply(plot_data,2,length))
plot_mean <- rbind(plot_mean,apply(plot_data,2,mean))
plot_sd <- rbind(plot_sd,apply(plot_data,2,sd))
plot_data <- C48parameters
plot_n <- rbind(plot_n,apply(plot_data,2,length))
plot_mean <- rbind(plot_mean,apply(plot_data,2,mean))
plot_sd <- rbind(plot_sd,apply(plot_data,2,sd))

pdf("FigR1.MethyRatioTransitionPerStage.Genebody.pdf",width=4,height=4)
par(mar=c(4,4,4,4))
plot(1,type="n",xaxt="n",bty="n",xlab="",ylab="Parameter value",ylim=c(0,1),xlim=c(1,3),main="Human Genebody")
box(lwd=2)
v1 = plot_mean[,1]
# v2 = v1 - qt(0.975, plot_n[,1]-1) * plot_sd[,1] / sqrt(plot_n[,1])
v2 = v1 - 1.96 * plot_sd[,1] / sqrt(plot_n[,1])
# v3 = v1 + qt(0.975, plot_n[,1]-1) * plot_sd[,1] / sqrt(plot_n[,1])
v3 = v1 + 1.96 * plot_sd[,1] / sqrt(plot_n[,1])
points(v1,lwd=3,type="l",col=cccol[1])
polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
v1 = plot_mean[,2]
# v2 = v1 - qt(0.975, plot_n[,2]-1) * plot_sd[,2] / sqrt(plot_n[,2])
v2 = v1 - 1.96 * plot_sd[,2] / sqrt(plot_n[,2])
# v3 = v1 + qt(0.975, plot_n[,2]-1) * plot_sd[,2] / sqrt(plot_n[,2])
v3 = v1 + 1.96 * plot_sd[,2] / sqrt(plot_n[,2])
points(v1,lwd=3,type="l",col=cccol[2])
polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
v1 = plot_mean[,3]
# v2 = v1 - qt(0.975, plot_n[,3]-1) * plot_sd[,3] / sqrt(plot_n[,3])
v2 = v1 - 1.96 * plot_sd[,3] / sqrt(plot_n[,3])
# v3 = v1 + qt(0.975, plot_n[,3]-1) * plot_sd[,3] / sqrt(plot_n[,3])
v3 = v1 + 1.96 * plot_sd[,3] / sqrt(plot_n[,3])
points(v1,lwd=3,type="l",col=cccol[3])
polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
axis(side=1,at=seq(3),labels=c("1-cell -> 2cell","2-cell -> 4-cell","4-cell -> 8-cell"))
legend("topleft",c("u","d","p"),bty="n",col=cccol,lwd=2)
dev.off()


# Enhancer
C12parameters <- Enhancer[seq(2),]
C24parameters <- Enhancer[seq(3,10),]
C48parameters <- Enhancer[seq(11,42),]

plot_n <- c()
plot_mean <- c()
plot_sd <- c()

plot_data <- C12parameters
plot_n <- rbind(plot_n,apply(plot_data,2,length))
plot_mean <- rbind(plot_mean,apply(plot_data,2,mean))
plot_sd <- rbind(plot_sd,apply(plot_data,2,sd))
plot_data <- C24parameters
plot_n <- rbind(plot_n,apply(plot_data,2,length))
plot_mean <- rbind(plot_mean,apply(plot_data,2,mean))
plot_sd <- rbind(plot_sd,apply(plot_data,2,sd))
plot_data <- C48parameters
plot_n <- rbind(plot_n,apply(plot_data,2,length))
plot_mean <- rbind(plot_mean,apply(plot_data,2,mean))
plot_sd <- rbind(plot_sd,apply(plot_data,2,sd))

pdf("FigR1.MethyRatioTransitionPerStage.Enhancer.pdf",width=4,height=4)
par(mar=c(4,4,4,4))
plot(1,type="n",xaxt="n",bty="n",xlab="",ylab="Parameter value",ylim=c(0,1),xlim=c(1,3),main="Human Enhancer")
box(lwd=2)
v1 = plot_mean[,1]
# v2 = v1 - qt(0.975, plot_n[,1]-1) * plot_sd[,1] / sqrt(plot_n[,1])
v2 = v1 - 1.96 * plot_sd[,1] / sqrt(plot_n[,1])
# v3 = v1 + qt(0.975, plot_n[,1]-1) * plot_sd[,1] / sqrt(plot_n[,1])
v3 = v1 + 1.96 * plot_sd[,1] / sqrt(plot_n[,1])
points(v1,lwd=3,type="l",col=cccol[1])
polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
v1 = plot_mean[,2]
# v2 = v1 - qt(0.975, plot_n[,2]-1) * plot_sd[,2] / sqrt(plot_n[,2])
v2 = v1 - 1.96 * plot_sd[,2] / sqrt(plot_n[,2])
# v3 = v1 + qt(0.975, plot_n[,2]-1) * plot_sd[,2] / sqrt(plot_n[,2])
v3 = v1 + 1.96 * plot_sd[,2] / sqrt(plot_n[,2])
points(v1,lwd=3,type="l",col=cccol[2])
polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
v1 = plot_mean[,3]
# v2 = v1 - qt(0.975, plot_n[,3]-1) * plot_sd[,3] / sqrt(plot_n[,3])
v2 = v1 - 1.96 * plot_sd[,3] / sqrt(plot_n[,3])
# v3 = v1 + qt(0.975, plot_n[,3]-1) * plot_sd[,3] / sqrt(plot_n[,3])
v3 = v1 + 1.96 * plot_sd[,3] / sqrt(plot_n[,3])
points(v1,lwd=3,type="l",col=cccol[3])
polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
axis(side=1,at=seq(3),labels=c("1-cell -> 2cell","2-cell -> 4-cell","4-cell -> 8-cell"))
legend("topleft",c("u","d","p"),bty="n",col=cccol,lwd=2)
dev.off()