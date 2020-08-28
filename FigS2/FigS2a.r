###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E5003380","#0F231F80","#1187CD80")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)

red_col <- c(cccol[1],cccol80[1],cccol50[1])
blue_col <- c(cccol[2],cccol80[2],cccol50[2])
yellow_col <- c(cccol[3],cccol80[3],cccol50[3])
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
para <- read.table("NewtonParameters_CellCycle_udp_eachStage.txt",row.names=1)
samples <- rownames(para)

dr_epibolyT24hpf_1OD <- grep("dr_epibolyT24hpf_1OD_10",samples)
dr_epibolyT24hpftet123_1OD <- grep("dr_epibolyT24hpftet123_1OD_10",samples)
dr_epibolyT24hpftet3_1OD <- grep("dr_epibolyT24hpftet3_1OD_10",samples)

# dr_epibolyT24hpf_1OD <- grep("dr_epibolyT24hpf_1OD_",samples)
# dr_epibolyT24hpftet123_1OD <- grep("dr_epibolyT24hpftet123_1OD_",samples)
# dr_epibolyT24hpftet3_1OD <- grep("dr_epibolyT24hpftet3_1OD_",samples)

# x <- seq(20)

###################################################################################################
################################                PLOT               ################################
###################################################################################################
# pdf("FigS2a.MethylationParamatersinDRdevelopment_CellCycle.pdf",width=6,height=8)
# par(mfrow=c(3,1),mar=c(4,4,2,2))
# plot(1,type="n",xaxt="n",bty="n",xlab="cell cycle",ylab="Parameter value (u)",main="",ylim=c(0.05,0.12),xlim=c(1,length(x)));box(lwd=2)
# axis(side=1,at=seq(length(x)),labels=x,las=2)
# points(para[dr_epibolyT24hpf_1OD,1],pch=16,type="p",lwd=2,col=red_col[1])
# points(para[dr_epibolyT24hpftet3_1OD,1],pch=2,type="p",lwd=2,col=red_col[2])
# points(para[dr_epibolyT24hpftet123_1OD,1],pch=8,type="p",lwd=2,col=red_col[3])
# abline(v=seq(19)+0.5,lty=2,col=cccol50[8])
# legend("bottomright",c("epiboly~24hpf","epiboly~24hpf Tet3-/-","epiboly~24hpf Tet1/2/3-/-"),col=red_col,pch=c(16,2,8),bty="n")
# plot(1,type="n",xaxt="n",bty="n",xlab="cell cycle",ylab="Parameter value (d)",main="",ylim=c(0.09,0.16),xlim=c(1,length(x)));box(lwd=2)
# axis(side=1,at=seq(length(x)),labels=x,las=2)
# points(para[dr_epibolyT24hpf_1OD,2],pch=16,type="p",lwd=2,col=blue_col[1])
# points(para[dr_epibolyT24hpftet3_1OD,2],pch=2,type="p",lwd=2,col=blue_col[2])
# points(para[dr_epibolyT24hpftet123_1OD,2],pch=8,type="p",lwd=2,col=blue_col[3])
# abline(v=seq(19)+0.5,lty=2,col=cccol50[8])
# legend("topright",c("epiboly~24hpf","epiboly~24hpf Tet3-/-","epiboly~24hpf Tet1/2/3-/-"),col=blue_col,pch=c(16,2,8),bty="n")
# plot(1,type="n",xaxt="n",bty="n",xlab="cell cycle",ylab="Parameter value (p)",main="",ylim=c(0.83,0.90),xlim=c(1,length(x)));box(lwd=2)
# axis(side=1,at=seq(length(x)),labels=x,las=2)
# points(para[dr_epibolyT24hpf_1OD,3],pch=16,type="p",lwd=2,col=yellow_col[1])
# points(para[dr_epibolyT24hpftet3_1OD,3],pch=2,type="p",lwd=2,col=yellow_col[2])
# points(para[dr_epibolyT24hpftet123_1OD,3],pch=8,type="p",lwd=2,col=yellow_col[3])
# abline(v=seq(19)+0.5,lty=2,col=cccol50[8])
# legend("bottomright",c("epiboly~24hpf","epiboly~24hpf Tet3-/-","epiboly~24hpf Tet1/2/3-/-"),col=yellow_col,pch=c(16,2,8),bty="n")
# dev.off()

pdf("FigS2a.MethylationParamatersinDRdevelopment_10CellCycle.pdf",width=6,height=3.2)
par(mar=c(6,6,2,2),mfrow=c(1,3))
plot_n <- c(length(para[dr_epibolyT24hpf_1OD,1]),length(para[dr_epibolyT24hpftet3_1OD,1]),length(para[dr_epibolyT24hpftet123_1OD,1]))
plot_mean <- c(mean(para[dr_epibolyT24hpf_1OD,1]),mean(para[dr_epibolyT24hpftet3_1OD,1]),mean(para[dr_epibolyT24hpftet123_1OD,1]))
plot_sd <- c(sd(para[dr_epibolyT24hpf_1OD,1]),sd(para[dr_epibolyT24hpftet3_1OD,1]),sd(para[dr_epibolyT24hpftet123_1OD,1]))
bp <- barplot(plot_mean,col=cccol80[1],border=cccol[1],ylab="parameter value (u)",width=1, space=c(0.8,0.2,0.2),xlim=c(0,5),names.arg=c("epiboly~24hpf","epiboly~24hpf Tet3-/-","epiboly~24hpf Tet1/2/3-/-"),las=2,ylim=c(0.084,0.1),xpd=F);box(lwd=2)
# beeswarm(at=bp,list(para[dr_epibolyT24hpf_1OD,1],para[dr_epibolyT24hpftet3_1OD,1],para[dr_epibolyT24hpftet123_1OD,1]),xaxt="n",spacing = 0.5,col = cccol80[1],pch = 16,method = "swarm",corral = "wrap",add=T)
# HalfErrorBar(bp,plot_mean, qt(0.975, plot_n-1) * plot_sd / sqrt(plot_n),col=cccol[1],lwd=1)
# HalfErrorBar(bp,plot_mean, 1.96*plot_sd/sqrt(plot_n),col=cccol[1],lwd=1)

plot_n <- c(length(para[dr_epibolyT24hpf_1OD,2]),length(para[dr_epibolyT24hpftet3_1OD,2]),length(para[dr_epibolyT24hpftet123_1OD,2]))
plot_mean <- c(mean(para[dr_epibolyT24hpf_1OD,2]),mean(para[dr_epibolyT24hpftet3_1OD,2]),mean(para[dr_epibolyT24hpftet123_1OD,2]))
plot_sd <- c(sd(para[dr_epibolyT24hpf_1OD,2]),sd(para[dr_epibolyT24hpftet3_1OD,2]),sd(para[dr_epibolyT24hpftet123_1OD,2]))
bp <- barplot(plot_mean,col=cccol80[2],border=cccol[2],ylab="parameter value (d)",width=1, space=c(0.8,0.2,0.2),xlim=c(0,5),names.arg=c("epiboly~24hpf","epiboly~24hpf Tet3-/-","epiboly~24hpf Tet1/2/3-/-"),las=2,ylim=c(0.084,0.1),xpd=F);box(lwd=2)
# beeswarm(at=bp,list(para[dr_epibolyT24hpf_1OD,2],para[dr_epibolyT24hpftet3_1OD,2],para[dr_epibolyT24hpftet123_1OD,2]),xaxt="n",spacing = 0.5,col = cccol80[2],pch = 16,method = "swarm",corral = "wrap",add=T)
# HalfErrorBar(bp,plot_mean, qt(0.975, plot_n-1) * plot_sd / sqrt(plot_n),col=cccol[2],lwd=1)
# HalfErrorBar(bp,plot_mean, 1.96*plot_sd/sqrt(plot_n),col=cccol[2],lwd=1)

plot_n <- c(length(para[dr_epibolyT24hpf_1OD,3]),length(para[dr_epibolyT24hpftet3_1OD,3]),length(para[dr_epibolyT24hpftet123_1OD,3]))
plot_mean <- c(mean(para[dr_epibolyT24hpf_1OD,3]),mean(para[dr_epibolyT24hpftet3_1OD,3]),mean(para[dr_epibolyT24hpftet123_1OD,3]))
plot_sd <- c(sd(para[dr_epibolyT24hpf_1OD,3]),sd(para[dr_epibolyT24hpftet3_1OD,3]),sd(para[dr_epibolyT24hpftet123_1OD,3]))
bp <- barplot(plot_mean,col=cccol80[3],border=cccol[3],ylab="parameter value (p)",width=1, space=c(0.8,0.2,0.2),xlim=c(0,5),names.arg=c("epiboly~24hpf","epiboly~24hpf Tet3-/-","epiboly~24hpf Tet1/2/3-/-"),las=2,ylim=c(0.884,0.9),xpd=F);box(lwd=2)
# beeswarm(at=bp,list(para[dr_epibolyT24hpf_1OD,3],para[dr_epibolyT24hpftet3_1OD,3],para[dr_epibolyT24hpftet123_1OD,3]),xaxt="n",spacing = 0.5,col = cccol80[3],pch = 16,method = "swarm",corral = "wrap",add=T)
# HalfErrorBar(bp,plot_mean, qt(0.975, plot_n-1) * plot_sd / sqrt(plot_n),col=cccol[3],lwd=1)
# HalfErrorBar(bp,plot_mean, 1.96*plot_sd/sqrt(plot_n),col=cccol[3],lwd=1)
dev.off()
