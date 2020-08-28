###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E5003380","#0F231F80","#1187CD80")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)

HalfErrorBar <- function(x, y, upper, lower=upper, length=0.03,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	arrows(x,y+upper, x, y, angle=90, code=1, length=length, ...)
}


###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
data <- read.table("GSE130735_mmEmbryo_WGBS_DNMTKO_parameters_100.txt",header=T)
samples <- rownames(data)
WT <- samples[grep("WT",samples)]
D1KO <- samples[grep("D1KO",samples)]
D3abKO <- samples[grep("D3abKO",samples)]


###################################################################################################
################################                PLOT               ################################
###################################################################################################
pdf("FigS2b.GSE130735_mmEmbryo_WGBS_DNMTKO_parameters.pdf",width=6,height=3.2)
par(mar=c(6,6,2,2),mfrow=c(1,3))
plot_n <- c(length(data[WT,1]),length(data[D1KO,1]),length(data[D3abKO,1]))
plot_mean <- c(mean(data[WT,1]),mean(data[D1KO,1]),mean(data[D3abKO,1]))
plot_sd <- c(sd(data[WT,1]),sd(data[D1KO,1]),sd(data[D3abKO,1]))
bp <- barplot(plot_mean,col=cccol80[1],border=cccol[1],ylab="parameter value (u)",width=1, space=c(0.8,0.2,0.2),xlim=c(0,5),names.arg=c("WT","Dnmt1-/-","DKO"),las=2,ylim=c(0,0.4));box(lwd=2)
beeswarm(at=bp,list(data[WT,1],data[D1KO,1],data[D3abKO,1]),xaxt="n",spacing = 0.5,col = cccol80[1],pch = 16,method = "swarm",corral = "wrap",add=T)
# HalfErrorBar(bp,plot_mean, qt(0.975, plot_n-1) * plot_sd / sqrt(plot_n),col=cccol[1],lwd=1)
# HalfErrorBar(bp,plot_mean, 1.96*plot_sd/sqrt(plot_n),col=cccol[1],lwd=1)

plot_n <- c(length(data[WT,2]),length(data[D1KO,2]),length(data[D3abKO,2]))
plot_mean <- c(mean(data[WT,2]),mean(data[D1KO,2]),mean(data[D3abKO,2]))
plot_sd <- c(sd(data[WT,2]),sd(data[D1KO,2]),sd(data[D3abKO,2]))
bp <- barplot(plot_mean,col=cccol80[2],border=cccol[2],ylab="parameter value (d)",width=1, space=c(0.8,0.2,0.2),xlim=c(0,5),names.arg=c("WT","Dnmt1-/-","DKO"),las=2,ylim=c(0,0.4));box(lwd=2)
beeswarm(at=bp,list(data[WT,2],data[D1KO,2],data[D3abKO,2]),xaxt="n",spacing = 0.5,col = cccol80[2],pch = 16,method = "swarm",corral = "wrap",add=T)
# HalfErrorBar(bp,plot_mean, qt(0.975, plot_n-1) * plot_sd / sqrt(plot_n),col=cccol[2],lwd=1)
# HalfErrorBar(bp,plot_mean, 1.96*plot_sd/sqrt(plot_n),col=cccol[2],lwd=1)

plot_n <- c(length(data[WT,3]),length(data[D1KO,3]),length(data[D3abKO,3]))
plot_mean <- c(mean(data[WT,3]),mean(data[D1KO,3]),mean(data[D3abKO,3]))
plot_sd <- c(sd(data[WT,3]),sd(data[D1KO,3]),sd(data[D3abKO,3]))
bp <- barplot(plot_mean,col=cccol80[3],border=cccol[3],ylab="parameter value (p)",width=1, space=c(0.8,0.2,0.2),xlim=c(0,5),names.arg=c("WT","Dnmt1-/-","DKO"),las=2,ylim=c(0,1));box(lwd=2)
beeswarm(at=bp,list(data[WT,3],data[D1KO,3],data[D3abKO,3]),xaxt="n",spacing = 0.5,col = cccol80[3],pch = 16,method = "swarm",corral = "wrap",add=T)
# HalfErrorBar(bp,plot_mean, qt(0.975, plot_n-1) * plot_sd / sqrt(plot_n),col=cccol[3],lwd=1)
# HalfErrorBar(bp,plot_mean, 1.96*plot_sd/sqrt(plot_n),col=cccol[3],lwd=1)
dev.off()
