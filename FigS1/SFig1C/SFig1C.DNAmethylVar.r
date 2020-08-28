###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
sys_argv <- commandArgs(T)
# sys_argv <- c("scBS_2C_11_2")
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)

Boxplot_List <- function(plot_list,MAIN,NAMES,YLAB){
	n <- length(plot_list)
	xpos <- 0:(n-1)+1.5
	p_value <- c()
	for (i in seq(n-1)){
		tryCatch({p_value <- c(p_value,t.test(plot_list[[i]],plot_list[[i+1]])$p.value)},
			error = function(err){p_value <- c(p_value,1)})
	}
	mark <- symnum(p_value, cutpoints=c(0,0.001,0.01,0.05,1), symbols=c("***","**","*","-"))
	bp <- boxplot(plot_list,at=0:(n-1)+1, xlim=c(0.5,n+0.5),outline=F,plot=F)
	ylim <- range(bp$stats,na.rm=T)
	dist <- (ylim[2]-ylim[1])/20
	ylim[2] <- ylim[2]+1*dist
	boxplot(plot_list,at=0:(n-1)+1,boxwex=0.6, xlim=c(0.5,n+0.5),ylim=ylim,col=cccol50[c(9,6,7,8)],border=cccol[c(9,6,7,8)],outline=F,names=NAMES,main=MAIN,lwd=2,ylab=YLAB,las=2)
	box(lwd=2)
	ypos_1 <- bp$stats[5,][1:n-1]
	ypos_2 <- bp$stats[5,][2:n]
	# legend("topleft",c("***:p<0.001","**:p<0.01","*:p<0.05"),text.col="red",bty="n")
	if (length(mark)!=0){
		for(i in 1:length(mark)){
			if(!is.na(mark[i]) & mark[i]!="?"){
				segments(xpos[i]-.4, ypos_1[i]+dist/2, xpos[i]-.4, max(ypos_1[i], ypos_2[i])+dist)
				segments(xpos[i]+.4, ypos_2[i]+dist/2, xpos[i]+.4, max(ypos_1[i], ypos_2[i])+dist)
				segments(xpos[i]-.4, max(ypos_1[i], ypos_2[i])+dist, xpos[i]-0.2, max(ypos_1[i], ypos_2[i])+dist)
				segments(xpos[i]+.4, max(ypos_1[i], ypos_2[i])+dist, xpos[i]+0.2, max(ypos_1[i], ypos_2[i])+dist)
				text(x=xpos[i], y=max(ypos_1[i], ypos_2[i])+dist, label=mark[i], col="black")
			}
		}
	}
	for (i in seq(n)){
		text(x=i, y=bp$stats[3,i], label=length(na.omit(plot_list[[i]])), col="white",cex=0.3)
	}
}

SFigBoxplot_List <- function(plot_list,MAIN,NAMES,YLAB){
	n <- length(plot_list)
	xpos <- 0:(n-1)+1.5
	p_value <- c()
	for (i in seq(n-1)){
		tryCatch({p_value <- c(p_value,t.test(plot_list[[i]],plot_list[[n]])$p.value)},
			error = function(err){p_value <- c(p_value,1)})
	}
	mark <- symnum(p_value, cutpoints=c(0,0.001,0.01,0.05,1), symbols=c("***","**","*","-"))
	bp <- boxplot(plot_list,at=0:(n-1)+1, xlim=c(0.5,n+0.5),outline=F,plot=F)
	ylim <- range(bp$stats,na.rm=T)
	dist <- (ylim[2]-ylim[1])/20
	ylim[2] <- ylim[2]+2*dist
	boxplot(plot_list,at=0:(n-1)+1,boxwex=0.6, xlim=c(0.5,n+0.5),ylim=ylim,col=cccol50[c(6,7,8)],border=cccol[c(6,7,8)],outline=F,names=NAMES,main=MAIN,lwd=2,ylab=YLAB,las=2)
	box(lwd=2)
	
	for (i in seq(n-1)){
		text(x=i, y=bp$stats[5,i]+dist, label=mark[i], col=cccol[8],cex=1)
	}

	# for (i in seq(n)){
	# 	text(x=i, y=bp$stats[3,i], label=length(na.omit(plot_list[[i]])), col="white",cex=0.3)
	# }
}
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
name <- sys_argv[1]
# each,var,std,mean,c_counts
PromoterRegion <- read.table(paste("DNAmethylVar_",name,"_PromoterRegion.txt",sep=""),row.names=1)
GenebodyRegion <- read.table(paste("DNAmethylVar_",name,"_GenebodyRegion.txt",sep=""),row.names=1)
RandomRegion <- read.table(paste("DNAmethylVar_",name,"_RandomRegion.txt",sep=""),row.names=1)

###################################################################################################
################################                PLOT               ################################
###################################################################################################
pdf(paste("SFig1C.",name,"_DNAmethylVariance.pdf",sep=""),width=3,height=4.2)
par(mgp=c(4,1,0),mar=c(8,6,4,2))
CpGCut <- 50
SFigBoxplot_List(list(PromoterRegion[which(PromoterRegion[,4] > CpGCut),1],GenebodyRegion[which(GenebodyRegion[,4] > CpGCut),1],RandomRegion[which(RandomRegion[,4] > CpGCut),1]),"",c("Promoter","Genebody","Random region"),"Methylation difference among the region")
dev.off()

# pdf(paste("SFig1C.",name,"_DNAmethylCV2.pdf",sep=""),width=3,height=4.2)
# par(mgp=c(4,1,0),mar=c(8,6,4,2))
# CpGCut <- 50
# Boxplot_List(list(PromoterRegion[which(PromoterRegion[,4] > CpGCut),1]/(PromoterRegion[which(PromoterRegion[,4] > CpGCut),3])^2,GenebodyRegion[which(GenebodyRegion[,4] > CpGCut),1]/(GenebodyRegion[which(GenebodyRegion[,4] > CpGCut),3])^2,CGIRegion[which(CGIRegion[,4] > CpGCut),1]/(CGIRegion[which(CGIRegion[,4] > CpGCut),3])^2,RandomRegion[which(RandomRegion[,4] > CpGCut),1]/(RandomRegion[which(RandomRegion[,4] > CpGCut),3])^2,CGI[which(CGI[,4] > CpGCut),1]/(CGI[which(CGI[,4] > CpGCut),3])^2),paste("CpG > ",CpGCut,"_Region_CV2",sep=""),c("PromoterRegion","GenebodyRegion","CGIRegion","RandomRegion","CGI"),"Squared Coefficient of Variation")
# dev.off()