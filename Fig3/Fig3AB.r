###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001330","#16557A30","#C7A60930","#87C23230","#64C0AB30","#A14C9430","#15A08C30","#8B7E7530","#1E7CAF30","#EA425F30","#46489A30","#E5003330","#0F231F30","#1187CD30")
cccolBR <- c("#16557A","#443F60","#712A46","#A0152C","#CE0013")
cccolBR30 <- c("#16557A30","#443F6030","#712A4630","#A0152C30","#CE001330")

library(ineq)
MethylationHeterogeneityObservation <- function(x){
	a1 <- length(which(x==0))
	a2 <- length(which(x==1/4))
	a3 <- length(which(x==1/2))
	a4 <- length(which(x==3/4))
	a5 <- length(which(x==1))
	total <- a1+a2+a3+a4+a5
	if (total >=2){
		# return(1-Gini(c(a1/total,a2/total,a3/total,a4/total,a5/total)))
		return(1-Gini(c(a1/total,a2/total,a3/total,a4/total,a5/total),corr=T))
	}else{
		return(NA)
	}
}

GiniIndexBoxplot_List <- function(plot_list,out_name){
	n <- 5
	xpos <- 0:(n-1)+1.5
	p_value <- c()
	for (i in seq(n-1)){
		tryCatch({p_value <- c(p_value,t.test(plot_list[[i]],plot_list[[i+1]])$p.value)},
			error = function(err){p_value <- c(p_value,1)})
	}
	mark <- symnum(p_value, cutpoints=c(0,0.001,0.01,0.05,1), symbols=c("***","**","*","-"))
	bp <- boxplot(plot_list,at=0:(n-1)+1, xlim=c(0.5,n+0.5),outline=F,plot=F)
	ylim <- c(0,1)
	dist <- (ylim[2]-ylim[1])/20
	ylim[2] <- ylim[2]
	boxplot(plot_list,at=0:(n-1)+1,boxwex=0.6,xlim=c(0.5,n+0.5),ylim=ylim,col=cccolBR30,border=cccolBR,outline=F,names=c("S1","S2","S3","S4","S5"),main=out_name,lwd=2,ylab="Heterogeneity")
	box(lwd=2)
	ypos_1 <- bp$stats[5,][1:n-1]
	ypos_2 <- bp$stats[5,][2:n]
	# legend("topleft",c("***:p<0.001","**:p<0.01","*:p<0.05"),text.col=cccol[1],bty="n")
	for(i in 1:length(mark)){
		if(!is.na(mark[i])){
			segments(xpos[i]-.45, ypos_1[i]+dist/2, xpos[i]-.45, max(ypos_1[i], ypos_2[i])+dist)
			segments(xpos[i]+.45, ypos_2[i]+dist/2, xpos[i]+.45, max(ypos_1[i], ypos_2[i])+dist)
			segments(xpos[i]-.45, max(ypos_1[i], ypos_2[i])+dist, xpos[i]-0.25, max(ypos_1[i], ypos_2[i])+dist)
			segments(xpos[i]+.45, max(ypos_1[i], ypos_2[i])+dist, xpos[i]+0.25, max(ypos_1[i], ypos_2[i])+dist)
			text(x=xpos[i], y=max(ypos_1[i], ypos_2[i])+dist, label=mark[i], col="red")
		}
	}
}
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
######## methylation class
hg_methyclass <- read.table("../../Data/hg_MethylationClass.txt",header=T,row.names=1)

hg_genes <- row.names(hg_methyclass)
hg_all_cells <- colnames(hg_methyclass) 

hg_GV_id <- grep("GV", hg_all_cells)
hg_MII_id <- grep("MII", hg_all_cells)
hg_Sp_id <- grep("Sp", hg_all_cells)
hg_LPN_id <- grep("LPN", hg_all_cells)
hg_2C_6_id <- grep("2C_6", hg_all_cells)
hg_4C_7_id <- grep("4C_7", hg_all_cells)
hg_4C_8_id <- grep("4C_8", hg_all_cells)
hg_8C_6_id <- grep("8C_6", hg_all_cells)
hg_8C_9_id <- grep("8C_9", hg_all_cells)
hg_Mor_id <- grep("Mor", hg_all_cells)
hg_ICM_id <- grep("ICM", hg_all_cells)
hg_TE_id <- grep("TE", hg_all_cells)
hg_BST_id <- grep("BST", hg_all_cells)

# hg_stages <- c("GV","MII","Sperm","EPN","MPN","LPN","2C","4C","8C","Morula","ICM","TE","Blastocyst")

###################################################################################################
################################                PLOT               ################################
###################################################################################################

# Observed heterogeneity value
hg_2C_6_obs <- apply(na.omit(hg_methyclass[,hg_2C_6_id]),1,MethylationHeterogeneityObservation)
hg_4C_7_obs <- apply(na.omit(hg_methyclass[,hg_4C_7_id]),1,MethylationHeterogeneityObservation)
# hg_4C_8_obs <- apply(na.omit(hg_methyclass[,hg_4C_8_id]),1,MethylationHeterogeneityObservation)
# hg_8C_6_obs <- apply(na.omit(hg_methyclass[,hg_8C_6_id]),1,MethylationHeterogeneityObservation)
hg_8C_9_obs <- apply(na.omit(hg_methyclass[,hg_8C_9_id]),1,MethylationHeterogeneityObservation)

selected_LPN <- "hg_LPN_8"
c1_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="0")]
c2_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="0.25")]
c3_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="0.5")]
c4_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="0.75")]
c5_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="1")]
# pdf("Fig2DEF.MethylationHeterogeneityObservation.pdf",width=7,height=4)
# par(mfrow=c(1,3))
# plot_list <- list(hg_2C_6_obs[c1_gene],hg_2C_6_obs[c2_gene],hg_2C_6_obs[c3_gene],hg_2C_6_obs[c4_gene],hg_2C_6_obs[c5_gene])
# GiniIndexBoxplot_List(plot_list,paste("Observed LPN->2C",sep=""))
# plot_list <- list(hg_4C_7_obs[c1_gene],hg_4C_7_obs[c2_gene],hg_4C_7_obs[c3_gene],hg_4C_7_obs[c4_gene],hg_4C_7_obs[c5_gene])
# GiniIndexBoxplot_List(plot_list,paste("Observed LPN->4C",sep=""))
# # plot_list <- list(hg_4C_8_obs[c1_gene],hg_4C_8_obs[c2_gene],hg_4C_8_obs[c3_gene],hg_4C_8_obs[c4_gene],hg_4C_8_obs[c5_gene])
# # GiniIndexBoxplot_List(plot_list,paste("Observed LPN->4C",sep=""))
# plot_list <- list(hg_8C_6_obs[c1_gene],hg_8C_6_obs[c2_gene],hg_8C_6_obs[c3_gene],hg_8C_6_obs[c4_gene],hg_8C_6_obs[c5_gene])
# GiniIndexBoxplot_List(plot_list,paste("Observed LPN->8C",sep=""))
# # plot_list <- list(hg_8C_9_obs[c1_gene],hg_8C_9_obs[c2_gene],hg_8C_9_obs[c3_gene],hg_8C_9_obs[c4_gene],hg_8C_9_obs[c5_gene])
# # GiniIndexBoxplot_List(plot_list,paste("Observed LPN->8C",sep=""))
# dev.off()

f_mean <- c(mean(hg_4C_7_obs[c1_gene],na.rm=T),mean(hg_4C_7_obs[c2_gene],na.rm=T),mean(hg_4C_7_obs[c3_gene],na.rm=T),mean(hg_4C_7_obs[c4_gene],na.rm=T),mean(hg_4C_7_obs[c5_gene],na.rm=T))
e_mean <- c(mean(hg_8C_9_obs[c1_gene],na.rm=T),mean(hg_8C_9_obs[c2_gene],na.rm=T),mean(hg_8C_9_obs[c3_gene],na.rm=T),mean(hg_8C_9_obs[c4_gene],na.rm=T),mean(hg_8C_9_obs[c5_gene],na.rm=T))

pdf("Fig3AB.MethylationHeterogeneityObservation.pdf",width=7,height=4)
par(mfrow=c(1,2))
plot_list <- list(hg_4C_7_obs[c1_gene],hg_4C_7_obs[c2_gene],hg_4C_7_obs[c3_gene],hg_4C_7_obs[c4_gene],hg_4C_7_obs[c5_gene])
GiniIndexBoxplot_List(plot_list,paste("Observed 4-cell",sep=""))
points(seq(5),f_mean,type="b",bg=cccolBR,pch=21,lwd=2)
plot_list <- list(hg_8C_9_obs[c1_gene],hg_8C_9_obs[c2_gene],hg_8C_9_obs[c3_gene],hg_8C_9_obs[c4_gene],hg_8C_9_obs[c5_gene])
GiniIndexBoxplot_List(plot_list,paste("Observed 8-cell",sep=""))
points(seq(5),e_mean,type="b",bg=cccolBR,pch=21,lwd=2)
dev.off()
