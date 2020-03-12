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
mm_methyclass <- read.table("../Data/mm_DNAMethylationLevel_Promoter.txt",header=T,row.names=1)
mm_genes <- rownames(mm_methyclass)
mm_all_cells <- colnames(mm_methyclass)
# selected 1cell :
#         "Zygote_1"
# selected 2cell :
#         "C2_e1_1","C2_e1_2"
# selected 4cell :
#         "C4_e6_1","C4_e6_2","C4_e6_3","C4_e6_4"
# selected 8cell :
#         "C8_e3_1","C8_e3_2","C8_e3_3","C8_e3_4","C8_e3_5","C8_e3_6","C8_e3_7","C8_e3_8"
mm_Sperm_id <- grep("Sperm", mm_all_cells)
mm_Oocyte_id <- grep("Oocyte", mm_all_cells)
mm_Zygote_id <- "Zygote_1"
mm_c2_id <- c("C2_e1_1","C2_e1_2")
mm_c4_id <- c("C4_e6_1","C4_e6_2","C4_e6_3","C4_e6_4")
mm_c8_id <- c("C8_e3_1","C8_e3_2","C8_e3_3","C8_e3_4","C8_e3_5","C8_e3_6","C8_e3_7","C8_e3_8")

###################################################################################################
################################                PLOT               ################################
###################################################################################################

# Observed heterogeneity value
mm_2C_obs <- apply(na.omit(mm_methyclass[,mm_c2_id]),1,MethylationHeterogeneityObservation)
mm_4C_obs <- apply(na.omit(mm_methyclass[,mm_c4_id]),1,MethylationHeterogeneityObservation)
mm_8C_obs <- apply(na.omit(mm_methyclass[,mm_c8_id]),1,MethylationHeterogeneityObservation)

selected_Zygote <- mm_Zygote_id
c1_gene <- mm_genes[which(mm_methyclass[,selected_Zygote]=="0")]
c2_gene <- mm_genes[which(mm_methyclass[,selected_Zygote]=="0.25")]
c3_gene <- mm_genes[which(mm_methyclass[,selected_Zygote]=="0.5")]
c4_gene <- mm_genes[which(mm_methyclass[,selected_Zygote]=="0.75")]
c5_gene <- mm_genes[which(mm_methyclass[,selected_Zygote]=="1")]

f_mean <- c(mean(mm_4C_obs[c1_gene],na.rm=T),mean(mm_4C_obs[c2_gene],na.rm=T),mean(mm_4C_obs[c3_gene],na.rm=T),mean(mm_4C_obs[c4_gene],na.rm=T),mean(mm_4C_obs[c5_gene],na.rm=T))
e_mean <- c(mean(mm_8C_obs[c1_gene],na.rm=T),mean(mm_8C_obs[c2_gene],na.rm=T),mean(mm_8C_obs[c3_gene],na.rm=T),mean(mm_8C_obs[c4_gene],na.rm=T),mean(mm_8C_obs[c5_gene],na.rm=T))

pdf("SFig4AB.MethylationHeterogeneityObservation.pdf",width=7,height=4)
par(mfrow=c(1,2))
plot_list <- list(mm_4C_obs[c1_gene],mm_4C_obs[c2_gene],mm_4C_obs[c3_gene],mm_4C_obs[c4_gene],mm_4C_obs[c5_gene])
GiniIndexBoxplot_List(plot_list,paste("Observed 4-cell",sep=""))
points(seq(5),f_mean,type="b",bg=cccolBR,pch=21,lwd=2)
plot_list <- list(mm_8C_obs[c1_gene],mm_8C_obs[c2_gene],mm_8C_obs[c3_gene],mm_8C_obs[c4_gene],mm_8C_obs[c5_gene])
GiniIndexBoxplot_List(plot_list,paste("Observed 8-cell",sep=""))
points(seq(5),e_mean,type="b",bg=cccolBR,pch=21,lwd=2)
dev.off()
