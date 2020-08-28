###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#64C0AB50","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
# sys_argv <- commandArgs(T)
# sample <- sys_argv[1]
sample <- "hg_LPN_8"

library("MethylTransition")

HalfErrorBar <- function(x, y, upper, lower=upper, length=0.05,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	for (i in seq(length(x))){
		if (y[i]>0){arrows(x[i],y[i]+upper[i], x[i], y[i], angle=90, code=1, length=length, ...)}
		else{arrows(x[i],y[i],x[i],y[i]-lower[i], angle=90, code=2, length=length, ...)}
	}
}

Boxplot_List <- function(plot_list,out_name,YLAB,XNAMES=c("LHG","MHG","HHG"),YLIM){
	n <- length(plot_list)
	xpos <- 0:(n-1)+1.5
	p_value <- c()
	for (i in seq(n-1)){
		p_value <- tryCatch(c(p_value,t.test(plot_list[[i]],plot_list[[i+1]])$p.value),error=function(a) c(p_value,-1))
	}
	mark <- symnum(p_value, cutpoints=c(-1,-0.1,0.001,0.01,0.05,1), symbols=c(NA,"***","**","*","-"))
	bp <- boxplot(plot_list,at=0:(n-1)+1, xlim=c(0.5,n+0.5),outline=F,plot=F)
	ylim <- range(bp$stats,na.rm=T)
	dist <- (ylim[2]-ylim[1])/20
	# ylim[2] <- ylim[2]+4*dist
	ylim[2] <- ylim[2]+1*dist
	boxplot(plot_list,at=0:(n-1)+1, xlim=c(0.5,n+0.5),ylim=YLIM,col=cccol50,border=cccol,outline=F,names=XNAMES,main=out_name,lwd=2,ylab=YLAB,las=2)
	box(lwd=2)
	ypos_1 <- bp$stats[5,][1:n-1]
	ypos_2 <- bp$stats[5,][2:n]
	# legend("topleft",c("***:p<0.001","**:p<0.01","*:p<0.05"),text.col=cccol[1],bty="n")
	for(i in 1:length(mark)){
		if(!is.na(mark[i])){
			segments(xpos[i]-.4, ypos_1[i]+dist/2, xpos[i]-.4, max(ypos_1[i], ypos_2[i])+dist)
			segments(xpos[i]+.4, ypos_2[i]+dist/2, xpos[i]+.4, max(ypos_1[i], ypos_2[i])+dist)
			segments(xpos[i]-.4, max(ypos_1[i], ypos_2[i])+dist, xpos[i]-0.2, max(ypos_1[i], ypos_2[i])+dist)
			segments(xpos[i]+.4, max(ypos_1[i], ypos_2[i])+dist, xpos[i]+0.2, max(ypos_1[i], ypos_2[i])+dist)
			text(x=xpos[i], y=max(ypos_1[i], ypos_2[i])+dist, label=mark[i], col="red")
		}
	}
}

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################

######## methylation class
hg_methyclass <- read.table("../Data/hg_MethylationClass.txt",header=T,row.names=1)

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


# Different Classes of genes 
hg_4C_7_c1_MHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c1_MHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c2_MHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c2_MHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c3_MHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c3_MHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c4_MHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c4_MHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c5_MHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c5_MHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c1_MHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c1_MHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c2_MHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c2_MHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c3_MHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c3_MHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c4_MHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c4_MHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c5_MHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c5_MHG.txt")[,1]),error = function(a) NA)

hg_4C_7_c1_LHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c1_LHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c2_LHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c2_LHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c3_LHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c3_LHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c4_LHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c4_LHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c5_LHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c5_LHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c1_LHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c1_LHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c2_LHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c2_LHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c3_LHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c3_LHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c4_LHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c4_LHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c5_LHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c5_LHG.txt")[,1]),error = function(a) NA)

hg_4C_7_c1_HHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c1_HHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c2_HHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c2_HHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c3_HHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c3_HHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c4_HHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c4_HHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c5_HHG <- tryCatch(as.vector(read.table("../Fig4/hg_4C_7_c5_HHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c1_HHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c1_HHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c2_HHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c2_HHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c3_HHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c3_HHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c4_HHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c4_HHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c5_HHG <- tryCatch(as.vector(read.table("../Fig4/hg_8C_9_c5_HHG.txt")[,1]),error = function(a) NA)

# CGratio
# CG_ratio <- read.table("../Data/hg19_Promoter.CGratio",header=F,row.names=1)
CG_ratio <- round(read.table("../Data/hg19.CGRatio.txt",header=F,row.names=1),5)
# GCcontent
GC_content <- round(read.table("../Data/hg19.GCContent.txt",header=F,row.names=1),5)
# DHSratio
DHS_ratio_4C_1 <- read.table("../Data/hg_4C_DHS_r1.PromoterSignal",header=F,row.names=1)
DHS_ratio_8C_1 <- read.table("../Data/hg_8C_DHS_r1.PromoterSignal",header=F,row.names=1)
ATAC_ratio_4C_1 <- read.table("../Data/GSE101571_hg_C4_1.ATAC.PromoterSignal",header=F,row.names=1)
ATAC_ratio_8C_1 <- read.table("../Data/GSE101571_hg_C8_1.ATAC.PromoterSignal",header=F,row.names=1)
# H3K9me3 PromoterSignal (ES)
H3K9me3_ES <- read.table("../Data/H3K9me3_r4.PromoterSignal",header=F,row.names=1)
# H3K27me3 PromoterSignal (8-Cell Cut&Run)
H3K27me3_8C <- read.table("../Data/hg_8C_2PN_H3K27me3_r1.PromoterSignal",header=F,row.names=1)
H3K27me3_4C <- read.table("../Data/hg_4C_3PN_H3K27me3_r1.PromoterSignal",header=F,row.names=1)
# H3K4me3 PromoterSignal (8-Cell Cut&Run)
H3K4me3_8C <- read.table("../Data/hg_8C_2PN_H3K4me3_r1.PromoterSignal",header=F,row.names=1)
H3K4me3_4C <- read.table("../Data/hg_4C_3PN_H3K4me3_r1.PromoterSignal",header=F,row.names=1)

hg19_GeneAgeClass <- read.table("../Data/hg19.GeneAgeClass",row.names=1)
###################################################################################################
################################                PLOT               ################################
###################################################################################################

pdf(paste("Fig4C.Boxplot_CGRatio_",sample,".pdf",sep=""),width=9,height=3.6)
par(mfrow=c(1,5))
# Boxplot_List(list(CG_ratio[hg_4C_7_c1_LHG,],CG_ratio[hg_4C_7_c1_MHG,],CG_ratio[hg_4C_7_c1_HHG,]),"S1","CpG ratio",YLIM=c(0.1,1.1))
# Boxplot_List(list(CG_ratio[hg_4C_7_c2_LHG,],CG_ratio[hg_4C_7_c2_MHG,],CG_ratio[hg_4C_7_c2_HHG,]),"S2","",YLIM=c(0.1,1.1))
# Boxplot_List(list(CG_ratio[hg_4C_7_c3_LHG,],CG_ratio[hg_4C_7_c3_MHG,],CG_ratio[hg_4C_7_c3_HHG,]),"S3","",YLIM=c(0.1,1.1))
# Boxplot_List(list(CG_ratio[hg_4C_7_c4_LHG,],CG_ratio[hg_4C_7_c4_MHG,],CG_ratio[hg_4C_7_c4_HHG,]),"S4","",YLIM=c(0.1,1.1))
# Boxplot_List(list(CG_ratio[hg_4C_7_c5_LHG,],CG_ratio[hg_4C_7_c5_MHG,],CG_ratio[hg_4C_7_c5_HHG,]),"S5","",YLIM=c(0.1,1.1))
Boxplot_List(list(CG_ratio[hg_8C_9_c1_LHG,],CG_ratio[hg_8C_9_c1_MHG,],CG_ratio[hg_8C_9_c1_HHG,]),"S1","CpG ratio",YLIM=c(0.1,1.1))
Boxplot_List(list(CG_ratio[hg_8C_9_c2_LHG,],CG_ratio[hg_8C_9_c2_MHG,],CG_ratio[hg_8C_9_c2_HHG,]),"S2","",YLIM=c(0.1,1.1))
Boxplot_List(list(CG_ratio[hg_8C_9_c3_LHG,],CG_ratio[hg_8C_9_c3_MHG,],CG_ratio[hg_8C_9_c3_HHG,]),"S3","",YLIM=c(0.1,1.1))
Boxplot_List(list(CG_ratio[hg_8C_9_c4_LHG,],CG_ratio[hg_8C_9_c4_MHG,],CG_ratio[hg_8C_9_c4_HHG,]),"S4","",YLIM=c(0.1,1.1))
Boxplot_List(list(CG_ratio[hg_8C_9_c5_LHG,],CG_ratio[hg_8C_9_c5_MHG,],CG_ratio[hg_8C_9_c5_HHG,]),"S5","",YLIM=c(0.1,1.1))
dev.off()

pdf(paste("FigS4B.Boxplot_ATACRatio_",sample,".pdf",sep=""),width=9,height=3.6)
par(mfrow=c(1,5))
# Boxplot_List(list(ATAC_ratio_4C_1[hg_4C_7_c1_LHG,],ATAC_ratio_4C_1[hg_4C_7_c1_MHG,],ATAC_ratio_4C_1[hg_4C_7_c1_HHG,]),"S1","Chromatin accessibility",YLIM=c(0,2))
# Boxplot_List(list(ATAC_ratio_4C_1[hg_4C_7_c2_LHG,],ATAC_ratio_4C_1[hg_4C_7_c2_MHG,],ATAC_ratio_4C_1[hg_4C_7_c2_HHG,]),"S2","",YLIM=c(0,2))
# Boxplot_List(list(ATAC_ratio_4C_1[hg_4C_7_c3_LHG,],ATAC_ratio_4C_1[hg_4C_7_c3_MHG,],ATAC_ratio_4C_1[hg_4C_7_c3_HHG,]),"S3","",YLIM=c(0,2))
# Boxplot_List(list(ATAC_ratio_4C_1[hg_4C_7_c4_LHG,],ATAC_ratio_4C_1[hg_4C_7_c4_MHG,],ATAC_ratio_4C_1[hg_4C_7_c4_HHG,]),"S4","",YLIM=c(0,2))
# Boxplot_List(list(ATAC_ratio_4C_1[hg_4C_7_c5_LHG,],ATAC_ratio_4C_1[hg_4C_7_c5_MHG,],ATAC_ratio_4C_1[hg_4C_7_c5_HHG,]),"S5","",YLIM=c(0,2))
Boxplot_List(list(ATAC_ratio_8C_1[hg_8C_9_c1_LHG,],ATAC_ratio_8C_1[hg_8C_9_c1_MHG,],ATAC_ratio_8C_1[hg_8C_9_c1_HHG,]),"S1","Chromatin accessibility",YLIM=c(0,2.5))
Boxplot_List(list(ATAC_ratio_8C_1[hg_8C_9_c2_LHG,],ATAC_ratio_8C_1[hg_8C_9_c2_MHG,],ATAC_ratio_8C_1[hg_8C_9_c2_HHG,]),"S2","",YLIM=c(0,2.5))
Boxplot_List(list(ATAC_ratio_8C_1[hg_8C_9_c3_LHG,],ATAC_ratio_8C_1[hg_8C_9_c3_MHG,],ATAC_ratio_8C_1[hg_8C_9_c3_HHG,]),"S3","",YLIM=c(0,2.5))
Boxplot_List(list(ATAC_ratio_8C_1[hg_8C_9_c4_LHG,],ATAC_ratio_8C_1[hg_8C_9_c4_MHG,],ATAC_ratio_8C_1[hg_8C_9_c4_HHG,]),"S4","",YLIM=c(0,2.5))
Boxplot_List(list(ATAC_ratio_8C_1[hg_8C_9_c5_LHG,],ATAC_ratio_8C_1[hg_8C_9_c5_MHG,],ATAC_ratio_8C_1[hg_8C_9_c5_HHG,]),"S5","",YLIM=c(0,2.5))
dev.off()

pdf(paste("FigS4C.Boxplot_H3K4me3_",sample,".pdf",sep=""),width=9,height=3.6)
par(mfrow=c(1,5))
# Boxplot_List(list(H3K4me3_4C[hg_4C_7_c1_LHG,],H3K4me3_4C[hg_4C_7_c1_MHG,],H3K4me3_4C[hg_4C_7_c1_HHG,]),"S1","H3K4me3 signal",YLIM=c(0,10))
# Boxplot_List(list(H3K4me3_4C[hg_4C_7_c2_LHG,],H3K4me3_4C[hg_4C_7_c2_MHG,],H3K4me3_4C[hg_4C_7_c2_HHG,]),"S2","",YLIM=c(0,10))
# Boxplot_List(list(H3K4me3_4C[hg_4C_7_c3_LHG,],H3K4me3_4C[hg_4C_7_c3_MHG,],H3K4me3_4C[hg_4C_7_c3_HHG,]),"S3","",YLIM=c(0,10))
# Boxplot_List(list(H3K4me3_4C[hg_4C_7_c4_LHG,],H3K4me3_4C[hg_4C_7_c4_MHG,],H3K4me3_4C[hg_4C_7_c4_HHG,]),"S4","",YLIM=c(0,10))
# Boxplot_List(list(H3K4me3_4C[hg_4C_7_c5_LHG,],H3K4me3_4C[hg_4C_7_c5_MHG,],H3K4me3_4C[hg_4C_7_c5_HHG,]),"S5","",YLIM=c(0,10))
Boxplot_List(list(H3K4me3_8C[hg_8C_9_c1_LHG,],H3K4me3_8C[hg_8C_9_c1_MHG,],H3K4me3_8C[hg_8C_9_c1_HHG,]),"S1","H3K4me3 signal",YLIM=c(0,60))
Boxplot_List(list(H3K4me3_8C[hg_8C_9_c2_LHG,],H3K4me3_8C[hg_8C_9_c2_MHG,],H3K4me3_8C[hg_8C_9_c2_HHG,]),"S2","",YLIM=c(0,60))
Boxplot_List(list(H3K4me3_8C[hg_8C_9_c3_LHG,],H3K4me3_8C[hg_8C_9_c3_MHG,],H3K4me3_8C[hg_8C_9_c3_HHG,]),"S3","",YLIM=c(0,60))
Boxplot_List(list(H3K4me3_8C[hg_8C_9_c4_LHG,],H3K4me3_8C[hg_8C_9_c4_MHG,],H3K4me3_8C[hg_8C_9_c4_HHG,]),"S4","",YLIM=c(0,60))
Boxplot_List(list(H3K4me3_8C[hg_8C_9_c5_LHG,],H3K4me3_8C[hg_8C_9_c5_MHG,],H3K4me3_8C[hg_8C_9_c5_HHG,]),"S5","",YLIM=c(0,60))
dev.off()