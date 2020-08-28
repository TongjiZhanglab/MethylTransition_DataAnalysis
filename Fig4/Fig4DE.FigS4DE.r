###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#64C0AB50","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")

library("MethylTransition")
library(ineq)
MethylationHeterogeneity <- function(x){
	library(ineq)
	a1 <- length(which(x<=1/8))
	a2 <- length(which(x>1/8 & x<=3/8))
	a3 <- length(which(x>3/8 & x<=5/8))
	a4 <- length(which(x>5/8 & x<=7/8))
	a5 <- length(which(x>7/8))
	total <- a1+a2+a3+a4+a5
	if (total >=2){
		# return(1-Gini(c(a1/total,a2/total,a3/total,a4/total,a5/total)))
		return(1-Gini(c(a1/total,a2/total,a3/total,a4/total,a5/total),corr=T))	
	}else{
		return(NA)
	}
}

MethylationTransMatrix <- function(u,p,d){
	#   0-0 0-1 1-0 1-1
	# 0-0 (1-u)(1-u)  (1-u)(1-p)d d(1-u)(1-p) 0
	# 0-1 u(1-u)  (1-p)(1-u)(1-d) d(u+p-up) 0
	# 1-0 u(1-u)  (u+p-up)d (1-d)(1-u)(1-p) 0
	# 1-1 uu  (u+p-up)(1-d) (1-d)(u+p-up) 1

	Q00_00 <- (1-u)*(1-u); Q01_00 <- (1-u)*(1-p)*d; Q10_00 <- d*(1-u)*(1-p) ; Q11_00 <- 0
	Q00_01 <- u*(1-u); Q01_01 <- (1-p)*(1-u)*(1-d); Q10_01 <- d*(u+p-u*p) ; Q11_01 <- 0
	Q00_10 <- u*(1-u); Q01_10 <- (u+p-u*p)*d; Q10_10 <- (1-d)*(1-u)*(1-p) ; Q11_10 <- 0
	Q00_11 <- u*u; Q01_11 <- (u+p-u*p)*(1-d); Q10_11 <- (1-d)*(u+p-u*p) ; Q11_11 <- 1
	# Q <- matrix(c(Q00_00,Q01_00,Q10_00,Q11_00,Q00_01,Q01_01,Q10_01,Q11_01,Q00_10,Q01_10,Q10_10,Q11_10,Q00_11,Q01_11,Q10_11,Q11_11),c(4,4),byrow = T)
	# print(Q)
	x00_00 <- Q00_00; x01_00 <- 1/2*(Q00_00+Q01_00); x10_00 <- 1/2*(Q00_00+Q10_00); x11_00 <- 1/2*(Q01_00+Q10_00)
	x00_01 <- Q00_01; x01_01 <- 1/2*(Q00_01+Q01_01); x10_01 <- 1/2*(Q00_01+Q10_01); x11_01 <- 1/2*(Q01_01+Q10_01)
	x00_10 <- Q00_10; x01_10 <- 1/2*(Q00_10+Q01_10); x10_10 <- 1/2*(Q00_10+Q10_10); x11_10 <- 1/2*(Q01_10+Q10_10)
	x00_11 <- Q00_11; x01_11 <- 1/2*(Q00_11+Q01_11); x10_11 <- 1/2*(Q00_11+Q10_11); x11_11 <- 1/2*(Q01_11+Q10_11)
	# x <- matrix(c(x00_00,x01_00,x10_00,x11_00,x00_01,x01_01,x10_01,x11_01,x00_10,x01_10,x10_10,x11_10,x00_11,x01_11,x10_11,x11_11),c(4,4),byrow = T)
	# print(x)
	tmp_matrix <- c(x00_00*x00_00,1/4*(x00_00*x01_00+x00_00*x10_00+x01_00*x00_00+x10_00*x00_00),1/6*(x00_00*x11_00+x01_00*x01_00+x01_00*x10_00+x10_00*x01_00+x10_00*x10_00+x11_00*x00_00),1/4*(x01_00*x11_00+x10_00*x11_00+x11_00*x01_00+x11_00*x10_00),x11_00*x11_00,x00_00*x00_01+x00_00*x00_10+x00_01*x00_00+x00_10*x00_00,1/4*(x00_00*x01_01+x00_00*x10_01+x01_00*x00_01+x10_00*x00_01+x00_00*x01_10+x00_00*x10_10+x01_00*x00_10+x10_00*x00_10+x00_01*x01_00+x00_01*x10_00+x01_01*x00_00+x10_01*x00_00+x00_10*x01_00+x00_10*x10_00+x01_10*x00_00+x10_10*x00_00),1/6*(x00_00*x11_01+x01_00*x01_01+x01_00*x10_01+x10_00*x01_01+x10_00*x10_01+x11_00*x00_01+x00_00*x11_10+x01_00*x01_10+x01_00*x10_10+x10_00*x01_10+x10_00*x10_10+x11_00*x00_10+x00_01*x11_00+x01_01*x01_00+x01_01*x10_00+x10_01*x01_00+x10_01*x10_00+x11_01*x00_00+x00_10*x11_00+x01_10*x01_00+x01_10*x10_00+x10_10*x01_00+x10_10*x10_00+x11_10*x00_00),1/4*(x01_00*x11_01+x10_00*x11_01+x11_00*x01_01+x11_00*x10_01+x01_00*x11_10+x10_00*x11_10+x11_00*x01_10+x11_00*x10_10+x01_01*x11_00+x10_01*x11_00+x11_01*x01_00+x11_01*x10_00+x01_10*x11_00+x10_10*x11_00+x11_10*x01_00+x11_10*x10_00),x11_00*x11_01+x11_00*x11_10+x11_01*x11_00+x11_10*x11_00,x00_00*x00_11+x00_01*x00_01+x00_01*x00_10+x00_10*x00_01+x00_10*x00_10+x00_11*x00_00,1/4*(x00_00*x01_11+x00_00*x10_11+x01_00*x00_11+x10_00*x00_11+x00_01*x01_01+x00_01*x10_01+x01_01*x00_01+x10_01*x00_01+x00_01*x01_10+x00_01*x10_10+x01_01*x00_10+x10_01*x00_10+x00_10*x01_01+x00_10*x10_01+x01_10*x00_01+x10_10*x00_01+x00_10*x01_10+x00_10*x10_10+x01_10*x00_10+x10_10*x00_10+x00_11*x01_00+x00_11*x10_00+x01_11*x00_00+x10_11*x00_00),1/6*(x00_00*x11_11+x01_00*x01_11+x01_00*x10_11+x10_00*x01_11+x10_00*x10_11+x11_00*x00_11+x00_01*x11_01+x01_01*x01_01+x01_01*x10_01+x10_01*x01_01+x10_01*x10_01+x11_01*x00_01+x00_01*x11_10+x01_01*x01_10+x01_01*x10_10+x10_01*x01_10+x10_01*x10_10+x11_01*x00_10+x00_10*x11_01+x01_10*x01_01+x01_10*x10_01+x10_10*x01_01+x10_10*x10_01+x11_10*x00_01+x00_10*x11_10+x01_10*x01_10+x01_10*x10_10+x10_10*x01_10+x10_10*x10_10+x11_10*x00_10+x00_11*x11_00+x01_11*x01_00+x01_11*x10_00+x10_11*x01_00+x10_11*x10_00+x11_11*x00_00),1/4*(x01_00*x11_11+x10_00*x11_11+x11_00*x01_11+x11_00*x10_11+x01_01*x11_01+x10_01*x11_01+x11_01*x01_01+x11_01*x10_01+x01_01*x11_10+x10_01*x11_10+x11_01*x01_10+x11_01*x10_10+x01_10*x11_01+x10_10*x11_01+x11_10*x01_01+x11_10*x10_01+x01_10*x11_10+x10_10*x11_10+x11_10*x01_10+x11_10*x10_10+x01_11*x11_00+x10_11*x11_00+x11_11*x01_00+x11_11*x10_00),x11_00*x11_11+x11_01*x11_01+x11_01*x11_10+x11_10*x11_01+x11_10*x11_10+x11_11*x11_00,x00_01*x00_11+x00_10*x00_11+x00_11*x00_01+x00_11*x00_10,1/4*(x00_01*x01_11+x00_01*x10_11+x01_01*x00_11+x10_01*x00_11+x00_10*x01_11+x00_10*x10_11+x01_10*x00_11+x10_10*x00_11+x00_11*x01_01+x00_11*x10_01+x01_11*x00_01+x10_11*x00_01+x00_11*x01_10+x00_11*x10_10+x01_11*x00_10+x10_11*x00_10),1/6*(x00_01*x11_11+x01_01*x01_11+x01_01*x10_11+x10_01*x01_11+x10_01*x10_11+x11_01*x00_11+x00_10*x11_11+x01_10*x01_11+x01_10*x10_11+x10_10*x01_11+x10_10*x10_11+x11_10*x00_11+x00_11*x11_01+x01_11*x01_01+x01_11*x10_01+x10_11*x01_01+x10_11*x10_01+x11_11*x00_01+x00_11*x11_10+x01_11*x01_10+x01_11*x10_10+x10_11*x01_10+x10_11*x10_10+x11_11*x00_10),1/4*(x01_01*x11_11+x10_01*x11_11+x11_01*x01_11+x11_01*x10_11+x01_10*x11_11+x10_10*x11_11+x11_10*x01_11+x11_10*x10_11+x01_11*x11_01+x10_11*x11_01+x11_11*x01_01+x11_11*x10_01+x01_11*x11_10+x10_11*x11_10+x11_11*x01_10+x11_11*x10_10),x11_01*x11_11+x11_10*x11_11+x11_11*x11_01+x11_11*x11_10,x00_11*x00_11,1/4*(x00_11*x01_11+x00_11*x10_11+x01_11*x00_11+x10_11*x00_11),1/6*(x00_11*x11_11+x01_11*x01_11+x01_11*x10_11+x10_11*x01_11+x10_11*x10_11+x11_11*x00_11),1/4*(x01_11*x11_11+x10_11*x11_11+x11_11*x01_11+x11_11*x10_11),x11_11*x11_11)
	return(matrix(tmp_matrix,c(5,5),byrow = T))
}

TransMatrixMultiplication <- function(u,p,d,right_matrix){
	left_matrix <- MethylationTransMatrix(u,p,d)
	return (left_matrix%*%right_matrix)
}

MethylationTransion <- function(para,total,start,seed){
	para_matrix_plus <- para * total
	para_1 <- c(rep(0,para_matrix_plus[1,1]),rep(0.25,para_matrix_plus[2,1]),rep(0.5,para_matrix_plus[3,1]),rep(0.75,para_matrix_plus[4,1]),rep(1,para_matrix_plus[5,1]))
	para_2 <- c(rep(0,para_matrix_plus[1,2]),rep(0.25,para_matrix_plus[2,2]),rep(0.5,para_matrix_plus[3,2]),rep(0.75,para_matrix_plus[4,2]),rep(1,para_matrix_plus[5,2]))
	para_3 <- c(rep(0,para_matrix_plus[1,3]),rep(0.25,para_matrix_plus[2,3]),rep(0.5,para_matrix_plus[3,3]),rep(0.75,para_matrix_plus[4,3]),rep(1,para_matrix_plus[5,3]))
	para_4 <- c(rep(0,para_matrix_plus[1,4]),rep(0.25,para_matrix_plus[2,4]),rep(0.5,para_matrix_plus[3,4]),rep(0.75,para_matrix_plus[4,4]),rep(1,para_matrix_plus[5,4]))
	para_5 <- c(rep(0,para_matrix_plus[1,5]),rep(0.25,para_matrix_plus[2,5]),rep(0.5,para_matrix_plus[3,5]),rep(0.75,para_matrix_plus[4,5]),rep(1,para_matrix_plus[5,5]))
	end <- c()
	for (each in seq(total)){
		if (start[each] == 0.00){set.seed(seed*each+log2(1.1));end <- c(end,sample(para_1,1))}
		if (start[each] == 0.25){set.seed(seed*each+log2(2.2));end <- c(end,sample(para_2,1))}
		if (start[each] == 0.50){set.seed(seed*each+log2(3.3));end <- c(end,sample(para_3,1))}
		if (start[each] == 0.75){set.seed(seed*each+log2(4.4));end <- c(end,sample(para_4,1))}
		if (start[each] == 1.00){set.seed(seed*each+log2(5.5));end <- c(end,sample(para_5,1))}
	}
	return(end)
}

MeanSquareError <- function(predict,observe){
	return(mean((predict-observe)^2,na.rm=T))
}


ConditionalParametersPrediction <- function(start,genes,out_name){
	start <- start[genes]
	predict_heter <- c()
	e_para_mat <- c()
	f_para_mat <- c()
	t_para_mat <- c()
	for (each_2C in hg_all_cells[hg_2C_6_id]){
		para_LPN2C <- tryCatch(ParameterEstimation(hg_methyclass[genes,selected_LPN],hg_methyclass[genes,each_2C],iter=200)$estimated_parameters,error=function(a) return(NA))
		if (is.na(para_LPN2C) == FALSE){
			t_para_mat <- rbind(t_para_mat,para_LPN2C)
			LPN2C_matrix <- MethylationTransMatrix(para_LPN2C[1],para_LPN2C[3],para_LPN2C[2])
			for (each_4C in hg_all_cells[hg_4C_7_id]){
				para_2C4C <- tryCatch(ParameterEstimation(hg_methyclass[genes,each_2C],hg_methyclass[genes,each_4C],iter=200)$estimated_parameters,error=function(a) return(NA))
				if (is.na(para_2C4C) == FALSE){
					f_para_mat <- rbind(f_para_mat,para_2C4C)
					LPN4C_matrix <- TransMatrixMultiplication(para_2C4C[1],para_2C4C[3],para_2C4C[2],LPN2C_matrix)
					for (each_8C in hg_all_cells[hg_8C_9_id]){
						para_4C8C <- tryCatch(ParameterEstimation(hg_methyclass[genes,each_4C],hg_methyclass[genes,each_8C],iter=200)$estimated_parameters,error=function(a) return(NA))
						if (is.na(para_4C8C) == FALSE){
							e_para_mat <- rbind(e_para_mat,para_4C8C)
							LPN8C_matrix <- TransMatrixMultiplication(para_4C8C[1],para_4C8C[3],para_4C8C[2],LPN4C_matrix)
							t_cell_1 <- MethylationTransion(LPN2C_matrix,length(genes),start,1.1);names(t_cell_1) <- genes
							t_cell_2 <- MethylationTransion(LPN2C_matrix,length(genes),start,2.1);names(t_cell_2) <- genes
							f_cell_1 <- MethylationTransion(LPN4C_matrix,length(genes),t_cell_1,1.2);names(f_cell_1) <- genes
							f_cell_2 <- MethylationTransion(LPN4C_matrix,length(genes),t_cell_1,2.2);names(f_cell_2) <- genes
							f_cell_3 <- MethylationTransion(LPN4C_matrix,length(genes),t_cell_2,3.2);names(f_cell_3) <- genes
							f_cell_4 <- MethylationTransion(LPN4C_matrix,length(genes),t_cell_2,4.2);names(f_cell_4) <- genes
							e_cell_1 <- MethylationTransion(LPN8C_matrix,length(genes),f_cell_1,1.3);names(e_cell_1) <- genes
							e_cell_2 <- MethylationTransion(LPN8C_matrix,length(genes),f_cell_1,2.3);names(e_cell_2) <- genes
							e_cell_3 <- MethylationTransion(LPN8C_matrix,length(genes),f_cell_2,3.3);names(e_cell_3) <- genes
							e_cell_4 <- MethylationTransion(LPN8C_matrix,length(genes),f_cell_2,4.3);names(e_cell_4) <- genes
							e_cell_5 <- MethylationTransion(LPN8C_matrix,length(genes),f_cell_3,5.3);names(e_cell_5) <- genes
							e_cell_6 <- MethylationTransion(LPN8C_matrix,length(genes),f_cell_3,6.3);names(e_cell_6) <- genes
							e_cell_7 <- MethylationTransion(LPN8C_matrix,length(genes),f_cell_4,7.3);names(e_cell_7) <- genes
							e_cell_8 <- MethylationTransion(LPN8C_matrix,length(genes),f_cell_4,8.3);names(e_cell_8) <- genes
							e_hetergeneity <- apply(cbind(e_cell_1,e_cell_2,e_cell_3,e_cell_4,e_cell_5,e_cell_6,e_cell_7,e_cell_8),1,MethylationHeterogeneity)
							f_hetergeneity <- apply(cbind(f_cell_1,f_cell_2,f_cell_3,f_cell_4),1,MethylationHeterogeneity)
							t_hetergeneity <- apply(cbind(t_cell_1,t_cell_2),1,MethylationHeterogeneity)
							predict_heter <- rbind(predict_heter, e_hetergeneity)
						}else{predict_heter <- rbind(predict_heter, rep(NA,length(genes)))}
					}
				}else{for (i in seq(8)){predict_heter <- rbind(predict_heter, rep(NA,length(genes)))}}
			}
		}else{for (i in seq(8*4)){predict_heter <- rbind(predict_heter, rep(NA,length(genes)))}}
	}
	write.table(t_para_mat,file=paste(out_name,"_1C2C_udp.txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(f_para_mat,file=paste(out_name,"_2C4C_udp.txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(e_para_mat,file=paste(out_name,"_4C8C_udp.txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	return(predict_heter)
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

selected_LPN <- "hg_LPN_8"
c1_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="0")]
c2_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="0.25")]
c3_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="0.5")]
c4_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="0.75")]
c5_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="1")]
all_genes <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)

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

# raw heterogeneity
# hg_2C_6_obs <- apply(na.omit(hg_methyclass[,hg_2C_6_id]),1,MethylationHeterogeneity)
hg_4C_7_obs <- apply(na.omit(hg_methyclass[,hg_4C_7_id]),1,MethylationHeterogeneity)
# hg_4C_8_obs <- apply(na.omit(hg_methyclass[,hg_4C_8_id]),1,MethylationHeterogeneity)
# hg_8C_6_obs <- apply(na.omit(hg_methyclass[,hg_8C_6_id]),1,MethylationHeterogeneity)
hg_8C_9_obs <- apply(na.omit(hg_methyclass[,hg_8C_9_id]),1,MethylationHeterogeneity)

# conditional gene list 
pdf("CGRatioDHSRatioHMSignal_Distribution.pdf",width=3.6,height=4)
hist(CG_ratio[,1],breaks=100,col=cccol[6],border=cccol[6],main="CGRatio");box(lwd=2)
abline(v=c(0.4,0.6),lty=2)
hist(ATAC_ratio_4C_1[,1],breaks=100,col=cccol[7],border=cccol[7],main="ATAC_ratio_4C_1");box(lwd=2)
abline(v=c(0.2,0.8),lty=2)
hist(ATAC_ratio_8C_1[,1],breaks=100,col=cccol[7],border=cccol[7],main="ATAC_ratio_8C_1");box(lwd=2)
abline(v=c(0.2,0.8),lty=2)
hist(log2(H3K9me3_ES[,1]+1),breaks=100,col=cccol[9],border=cccol[9],main="H3K9me3_ES");box(lwd=2)
abline(v=c(1,3),lty=2)
hist(log2(H3K27me3_8C[,1]+1),breaks=20,col=cccol[8],border=cccol[8],main="H3K27me3_8C");box(lwd=2)
abline(v=c(0.5,1.5),lty=2)
hist(log2(H3K27me3_4C[,1]+1),breaks=20,col=cccol[8],border=cccol[8],main="H3K27me3_4C");box(lwd=2)
abline(v=c(0.5,1.5),lty=2)
hist(log2(H3K4me3_8C[,1]+1),breaks=20,col=cccol[10],border=cccol[10],main="H3K4me3_8C");box(lwd=2)
abline(v=c(2,4),lty=2)
hist(log2(H3K4me3_4C[,1]+1),breaks=20,col=cccol[10],border=cccol[10],main="H3K4me3_4C");box(lwd=2)
abline(v=c(2,4),lty=2)
dev.off()

HCGratio_cut <- 0.6
LCGratio_cut <- 0.4
HATACratio_cut <- 0.8
LATACratio_cut <- 0.2
HK9_cut <- 2^3-1
LK9_cut <- 2^1-1
HK27_cut <- 2^1.5-1
LK27_cut <- 2^0.5-1
HK4_cut <- 2^4-1
LK4_cut <- 2^2-1

HCGratio <- intersect(all_genes,rownames(CG_ratio)[which(CG_ratio[,1] >= HCGratio_cut)])
MCGratio <- intersect(all_genes,rownames(CG_ratio)[which((CG_ratio[,1] >= LCGratio_cut) & (CG_ratio[,1] < HCGratio_cut))])
LCGratio <- intersect(all_genes,rownames(CG_ratio)[which(CG_ratio[,1] < LCGratio_cut)])
HATACratio <- intersect(all_genes,rownames(ATAC_ratio_8C_1)[which(ATAC_ratio_8C_1[,1] >= HATACratio_cut)])
MATACratio <- intersect(all_genes,rownames(ATAC_ratio_8C_1)[which((ATAC_ratio_8C_1[,1] >= LATACratio_cut) & (ATAC_ratio_8C_1[,1] < HATACratio_cut))])
LATACratio <- intersect(all_genes,rownames(ATAC_ratio_8C_1)[which(ATAC_ratio_8C_1[,1] < LATACratio_cut)])
HK4 <- intersect(all_genes,rownames(H3K4me3_8C)[which(H3K4me3_8C[,1] >= HK4_cut)])
MK4 <- intersect(all_genes,rownames(H3K4me3_8C)[which((H3K4me3_8C[,1] >= LK4_cut) & (H3K4me3_8C[,1] < HK4_cut))])
LK4 <- intersect(all_genes,rownames(H3K4me3_8C)[which(H3K4me3_8C[,1] < LK4_cut)])

######## prediction heterogeneity for all promoters
hg_parameters <- read.table("../Data/NewtonParameters_udp.txt",header=F,row.names=1)

start <- hg_methyclass[c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene),selected_LPN]
names(start) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
total <- length(start)

e_hetergeneity_matrix <- c()
for (each_2C in hg_all_cells[hg_2C_6_id]){
	para_LPN2C <- hg_parameters[paste(selected_LPN,each_2C,sep="__"),]
	LPN2C_matrix <- MethylationTransMatrix(para_LPN2C[,1],para_LPN2C[,3],para_LPN2C[,2])
	for (each_4C in hg_all_cells[hg_4C_7_id]){
		para_2C4C <- hg_parameters[paste(each_2C,each_4C,sep="__"),]
		LPN4C_matrix <- TransMatrixMultiplication(para_2C4C[,1],para_2C4C[,3],para_2C4C[,2],LPN2C_matrix)
		for (each_8C in hg_all_cells[hg_8C_9_id]){
			para_4C8C <- hg_parameters[paste(each_4C,each_8C,sep="__"),]
			LPN8C_matrix <- TransMatrixMultiplication(para_4C8C[,1],para_4C8C[,3],para_4C8C[,2],LPN4C_matrix)
			t_cell_1 <- MethylationTransion(LPN2C_matrix,total,start,1.1);names(t_cell_1) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
			t_cell_2 <- MethylationTransion(LPN2C_matrix,total,start,2.1);names(t_cell_2) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
			f_cell_1 <- MethylationTransion(LPN4C_matrix,total,t_cell_1,1.2);names(f_cell_1) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
			f_cell_2 <- MethylationTransion(LPN4C_matrix,total,t_cell_1,2.2);names(f_cell_2) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
			f_cell_3 <- MethylationTransion(LPN4C_matrix,total,t_cell_2,3.2);names(f_cell_3) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
			f_cell_4 <- MethylationTransion(LPN4C_matrix,total,t_cell_2,4.2);names(f_cell_4) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
			e_cell_1 <- MethylationTransion(LPN8C_matrix,total,f_cell_1,1.3);names(e_cell_1) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
			e_cell_2 <- MethylationTransion(LPN8C_matrix,total,f_cell_1,2.3);names(e_cell_2) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
			e_cell_3 <- MethylationTransion(LPN8C_matrix,total,f_cell_2,3.3);names(e_cell_3) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
			e_cell_4 <- MethylationTransion(LPN8C_matrix,total,f_cell_2,4.3);names(e_cell_4) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
			e_cell_5 <- MethylationTransion(LPN8C_matrix,total,f_cell_3,5.3);names(e_cell_5) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
			e_cell_6 <- MethylationTransion(LPN8C_matrix,total,f_cell_3,6.3);names(e_cell_6) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
			e_cell_7 <- MethylationTransion(LPN8C_matrix,total,f_cell_4,7.3);names(e_cell_7) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
			e_cell_8 <- MethylationTransion(LPN8C_matrix,total,f_cell_4,8.3);names(e_cell_8) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
			e_hetergeneity <- apply(cbind(e_cell_1,e_cell_2,e_cell_3,e_cell_4,e_cell_5,e_cell_6,e_cell_7,e_cell_8),1,MethylationHeterogeneity)
			f_hetergeneity <- apply(cbind(f_cell_1,f_cell_2,f_cell_3,f_cell_4),1,MethylationHeterogeneity)
			t_hetergeneity <- apply(cbind(t_cell_1,t_cell_2),1,MethylationHeterogeneity)
			e_hetergeneity_matrix <- cbind(e_hetergeneity_matrix,e_hetergeneity)
		}
	}
}
rownames(e_hetergeneity_matrix) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)

RawParaMSE_CGratio <- c()
RawParaMSE_ATACratio <- c()
RawParaMSE_K4 <- c()
for (i in seq(ncol(e_hetergeneity_matrix))){
	RawParaMSE_CGratio <- c(RawParaMSE_CGratio,MeanSquareError(e_hetergeneity_matrix[c(HCGratio,MCGratio,LCGratio),i],hg_8C_9_obs[c(HCGratio,MCGratio,LCGratio)]))
	RawParaMSE_ATACratio <- c(RawParaMSE_ATACratio,MeanSquareError(e_hetergeneity_matrix[c(HATACratio,MATACratio,LATACratio),i],hg_8C_9_obs[c(HATACratio,MATACratio,LATACratio)]))
	RawParaMSE_K4 <- c(RawParaMSE_K4,MeanSquareError(e_hetergeneity_matrix[c(HK4,MK4,LK4),i],hg_8C_9_obs[c(HK4,MK4,LK4)]))
}

######## prediction heterogeneity for conditional promoters
predict_heter_CGratio <- c()
predict_heter_ATACratio <- c()
predict_heter_K4 <- c()

predict_heter_CGratio <- cbind(predict_heter_CGratio,ConditionalParametersPrediction(start,HCGratio,"HCGratio"))
predict_heter_CGratio <- cbind(predict_heter_CGratio,ConditionalParametersPrediction(start,MCGratio,"MCGratio"))
predict_heter_CGratio <- cbind(predict_heter_CGratio,ConditionalParametersPrediction(start,LCGratio,"LCGratio"))
predict_heter_ATACratio <- cbind(predict_heter_ATACratio,ConditionalParametersPrediction(start,HATACratio,"HATACratio"))
predict_heter_ATACratio <- cbind(predict_heter_ATACratio,ConditionalParametersPrediction(start,MATACratio,"MATACratio"))
predict_heter_ATACratio <- cbind(predict_heter_ATACratio,ConditionalParametersPrediction(start,LATACratio,"LATACratio"))
predict_heter_K4 <- cbind(predict_heter_K4,ConditionalParametersPrediction(start,HK4,"HK4"))
predict_heter_K4 <- cbind(predict_heter_K4,ConditionalParametersPrediction(start,MK4,"MK4"))
predict_heter_K4 <- cbind(predict_heter_K4,ConditionalParametersPrediction(start,LK4,"LK4"))

ConParaMSE_CGratio <- apply(predict_heter_CGratio,1,MeanSquareError,hg_8C_9_obs[c(HCGratio,MCGratio,LCGratio)])
ConParaMSE_ATACratio <- apply(predict_heter_ATACratio,1,MeanSquareError,hg_8C_9_obs[c(HATACratio,MATACratio,LATACratio)])
ConParaMSE_K4 <- apply(predict_heter_K4,1,MeanSquareError,hg_8C_9_obs[c(HK4,MK4,LK4)])

# MSE <- list("ConParaMSE_CGratio"=ConParaMSE_CGratio,"ConParaMSE_ATACratio"=ConParaMSE_ATACratio,"ConParaMSE_K4"=ConParaMSE_K4,"RawParaMSE_CGratio"=RawParaMSE_CGratio,"RawParaMSE_ATACratio"=RawParaMSE_ATACratio,"RawParaMSE_K4"=RawParaMSE_K4)
# save(MSE,file="MSE.Rdata")
load("MSE.Rdata")
RawParaMSE <- MSE[["RawParaMSE_CGratio"]]
ConParaMSE_CGratio <- MSE[["ConParaMSE_CGratio"]]
ConParaMSE_ATACratio <- MSE[["ConParaMSE_ATACratio"]]
ConParaMSE_K4 <- MSE[["ConParaMSE_K4"]]

pdf("Fig4D.ConditionalParameters.pdf",width=5.2,height=5)
par(mar=c(8,6,4,4))
HCGratio_1C2C_para <- read.table("HCGratio_1C2C_udp.txt",header=T);MCGratio_1C2C_para <- read.table("MCGratio_1C2C_udp.txt",header=T);LCGratio_1C2C_para <- read.table("LCGratio_1C2C_udp.txt",header=T)
HCGratio_2C4C_para <- read.table("HCGratio_2C4C_udp.txt",header=T);MCGratio_2C4C_para <- read.table("MCGratio_2C4C_udp.txt",header=T);LCGratio_2C4C_para <- read.table("LCGratio_2C4C_udp.txt",header=T)
HCGratio_4C8C_para <- read.table("HCGratio_4C8C_udp.txt",header=T);MCGratio_4C8C_para <- read.table("MCGratio_4C8C_udp.txt",header=T);LCGratio_4C8C_para <- read.table("LCGratio_4C8C_udp.txt",header=T)
plot(1,type="n",xaxt="n",bty="o",xlab="",ylab="Parameter value",ylim=c(0,1),xlim=c(1,3),main="Promoters stratified by CpG ratio");box(lwd=2)
n <- rbind(apply(HCGratio_1C2C_para,2,length),apply(HCGratio_2C4C_para,2,length),apply(HCGratio_4C8C_para,2,length),
           apply(MCGratio_1C2C_para,2,length),apply(MCGratio_2C4C_para,2,length),apply(MCGratio_4C8C_para,2,length),
           apply(LCGratio_1C2C_para,2,length),apply(LCGratio_2C4C_para,2,length),apply(LCGratio_4C8C_para,2,length))
mean <- rbind(apply(HCGratio_1C2C_para,2,mean),apply(HCGratio_2C4C_para,2,mean),apply(HCGratio_4C8C_para,2,mean),
              apply(MCGratio_1C2C_para,2,mean),apply(MCGratio_2C4C_para,2,mean),apply(MCGratio_4C8C_para,2,mean),
              apply(LCGratio_1C2C_para,2,mean),apply(LCGratio_2C4C_para,2,mean),apply(LCGratio_4C8C_para,2,mean))
sd <- rbind(apply(HCGratio_1C2C_para,2,sd),apply(HCGratio_2C4C_para,2,sd),apply(HCGratio_4C8C_para,2,sd),
            apply(MCGratio_1C2C_para,2,sd),apply(MCGratio_2C4C_para,2,sd),apply(MCGratio_4C8C_para,2,sd),
            apply(LCGratio_1C2C_para,2,sd),apply(LCGratio_2C4C_para,2,sd),apply(LCGratio_4C8C_para,2,sd))
points(mean[c(1,2,3),1],type="b",col=cccol[1],pch=19,lty=1,lwd=2); v1=mean[c(1,2,3),1]; v2=v1-1.96*sd[c(1,2,3),1]/sqrt(n[c(1,2,3),1]); v3=v1+1.96*sd[c(1,2,3),1]/sqrt(n[c(1,2,3),1]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(1,2,3),2],type="b",col=cccol[2],pch=19,lty=1,lwd=2); v1=mean[c(1,2,3),2]; v2=v1-1.96*sd[c(1,2,3),2]/sqrt(n[c(1,2,3),2]); v3=v1+1.96*sd[c(1,2,3),2]/sqrt(n[c(1,2,3),2]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(1,2,3),3],type="b",col=cccol[3],pch=19,lty=1,lwd=2); v1=mean[c(1,2,3),3]; v2=v1-1.96*sd[c(1,2,3),3]/sqrt(n[c(1,2,3),3]); v3=v1+1.96*sd[c(1,2,3),3]/sqrt(n[c(1,2,3),3]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(4,5,6),1],type="b",col=cccol[1],pch=15,lty=2,lwd=2); v1=mean[c(4,5,6),1]; v2=v1-1.96*sd[c(4,5,6),1]/sqrt(n[c(4,5,6),1]); v3=v1+1.96*sd[c(4,5,6),1]/sqrt(n[c(4,5,6),1]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(4,5,6),2],type="b",col=cccol[2],pch=15,lty=2,lwd=2); v1=mean[c(4,5,6),2]; v2=v1-1.96*sd[c(4,5,6),2]/sqrt(n[c(4,5,6),2]); v3=v1+1.96*sd[c(4,5,6),2]/sqrt(n[c(4,5,6),2]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(4,5,6),3],type="b",col=cccol[3],pch=15,lty=2,lwd=2); v1=mean[c(4,5,6),3]; v2=v1-1.96*sd[c(4,5,6),3]/sqrt(n[c(4,5,6),3]); v3=v1+1.96*sd[c(4,5,6),3]/sqrt(n[c(4,5,6),3]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(7,8,9),1],type="b",col=cccol[1],pch=17,lty=3,lwd=2); v1=mean[c(7,8,9),1]; v2=v1-1.96*sd[c(7,8,9),1]/sqrt(n[c(7,8,9),1]); v3=v1+1.96*sd[c(7,8,9),1]/sqrt(n[c(7,8,9),1]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(7,8,9),2],type="b",col=cccol[2],pch=17,lty=3,lwd=2); v1=mean[c(7,8,9),2]; v2=v1-1.96*sd[c(7,8,9),2]/sqrt(n[c(7,8,9),2]); v3=v1+1.96*sd[c(7,8,9),2]/sqrt(n[c(7,8,9),2]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(7,8,9),3],type="b",col=cccol[3],pch=17,lty=3,lwd=2); v1=mean[c(7,8,9),3]; v2=v1-1.96*sd[c(7,8,9),3]/sqrt(n[c(7,8,9),3]); v3=v1+1.96*sd[c(7,8,9),3]/sqrt(n[c(7,8,9),3]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
dev.off()

pdf("FigS4DE.ConditionalParameters.ATAC_K4.pdf",width=5,height=5)
par(mar=c(8,6,4,4))
HATACratio_1C2C_para <- read.table("HATACratio_1C2C_udp.txt",header=T);MATACratio_1C2C_para <- read.table("MATACratio_1C2C_udp.txt",header=T);LATACratio_1C2C_para <- read.table("LATACratio_1C2C_udp.txt",header=T)
HATACratio_2C4C_para <- read.table("HATACratio_2C4C_udp.txt",header=T);MATACratio_2C4C_para <- read.table("MATACratio_2C4C_udp.txt",header=T);LATACratio_2C4C_para <- read.table("LATACratio_2C4C_udp.txt",header=T)
HATACratio_4C8C_para <- read.table("HATACratio_4C8C_udp.txt",header=T);MATACratio_4C8C_para <- read.table("MATACratio_4C8C_udp.txt",header=T);LATACratio_4C8C_para <- read.table("LATACratio_4C8C_udp.txt",header=T)
HATACratio_mean <- cbind(apply(HATACratio_1C2C_para,2,mean),apply(HATACratio_2C4C_para,2,mean),apply(HATACratio_4C8C_para,2,mean))
MATACratio_mean <- cbind(apply(MATACratio_1C2C_para,2,mean),apply(MATACratio_2C4C_para,2,mean),apply(MATACratio_4C8C_para,2,mean))
LATACratio_mean <- cbind(apply(LATACratio_1C2C_para,2,mean),apply(LATACratio_2C4C_para,2,mean),apply(LATACratio_4C8C_para,2,mean))
plot(1,type="n",xaxt="n",bty="o",xlab="",ylab="Parameter value",ylim=c(0,1),xlim=c(1,3),main="Promoters stratified by chromatin accessibility");box(lwd=2)
n <- rbind(apply(HATACratio_1C2C_para,2,length),apply(HATACratio_2C4C_para,2,length),apply(HATACratio_4C8C_para,2,length),
           apply(MATACratio_1C2C_para,2,length),apply(MATACratio_2C4C_para,2,length),apply(MATACratio_4C8C_para,2,length),
           apply(LATACratio_1C2C_para,2,length),apply(LATACratio_2C4C_para,2,length),apply(LATACratio_4C8C_para,2,length))
mean <- rbind(apply(HATACratio_1C2C_para,2,mean),apply(HATACratio_2C4C_para,2,mean),apply(HATACratio_4C8C_para,2,mean),
              apply(MATACratio_1C2C_para,2,mean),apply(MATACratio_2C4C_para,2,mean),apply(MATACratio_4C8C_para,2,mean),
              apply(LATACratio_1C2C_para,2,mean),apply(LATACratio_2C4C_para,2,mean),apply(LATACratio_4C8C_para,2,mean))
sd <- rbind(apply(HATACratio_1C2C_para,2,sd),apply(HATACratio_2C4C_para,2,sd),apply(HATACratio_4C8C_para,2,sd),
            apply(MATACratio_1C2C_para,2,sd),apply(MATACratio_2C4C_para,2,sd),apply(MATACratio_4C8C_para,2,sd),
            apply(LATACratio_1C2C_para,2,sd),apply(LATACratio_2C4C_para,2,sd),apply(LATACratio_4C8C_para,2,sd))
points(mean[c(1,2,3),1],type="b",col=cccol[1],pch=19,lty=1,lwd=2); v1=mean[c(1,2,3),1]; v2=v1-1.96*sd[c(1,2,3),1]/sqrt(n[c(1,2,3),1]); v3=v1+1.96*sd[c(1,2,3),1]/sqrt(n[c(1,2,3),1]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(1,2,3),2],type="b",col=cccol[2],pch=19,lty=1,lwd=2); v1=mean[c(1,2,3),2]; v2=v1-1.96*sd[c(1,2,3),2]/sqrt(n[c(1,2,3),2]); v3=v1+1.96*sd[c(1,2,3),2]/sqrt(n[c(1,2,3),2]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(1,2,3),3],type="b",col=cccol[3],pch=19,lty=1,lwd=2); v1=mean[c(1,2,3),3]; v2=v1-1.96*sd[c(1,2,3),3]/sqrt(n[c(1,2,3),3]); v3=v1+1.96*sd[c(1,2,3),3]/sqrt(n[c(1,2,3),3]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(4,5,6),1],type="b",col=cccol[1],pch=15,lty=2,lwd=2); v1=mean[c(4,5,6),1]; v2=v1-1.96*sd[c(4,5,6),1]/sqrt(n[c(4,5,6),1]); v3=v1+1.96*sd[c(4,5,6),1]/sqrt(n[c(4,5,6),1]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(4,5,6),2],type="b",col=cccol[2],pch=15,lty=2,lwd=2); v1=mean[c(4,5,6),2]; v2=v1-1.96*sd[c(4,5,6),2]/sqrt(n[c(4,5,6),2]); v3=v1+1.96*sd[c(4,5,6),2]/sqrt(n[c(4,5,6),2]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(4,5,6),3],type="b",col=cccol[3],pch=15,lty=2,lwd=2); v1=mean[c(4,5,6),3]; v2=v1-1.96*sd[c(4,5,6),3]/sqrt(n[c(4,5,6),3]); v3=v1+1.96*sd[c(4,5,6),3]/sqrt(n[c(4,5,6),3]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(7,8,9),1],type="b",col=cccol[1],pch=17,lty=3,lwd=2); v1=mean[c(7,8,9),1]; v2=v1-1.96*sd[c(7,8,9),1]/sqrt(n[c(7,8,9),1]); v3=v1+1.96*sd[c(7,8,9),1]/sqrt(n[c(7,8,9),1]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(7,8,9),2],type="b",col=cccol[2],pch=17,lty=3,lwd=2); v1=mean[c(7,8,9),2]; v2=v1-1.96*sd[c(7,8,9),2]/sqrt(n[c(7,8,9),2]); v3=v1+1.96*sd[c(7,8,9),2]/sqrt(n[c(7,8,9),2]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(7,8,9),3],type="b",col=cccol[3],pch=17,lty=3,lwd=2); v1=mean[c(7,8,9),3]; v2=v1-1.96*sd[c(7,8,9),3]/sqrt(n[c(7,8,9),3]); v3=v1+1.96*sd[c(7,8,9),3]/sqrt(n[c(7,8,9),3]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)

HK4_1C2C_para <- read.table("HK4_1C2C_udp.txt",header=T);MK4_1C2C_para <- read.table("MK4_1C2C_udp.txt",header=T);LK4_1C2C_para <- read.table("LK4_1C2C_udp.txt",header=T)
HK4_2C4C_para <- read.table("HK4_2C4C_udp.txt",header=T);MK4_2C4C_para <- read.table("MK4_2C4C_udp.txt",header=T);LK4_2C4C_para <- read.table("LK4_2C4C_udp.txt",header=T)
HK4_4C8C_para <- read.table("HK4_4C8C_udp.txt",header=T);MK4_4C8C_para <- read.table("MK4_4C8C_udp.txt",header=T);LK4_4C8C_para <- read.table("LK4_4C8C_udp.txt",header=T)
HK4_mean <- cbind(apply(HK4_1C2C_para,2,mean),apply(HK4_2C4C_para,2,mean),apply(HK4_4C8C_para,2,mean))
MK4_mean <- cbind(apply(MK4_1C2C_para,2,mean),apply(MK4_2C4C_para,2,mean),apply(MK4_4C8C_para,2,mean))
LK4_mean <- cbind(apply(LK4_1C2C_para,2,mean),apply(LK4_2C4C_para,2,mean),apply(LK4_4C8C_para,2,mean))
plot(1,type="n",xaxt="n",bty="o",xlab="",ylab="Parameter value",ylim=c(0,1),xlim=c(1,3),main="Promoters stratified by H3K4me3 signal");box(lwd=2)
n <- rbind(apply(HK4_1C2C_para,2,length),apply(HK4_2C4C_para,2,length),apply(HK4_4C8C_para,2,length),
           apply(MK4_1C2C_para,2,length),apply(MK4_2C4C_para,2,length),apply(MK4_4C8C_para,2,length),
           apply(LK4_1C2C_para,2,length),apply(LK4_2C4C_para,2,length),apply(LK4_4C8C_para,2,length))
mean <- rbind(apply(HK4_1C2C_para,2,mean),apply(HK4_2C4C_para,2,mean),apply(HK4_4C8C_para,2,mean),
              apply(MK4_1C2C_para,2,mean),apply(MK4_2C4C_para,2,mean),apply(MK4_4C8C_para,2,mean),
              apply(LK4_1C2C_para,2,mean),apply(LK4_2C4C_para,2,mean),apply(LK4_4C8C_para,2,mean))
sd <- rbind(apply(HK4_1C2C_para,2,sd),apply(HK4_2C4C_para,2,sd),apply(HK4_4C8C_para,2,sd),
            apply(MK4_1C2C_para,2,sd),apply(MK4_2C4C_para,2,sd),apply(MK4_4C8C_para,2,sd),
            apply(LK4_1C2C_para,2,sd),apply(LK4_2C4C_para,2,sd),apply(LK4_4C8C_para,2,sd))
points(mean[c(1,2,3),1],type="b",col=cccol[1],pch=19,lty=1,lwd=2); v1=mean[c(1,2,3),1]; v2=v1-1.96*sd[c(1,2,3),1]/sqrt(n[c(1,2,3),1]); v3=v1+1.96*sd[c(1,2,3),1]/sqrt(n[c(1,2,3),1]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(1,2,3),2],type="b",col=cccol[2],pch=19,lty=1,lwd=2); v1=mean[c(1,2,3),2]; v2=v1-1.96*sd[c(1,2,3),2]/sqrt(n[c(1,2,3),2]); v3=v1+1.96*sd[c(1,2,3),2]/sqrt(n[c(1,2,3),2]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(1,2,3),3],type="b",col=cccol[3],pch=19,lty=1,lwd=2); v1=mean[c(1,2,3),3]; v2=v1-1.96*sd[c(1,2,3),3]/sqrt(n[c(1,2,3),3]); v3=v1+1.96*sd[c(1,2,3),3]/sqrt(n[c(1,2,3),3]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(4,5,6),1],type="b",col=cccol[1],pch=15,lty=2,lwd=2); v1=mean[c(4,5,6),1]; v2=v1-1.96*sd[c(4,5,6),1]/sqrt(n[c(4,5,6),1]); v3=v1+1.96*sd[c(4,5,6),1]/sqrt(n[c(4,5,6),1]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(4,5,6),2],type="b",col=cccol[2],pch=15,lty=2,lwd=2); v1=mean[c(4,5,6),2]; v2=v1-1.96*sd[c(4,5,6),2]/sqrt(n[c(4,5,6),2]); v3=v1+1.96*sd[c(4,5,6),2]/sqrt(n[c(4,5,6),2]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(4,5,6),3],type="b",col=cccol[3],pch=15,lty=2,lwd=2); v1=mean[c(4,5,6),3]; v2=v1-1.96*sd[c(4,5,6),3]/sqrt(n[c(4,5,6),3]); v3=v1+1.96*sd[c(4,5,6),3]/sqrt(n[c(4,5,6),3]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(7,8,9),1],type="b",col=cccol[1],pch=17,lty=3,lwd=2); v1=mean[c(7,8,9),1]; v2=v1-1.96*sd[c(7,8,9),1]/sqrt(n[c(7,8,9),1]); v3=v1+1.96*sd[c(7,8,9),1]/sqrt(n[c(7,8,9),1]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(7,8,9),2],type="b",col=cccol[2],pch=17,lty=3,lwd=2); v1=mean[c(7,8,9),2]; v2=v1-1.96*sd[c(7,8,9),2]/sqrt(n[c(7,8,9),2]); v3=v1+1.96*sd[c(7,8,9),2]/sqrt(n[c(7,8,9),2]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
points(mean[c(7,8,9),3],type="b",col=cccol[3],pch=17,lty=3,lwd=2); v1=mean[c(7,8,9),3]; v2=v1-1.96*sd[c(7,8,9),3]/sqrt(n[c(7,8,9),3]); v3=v1+1.96*sd[c(7,8,9),3]/sqrt(n[c(7,8,9),3]); polygon(c(1,1:length(v1),length(v1):2),c(v2[1],v3,v2[length(v1):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
dev.off()


CoffCols <- colorRampPalette(c(cccol[1],cccol[9],cccol[seq(2,8)]), bias=1)(7)
CoffCols50 <- paste(CoffCols,"50",sep="")
pdf("Fig4E.ConditionalPredictionMSE.pdf",width=4.3,height=5)
MSE_CGratio <- list(RawParaMSE,ConParaMSE_CGratio,ConParaMSE_ATACratio,ConParaMSE_K4)
par(mar=c(8,6,4,4))
n <- length(MSE_CGratio)
# xpos <- 0:(n-1)+1.5
# p_value <- c()
# for (i in seq(n-1)){
# 	p_value <- tryCatch(c(p_value,t.test(MSE_CGratio[[1]],MSE_CGratio[[i+1]])$p.value),error=function(a) c(p_value,-1))
# 	p_value <- c(p_value,-1)
# }
# mark <- symnum(p_value, cutpoints=c(-1,-0.1,0.001,0.01,0.05,1), symbols=c(NA,"***","**","*","-"))
bp <- boxplot(MSE_CGratio,at=0:(n-1)+1, xlim=c(0.5,n+0.5),outline=F,plot=F)
ylim <- range(bp$stats,na.rm=T)
dist <- (ylim[2]-ylim[1])/20
ylim[2] <- ylim[2]+3*dist
boxplot(MSE_CGratio,at=0:(n-1)+1,boxwex=0.6, xlim=c(0.5,n+0.5),ylim=ylim,col=c(cccol50[8],CoffCols50[c(2,4,3)]),border=c(cccol[8],CoffCols[c(2,4,3)]),outline=F,names=c("Original", "CpG ratio adjusted", "Chromatin accessibility adjusted", "H3K4me3 signal adjusted"),main="",lwd=2,ylab="Mean squared error",las=2)
box(lwd=2)
# ypos_1 <- bp$stats[5,][1:n-1]
# ypos_2 <- bp$stats[5,][2:n]
# legend("topleft",c("***:p<0.001","**:p<0.01","*:p<0.05"),text.col=cccol[1],bty="n")
# for(i in 1:length(mark)){
# 	if(!is.na(mark[i])){
# 		segments(xpos[i]-.45, ypos_1[i]+dist/2, xpos[i]-.45, max(ypos_1[i], ypos_2[i])+dist)
# 		segments(xpos[i]+.45, ypos_2[i]+dist/2, xpos[i]+.45, max(ypos_1[i], ypos_2[i])+dist)
# 		segments(xpos[i]-.45, max(ypos_1[i], ypos_2[i])+dist, xpos[i]-0.3, max(ypos_1[i], ypos_2[i])+dist)
# 		segments(xpos[i]+.45, max(ypos_1[i], ypos_2[i])+dist, xpos[i]+0.3, max(ypos_1[i], ypos_2[i])+dist)
# 		text(x=xpos[i], y=max(ypos_1[i], ypos_2[i])+dist, label=mark[i], col="black")
# 	}
# }
dev.off()
