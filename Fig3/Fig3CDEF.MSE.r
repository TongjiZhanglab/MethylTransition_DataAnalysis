###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)

MeanSquareError <- function(predict,observe){
	return(mean((predict-observe)^2))
}

cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001330","#16557A30","#C7A60930","#87C23230","#64C0AB30","#A14C9430","#15A08C30","#8B7E7530","#1E7CAF30","#EA425F30","#46489A30","#E5003330","#0F231F30","#1187CD30")
cccolBR <- c("#16557A","#443F60","#712A46","#A0152C","#CE0013")
cccolBR30 <- c("#16557A30","#443F6030","#712A4630","#A0152C30","#CE001330")
library(ineq)

MethylationClassSummary <- function(x){
	a1 <- length(which(x==0))
	a2 <- length(which(x==1/4))
	a3 <- length(which(x==1/2))
	a4 <- length(which(x==3/4))
	a5 <- length(which(x==1))
	total <- a1+a2+a3+a4+a5
	return(c(a1/total,a2/total,a3/total,a4/total,a5/total))
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

MethylationHeterogeneity <- function(x){
	library(ineq)
	a1 <- length(which(x<=1/8))
	a2 <- length(which(x>1/8 & x<=3/8))
	a3 <- length(which(x>3/8 & x<=5/8))
	a4 <- length(which(x>5/8 & x<=7/8))
	a5 <- length(which(x>7/8))
	total <- a1+a2+a3+a4+a5
	if (total >=2){
		return(1-Gini(c(a1/total,a2/total,a3/total,a4/total,a5/total),corr=T))
	}else{
		return(NA)
	}
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
######## prediction parameter
hg_parameters <- read.table("../Data/NewtonParameters_udp.txt",header=F,row.names=1)

selected_LPN <- "hg_LPN_8"
c1_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="0")]
c2_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="0.25")]
c3_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="0.5")]
c4_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="0.75")]
c5_gene <- hg_genes[which(hg_methyclass[,selected_LPN]=="1")]
start <- hg_methyclass[c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene),selected_LPN]
names(start) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
total <- length(start)

e_para_mat <- c();e_name <- c()
f_para_mat <- c();f_name <- c()
t_para_mat <- c();t_name <- c()
tmp_LPN2C <- c();tmp_2C4C <- c();tmp_4C8C <- c()

for (each_2C in hg_all_cells[hg_2C_6_id]){
	para_LPN2C <- hg_parameters[paste(selected_LPN,each_2C,sep="__"),]
	LPN2C_matrix <- MethylationTransMatrix(para_LPN2C[,1],para_LPN2C[,3],para_LPN2C[,2])
	tmp_LPN2C <- cbind(tmp_LPN2C,as.vector(MethylationTransMatrix(para_LPN2C[,1],para_LPN2C[,3],para_LPN2C[,2])))
	t_para_mat <- cbind(t_para_mat,as.vector(LPN2C_matrix))
	t_name <- c(t_name,each_2C)
	for (each_4C in hg_all_cells[hg_4C_7_id]){
		para_2C4C <- hg_parameters[paste(each_2C,each_4C,sep="__"),]
		tmp_2C4C <- cbind(tmp_2C4C,as.vector(MethylationTransMatrix(para_2C4C[,1],para_2C4C[,3],para_2C4C[,2])))
		LPN4C_matrix <- TransMatrixMultiplication(para_2C4C[,1],para_2C4C[,3],para_2C4C[,2],LPN2C_matrix)
		f_para_mat <- cbind(f_para_mat,as.vector(LPN4C_matrix))
		f_name <- c(f_name,paste(each_2C,each_4C,sep="__"))
		for (each_8C in hg_all_cells[hg_8C_9_id]){
			para_4C8C <- hg_parameters[paste(each_4C,each_8C,sep="__"),]
			tmp_4C8C <- cbind(tmp_4C8C,as.vector(MethylationTransMatrix(para_4C8C[,1],para_4C8C[,3],para_4C8C[,2])))
			LPN8C_matrix <- TransMatrixMultiplication(para_4C8C[,1],para_4C8C[,3],para_4C8C[,2],LPN4C_matrix)
			e_para_mat <- cbind(e_para_mat,as.vector(LPN8C_matrix))
			e_name <- c(e_name,paste(each_2C,each_4C,each_8C,sep="__"))
		}
	}
}
t_cell_para <- matrix(apply(tmp_LPN2C,1,mean),c(5,5))
f_cell_para <- matrix(apply(tmp_2C4C,1,mean),c(5,5))
e_cell_para <- matrix(apply(tmp_4C8C,1,mean),c(5,5))

t_cell_1 <- MethylationTransion(t_cell_para,total,start,1.1)
t_cell_2 <- MethylationTransion(t_cell_para,total,start,2.1)
f_cell_1 <- MethylationTransion(f_cell_para,total,t_cell_1,1.2)
f_cell_2 <- MethylationTransion(f_cell_para,total,t_cell_1,2.2)
f_cell_3 <- MethylationTransion(f_cell_para,total,t_cell_2,3.2)
f_cell_4 <- MethylationTransion(f_cell_para,total,t_cell_2,4.2)
e_cell_1 <- MethylationTransion(e_cell_para,total,f_cell_1,1.3)
e_cell_2 <- MethylationTransion(e_cell_para,total,f_cell_1,2.3)
e_cell_3 <- MethylationTransion(e_cell_para,total,f_cell_2,3.3)
e_cell_4 <- MethylationTransion(e_cell_para,total,f_cell_2,4.3)
e_cell_5 <- MethylationTransion(e_cell_para,total,f_cell_3,5.3)
e_cell_6 <- MethylationTransion(e_cell_para,total,f_cell_3,6.3)
e_cell_7 <- MethylationTransion(e_cell_para,total,f_cell_4,7.3)
e_cell_8 <- MethylationTransion(e_cell_para,total,f_cell_4,8.3)
e_hetergeneity <- apply(cbind(e_cell_1,e_cell_2,e_cell_3,e_cell_4,e_cell_5,e_cell_6,e_cell_7,e_cell_8),1,MethylationHeterogeneity)
f_hetergeneity <- apply(cbind(f_cell_1,f_cell_2,f_cell_3,f_cell_4),1,MethylationHeterogeneity)
t_hetergeneity <- apply(cbind(t_cell_1,t_cell_2),1,MethylationHeterogeneity)

names(e_hetergeneity) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
names(f_hetergeneity) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)
names(t_hetergeneity) <- c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene)


hg_4C_7_obs <- apply(na.omit(hg_methyclass[,hg_4C_7_id]),1,MethylationHeterogeneityObservation)
hg_8C_9_obs <- apply(na.omit(hg_methyclass[,hg_8C_9_id]),1,MethylationHeterogeneityObservation)


f_gene <- intersect(names(f_hetergeneity),names(hg_4C_7_obs))
f_mse <- MeanSquareError(f_hetergeneity[f_gene],hg_4C_7_obs[f_gene])
f_mse
# [1] 0.01816356
e_gene <- intersect(names(e_hetergeneity),names(hg_8C_9_obs))
e_mse <- MeanSquareError(e_hetergeneity[e_gene],hg_8C_9_obs[e_gene])
e_mse
# [1] 0.03868387
###################################################################################################
################################                PLOT               ################################
###################################################################################################
