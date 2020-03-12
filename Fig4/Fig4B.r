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
hg_4C_7_c1_MHG <- tryCatch(as.vector(read.table("hg_4C_7_c1_MHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c2_MHG <- tryCatch(as.vector(read.table("hg_4C_7_c2_MHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c3_MHG <- tryCatch(as.vector(read.table("hg_4C_7_c3_MHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c4_MHG <- tryCatch(as.vector(read.table("hg_4C_7_c4_MHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c5_MHG <- tryCatch(as.vector(read.table("hg_4C_7_c5_MHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c1_MHG <- tryCatch(as.vector(read.table("hg_8C_9_c1_MHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c2_MHG <- tryCatch(as.vector(read.table("hg_8C_9_c2_MHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c3_MHG <- tryCatch(as.vector(read.table("hg_8C_9_c3_MHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c4_MHG <- tryCatch(as.vector(read.table("hg_8C_9_c4_MHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c5_MHG <- tryCatch(as.vector(read.table("hg_8C_9_c5_MHG.txt")[,1]),error = function(a) NA)

hg_4C_7_c1_LHG <- tryCatch(as.vector(read.table("hg_4C_7_c1_LHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c2_LHG <- tryCatch(as.vector(read.table("hg_4C_7_c2_LHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c3_LHG <- tryCatch(as.vector(read.table("hg_4C_7_c3_LHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c4_LHG <- tryCatch(as.vector(read.table("hg_4C_7_c4_LHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c5_LHG <- tryCatch(as.vector(read.table("hg_4C_7_c5_LHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c1_LHG <- tryCatch(as.vector(read.table("hg_8C_9_c1_LHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c2_LHG <- tryCatch(as.vector(read.table("hg_8C_9_c2_LHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c3_LHG <- tryCatch(as.vector(read.table("hg_8C_9_c3_LHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c4_LHG <- tryCatch(as.vector(read.table("hg_8C_9_c4_LHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c5_LHG <- tryCatch(as.vector(read.table("hg_8C_9_c5_LHG.txt")[,1]),error = function(a) NA)

hg_4C_7_c1_HHG <- tryCatch(as.vector(read.table("hg_4C_7_c1_HHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c2_HHG <- tryCatch(as.vector(read.table("hg_4C_7_c2_HHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c3_HHG <- tryCatch(as.vector(read.table("hg_4C_7_c3_HHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c4_HHG <- tryCatch(as.vector(read.table("hg_4C_7_c4_HHG.txt")[,1]),error = function(a) NA)
hg_4C_7_c5_HHG <- tryCatch(as.vector(read.table("hg_4C_7_c5_HHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c1_HHG <- tryCatch(as.vector(read.table("hg_8C_9_c1_HHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c2_HHG <- tryCatch(as.vector(read.table("hg_8C_9_c2_HHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c3_HHG <- tryCatch(as.vector(read.table("hg_8C_9_c3_HHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c4_HHG <- tryCatch(as.vector(read.table("hg_8C_9_c4_HHG.txt")[,1]),error = function(a) NA)
hg_8C_9_c5_HHG <- tryCatch(as.vector(read.table("hg_8C_9_c5_HHG.txt")[,1]),error = function(a) NA)

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
#########################     prepare input for RandomForest      #################################
###################################################################################################
# install.packages("glmnet",dependencies=TRUE,repos="https://mirrors.tongji.edu.cn/CRAN/")
library(glmnet)

# raw heterogeneity
hg_2C_6_obs <- apply(na.omit(hg_methyclass[,hg_2C_6_id]),1,MethylationHeterogeneityObservation)
hg_4C_7_obs <- apply(na.omit(hg_methyclass[,hg_4C_7_id]),1,MethylationHeterogeneityObservation)
hg_4C_8_obs <- apply(na.omit(hg_methyclass[,hg_4C_8_id]),1,MethylationHeterogeneityObservation)
hg_8C_6_obs <- apply(na.omit(hg_methyclass[,hg_8C_6_id]),1,MethylationHeterogeneityObservation)
hg_8C_9_obs <- apply(na.omit(hg_methyclass[,hg_8C_9_id]),1,MethylationHeterogeneityObservation)

hg_4C_Coff <- c()
genes <- intersect(intersect(names(hg_4C_7_obs),rownames(hg19_GeneAgeClass)),c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene))
HeteroLabels <- rep(0,length(genes));names(HeteroLabels) <- genes
HeteroLabels[intersect(genes,c(hg_4C_7_c1_LHG,hg_4C_7_c2_LHG,hg_4C_7_c3_LHG,hg_4C_7_c4_LHG,hg_4C_7_c5_LHG))] <- 1
HeteroLabels[intersect(genes,c(hg_4C_7_c1_MHG,hg_4C_7_c2_MHG,hg_4C_7_c3_MHG,hg_4C_7_c4_MHG,hg_4C_7_c5_MHG))] <- 2
HeteroLabels[intersect(genes,c(hg_4C_7_c1_HHG,hg_4C_7_c2_HHG,hg_4C_7_c3_HHG,hg_4C_7_c4_HHG,hg_4C_7_c5_HHG))] <- 3
write.table(cbind(HeteroLabels[genes]),file="Fig4C.hg_4C_7_RandomForestLabel.txt",sep="\t",quote=F,col.names=F,row.names=T) # label_4C
feature_4C <- cbind(CG_ratio[genes,],ATAC_ratio_4C_1[genes,],H3K27me3_4C[genes,],H3K4me3_4C[genes,])
colnames(feature_4C) <- c("CGRatio","ATAC","H3K27me3","H3K4me3");rownames(feature_4C) <- genes
write.table(feature_4C,file="Fig4C.hg_4C_7_RandomForestFeature.txt",sep="\t",quote=F,col.names=T,row.names=T) # feature_4C
# cv.glmmod <- cv.glmnet(x=feature_4C,y=hg_4C_7_obs[genes],alpha=1,family="gaussian")
# write.table(coef(cv.glmmod,s="lambda.min")[,1],file="Fig4C.hg_4C_7_LassoCofficient.txt",sep="\t",quote=F,col.names=F,row.names=T)
# hg_4C_Coff <- rbind(hg_4C_Coff,coef(cv.glmmod,s="lambda.min")[2:5,1])


hg_8C_Coff <- c()
genes <- intersect(intersect(names(hg_8C_9_obs),rownames(hg19_GeneAgeClass)),c(c1_gene,c2_gene,c3_gene,c4_gene,c5_gene))
HeteroLabels <- rep(0,length(genes));names(HeteroLabels) <- genes
HeteroLabels[intersect(genes,c(hg_8C_9_c1_LHG,hg_8C_9_c2_LHG,hg_8C_9_c3_LHG,hg_8C_9_c4_LHG,hg_8C_9_c5_LHG))] <- 1
HeteroLabels[intersect(genes,c(hg_8C_9_c1_MHG,hg_8C_9_c2_MHG,hg_8C_9_c3_MHG,hg_8C_9_c4_MHG,hg_8C_9_c5_MHG))] <- 2
HeteroLabels[intersect(genes,c(hg_8C_9_c1_HHG,hg_8C_9_c2_HHG,hg_8C_9_c3_HHG,hg_8C_9_c4_HHG,hg_8C_9_c5_HHG))] <- 3
write.table(cbind(HeteroLabels[genes]),file="Fig4C.hg_8C_9_RandomForestLabel.txt",sep="\t",quote=F,col.names=F,row.names=T) # label_8C
feature_8C <- cbind(CG_ratio[genes,],ATAC_ratio_8C_1[genes,],H3K27me3_8C[genes,],H3K4me3_8C[genes,])
colnames(feature_8C) <- c("CGRatio","ATAC","H3K27me3","H3K4me3");rownames(feature_8C) <- genes
write.table(feature_8C,file="Fig4C.hg_8C_9_RandomForestFeature.txt",sep="\t",quote=F,col.names=T,row.names=T) # feature_8C
# cv.glmmod <- cv.glmnet(x=feature_8C,y=hg_8C_9_obs[genes],alpha=1,family="gaussian")
# write.table(coef(cv.glmmod,s="lambda.min")[,1],file="Fig4C.hg_8C_9_LassoCofficient.txt",sep="\t",quote=F,col.names=F,row.names=T)
# hg_8C_Coff <- rbind(hg_8C_Coff,coef(cv.glmmod,s="lambda.min")[2:5,1])


CoffCols <- colorRampPalette(c(cccol[1],cccol[9],cccol[seq(2,8)]), bias=1)(7)
# pdf("Fig4C.LassoCofficient.pdf",width=4,height=6)
# par(mar=c(10,6,4,4))
# barplot(hg_4C_Coff[1,],main="4-cell",ylab="Lasso cofficient",names.arg=c("CGRatio","ATAC","H3K27me3","H3K4me3"),bty="l",col=CoffCols[2:5],border=CoffCols[2:5],ylim=c(-0.35,0.35),xpd = FALSE,las=2);box(lwd=2)
# barplot(hg_8C_Coff[1,],main="8-cell",ylab="Lasso cofficient",names.arg=c("CGRatio","ATAC","H3K27me3","H3K4me3"),bty="l",col=CoffCols[2:5],border=CoffCols[2:5],ylim=c(-0.55,0.35),xpd = FALSE,las=2);box(lwd=2)
# dev.off()


# random forest
hg_4C_Imt <- read.table("../Data/Fig4B.hg_4C_7_RandomForestImportance.txt")[,1];names(hg_4C_Imt) <- c("CGRatio","ATAC","H3K27me3","H3K4me3")
hg_8C_Imt <- read.table("../Data/Fig4B.hg_8C_9_RandomForestImportance.txt")[,1];names(hg_8C_Imt) <- c("CGRatio","ATAC","H3K27me3","H3K4me3")
pdf("Fig4B.RandomForestImportance.pdf",width=3,height=5)
par(mar=c(8,6,4,4))
barplot(sort(hg_4C_Imt,decreasing=T),main="4-cell",ylab="RandomForest importance",bty="l",col=CoffCols[2:5],border=CoffCols[2:5],ylim=c(0,0.6),xpd = FALSE,las=2);box(lwd=2)
barplot(sort(hg_8C_Imt,decreasing=T),main="8-cell",ylab="RandomForest importance",bty="l",col=CoffCols[2:5],border=CoffCols[2:5],ylim=c(0,0.6),xpd = FALSE,las=2);box(lwd=2)
dev.off()
