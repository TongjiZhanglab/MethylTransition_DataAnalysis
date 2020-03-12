###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# hg_methyclass <- read.table("../Data/hg_MethylationClass.txt",header=T,row.names=1)
# selected_LPN <- "hg_LPN_8"
# MethylationClassSummary <- function(x){
# 	a1 <- length(which(x==0))
# 	a2 <- length(which(x==1/4))
# 	a3 <- length(which(x==1/2))
# 	a4 <- length(which(x==3/4))
# 	a5 <- length(which(x==1))
# 	total <- a1+a2+a3+a4+a5
# 	return(c(a1/total,a2/total,a3/total,a4/total,a5/total))
# }
# round(MethylationClassSummary(hg_methyclass[,selected_LPN]),3)

sys_argv <- c("0.686,0.138,0.053,0.063,0.060","Fig3F","0.01")
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
# install.packages("rgl")
# install.packages("plot3D")

# library(rgl)
library(plot3D)
# library(scales)

MethylationTransMatrix <- function(u,d,p){
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
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
start <- as.matrix(as.numeric(strsplit(sys_argv[1],",")[[1]]))
out_name <- sys_argv[2]
step <- as.numeric(sys_argv[3])
total <- 1000
start <- start*total
total <- sum(start)

u_range_start <- 0
u_range_end <- 1
p_range_start <- 0
p_range_end <- 1
d_range_start <- 0
d_range_end <- 1

paramaters <- c()
for (u in c(0.03)){
	for (p in seq(p_range_start,p_range_end,step)){
		for (d in seq(d_range_start,d_range_end,step)){
			u <- round(u,3)
			p <- round(p,3)
			d <- round(d,3)
			t <- MethylationTransMatrix(u,d,p)
			t <- t*total
			t1 <- c(rep(0,t[1,1]),rep(0.25,t[2,1]),rep(0.5,t[3,1]),rep(0.75,t[4,1]),rep(1,t[5,1]))
			t2 <- c(rep(0,t[1,2]),rep(0.25,t[2,2]),rep(0.5,t[3,2]),rep(0.75,t[4,2]),rep(1,t[5,2]))
			t3 <- c(rep(0,t[1,3]),rep(0.25,t[2,3]),rep(0.5,t[3,3]),rep(0.75,t[4,3]),rep(1,t[5,3]))
			t4 <- c(rep(0,t[1,4]),rep(0.25,t[2,4]),rep(0.5,t[3,4]),rep(0.75,t[4,4]),rep(1,t[5,4]))
			t5 <- c(rep(0,t[1,5]),rep(0.25,t[2,5]),rep(0.5,t[3,5]),rep(0.75,t[4,5]),rep(1,t[5,5]))
			
			t_cell_1 <- c();t_cell_2 <- c()
			f_cell_1 <- c();f_cell_2 <- c();f_cell_3 <- c();f_cell_4 <- c()
			e_cell_1 <- c();e_cell_2 <- c();e_cell_3 <- c();e_cell_4 <- c();e_cell_5 <- c();e_cell_6 <- c();e_cell_7 <- c();e_cell_8 <- c()
			for (each in seq(start[1,])){
				t_cell_1 <- c(t_cell_1,sample(t1,1));t_cell_2 <- c(t_cell_2,sample(t1,1))
			}
			for (each in seq(start[2,])){
				t_cell_1 <- c(t_cell_1,sample(t2,1));t_cell_2 <- c(t_cell_2,sample(t2,1))
			}
			for (each in seq(start[3,])){
				t_cell_1 <- c(t_cell_1,sample(t3,1));t_cell_2 <- c(t_cell_2,sample(t3,1))
			}
			for (each in seq(start[4,])){
				t_cell_1 <- c(t_cell_1,sample(t4,1));t_cell_2 <- c(t_cell_2,sample(t4,1))
			}
			for (each in seq(start[5,])){
				t_cell_1 <- c(t_cell_1,sample(t5,1));t_cell_2 <- c(t_cell_2,sample(t5,1))
			}
			# 2cell to 4cell
			for (each in seq(total)){
				# 2cell_1 -> 4cell_1+4cell_2
				if (t_cell_1[each] == 0.00){f_cell_1 <- c(f_cell_1,sample(t1,1));f_cell_2 <- c(f_cell_2,sample(t1,1))}
				if (t_cell_1[each] == 0.25){f_cell_1 <- c(f_cell_1,sample(t2,1));f_cell_2 <- c(f_cell_2,sample(t2,1))}
				if (t_cell_1[each] == 0.50){f_cell_1 <- c(f_cell_1,sample(t3,1));f_cell_2 <- c(f_cell_2,sample(t3,1))}
				if (t_cell_1[each] == 0.75){f_cell_1 <- c(f_cell_1,sample(t4,1));f_cell_2 <- c(f_cell_2,sample(t4,1))}
				if (t_cell_1[each] == 1.00){f_cell_1 <- c(f_cell_1,sample(t5,1));f_cell_2 <- c(f_cell_2,sample(t5,1))}
				# 2cell_2 -> 4cell_3+4cell_4
				if (t_cell_2[each] == 0.00){f_cell_3 <- c(f_cell_3,sample(t1,1));f_cell_4 <- c(f_cell_4,sample(t1,1))}
				if (t_cell_2[each] == 0.25){f_cell_3 <- c(f_cell_3,sample(t2,1));f_cell_4 <- c(f_cell_4,sample(t2,1))}
				if (t_cell_2[each] == 0.50){f_cell_3 <- c(f_cell_3,sample(t3,1));f_cell_4 <- c(f_cell_4,sample(t3,1))}
				if (t_cell_2[each] == 0.75){f_cell_3 <- c(f_cell_3,sample(t4,1));f_cell_4 <- c(f_cell_4,sample(t4,1))}
				if (t_cell_2[each] == 1.00){f_cell_3 <- c(f_cell_3,sample(t5,1));f_cell_4 <- c(f_cell_4,sample(t5,1))}
			}
			for (each in seq(total)){
				# 4cell_1 -> 8cell_1+8cell_2
				if (f_cell_1[each] == 0.00){e_cell_1 <- c(e_cell_1,sample(t1,1));e_cell_2 <- c(e_cell_2,sample(t1,1))}
				if (f_cell_1[each] == 0.25){e_cell_1 <- c(e_cell_1,sample(t2,1));e_cell_2 <- c(e_cell_2,sample(t2,1))}
				if (f_cell_1[each] == 0.50){e_cell_1 <- c(e_cell_1,sample(t3,1));e_cell_2 <- c(e_cell_2,sample(t3,1))}
				if (f_cell_1[each] == 0.75){e_cell_1 <- c(e_cell_1,sample(t4,1));e_cell_2 <- c(e_cell_2,sample(t4,1))}
				if (f_cell_1[each] == 1.00){e_cell_1 <- c(e_cell_1,sample(t5,1));e_cell_2 <- c(e_cell_2,sample(t5,1))}
				# 4cell_2 -> 8cell_3+8cell_4
				if (f_cell_2[each] == 0.00){e_cell_3 <- c(e_cell_3,sample(t1,1));e_cell_4 <- c(e_cell_4,sample(t1,1))}
				if (f_cell_2[each] == 0.25){e_cell_3 <- c(e_cell_3,sample(t2,1));e_cell_4 <- c(e_cell_4,sample(t2,1))}
				if (f_cell_2[each] == 0.50){e_cell_3 <- c(e_cell_3,sample(t3,1));e_cell_4 <- c(e_cell_4,sample(t3,1))}
				if (f_cell_2[each] == 0.75){e_cell_3 <- c(e_cell_3,sample(t4,1));e_cell_4 <- c(e_cell_4,sample(t4,1))}
				if (f_cell_2[each] == 1.00){e_cell_3 <- c(e_cell_3,sample(t5,1));e_cell_4 <- c(e_cell_4,sample(t5,1))}
				# 4cell_3 -> 8cell_5+8cell_6
				if (f_cell_3[each] == 0.00){e_cell_5 <- c(e_cell_5,sample(t1,1));e_cell_6 <- c(e_cell_6,sample(t1,1))}
				if (f_cell_3[each] == 0.25){e_cell_5 <- c(e_cell_5,sample(t2,1));e_cell_6 <- c(e_cell_6,sample(t2,1))}
				if (f_cell_3[each] == 0.50){e_cell_5 <- c(e_cell_5,sample(t3,1));e_cell_6 <- c(e_cell_6,sample(t3,1))}
				if (f_cell_3[each] == 0.75){e_cell_5 <- c(e_cell_5,sample(t4,1));e_cell_6 <- c(e_cell_6,sample(t4,1))}
				if (f_cell_3[each] == 1.00){e_cell_5 <- c(e_cell_5,sample(t5,1));e_cell_6 <- c(e_cell_6,sample(t5,1))}
				# 4cell_4 -> 8cell_7+8cell_8
				if (f_cell_4[each] == 0.00){e_cell_7 <- c(e_cell_7,sample(t1,1));e_cell_8 <- c(e_cell_8,sample(t1,1))}
				if (f_cell_4[each] == 0.25){e_cell_7 <- c(e_cell_7,sample(t2,1));e_cell_8 <- c(e_cell_8,sample(t2,1))}
				if (f_cell_4[each] == 0.50){e_cell_7 <- c(e_cell_7,sample(t3,1));e_cell_8 <- c(e_cell_8,sample(t3,1))}
				if (f_cell_4[each] == 0.75){e_cell_7 <- c(e_cell_7,sample(t4,1));e_cell_8 <- c(e_cell_8,sample(t4,1))}
				if (f_cell_4[each] == 1.00){e_cell_7 <- c(e_cell_7,sample(t5,1));e_cell_8 <- c(e_cell_8,sample(t5,1))}
			}
			e_cell <- cbind(e_cell_1,e_cell_2,e_cell_3,e_cell_4,e_cell_5,e_cell_6,e_cell_7,e_cell_8)
			hetergeneity <- apply(e_cell,1,MethylationHeterogeneity)
			paramaters <- rbind(paramaters,c(u,p,d,mean(hetergeneity[1:start[1,]]),mean(hetergeneity[(start[1,]+1):(start[1,]+start[2,])]),mean(hetergeneity[(start[1,]+start[2,]+1):(start[1,]+start[2,]+start[3,])]),mean(hetergeneity[(start[1,]+start[2,]+start[3,]+1):(start[1,]+start[2,]+start[3,]+start[4,])]),mean(hetergeneity[(start[1,]+start[2,]+start[3,]+start[4,]+1):total]),mean(hetergeneity)))
		}
	}
}
write.table(paramaters,file=paste(out_name,".paramaters.txt",sep=""),quote=F,sep="\t",col.names=F,row.names=F)
paramaters <- read.table(paste(out_name,".paramaters.txt",sep=""))


# tmp_matrix <- paramaters[which(paramaters[,4]<0.2 & paramaters[,8]<0.2),]

summary(paramaters[which(paramaters[,2]>0.94 & paramaters[,3]<0.05),c(4,8)])

###################################################################################################
################################                PLOT               ################################
###################################################################################################
u <- paramaters[,1]
p <- paramaters[,2]
d <- paramaters[,3]
e1 <- paramaters[,4]
e2 <- paramaters[,5]
e3 <- paramaters[,6]
e4 <- paramaters[,7]
e5 <- paramaters[,8]
e <- paramaters[,9]

zmax <- max(paramaters[,c(4,5,6,7,8)]);zmin <- min(paramaters[,c(4,5,6,7,8)]);
ColorRamp <- colorRampPalette(c("#1b315e","#008792","#afdfe4","white","orange","red","#CE0013"), bias=1)(1000*(zmax-zmin))
pdf("Fig3F.5ClassupHeterogeneity3Dsurf.pdf",width=3,height=12.5)
par(mfrow=c(5,1))
for(each_h in seq(4,8)){
	tmp_u <- 0.03
	tmp_matrix <- matrix(paramaters[which(paramaters[,1]==tmp_u),each_h],c(1/step+1,1/step+1),byrow=T)
	tmp_x <- matrix(paramaters[which(paramaters[,1]==tmp_u),2],c(1/step+1,1/step+1),byrow=T)
	tmp_y <- matrix(paramaters[which(paramaters[,1]==tmp_u),3],c(1/step+1,1/step+1),byrow=T)
	surf3D(tmp_x, tmp_y, tmp_matrix, ticktype = "simple", pch = 16,type="l", bty = "b2",lighting = F, ltheta = 45, theta = 45, xlim = c(0,1), ylim = c(0,1), zlim = c(0,1), col=ColorRamp[round((min(tmp_matrix)-zmin)*1000):round((max(tmp_matrix)-zmin)*1000)],NAcol = "grey", shade = 0.05,xlab = "p", ylab = "d", zlab = "Heterogeneity", main=paste("S",each_h-3,sep=""),colkey=F)
}
dev.off()
pdf("Fig3F.5ClassupHeterogeneity3Dsurf_legend.pdf",width=2.5,height=2.5)
par(mar=c(2,6,2,6))
ColorLevels <- seq(to=zmax,from=zmin, length=1000)   #number sequence
image(1,ColorLevels,t(matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1)),col=t(ColorRamp), xlab="",ylab="",xaxt="n",yaxt="n",useRaster=F);box(lwd=2)
axis(side=4,c(zmin,round((zmax+zmin)/2,2),zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,2),round(zmax,2)),cex.axis=1,las=2)
dev.off()

