###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001330","#16557A30","#C7A60930","#87C23250","#64C0AB50","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")


error.bar <- function(x, y, upper, lower=upper, length=0.02,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
trans_para <- read.table("../../Data/NewtonParameters_udp.txt",header=F,row.names=1)
trans_para <- na.omit(trans_para)
samples <- row.names(trans_para)

hg_LPN_2C <- samples[grep("LPN",samples)]
hg_2C_4C <- samples[grep("2C",samples)][grep("4C",samples[grep("2C",samples)])]
hg_4C_8C <- samples[grep("4C",samples)][grep("8C",samples[grep("4C",samples)])]

###################################################################################################
################################                PLOT               ################################
###################################################################################################

plot_n <- c()
plot_mean <- c()
plot_sd <- c()

plot_data <- trans_para[hg_LPN_2C,1:3]
plot_n <- rbind(plot_n,apply(plot_data,2,length))
plot_mean <- rbind(plot_mean,apply(plot_data,2,mean))
plot_sd <- rbind(plot_sd,apply(plot_data,2,sd))
plot_data <- trans_para[hg_2C_4C,1:3]
plot_n <- rbind(plot_n,apply(plot_data,2,length))
plot_mean <- rbind(plot_mean,apply(plot_data,2,mean))
plot_sd <- rbind(plot_sd,apply(plot_data,2,sd))
plot_data <- trans_para[hg_4C_8C,1:3]
plot_n <- rbind(plot_n,apply(plot_data,2,length))
plot_mean <- rbind(plot_mean,apply(plot_data,2,mean))
plot_sd <- rbind(plot_sd,apply(plot_data,2,sd))

for (selected_LPN in c("hg_LPN_11","hg_LPN_12","hg_LPN_8","hg_LPN_9")){
	u_matrix <- c()
	d_matrix <- c()
	p_matrix <- c()
	for (x in seq(length(hg_LPN_2C))){
		for (y in seq(length(hg_2C_4C))){
			for (z in seq(length(hg_4C_8C))){
				tmp_1 <- strsplit(hg_LPN_2C[x],"__")[[1]][2]
				tmp_2 <- strsplit(hg_2C_4C[y],"__")[[1]][1]
				tmp_3 <- strsplit(hg_2C_4C[y],"__")[[1]][2]
				tmp_4 <- strsplit(hg_4C_8C[z],"__")[[1]][1]
				if (tmp_1==tmp_2 & tmp_3==tmp_4 & grepl(selected_LPN,hg_LPN_2C[x]) & grepl("2C_6",hg_LPN_2C[x]) & grepl("4C_7",hg_2C_4C[y]) & grepl("8C_9",hg_4C_8C[z])){
					u_matrix <- rbind(u_matrix,trans_para[c(hg_LPN_2C[x],hg_2C_4C[y],hg_4C_8C[z]),1])
					d_matrix <- rbind(d_matrix,trans_para[c(hg_LPN_2C[x],hg_2C_4C[y],hg_4C_8C[z]),2])
					p_matrix <- rbind(p_matrix,trans_para[c(hg_LPN_2C[x],hg_2C_4C[y],hg_4C_8C[z]),3])
				}
			}
		}
	}
	# pdf(paste("Fig3C.MethyRatioTransitionPerStage_",selected_LPN,".pdf",sep=""),width=5.5,height=4)
	pdf(paste("Fig3C.MethyRatioTransitionPerStage_",selected_LPN,".pdf",sep=""),width=5,height=7)
	par(mar=c(6,4,4,2))
	plot(1,type="n",xaxt="n",bty="n",xlab="",ylab="Parameter value",ylim=c(0,1),xlim=c(1,nrow(plot_mean)),main="Human")
	box(lwd=2)
	v1 = apply(u_matrix,2,mean)
	v2 = apply(u_matrix,2,quantile,0.1)
	v3 = apply(u_matrix,2,quantile,0.9)
	points(v1,lwd=3,type="l",col=cccol[1])
	polygon(c(1,1:ncol(u_matrix),ncol(u_matrix):2),c(v2[1],v3,v2[ncol(u_matrix):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
	v1 = apply(d_matrix,2,mean)
	v2 = apply(d_matrix,2,quantile,0.1)
	v3 = apply(d_matrix,2,quantile,0.9)
	points(v1,lwd=3,type="l",col=cccol[2])
	polygon(c(1,1:ncol(d_matrix),ncol(d_matrix):2),c(v2[1],v3,v2[ncol(d_matrix):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
	v1 = apply(p_matrix,2,mean)
	v2 = apply(p_matrix,2,quantile,0.1)
	v3 = apply(p_matrix,2,quantile,0.9)
	points(v1,lwd=3,type="l",col=cccol[3])
	polygon(c(1,1:ncol(p_matrix),ncol(p_matrix):2),c(v2[1],v3,v2[ncol(p_matrix):2]),col=adjustcolor("grey", alpha.f = 0.3),border=NA)
	axis(side=1,at=seq(nrow(plot_mean)),labels=c("1-cell -> 2cell","2-cell -> 4-cell","4-cell -> 8-cell"))
	legend("topleft",c("u","d","p"),bty="n",col=cccol,lwd=2)
	dev.off()
}