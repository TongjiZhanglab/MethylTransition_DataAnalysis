###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol01 <- c("#CE001301","#16557A01","#C7A60901","#87C23201","#00879201","#A14C9401","#15A08C01","#8B7E7501","#1E7CAF01","#EA425F01","#46489A01","#E5003301","#0F231F01","#1187CD01")
cccol05 <- c("#CE001305","#16557A05","#C7A60905","#87C23205","#00879205","#A14C9405","#15A08C05","#8B7E7505","#1E7CAF05","#EA425F05","#46489A05","#E5003305","#0F231F05","#1187CD05")
cccol30 <- c("#CE001330","#16557A30","#C7A60930","#87C23230","#00879230","#A14C9430","#15A08C30","#8B7E7530","#1E7CAF30","#EA425F30","#46489A30","#E5003330","#0F231F30","#1187CD30")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E8003380","#0F231F80","#1187CD80")
options(scipen = 200)
library(MethylTransition)
library(beeswarm)

ErrorBar <- function(x, y, upper, lower=upper, length=0.05,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
HalfErrorBar <- function(x, y, upper, lower=upper, length=0.03,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	arrows(x,y+upper, x, y, angle=90, code=1, length=length, ...)
}
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
hg_methyclass <- read.table("../Data/hg_MethylationClass.txt",header=T,row.names=1)
hg_genes <- rownames(hg_methyclass)
hg_all_cells <- colnames(hg_methyclass)

hg_1C <- as.vector(read.table("Fig2A.1C_cells.txt")[,1])
hg_2C <- as.vector(read.table("Fig2A.2C_cells.txt")[,1])
hg_4C <- as.vector(read.table("Fig2A.4C_cells.txt")[,1])
hg_8C <- as.vector(read.table("Fig2A.8C_cells.txt")[,1])
parameters <- c()
rownames <- c()
for (each_1C in hg_1C){
	for (each_2C in hg_2C){
		tmp_para <- tryCatch(ParameterEstimation(hg_methyclass[,each_1C],hg_methyclass[,each_2C],iter=100)$estimated_parameters,error = function(a) c(NA,NA,NA))
		parameters <- rbind(parameters,tmp_para)
		rownames <- c(rownames,paste(each_1C,each_2C,sep="__"))
	}
}
for (each_2C in hg_2C){
	for (each_4C in hg_4C){
		tmp_para <- tryCatch(ParameterEstimation(hg_methyclass[,each_2C],hg_methyclass[,each_4C],iter=100)$estimated_parameters,error = function(a) c(NA,NA,NA))
		parameters <- rbind(parameters,tmp_para)
		rownames <- c(rownames,paste(each_2C,each_4C,sep="__"))
	}
}
for (each_4C in hg_4C){
	for (each_8C in hg_8C){
		tmp_para <- tryCatch(ParameterEstimation(hg_methyclass[,each_4C],hg_methyclass[,each_8C],iter=100)$estimated_parameters,error = function(a) c(NA,NA,NA))
		parameters <- rbind(parameters,tmp_para)
		rownames <- c(rownames,paste(each_4C,each_8C,sep="__"))
	}
}
rownames(parameters) <- rownames
write.table(file="hg_NewtonParameters_udp.2A.txt",parameters,quote=F,sep="\t",col.names=F,row.names=T)
trans_para <- parameters
###################################################################################################
################################                PLOT               ################################
###################################################################################################
trans_para <- read.table("hg_NewtonParameters_udp.2A.txt",header=F,row.names=1)
trans_para <- na.omit(trans_para)
samples <- row.names(trans_para)

hg_LPN_2C <- samples[grep("LPN",samples)]
hg_2C_4C <- samples[grep("2C",samples)][grep("4C",samples[grep("2C",samples)])]
hg_4C_8C <- samples[grep("4C",samples)][grep("8C",samples[grep("4C",samples)])]


pdf("Fig2A.MethyRatioTransitionPerStage.pdf",width=8,height=4)
par(mar=c(2,6,4,4),mfrow=c(1,3))
boxplot(trans_para[hg_LPN_2C,],xaxt="n",col = cccol30,boxwex=0.6,border=cccol,outline=F,main="",xlab="",ylab="Parameter Distance",ylim=c(0,1),las=2,lwd=2);box(lwd=2)
axis(at=seq(3),side=1,labels=c("u","d","p"))
beeswarm(trans_para[hg_LPN_2C,],xaxt="n",spacing = 1,col = cccol50,pch = 16,method = "swarm",corral = "wrap",add=TRUE)
boxplot(trans_para[hg_2C_4C,],xaxt="n",col = cccol30,boxwex=0.6,border=cccol,outline=F,main="",xlab="",ylab="Parameter Distance",ylim=c(0,1),las=2,lwd=2);box(lwd=2)
axis(at=seq(3),side=1,labels=c("u","d","p"))
beeswarm(trans_para[hg_2C_4C,],xaxt="n",spacing = 1,col = cccol50,pch = 16,method = "swarm",corral = "wrap",add=TRUE)
boxplot(trans_para[hg_4C_8C,],xaxt="n",col = cccol30,boxwex=0.6,border=cccol,outline=F,main="",xlab="",ylab="Parameter Distance",ylim=c(0,1),las=2,lwd=2);box(lwd=2)
axis(at=seq(3),side=1,labels=c("u","d","p"))
beeswarm(trans_para[hg_4C_8C,],xaxt="n",spacing = 1,col = cccol50,pch = 16,method = "swarm",corral = "wrap",add=TRUE)
dev.off()
