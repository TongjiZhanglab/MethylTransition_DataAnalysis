###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#64C0AB50","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")

library(ineq)
MethylationHeterogeneityObservation <- function(x){
	a1 <- length(which(x==0))
	a2 <- length(which(x==1/4))
	a3 <- length(which(x==1/2))
	a4 <- length(which(x==3/4))
	a5 <- length(which(x==1))
	total <- a1+a2+a3+a4+a5
	if (total >=2){
		return(1-Gini(c(a1/total,a2/total,a3/total,a4/total,a5/total),corr=T))
	}else{
		return(NA)
	}
}

Boxplot_List <- function(plot_list,MAIN,NAMES,YLAB,COL,COL50){
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
	boxplot(plot_list,at=0:(n-1)+1,boxwex=0.6, xlim=c(0.5,n+0.5),ylim=ylim,col=COL50,border=COL,outline=F,names=NAMES,main=MAIN,lwd=2,ylab=YLAB,las=2)
	box(lwd=2)
	ypos_1 <- bp$stats[5,][1:n-1]
	ypos_2 <- bp$stats[5,][2:n]
	# legend("topleft",c("***:p<0.001","**:p<0.01","*:p<0.05"),text.col="red",bty="n")
	if (length(mark)!=0){
		for(i in 1:length(mark)){
			if(!is.na(mark[i]) & mark[i]!="?"){
				segments(xpos[i]-.45, ypos_1[i]+dist/2, xpos[i]-.45, max(ypos_1[i], ypos_2[i])+dist)
				segments(xpos[i]+.45, ypos_2[i]+dist/2, xpos[i]+.45, max(ypos_1[i], ypos_2[i])+dist)
				segments(xpos[i]-.45, max(ypos_1[i], ypos_2[i])+dist, xpos[i]-0.3, max(ypos_1[i], ypos_2[i])+dist)
				segments(xpos[i]+.45, max(ypos_1[i], ypos_2[i])+dist, xpos[i]+0.3, max(ypos_1[i], ypos_2[i])+dist)
				text(x=xpos[i], y=max(ypos_1[i], ypos_2[i])+dist, label=mark[i], col="black")
			}
		}
	}
	# for (i in seq(n)){
	# 	text(x=i, y=bp$stats[3,i], label=length(na.omit(plot_list[[i]])), col="white",cex=0.3)
	# }
}
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################

########## DNA methylation ##########
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

hg_2C_6_Gini <- apply(na.omit(hg_methyclass[,hg_2C_6_id]),1,MethylationHeterogeneityObservation)
hg_4C_7_Gini <- apply(na.omit(hg_methyclass[,hg_4C_7_id]),1,MethylationHeterogeneityObservation)
hg_4C_8_Gini <- apply(na.omit(hg_methyclass[,hg_4C_8_id]),1,MethylationHeterogeneityObservation)
hg_8C_6_Gini <- apply(na.omit(hg_methyclass[,hg_8C_6_id]),1,MethylationHeterogeneityObservation)
hg_8C_9_Gini <- apply(na.omit(hg_methyclass[,hg_8C_9_id]),1,MethylationHeterogeneityObservation)

hg_2C_6_mean <- apply(na.omit(hg_methyclass[,hg_2C_6_id]),1,mean)
hg_4C_7_mean <- apply(na.omit(hg_methyclass[,hg_4C_7_id]),1,mean)
hg_4C_8_mean <- apply(na.omit(hg_methyclass[,hg_4C_8_id]),1,mean)
hg_8C_6_mean <- apply(na.omit(hg_methyclass[,hg_8C_6_id]),1,mean)
hg_8C_9_mean <- apply(na.omit(hg_methyclass[,hg_8C_9_id]),1,mean)

########## gene expression ##########
data <- read.table("../Data/nsmb.2660-S2.txt",header=T,row.names=1)
Oocyte <- 1:3; Zygote <- 4:6; cell2 <- 7:12; cell4 <- 13:24; cell8 <- 25:44; Morula <- 45:60; 
MTE <- c(64,66,67,69,72,76:79);
PTE <- c(61:63,65,68,70,71,81,82); 
PE <- c(84:90);
EPI <- c(73:75,80,83);
hESC0 <- 91:98; hESC10 <- 99:124
blastocyst <- 61:90;

avg <- cbind(apply(data[,Oocyte],1,mean),apply(data[,Zygote],1,mean),apply(data[,cell2],1,mean),apply(data[,cell4],1,mean),apply(data[,cell8],1,mean),apply(data[,Morula],1,mean),apply(data[,MTE],1,mean),apply(data[,PTE],1,mean),apply(data[,PE],1,mean),apply(data[,EPI],1,mean),apply(data[,hESC0],1,mean),apply(data[,hESC10],1,mean))
time_point <- c("Oocyte","Zygote","X2cell","X4cell","X8cell","Morula","MTE","PTE","PE","EPI","hESC0","hESC10")
dev_labels <- c("Oocyte","Zygote","2cell","4cell","8cell","Morula","MTE","PTE","PE","EPI","hESC0","hESC10")
colnames(avg) <- time_point
development_path <- time_point
dData <- log2(avg+1)

# data_log <- log2(data[,1:124]+1)

E1_4C <- grep("4cell_e1",colnames(data))
E2_4C <- grep("4cell_e2",colnames(data))
E3_4C <- grep("4cell_e3",colnames(data))
E2_8C <- grep("8cell_e2",colnames(data))
E3_8C <- grep("8cell_e3",colnames(data))
E1_Morula <- grep("Morulae_1",colnames(data))
E2_Morula <- grep("Morulae_2",colnames(data))
E1_blastocyst <- grep("Late_blastocyst_1",colnames(data)) # MTE & PTE
E2_blastocyst <- grep("Late_blastocyst_2",colnames(data)) # MTE & PTE & EPI
E3_blastocyst <- grep("Late_blastocyst_3",colnames(data)) # PE & EPI(83)

mean_matrix <- cbind(apply(data[,E1_4C],1,mean,na.rm=T),
					 apply(data[,E2_4C],1,mean,na.rm=T),
					 apply(data[,E3_4C],1,mean,na.rm=T),
					 apply(data[,E2_8C],1,mean,na.rm=T),
					 apply(data[,E3_8C],1,mean,na.rm=T),
					 apply(data[,E1_Morula],1,mean,na.rm=T),
					 apply(data[,E2_Morula],1,mean,na.rm=T),
					 apply(data[,E1_blastocyst],1,mean,na.rm=T),
					 apply(data[,E2_blastocyst],1,mean,na.rm=T),
					 apply(data[,E3_blastocyst],1,mean,na.rm=T))
var_matrix <- cbind(apply(data[,E1_4C],1,var,na.rm=T),
					 apply(data[,E2_4C],1,var,na.rm=T),
					 apply(data[,E3_4C],1,var,na.rm=T),
					 apply(data[,E2_8C],1,var,na.rm=T),
					 apply(data[,E3_8C],1,var,na.rm=T),
					 apply(data[,E1_Morula],1,var,na.rm=T),
					 apply(data[,E2_Morula],1,var,na.rm=T),
					 apply(data[,E1_blastocyst],1,var,na.rm=T),
					 apply(data[,E2_blastocyst],1,var,na.rm=T),
					 apply(data[,E3_blastocyst],1,var,na.rm=T))
CV2 <- var_matrix/(mean_matrix)^2
embryo_labels <- c("4C_e1","4C_e2","4C_e3","8C_e2","8C_e3","Morula_e1","Morula_e2","blastocyst_e1","blastocyst_e2","blastocyst_e3")
colnames(CV2) <- embryo_labels
colnames(var_matrix) <- embryo_labels
colnames(mean_matrix) <- embryo_labels

# ########## typical genes ##########
# The specific gene is predicted from /mnt/Storage/home/zhaochengchen/Work/1.scMethylome/data/GSE71318_hgICMTE_RNAseq/DESeq2
IT_fpkm <- read.table("../Data/GSE71318_hgICMTEFPKM.txt",header=T,row.names=1)
IT_mean <- cbind(apply(IT_fpkm[,c("ICM_1_FPKM","ICM_2_FPKM","ICM_3_FPKM")],1,mean,na.rm=T),apply(IT_fpkm[,c("TE_1_FPKM","TE_2_FPKM","TE_3_FPKM")],1,mean,na.rm=T))
# IT_mean <- log2(IT_mean+1)
ICM_specific <- row.names(read.table("../Data/GSE71318_ICM_DESeqGenes.txt",header=T,row.names=1))
TE_specific <- row.names(read.table("../Data/GSE71318_TE_DESeqGenes.txt",header=T,row.names=1))

ICM_specific <- intersect(intersect(ICM_specific,rownames(CV2)),names(hg_8C_9_Gini))
TE_specific <- intersect(intersect(TE_specific,rownames(CV2)),names(hg_8C_9_Gini))
CFGene <- c(ICM_specific,TE_specific)
print(length(ICM_specific))
print(length(TE_specific))
print(length(CFGene))

all_expressed_genes <- row.names(IT_mean)[which(IT_mean[,1]>1 & IT_mean[,2]>1)]
print(length(all_expressed_genes))
others <- setdiff(row.names(IT_mean),CFGene)
expressed_genes <- others[which(IT_mean[others,1]>1 & IT_mean[others,2]>1)]

###################################################################################################
################################                PLOT               ################################
###################################################################################################
# plot_list <- list(hg_8C_9_Gini[CFGene],hg_8C_9_Gini[setdiff(hg_genes,CFGene)])
# Boxplot_List(plot_list,"",c("Cell fate determination related genes","Others"),"DNA methylation heterogeneity",COL,COL50)

# a <- table(hg_8C_9_Gini[CFGene])
# b <- table(hg_8C_9_Gini[setdiff(hg_genes,CFGene)])
# a/sum(a)*100
# b/sum(b)*100

# 1-Gini(c(4,4,0,0,0)*0.125,corr=T)

# 1-Gini(c(8,0,0,0,0)*0.125,corr=T)
# 1-Gini(c(7,1,0,0,0)*0.125,corr=T)
# 1-Gini(c(6,2,0,0,0)*0.125,corr=T)
# 1-Gini(c(6,1,1,0,0)*0.125,corr=T)
# 1-Gini(c(5,3,0,0,0)*0.125,corr=T)
# 1-Gini(c(5,2,1,0,0)*0.125,corr=T)
# 1-Gini(c(5,1,1,1,0)*0.125,corr=T)
# 1-Gini(c(4,4,0,0,0)*0.125,corr=T)
# 1-Gini(c(4,3,1,0,0)*0.125,corr=T)
# 1-Gini(c(4,2,2,0,0)*0.125,corr=T)
# 1-Gini(c(4,2,1,1,0)*0.125,corr=T)
# 1-Gini(c(4,1,1,1,1)*0.125,corr=T)
# # 1-Gini(c(3,5,0,0,0)*0.125,corr=T)
# # 1-Gini(c(3,4,1,0,0)*0.125,corr=T)
# 1-Gini(c(3,3,2,0,0)*0.125,corr=T)
# 1-Gini(c(3,2,1,2,0)*0.125,corr=T)
# 1-Gini(c(3,3,1,1,0)*0.125,corr=T)
# 1-Gini(c(3,2,1,1,1)*0.125,corr=T)
# # 1-Gini(c(2,6,0,0,0)*0.125,corr=T)
# # 1-Gini(c(2,5,1,0,0)*0.125,corr=T)
# # 1-Gini(c(2,4,1,1,0)*0.125,corr=T)
# # 1-Gini(c(2,3,2,1,0)*0.125,corr=T)

pdf("Fig5D.MethylHeterogeneityofCFGenes.Boxplot.pdf",width=2.5,height=4.5)
par(mar=c(6,4,4,2))
COL <- c("#2874a6",cccol[8])
COL50 <- paste(COL,"50",sep="")
plot_list <- list(hg_8C_9_Gini[CFGene],hg_8C_9_Gini[expressed_genes])
Boxplot_List(plot_list,"",c("Cell fate determination related genes","Others"),"DNA methylation heterogeneity",COL,COL50)
dev.off()


pdf("Fig5E.ExpreHeterogeneityofCFGenes.Boxplot.pdf",width=2.5,height=4.5)
par(mar=c(6,4,4,2))
COL <- c(cccol[7],cccol[8])
COL50 <- paste(COL,"50",sep="")
plot_list <- list(CV2[CFGene,"8C_e2"],CV2[intersect(expressed_genes,row.names(CV2)),"8C_e2"])
Boxplot_List(plot_list,"",c("Cell fate determination related genes","Others"),"Expression heterogeneity",COL,COL50)
dev.off()

# # pdf("Fig5F.MH_EH_distribution.pdf",width=4,height=4)
# # hist(hg_8C_9_Gini,breaks=50,col=cccol[1],border=cccol[1],main="MH_C8")
# # abline(v=quantile(hg_8C_9_Gini,c(0.5,0.9),na.rm=T),lty=2)
# # hist(CV2[,"8C_e2"],breaks=50,col=cccol[2],border=cccol[2],main="EH_C8")
# # abline(v=quantile(CV2[,"8C_e2"],c(0.5,0.9),na.rm=T),lty=2)
# # dev.off()

# # quantile(hg_8C_9_Gini,seq(0,1,0.1),na.rm=T)
# # quantile(CV2[,"8C_e2"],seq(0,1,0.1),na.rm=T)         

# m_cut <- 0.3
# e_cut <- 0.5
# predict_genes <- intersect(CFGene[hg_8C_9_Gini[CFGene]>m_cut],CFGene[CV2[CFGene,"8C_e2"]>e_cut])

# length(predict_genes)
# write.table(cbind(predict_genes,hg_8C_9_Gini[predict_genes],round(CV2[predict_genes,"8C_e2"],4)),file="STable2.MethylHeterVSExpreHeter.PredictedGenes.txt",col.names = F,row.names = F,sep="\t",quote=F)

# MH_mean_C8 <- hg_8C_9_Gini[CFGene]
# EH_mean_C8 <- CV2[CFGene,"8C_e2"]
# HH <- intersect(CFGene[hg_8C_9_Gini[CFGene]>m_cut],CFGene[CV2[CFGene,"8C_e2"]>e_cut])
# HL <- intersect(CFGene[hg_8C_9_Gini[CFGene]>m_cut],CFGene[CV2[CFGene,"8C_e2"]<=e_cut])
# LH <- intersect(CFGene[hg_8C_9_Gini[CFGene]<=m_cut],CFGene[CV2[CFGene,"8C_e2"]>e_cut])
# LL <- intersect(CFGene[hg_8C_9_Gini[CFGene]<=m_cut],CFGene[CV2[CFGene,"8C_e2"]<=e_cut])
# boundary <- c(length(HH),length(c(HH,HL)),length(c(HH,HL,LH)))+0.5

# marker <- c("GDF3","CUBN","IGF1","GCM1")

# pdf("Fig5F.MethylHeterVSExpreHeter.Heatmap.pdf",width=3,height=8)
# layout(matrix(c(rep(1,16),rep(2,3),rep(3,16),rep(4,3)),nrow=19,ncol=2,byrow=F))
# par(oma=c(0.5,0.5,0.5,0.5),mar=c(2,2,2,2))
# mty_matrix <- c(MH_mean_C8[c(HH,HL,LH,LL)])
# zmax <- max(mty_matrix,na.rm=T)
# zmin <- min(mty_matrix,na.rm=T)
# ColorRamp <- colorRampPalette(c("white","#d6eaf8","#3498db","#2874a6","#1b4f72"), bias=1)(10000)   #color list
# ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence
# image(1, 1:length(mty_matrix[c(HH,HL,LH,LL)]), t(mty_matrix[c(HH,HL,LH,LL)]), xaxt="n", yaxt="n", col=ColorRamp, xlab="", ylab="");box(lwd=2)
# axis(side=2,at=c(which(c(HH,HL,LH,LL)==marker[1]),which(c(HH,HL,LH,LL)==marker[2]),which(c(HH,HL,LH,LL)==marker[3]),which(c(HH,HL,LH,LL)==marker[4]))+0.5,labels=marker,cex.axis=0.3,las=2)
# # axis(side=3,1,labels=c("8-cell"),cex.axis=1)
# abline(h=boundary,lty=1,lwd=2)
# image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
# axis(side=1,c(zmin,(zmin+zmax)/2,zmax),labels=c(round(zmin,2),round((zmin+zmax)/2,2),round(zmax,2)))

# exp_matrix <- c(EH_mean_C8[c(HH,HL,LH,LL)])
# zmax <- quantile(exp_matrix,0.99,na.rm=T)
# zmin <- quantile(exp_matrix,0.01,na.rm=T)
# exp_matrix[exp_matrix<zmin] <- zmin
# exp_matrix[exp_matrix>zmax] <- zmax
# ColorRamp <- colorRampPalette(c("white","#d5f5e3","#2ecc71","#239b56","#186a3b"), bias=1)(10000)   #color list
# ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence
# image(1, 1:length(exp_matrix[c(HH,HL,LH,LL)]), t(exp_matrix[c(HH,HL,LH,LL)]), xaxt="n", yaxt="n", col=ColorRamp, xlab="", ylab="");box(lwd=2)
# # axis(side=3,1,labels=c("8-cell"),cex.axis=1)
# abline(h=boundary,lty=1,lwd=2)
# image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
# axis(side=1,c(zmin,(zmin+zmax)/2,zmax),labels=c(round(zmin,2),round((zmin+zmax)/2,2),round(zmax,2)))
# dev.off()
