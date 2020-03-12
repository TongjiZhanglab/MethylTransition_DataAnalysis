###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#64C0AB50","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
library(gplots)
# library(beeswarm)

cccol <- colorRampPalette(c(cccol[6],cccol[7]), bias=1)(5)
cccol50 <- paste(cccol,"50",sep="")

library(ineq)
MethylationHetergenityObservation <- function(x){
	a1 <- length(which(x==0))
	a2 <- length(which(x==1/4))
	a3 <- length(which(x==1/2))
	a4 <- length(which(x==3/4))
	a5 <- length(which(x==1))
	total <- a1+a2+a3+a4+a5
	return(1-Gini(c(a1/total,a2/total,a3/total,a4/total,a5/total),corr=T))
}

Boxplot_List <- function(plot_list,MAIN,NAMES,YLAB){
	n <- length(plot_list)
	xpos <- 0:(n-1)+1.5
	p_value <- c()
	for (i in seq(n-1)){
		p_value <- c(p_value,t.test(plot_list[[i]],plot_list[[i+1]])$p.value)
	}
	mark <- symnum(p_value, cutpoints=c(0,0.001,0.01,0.05,1), symbols=c("***","**","*","-"))
	bp <- boxplot(plot_list,at=0:(n-1)+1, xlim=c(0.5,n+0.5),outline=F,plot=F)
	ylim <- range(bp$stats,na.rm=T)
	dist <- (ylim[2]-ylim[1])/20
	ylim[2] <- ylim[2]+1*dist
	boxplot(plot_list,at=0:(n-1)+1, xlim=c(0.5,n+0.5),ylim=c(0,13),col=cccol50,border=cccol,outline=F,names=NAMES,main=MAIN,lwd=2,ylab=YLAB)
	box(lwd=2)
	ypos_1 <- bp$stats[5,][1:n-1]
	ypos_2 <- bp$stats[5,][2:n]
	for(i in 1:length(mark)){
		if(!is.na(mark[i])){
			segments(xpos[i]-.4, ypos_1[i]+dist/2, xpos[i]-.4, max(ypos_1[i], ypos_2[i])+dist)
			segments(xpos[i]+.4, ypos_2[i]+dist/2, xpos[i]+.4, max(ypos_1[i], ypos_2[i])+dist)
			segments(xpos[i]-.4, max(ypos_1[i], ypos_2[i])+dist, xpos[i]-0.2, max(ypos_1[i], ypos_2[i])+dist)
			segments(xpos[i]+.4, max(ypos_1[i], ypos_2[i])+dist, xpos[i]+0.2, max(ypos_1[i], ypos_2[i])+dist)
			text(x=xpos[i], y=max(ypos_1[i], ypos_2[i])+dist, label=mark[i], col=cccol[1])
		}
	}
}

Boxplot_List_2 <- function(plot_list,MAIN,NAMES,XLAB,YLAB){
	cccol <- colorRampPalette(c("#A14C94","#1E7CAF"), bias=1)(10)
	cccol50 <- paste(cccol,"50",sep="")
	n <- length(plot_list)
	xpos <- 0:(n-1)+1.5
	bp <- boxplot(plot_list,at=0:(n-1)+1, xlim=c(0.5,n+0.5),outline=F,plot=F)
	ylim <- range(bp$stats,na.rm=T)
	dist <- (ylim[2]-ylim[1])/20
	ylim[2] <- ylim[2]+3*dist
	boxplot(plot_list,at=0:(n-1)+1,boxwex=0.7,xlim=c(0.5,n+0.5),ylim=ylim,col=cccol50,border=cccol,outline=F,names=NAMES,main=MAIN,lwd=2,ylab=YLAB,xlab=XLAB,las=2)
	box(lwd=2)
}
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
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


sample <- "hg_LPN_8"

hg_methyclass <- read.table("../Data/hg_MethylationClass.txt",header=T,row.names=1)
hg_genes <- row.names(hg_methyclass)
hg_all_cells <- colnames(hg_methyclass)
hg_GV_id <- grep("GV", hg_all_cells)
hg_MII_id <- grep("MII", hg_all_cells)
hg_Sp_id <- grep("Sp", hg_all_cells)
hg_EPN_id <- grep("EPN", hg_all_cells)
hg_MPN_id <- grep("MPN", hg_all_cells)
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

c1_gene <- intersect(hg_genes[which(hg_methyclass[,sample]=="0")],rownames(CV2))
c2_gene <- intersect(hg_genes[which(hg_methyclass[,sample]=="0.25")],rownames(CV2))
c3_gene <- intersect(hg_genes[which(hg_methyclass[,sample]=="0.5")],rownames(CV2))
c4_gene <- intersect(hg_genes[which(hg_methyclass[,sample]=="0.75")],rownames(CV2))
c5_gene <- intersect(hg_genes[which(hg_methyclass[,sample]=="1")],rownames(CV2))

hg_2C_6_Gini <- apply(hg_methyclass[,hg_2C_6_id],1,MethylationHetergenityObservation)
hg_4C_7_Gini <- apply(hg_methyclass[,hg_4C_7_id],1,MethylationHetergenityObservation)
hg_4C_8_Gini <- apply(hg_methyclass[,hg_4C_8_id],1,MethylationHetergenityObservation)
hg_8C_6_Gini <- apply(hg_methyclass[,hg_8C_6_id],1,MethylationHetergenityObservation)
hg_8C_9_Gini <- apply(hg_methyclass[,hg_8C_9_id],1,MethylationHetergenityObservation)

###################################################################################################
################################                PLOT               ################################
###################################################################################################

pdf("Fig5A.MethyHeteroVSGeneHetero.pdf",width=4,height=5)
par(mar=c(6,4,4,2))
mh_0 <- intersect(rownames(CV2),names(which(hg_8C_9_Gini <= 0.1)))
mh_1 <- intersect(rownames(CV2),names(which(hg_8C_9_Gini <= 0.2 & hg_8C_9_Gini > 0.1)))
mh_2 <- intersect(rownames(CV2),names(which(hg_8C_9_Gini <= 0.3 & hg_8C_9_Gini > 0.2)))
mh_3 <- intersect(rownames(CV2),names(which(hg_8C_9_Gini <= 0.4 & hg_8C_9_Gini > 0.3)))
mh_4 <- intersect(rownames(CV2),names(which(hg_8C_9_Gini <= 0.5 & hg_8C_9_Gini > 0.4)))
mh_5 <- intersect(rownames(CV2),names(which(hg_8C_9_Gini <= 0.6 & hg_8C_9_Gini > 0.5)))
mh_6 <- intersect(rownames(CV2),names(which(hg_8C_9_Gini <= 0.7 & hg_8C_9_Gini > 0.6)))
mh_7 <- intersect(rownames(CV2),names(which(hg_8C_9_Gini <= 0.8 & hg_8C_9_Gini > 0.7)))
mh_8 <- intersect(rownames(CV2),names(which(hg_8C_9_Gini > 0.8)))

plot_list <- list(CV2[mh_0,"8C_e2"],CV2[mh_1,"8C_e2"],CV2[mh_2,"8C_e2"],CV2[mh_3,"8C_e2"],CV2[mh_4,"8C_e2"],CV2[mh_5,"8C_e2"],CV2[mh_6,"8C_e2"],CV2[mh_7,"8C_e2"],CV2[mh_8,"8C_e2"])
Boxplot_List_2(plot_list,"",c("0~0.1","0.1~0.2","0.2~0.3","0.3~0.4","0.4~0.5","0.5~0.6","0.6~0.7","0.7~0.8","0.8~1"),"DNA methylation heterogeneity","Expression heterogeneity")
legend("topright",bty="n",c(paste("cor = ",round(cor(as.numeric(CV2[intersect(rownames(CV2),names(hg_8C_9_Gini)),"8C_e2"]),as.numeric(hg_8C_9_Gini[intersect(rownames(CV2),names(hg_8C_9_Gini))]),use="pair"),2),sep="")))

plot(hg_8C_9_Gini[intersect(rownames(CV2),names(hg_8C_9_Gini))],CV2[intersect(rownames(CV2),names(hg_8C_9_Gini)),"8C_e2"],pch=".")
dev.off()

