###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)

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

Boxplot_List <- function(plot_list,MAIN,NAMES,YLAB,COL){
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
	boxplot(plot_list,at=0:(n-1)+1,boxwex=0.5, xlim=c(0.5,n+0.5),ylim=ylim,col=paste(COL,"50",sep=""),border=COL,outline=F,names=NAMES,main=MAIN,lwd=2,ylab=YLAB,las=2)
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
				text(x=xpos[i], y=max(ypos_1[i], ypos_2[i])+dist, label=mark[i], col=cccol[1])
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
######## methylation in Promoter
mm_methyclass <- read.table("../Data/mm_DNAMethylationLevel_Promoter.txt",header=T,row.names=1)
mm_genes <- rownames(mm_methyclass)
mm_all_cells <- colnames(mm_methyclass)
# selected 1cell :
# 				"Zygote_1"
# selected 2cell :
# 				"C2_e1_1","C2_e1_2"
# selected 4cell :
# 				"C4_e6_1","C4_e6_2","C4_e6_3","C4_e6_4"
# selected 8cell :
# 				"C8_e3_1","C8_e3_2","C8_e3_3","C8_e3_4","C8_e3_5","C8_e3_6","C8_e3_7","C8_e3_8"
# selected Morular :
# 				"Morular_e1_1","Morular_e1_2","Morular_e1_3","Morular_e1_4","Morular_e1_5","Morular_e1_6"
mm_Sperm_id <- mm_all_cells[grep("Sperm", mm_all_cells)]
mm_Oocyte_id <- mm_all_cells[grep("Oocyte", mm_all_cells)]
mm_Zygote_id <- mm_all_cells[grep("Zygote", mm_all_cells)]
# mm_Zygote_id <- "Zygote_1"
mm_c2_id <- c("C2_e1_1","C2_e1_2")
mm_c4_id <- c("C4_e6_1","C4_e6_2","C4_e6_3","C4_e6_4")
mm_c8_id <- c("C8_e3_1","C8_e3_2","C8_e3_3","C8_e3_4","C8_e3_5","C8_e3_6","C8_e3_7","C8_e3_8")
mm_ml_id <- c("Morular_e1_1","Morular_e1_2","Morular_e1_3","Morular_e1_4","Morular_e1_5","Morular_e1_6")
# # DNA methyaltion heterogeneity
# mm_2C_1_Gini_Promoter <- apply(mm_methyclass[,mm_c2_id],1,MethylationHeterogeneity)
# mm_4C_6_Gini_Promoter <- apply(mm_methyclass[,mm_c4_id],1,MethylationHeterogeneity)
# mm_8C_3_Gini_Promoter <- apply(mm_methyclass[,mm_c8_id],1,MethylationHeterogeneity)
# mm_ml_1_Gini_Promoter <- apply(mm_methyclass[,mm_ml_id],1,MethylationHeterogeneity)

# Average DNA methyaltion in zygotes
mm_zygote_avg_Promoter <- apply(mm_methyclass[,mm_Zygote_id],1,mean,na.rm=T)
# mm_zygote_avg_Promoter <- mm_methyclass[,mm_Zygote_id];names(mm_zygote_avg_Promoter) <- mm_genes

mm_c4_avg_Promoter <- apply(mm_methyclass[,mm_c4_id],1,mean,na.rm=T)
mm_c8_avg_Promoter <- apply(mm_methyclass[,mm_c8_id],1,mean,na.rm=T)


# ######## methylation in Enhancer
# mm_methyclass <- read.table("../Data/mm_DNAMethylationLevel_Enhancer.txt",header=T,row.names=1)
# mm_all_cells <- colnames(mm_methyclass)
# # selected 1cell :
# # 				"Zygote_1"
# # selected 2cell :
# # 				"C2_e1_1","C2_e1_2"
# # selected 4cell :
# # 				"C4_e6_1","C4_e6_2","C4_e6_3","C4_e6_4"
# # selected 8cell :
# # 				"C8_e3_1","C8_e3_2","C8_e3_3","C8_e3_4","C8_e3_5","C8_e3_6","C8_e3_7","C8_e3_8"
# # selected Morular :
# # 				"Morular_e1_1","Morular_e1_2","Morular_e1_3","Morular_e1_4","Morular_e1_5","Morular_e1_6"
# mm_Sperm_id <- mm_all_cells[grep("Sperm", mm_all_cells)]
# mm_Oocyte_id <- mm_all_cells[grep("Oocyte", mm_all_cells)]
# mm_Zygote_id <- mm_all_cells[grep("Zygote", mm_all_cells)]
# mm_c2_id <- c("C2_e1_1","C2_e1_2")
# mm_c4_id <- c("C4_e6_1","C4_e6_2","C4_e6_3","C4_e6_4")
# mm_c8_id <- c("C8_e3_1","C8_e3_2","C8_e3_3","C8_e3_4","C8_e3_5","C8_e3_6","C8_e3_7","C8_e3_8")
# mm_ml_id <- c("Morular_e1_1","Morular_e1_2","Morular_e1_3","Morular_e1_4","Morular_e1_5","Morular_e1_6")
# 	# DNA methyaltion heterogeneity
# mm_2C_1_Gini_Enhancer <- apply(mm_methyclass[,mm_c2_id],1,MethylationHeterogeneity)
# mm_4C_6_Gini_Enhancer <- apply(mm_methyclass[,mm_c4_id],1,MethylationHeterogeneity)
# mm_8C_3_Gini_Enhancer <- apply(mm_methyclass[,mm_c8_id],1,MethylationHeterogeneity)
# mm_ml_1_Gini_Enhancer <- apply(mm_methyclass[,mm_ml_id],1,MethylationHeterogeneity)
# 	# Average DNA methyaltion in zygotes
# mm_zygote_avg_Enhancer <- apply(mm_methyclass[,mm_Zygote_id],1,mean,na.rm=T)


######## gene list
KO_8C <- as.vector(read.table("8C.KO_ExpHeterogeneityGenes.txt")[,1])
WT_8C <- as.vector(read.table("8C.WT_ExpHeterogeneityGenes.txt")[,1])

# KO_4C <- as.vector(read.table("4C.KO_ExpHeterogeneityGenes.txt")[,1])
# WT_4C <- as.vector(read.table("4C.WT_ExpHeterogeneityGenes.txt")[,1])

###################################################################################################
################################                PLOT               ################################
###################################################################################################
pdf("Fig5C.ExpHeterChange_DNAmethylationBoxplot_Promoter.pdf",width=2.5,height=4.5)
par(mar=c(6,4,4,2))
plot_list <- list(mm_zygote_avg_Promoter[c(KO_8C)],mm_zygote_avg_Promoter[setdiff(mm_genes,KO_8C)])
Boxplot_List(plot_list,"",c("Dppa3-KO","Others"),"DNA methylation state",c(cccol[1],cccol[8]))
dev.off()

# pdf("Fig5C.ExpHeterChange_DNAmethylationBoxplot_Promoter.48C.pdf",width=5,height=4.5)
# par(mfrow=c(1,2),mar=c(6,4,4,2))
# plot_list <- list(mm_c4_avg_Promoter[KO_8C],mm_c4_avg_Promoter[random],mm_c4_avg_Promoter[WT_8C])
# Boxplot_List(plot_list,"4-cell",c("Dppa3-KO","Control","WT"),"DNA methylation state",c(cccol[1],cccol[8],cccol[2]))
# plot_list <- list(mm_c8_avg_Promoter[KO_8C],mm_c8_avg_Promoter[random],mm_c8_avg_Promoter[WT_8C])
# Boxplot_List(plot_list,"8-cell",c("Dppa3-KO","Control","WT"),"DNA methylation state",c(cccol[1],cccol[8],cccol[2]))
# dev.off()