###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E8003380","#0F231F80","#1187CD80")

cccolBR <- c("#16557A","#443F60","#712A46","#A0152C","#CE0013")
cccolBR30 <- c("#16557A30","#443F6030","#712A4630","#A0152C30","#CE001330")
cccolBR50 <- c("#16557A50","#443F6050","#712A4650","#A0152C50","#CE001350")
cccolBR80 <- c("#16557A80","#443F6080","#712A4680","#A0152C80","#CE001380")

options(scipen = 200)

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
######## methylation
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
mm_Sperm_id <- grep("Sperm", mm_all_cells)
mm_Oocyte_id <- grep("Oocyte", mm_all_cells)
mm_Zygote_id <- "Zygote_1"
mm_c2_id <- c("C2_e1_1","C2_e1_2")
mm_c4_id <- c("C4_e6_1","C4_e6_2","C4_e6_3","C4_e6_4")
mm_c8_id <- c("C8_e3_1","C8_e3_2","C8_e3_3","C8_e3_4","C8_e3_5","C8_e3_6","C8_e3_7","C8_e3_8")
c1_gene <- mm_genes[which(mm_methyclass[,mm_Zygote_id]<1/8)]
c2_gene <- mm_genes[which(mm_methyclass[,mm_Zygote_id]>1/8 & mm_methyclass[,mm_Zygote_id]<=3/8)]
c3_gene <- mm_genes[which(mm_methyclass[,mm_Zygote_id]>3/8 & mm_methyclass[,mm_Zygote_id]<=5/8)]
c4_gene <- mm_genes[which(mm_methyclass[,mm_Zygote_id]>5/8 & mm_methyclass[,mm_Zygote_id]<=7/8)]
c5_gene <- mm_genes[which(mm_methyclass[,mm_Zygote_id]>7/8)]

######## expression
GSE45719_exp <- read.table("../Data/GSE45719_mmEmbryoFPKM.txt",header=T,row.names=1)
GSE65160_exp <- read.table("../Data/GSE65160_mmEmbryoFPKM.txt",header=T,row.names=1)
inhouse_exp <- read.table("../Data/inhouse_mmEmbryoFPKM_4C8C.txt",header=T,row.names=1)

# selected embyros
colnames(GSE45719_exp) <- paste("GSE45719",colnames(GSE45719_exp),sep="_")
colnames(GSE65160_exp) <- paste("GSE65160",colnames(GSE65160_exp),sep="_")
colnames(inhouse_exp) <- paste("inhouse",colnames(inhouse_exp),sep="_")
GSE45719_4C_e2 <- colnames(GSE45719_exp)[grep("c4_e2",colnames(GSE45719_exp))]
GSE45719_4C_e4 <- colnames(GSE45719_exp)[grep("c4_e4",colnames(GSE45719_exp))]
GSE45719_8C_e1 <- colnames(GSE45719_exp)[grep("c8_e1",colnames(GSE45719_exp))]
GSE45719_8C_e2 <- colnames(GSE45719_exp)[grep("c8_e2",colnames(GSE45719_exp))]
GSE45719_8C_e5 <- colnames(GSE45719_exp)[grep("c8_e5",colnames(GSE45719_exp))]
GSE45719_8C_e8 <- colnames(GSE45719_exp)[grep("c8_e8",colnames(GSE45719_exp))]

GSE65160_4C_e1 <- colnames(GSE65160_exp)[grep("c4_e1",colnames(GSE65160_exp))]
GSE65160_4C_e2 <- colnames(GSE65160_exp)[grep("c4_e2",colnames(GSE65160_exp))]
GSE65160_4C_e3 <- colnames(GSE65160_exp)[grep("c4_e3",colnames(GSE65160_exp))]
GSE65160_8C_e1 <- colnames(GSE65160_exp)[grep("c8_e1",colnames(GSE65160_exp))]
GSE65160_8C_e2 <- colnames(GSE65160_exp)[grep("c8_e2",colnames(GSE65160_exp))]
GSE65160_8C_e3 <- colnames(GSE65160_exp)[grep("c8_e3",colnames(GSE65160_exp))]

inhouse_4C_e1 <- colnames(inhouse_exp)[grep("C4_e1",colnames(inhouse_exp))]
inhouse_4C_e2 <- colnames(inhouse_exp)[grep("C4_e2",colnames(inhouse_exp))]
inhouse_4C_e3 <- colnames(inhouse_exp)[grep("C4_e3",colnames(inhouse_exp))]
inhouse_8C_e1 <- colnames(inhouse_exp)[grep("C8_e1",colnames(inhouse_exp))]
inhouse_8C_e2 <- colnames(inhouse_exp)[grep("C8_e2",colnames(inhouse_exp))]

GSE65160_exp_merge <- c()
for (i in seq(4)){GSE65160_exp_merge <- cbind(GSE65160_exp_merge,apply(GSE65160_exp[,GSE65160_4C_e1[(2*i-1):(2*i)]],1,mean))}
for (i in seq(4)){GSE65160_exp_merge <- cbind(GSE65160_exp_merge,apply(GSE65160_exp[,GSE65160_4C_e2[(2*i-1):(2*i)]],1,mean))}
for (i in seq(4)){GSE65160_exp_merge <- cbind(GSE65160_exp_merge,apply(GSE65160_exp[,GSE65160_4C_e3[(2*i-1):(2*i)]],1,mean))}
for (i in seq(8)){GSE65160_exp_merge <- cbind(GSE65160_exp_merge,apply(GSE65160_exp[,GSE65160_8C_e1[(2*i-1):(2*i)]],1,mean))}
for (i in seq(8)){GSE65160_exp_merge <- cbind(GSE65160_exp_merge,apply(GSE65160_exp[,GSE65160_8C_e2[(2*i-1):(2*i)]],1,mean))}
for (i in seq(8)){GSE65160_exp_merge <- cbind(GSE65160_exp_merge,apply(GSE65160_exp[,GSE65160_8C_e3[(2*i-1):(2*i)]],1,mean))}
colnames(GSE65160_exp_merge) <- paste("GSE65160",c("C4_e1_1","C4_e1_2","C4_e1_3","C4_e1_4","C4_e2_1","C4_e2_2","C4_e2_3","C4_e2_4","C4_e3_1","C4_e3_2","C4_e3_3","C4_e3_4","C8_e1_1","C8_e1_2","C8_e1_3","C8_e1_4","C8_e1_5","C8_e1_6","C8_e1_7","C8_e1_8","C8_e2_1","C8_e2_2","C8_e2_3","C8_e2_4","C8_e2_5","C8_e2_6","C8_e2_7","C8_e2_8","C8_e3_1","C8_e3_2","C8_e3_3","C8_e3_4","C8_e3_5","C8_e3_6","C8_e3_7","C8_e3_8"),"FPKM",sep="_")
GSE65160_exp <- GSE65160_exp_merge
GSE65160_4C_e1 <- colnames(GSE65160_exp)[grep("C4_e1",colnames(GSE65160_exp))]
GSE65160_4C_e2 <- colnames(GSE65160_exp)[grep("C4_e2",colnames(GSE65160_exp))]
GSE65160_4C_e3 <- colnames(GSE65160_exp)[grep("C4_e3",colnames(GSE65160_exp))]
GSE65160_8C_e1 <- colnames(GSE65160_exp)[grep("C8_e1",colnames(GSE65160_exp))]
GSE65160_8C_e2 <- colnames(GSE65160_exp)[grep("C8_e2",colnames(GSE65160_exp))]
GSE65160_8C_e3 <- colnames(GSE65160_exp)[grep("C8_e3",colnames(GSE65160_exp))]

GSE45719_exp <- GSE45719_exp[,c(GSE45719_4C_e2,GSE45719_4C_e4,GSE45719_8C_e1,GSE45719_8C_e2,GSE45719_8C_e5,GSE45719_8C_e8)]

genes <- intersect(intersect(rownames(GSE45719_exp),rownames(GSE65160_exp)),rownames(inhouse_exp))
exp <- cbind(GSE45719_exp[genes,],GSE65160_exp[genes,],inhouse_exp[genes,])
log_exp <- log2(exp+1)

# # normalize
# library(edgeR)
# batch <- as.factor(c(rep(1,ncol(GSE45719_exp)),rep(2,ncol(GSE65160_exp))))
# rm_WT <- removeBatchEffect(cbind(GSE45719_exp[genes,],GSE65160_exp[genes,]),batch=batch)
# rm_WT[rm_WT < 0] <- 0
# rm_exp <- cbind(rm_WT,inhouse_exp[genes,])
# rm_log_exp <- log2(rm_exp+1)

ConfirmedCF <- c("Carm1","Cdx2","Dab2","Eed","Eset","Fgf4","Fgfr2","Gata3","Gata6","Kdm6b","Klf5","Lats2","Lrp2","Nanog","Pou5f1","Setdb1","Smarca4","Smarcc1","Sox17","Sox2","Suv39h1","Taz","Tead4","Tfap2a","Tfap2c","Yap1")
markers <- c("Neat1","LincGET","Carm1","Cdx2","Dab2","Eed","Eset","Fgf4","Fgfr2","Gata3","Gata6","Kdm6b","Klf5","Lats2","Lrp2","Nanog","Pou5f1","Setdb1","Smarca4","Smarcc1","Sox17","Sox2","Suv39h1","Taz","Tead4","Tfap2a","Tfap2c","Yap1")


mean_matrix <- cbind(apply(exp[,GSE45719_4C_e2],1,mean,na.rm=T),
					 apply(exp[,GSE45719_4C_e4],1,mean,na.rm=T),
					 apply(exp[,GSE65160_4C_e1],1,mean,na.rm=T),
					 apply(exp[,GSE65160_4C_e2],1,mean,na.rm=T),
					 apply(exp[,GSE65160_4C_e3],1,mean,na.rm=T),
					 apply(exp[,inhouse_4C_e1],1,mean,na.rm=T),
					 apply(exp[,inhouse_4C_e2],1,mean,na.rm=T),
					 apply(exp[,inhouse_4C_e3],1,mean,na.rm=T))
var_matrix <- cbind(apply(exp[,GSE45719_4C_e2],1,var,na.rm=T),
					 apply(exp[,GSE45719_4C_e4],1,var,na.rm=T),
					 apply(exp[,GSE65160_4C_e1],1,var,na.rm=T),
					 apply(exp[,GSE65160_4C_e2],1,var,na.rm=T),
					 apply(exp[,GSE65160_4C_e3],1,var,na.rm=T),
					 apply(exp[,inhouse_4C_e1],1,var,na.rm=T),
					 apply(exp[,inhouse_4C_e2],1,var,na.rm=T),
					 apply(exp[,inhouse_4C_e3],1,var,na.rm=T))
CV2_4C <- var_matrix/(mean_matrix)^2
CV2_4C[mean_matrix == 0] <- 0
embryo_labels <- c("GSE45719_4C_e2","GSE45719_4C_e4","GSE65160_4C_e1","GSE65160_4C_e2","GSE65160_4C_e3","inhouse_4C_e1","inhouse_4C_e2","inhouse_4C_e3")
colnames(CV2_4C) <- embryo_labels
colnames(var_matrix) <- embryo_labels
colnames(mean_matrix) <- embryo_labels

mean_matrix <- cbind(apply(exp[,GSE45719_8C_e1],1,mean,na.rm=T),
					 apply(exp[,GSE45719_8C_e2],1,mean,na.rm=T),
					 apply(exp[,GSE45719_8C_e5],1,mean,na.rm=T),
					 apply(exp[,GSE45719_8C_e8],1,mean,na.rm=T),
					 apply(exp[,GSE65160_8C_e1],1,mean,na.rm=T),
					 apply(exp[,GSE65160_8C_e2],1,mean,na.rm=T),
					 apply(exp[,GSE65160_8C_e3],1,mean,na.rm=T),
					 apply(exp[,inhouse_8C_e1],1,mean,na.rm=T),
					 apply(exp[,inhouse_8C_e2],1,mean,na.rm=T))
var_matrix <- cbind(apply(exp[,GSE45719_8C_e1],1,var,na.rm=T),
					 apply(exp[,GSE45719_8C_e2],1,var,na.rm=T),
					 apply(exp[,GSE45719_8C_e5],1,var,na.rm=T),
					 apply(exp[,GSE45719_8C_e8],1,var,na.rm=T),
					 apply(exp[,GSE65160_8C_e1],1,var,na.rm=T),
					 apply(exp[,GSE65160_8C_e2],1,var,na.rm=T),
					 apply(exp[,GSE65160_8C_e3],1,var,na.rm=T),
					 apply(exp[,inhouse_8C_e1],1,var,na.rm=T),
					 apply(exp[,inhouse_8C_e2],1,var,na.rm=T))
CV2_8C <- var_matrix/(mean_matrix)^2
CV2_8C[mean_matrix == 0] <- 0
embryo_labels <- c("GSE45719_8C_e1","GSE45719_8C_e2","GSE45719_8C_e5","GSE45719_8C_e8","GSE65160_8C_e1","GSE65160_8C_e2","GSE65160_8C_e3","inhouse_8C_e1","inhouse_8C_e2")
colnames(CV2_8C) <- embryo_labels
colnames(var_matrix) <- embryo_labels
colnames(mean_matrix) <- embryo_labels

###################################################################################################
################################                PLOT               ################################
###################################################################################################
# red_col <- colorRampPalette(c(cccol[1],cccol[3]), bias=1)(6+2) # 4-cell
# blue_col <- colorRampPalette(c(cccol[2],cccol[4]), bias=1)(7+1) # 8-cell

# 4-cell
WT <- CV2_4C[,c("GSE45719_4C_e2","GSE45719_4C_e4","GSE65160_4C_e1","GSE65160_4C_e2","GSE65160_4C_e3")]
KO <- CV2_4C[,c("inhouse_4C_e1","inhouse_4C_e2","inhouse_4C_e3")]
log10pvalue <- c()
log2fc <- c()
for (each_gene in row.names(CV2_4C)){
	if (length(table(WT[each_gene,]))==1 & length(table(KO[each_gene,]))==1 ){
		if (table(KO[each_gene,])!=table(WT[each_gene,])){log10pvalue <- c(log10pvalue, -1)} else {log10pvalue <- c(log10pvalue, 0)}
	}else{
	log10pvalue <- c(log10pvalue, -log10(t.test(KO[each_gene,],WT[each_gene,])$p.value))
	}
	log2fc <- c(log2fc, log2(mean(KO[each_gene,])/mean(WT[each_gene,])))
}
log10pvalue[log10pvalue > 7] <- 7
log10pvalue[log10pvalue == -1] <- 7
names(log10pvalue) <- rownames(CV2_4C)
names(log2fc) <- rownames(CV2_4C)

# pdf("Fig5B.4C.ExpHeterogeneityChangeVolcanoPlot.pdf",width=4,height=4.5)
# par(mar=c(6,4,4,2))
# plot(log2fc,log10pvalue,pch=".",main="Expression heterogeneity change",xlim=c(-10,10),xlab="log2(KO/WT)",ylab="-log10(p-value)");box(lwd=2)
# abline(h=-log10(0.01), lty=2)
# abline(v=c(0), lty=2)
# KO_genes <- rownames(CV2_4C)[which(log2fc > 0 & log10pvalue >= -log10(0.01))]
# WT_genes <- rownames(CV2_4C)[which(log2fc < 0 & log10pvalue >= -log10(0.01))]
# points(log2fc[KO_genes],log10pvalue[KO_genes],pch=16,col=cccol[1],cex=0.3)
# # points(log2fc[WT_genes],log10pvalue[WT_genes],pch=16,col=cccol[2],cex=0.3)
# points(log2fc[intersect(KO_genes,markers)],log10pvalue[intersect(KO_genes,markers)],pch="*",col=cccol[2],cex=1)
# text(log2fc[intersect(KO_genes,markers)],log10pvalue[intersect(KO_genes,markers)],labels=intersect(KO_genes,markers),col=cccol[2],cex=0.5)
# dev.off()
# write.table(KO_genes,file="4C.KO_ExpHeterogeneityGenes.txt",col.names=F,row.names=F,quote=F)
# write.table(WT_genes,file="4C.WT_ExpHeterogeneityGenes.txt",col.names=F,row.names=F,quote=F)

# 8-cell
WT <- CV2_8C[,c("GSE45719_8C_e1","GSE45719_8C_e2","GSE45719_8C_e5","GSE45719_8C_e8","GSE65160_8C_e1","GSE65160_8C_e2","GSE65160_8C_e3")]
KO <- CV2_8C[,c("inhouse_8C_e1","inhouse_8C_e2")]
log10pvalue <- c()
log2fc <- c()
for (each_gene in row.names(CV2_8C)){
	if (length(table(WT[each_gene,]))==1 & length(table(KO[each_gene,]))==1 ){
		if (table(KO[each_gene,])!=table(WT[each_gene,])){log10pvalue <- c(log10pvalue, -1)} else {log10pvalue <- c(log10pvalue, 0)}
	}else{
	log10pvalue <- c(log10pvalue, -log10(t.test(KO[each_gene,],WT[each_gene,])$p.value))
	}
	log2fc <- c(log2fc, log2(mean(KO[each_gene,])/mean(WT[each_gene,])))
}
log10pvalue[log10pvalue > 7] <- 7
log10pvalue[log10pvalue == -1] <- 7
names(log10pvalue) <- rownames(CV2_8C)
names(log2fc) <- rownames(CV2_8C)

pdf("Fig5B.8C.ExpHeterogeneityChangeVolcanoPlot.pdf",width=4,height=4.5)
par(mar=c(6,4,4,2))
plot(log2fc,log10pvalue,pch=".",main="Expression heterogeneity change",xlim=c(-10,10),xlab="log2(KO/WT)",ylab="-log10(p-value)");box(lwd=2)
abline(h=-log10(0.001), lty=2)
abline(v=c(0), lty=2)
KO_genes <- rownames(CV2_8C)[which(log2fc > 0 & log10pvalue >= -log10(0.001))]
WT_genes <- rownames(CV2_8C)[which(log2fc < 0 & log10pvalue >= -log10(0.001))]
points(log2fc[KO_genes],log10pvalue[KO_genes],pch=16,col=cccol[1],cex=0.3)
points(log2fc[WT_genes],log10pvalue[WT_genes],pch=16,col=cccol[2],cex=0.3)
points(log2fc[intersect(KO_genes,markers)],log10pvalue[intersect(KO_genes,markers)],pch="*",col=cccol[1],cex=1)
text(log2fc[intersect(KO_genes,markers)],log10pvalue[intersect(KO_genes,markers)],labels=intersect(KO_genes,markers),col=cccol[4],cex=0.5)
text(c(8,-8),c(6.8,6.8),labels=paste("n = ",c(length(KO_genes),length(WT_genes)),sep=""),col=cccol,cex=0.8)
dev.off()

write.table(KO_genes,file="8C.KO_ExpHeterogeneityGenes.txt",col.names=F,row.names=F,quote=F)
write.table(WT_genes,file="8C.WT_ExpHeterogeneityGenes.txt",col.names=F,row.names=F,quote=F)


length(KO_genes)/length(WT_genes)

