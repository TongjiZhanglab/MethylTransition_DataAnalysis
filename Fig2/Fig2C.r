###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)

library("MethylTransition")

MethylationClassClass2ClassSummary <- function(x,y){
    cutoff_1 <- 1/8
    cutoff_2 <- 3/8
    cutoff_3 <- 5/8
    cutoff_4 <- 7/8
    a11 <- sum(x<=cutoff_1 & y<=cutoff_1,na.rm=T)
    a12 <- sum(x<=cutoff_1 & y>cutoff_1&y<=cutoff_2,na.rm=T)
    a13 <- sum(x<=cutoff_1 & y>cutoff_2&y<=cutoff_3,na.rm=T)
    a14 <- sum(x<=cutoff_1 & y>cutoff_3&y<=cutoff_4,na.rm=T)
    a15 <- sum(x<=cutoff_1 & y>cutoff_4,na.rm=T)
    a21 <- sum(x>cutoff_1&x<=cutoff_2 & y<=cutoff_1,na.rm=T)
    a22 <- sum(x>cutoff_1&x<=cutoff_2 & y>cutoff_1&y<=cutoff_2,na.rm=T)
    a23 <- sum(x>cutoff_1&x<=cutoff_2 & y>cutoff_2&y<=cutoff_3,na.rm=T)
    a24 <- sum(x>cutoff_1&x<=cutoff_2 & y>cutoff_3&y<=cutoff_4,na.rm=T)
    a25 <- sum(x>cutoff_1&x<=cutoff_2 & y>cutoff_4,na.rm=T)
    a31 <- sum(x>cutoff_2&x<=cutoff_3 & y<=cutoff_1,na.rm=T)
    a32 <- sum(x>cutoff_2&x<=cutoff_3 & y>cutoff_1&y<=cutoff_2,na.rm=T)
    a33 <- sum(x>cutoff_2&x<=cutoff_3 & y>cutoff_2&y<=cutoff_3,na.rm=T)
    a34 <- sum(x>cutoff_2&x<=cutoff_3 & y>cutoff_3&y<=cutoff_4,na.rm=T)
    a35 <- sum(x>cutoff_2&x<=cutoff_3 & y>cutoff_4,na.rm=T)
    a41 <- sum(x>cutoff_3&x<=cutoff_4 & y<=cutoff_1,na.rm=T)
    a42 <- sum(x>cutoff_3&x<=cutoff_4 & y>cutoff_1&y<=cutoff_2,na.rm=T)
    a43 <- sum(x>cutoff_3&x<=cutoff_4 & y>cutoff_2&y<=cutoff_3,na.rm=T)
    a44 <- sum(x>cutoff_3&x<=cutoff_4 & y>cutoff_3&y<=cutoff_4,na.rm=T)
    a45 <- sum(x>cutoff_3&x<=cutoff_4 & y>cutoff_4,na.rm=T)
    a51 <- sum(x>cutoff_4 & y<=cutoff_1,na.rm=T)
    a52 <- sum(x>cutoff_4 & y>cutoff_1&y<=cutoff_2,na.rm=T)
    a53 <- sum(x>cutoff_4 & y>cutoff_2&y<=cutoff_3,na.rm=T)
    a54 <- sum(x>cutoff_4 & y>cutoff_3&y<=cutoff_4,na.rm=T)
    a55 <- sum(x>cutoff_4 & y>cutoff_4,na.rm=T)
    a1_total <- a11+a12+a13+a14+a15
    a2_total <- a21+a22+a23+a24+a25
    a3_total <- a31+a32+a33+a34+a35
    a4_total <- a41+a42+a43+a44+a45
    a5_total <- a51+a52+a53+a54+a55
    # return(c(a11/a1_total,a12/a1_total,a13/a1_total,a14/a1_total,a15/a1_total,a21/a2_total,a22/a2_total,a23/a2_total,a24/a2_total,a25/a2_total,a31/a3_total,a32/a3_total,a33/a3_total,a34/a3_total,a35/a3_total,a41/a4_total,a42/a4_total,a43/a4_total,a44/a4_total,a45/a4_total,a51/a5_total,a52/a5_total,a53/a5_total,a54/a5_total,a55/a5_total))
    tmp_matrix <- matrix(c(a11/a1_total,a12/a1_total,a13/a1_total,a14/a1_total,a15/a1_total,a21/a2_total,a22/a2_total,a23/a2_total,a24/a2_total,a25/a2_total,a31/a3_total,a32/a3_total,a33/a3_total,a34/a3_total,a35/a3_total,a41/a4_total,a42/a4_total,a43/a4_total,a44/a4_total,a45/a4_total,a51/a5_total,a52/a5_total,a53/a5_total,a54/a5_total,a55/a5_total),c(5,5),byrow = F)
    return(as.vector(t(tmp_matrix)))
}

MethylationClassSummary <- function(x){
  cutoff_1 <- 1/8
  cutoff_2 <- 3/8
  cutoff_3 <- 5/8
  cutoff_4 <- 7/8
  a <- sum(x<=cutoff_1,na.rm=T)
  b <- sum(x>cutoff_1 & x<=cutoff_2,na.rm=T)
  c <- sum(x>cutoff_2 & x<=cutoff_3,na.rm=T)
  d <- sum(x>cutoff_3 & x<=cutoff_4,na.rm=T)
  e <- sum(x>cutoff_4,na.rm=T)
  total <- a+b+c+d+e
  return(c(a/total,b/total,c/total,d/total,e/total))
}

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
# selected 1cell : "hg_LPN_8"
# selected 2cell : "hg_2C_6_1"
######## methylation class
hg_methyclass <- read.table("../Data/hg_MethylationClass.txt",header=T,row.names=1)
hg_genes <- row.names(hg_methyclass)
hg_all_cells <- colnames(hg_methyclass)
selected_1C <- "hg_LPN_8"
selected_2C <- "hg_2C_6_1"

initial_classes <- round(MethylationClassSummary(hg_methyclass[,selected_1C]),3)
terminal_classes <- round(MethylationClassSummary(hg_methyclass[,selected_2C]),3)

decreasing <- MethylationClassClass2ClassSummary(hg_methyclass[,selected_1C],hg_methyclass[,selected_2C])
equilibrium <- c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1)
tmpa <- matrix(decreasing,c(5,5),byrow=T)

# u_max <- 0.5
# S1_limit <- u_max^4;S2_limit <- u_max^3;S3_limit <- u_max^2;S4_limit <- u_max^1;S5_limit <- u_max^0;
# S1_limit <- min(S1_limit,tmpa[1,1]);S2_limit <- min(S2_limit,tmpa[1,2]);S3_limit <- min(S3_limit,tmpa[1,3]);S4_limit <- min(S4_limit,tmpa[1,4]);S5_limit <- min(S5_limit,tmpa[1,5]);
# tmpb1 <- rev(tmpa[,1]);tmpb1[seq(4)] <- tmpb1[seq(4)] + (tmpb1[5] - S1_limit)/4; tmpb1[5] <- S1_limit
# tmpb2 <- rev(tmpa[,2]);tmpb2[seq(4)] <- tmpb2[seq(4)] + (tmpb2[5] - S2_limit)/4; tmpb2[5] <- S2_limit
# tmpb3 <- rev(tmpa[,3]);tmpb3[seq(4)] <- tmpb3[seq(4)] + (tmpb3[5] - S3_limit)/4; tmpb3[5] <- S3_limit
# tmpb4 <- rev(tmpa[,4]);tmpb4[seq(4)] <- tmpb4[seq(4)] + (tmpb4[5] - S4_limit)/4; tmpb4[5] <- S4_limit
# tmpb5 <- rev(tmpa[,5]);tmpb5[seq(4)] <- tmpb5[seq(4)] + (tmpb5[5] - S5_limit)/4; tmpb5[5] <- S5_limit
# increasing <- round(c(t(cbind(tmpb1,tmpb2,tmpb3,tmpb4,tmpb5))),5)

increasing <- c(tmpa/apply(tmpa,1,sum))
# matrix(decreasing,c(5,5),byrow=T)
# matrix(increasing,c(5,5),byrow=T)
###################################################################################################
################################                PLOT               ################################
###################################################################################################
parameters <- c()
parameters <- rbind(parameters,TMParameterEstimation(matrix(decreasing,c(5,5),byrow=T), iter = 50, cell_cycle = 1)$estimated_parameters)
methyl_ratio <- c()
methyl_ratio <- c(methyl_ratio,sum((matrix(decreasing,c(5,5),byrow=T)%*%initial_classes)*seq(0,1,0.25)))
transitionmatrix <- c()
transitionmatrix <- rbind(transitionmatrix,round(decreasing,5))

delta_de <- (equilibrium-decreasing)/10
for (i in seq(10)){
  set.seed(1)
  tmp_transition_matrix <- matrix(round(decreasing + delta_de*i,5),c(5,5),byrow=T)
  methyl_ratio <- c(methyl_ratio,sum((tmp_transition_matrix%*%initial_classes)*seq(0,1,0.25)))
  transitionmatrix <- rbind(transitionmatrix,c(t(tmp_transition_matrix)))
  parameters <- rbind(parameters,TMParameterEstimation(tmp_transition_matrix, iter = 100+i*20, cell_cycle = 1)$estimated_parameters)
}

delta_ei <- (increasing-equilibrium)/10
for (i in seq(10)){
  set.seed(1)
  tmp_transition_matrix <- matrix(round(equilibrium + delta_ei*i,5),c(5,5),byrow=T)
  methyl_ratio <- c(methyl_ratio,sum((tmp_transition_matrix%*%initial_classes)*seq(0,1,0.25)))
  transitionmatrix <- rbind(transitionmatrix,c(t(tmp_transition_matrix)))
  parameters <- rbind(parameters,TMParameterEstimation(tmp_transition_matrix, iter = 100+i*20, cell_cycle = 1)$estimated_parameters)
}

write.table(parameters,file="Fig2C.parameters.txt",col.names=F,row.names=F,quote=F,sep="\t")
write.table(methyl_ratio,file="Fig2C.methylratio.txt",col.names=F,row.names=F,quote=F,sep="\t")
write.table(transitionmatrix,file="Fig2C.transitionmatrix.txt",col.names=F,row.names=F,quote=F,sep="\t")

parameters <- read.table("Fig2C.parameters.txt")

pdf("Fig2C.StabilityDisturbing.pdf",width=6,height=5)
plot(1,type="n",xaxt="n",yaxt="n",xlab="",ylab="Parameter value",ylim=c(0,1),xlim=c(1,21));box(lwd=2)
abline(v=11,lty=2,lwd=2)
axis(side=1,at=c(1,11,21),labels=c("decreasing","equilibrium","increasing"))
axis(side=2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
points(seq(21),parameters[,1],lwd=3,col=cccol[1],type="b",pch=16)
points(seq(21),parameters[,2],lwd=3,col=cccol[2],type="b",pch=16)
points(seq(21),parameters[,3],lwd=3,col=cccol[3],type="b",pch=16)
legend("topleft",c("u","d","p"),col=cccol,pch=16,bty="n")
dev.off()

pdf("Fig2C.StabilityDisturbing_Legend.pdf",width=6,height=3)
ColorRamp <- colorRampPalette(c(cccol[2],cccol[1]), bias=1)(1000)
plot_vector <- methyl_ratio-methyl_ratio[10]
barplot(plot_vector,xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(-0.1,0.4),col=ColorRamp[round(methyl_ratio*1000)],space=0,border="NA")
dev.off()
