###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)

library("MethylTransition")

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
initial_classes <- round(c(0.69081809,0.11193227,0.06824454,0.05159268,0.07741241),3)
# terminal_classes <- round(c(0.72506497,0.11394390,0.05968277,0.03701048,0.06429788),3)

decreasing <- c(0.85689778,0.60611431,0.39070830,0.37898089,0.35636969,0.08355969,0.18874612,0.24219345,0.19639066,0.14987510,0.02751966,0.09525919,0.15232292,0.16242038,0.14737719,0.01365261,0.05228179,0.11576542,0.10934183,0.09408826,0.01837026,0.05759858,0.09900990,0.15286624,0.25228976)
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

delta_de <- (equilibrium-decreasing)/9
for (i in seq(9)){
  set.seed(1)
  tmp_transition_matrix <- matrix(round(decreasing + delta_de*i,5),c(5,5),byrow=T)
  # methyl_ratio <- c(methyl_ratio,sum((tmp_transition_matrix%*%initial_classes)*seq(0,1,0.25)))
  transitionmatrix <- rbind(transitionmatrix,c(t(tmp_transition_matrix)))
  # parameters <- rbind(parameters,TMParameterEstimation(tmp_transition_matrix, iter = 100+i*20, cell_cycle = 1)$estimated_parameters)
}

delta_ei <- (increasing-equilibrium)/9
for (i in seq(9)){
  set.seed(1)
  tmp_transition_matrix <- matrix(round(equilibrium + delta_ei*i,5),c(5,5),byrow=T)
  # methyl_ratio <- c(methyl_ratio,sum((tmp_transition_matrix%*%initial_classes)*seq(0,1,0.25)))
  transitionmatrix <- rbind(transitionmatrix,c(t(tmp_transition_matrix)))
  # parameters <- rbind(parameters,TMParameterEstimation(tmp_transition_matrix, iter = 100+i*20, cell_cycle = 1)$estimated_parameters)
}

write.table(parameters,file="Fig2B.parameters.txt",col.names=F,row.names=F,quote=F,sep="\t")
write.table(methyl_ratio,file="Fig2B.methylratio.txt",col.names=F,row.names=F,quote=F,sep="\t")
write.table(transitionmatrix,file="Fig2B.transitionmatrix.txt",col.names=F,row.names=F,quote=F,sep="\t")

parameters <- read.table("Fig2B.parameters.txt")

pdf("Fig2B.StabilityDisturbing.pdf",width=6,height=5)
plot(1,type="n",xaxt="n",yaxt="n",xlab="",ylab="Parameter value",ylim=c(0,1),xlim=c(1,19));box(lwd=2)
abline(v=10,lty=2,lwd=2)
axis(side=1,at=c(1,10,19),labels=c("decreasing","equilibrium","increasing"))
axis(side=2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
points(seq(19),parameters[,1],lwd=3,col=cccol[1],type="b",pch=16)
points(seq(19),parameters[,2],lwd=3,col=cccol[2],type="b",pch=16)
points(seq(19),parameters[,3],lwd=3,col=cccol[3],type="b",pch=16)
legend("topleft",c("u","d","p"),col=cccol,pch=16,bty="n")
dev.off()

pdf("Fig2B.StabilityDisturbing_Legend.pdf",width=6,height=3)
ColorRamp <- colorRampPalette(c(cccol[2],cccol[1]), bias=1)(1000)
plot_vector <- methyl_ratio-methyl_ratio[10]
barplot(plot_vector,xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(-0.1,0.4),col=ColorRamp[round(methyl_ratio*1000)],space=0,border="NA")
dev.off()
