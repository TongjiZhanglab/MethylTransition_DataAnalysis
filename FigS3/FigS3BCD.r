###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
cccol05 <- c("#CE001305","#16557A05","#C7A60905","#87C23205","#00879205","#A14C9405","#15A08C05","#8B7E7505","#1E7CAF05","#EA425F05","#46489A05","#E5003305","#0F231F05","#1187CD05")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E8003380","#0F231F80","#1187CD80")

library(sm)
vioplot <- function(x,...,range=1.5,h=NULL,ylim=NULL,names=NULL, horizontal=FALSE,
  col="magenta", border="black", lty=1, lwd=1, rectCol="black", colMed="white", pchMed=19, at, add=FALSE, wex=1,
  drawRect=TRUE)
{
    # process multiple datas
    datas <- list(x,...)
    n <- length(datas)

    if(missing(at)) at <- 1:n

    # pass 1
    #
    # - calculate base range
    # - estimate density
    #

    # setup parameters for density estimation
    upper  <- vector(mode="numeric",length=n)
    lower  <- vector(mode="numeric",length=n)
    q1     <- vector(mode="numeric",length=n)
    q3     <- vector(mode="numeric",length=n)
    med    <- vector(mode="numeric",length=n)
    base   <- vector(mode="list",length=n)
    height <- vector(mode="list",length=n)
    baserange <- c(Inf,-Inf)

    # global args for sm.density function-call
    args <- list(display="none")

    if (!(is.null(h)))
        args <- c(args, h=h)

    for(i in 1:n) {
        data<-datas[[i]]

        # calculate plot parameters
        #   1- and 3-quantile, median, IQR, upper- and lower-adjacent
        data.min <- min(data)
        data.max <- max(data)
        q1[i]<-quantile(data,0.25,na.rm=TRUE)
        q3[i]<-quantile(data,0.75,na.rm=TRUE)
        med[i]<-median(data)
        iqd <- q3[i]-q1[i]
        upper[i] <- min( q3[i] + range*iqd, data.max )
        lower[i] <- max( q1[i] - range*iqd, data.min )

        #   strategy:
        #       xmin = min(lower, data.min))
        #       ymax = max(upper, data.max))
        #

        est.xlim <- c( min(lower[i], data.min), max(upper[i], data.max) )

        # estimate density curve
        smout <- do.call("sm.density", c( list(data, xlim=est.xlim), args ) )

        # calculate stretch factor
        #
        #  the plots density heights is defined in range 0.0 ... 0.5
        #  we scale maximum estimated point to 0.4 per data
        #
        hscale <- 0.4/max(smout$estimate) * wex

        # add density curve x,y pair to lists
        base[[i]]   <- smout$eval.points
        height[[i]] <- smout$estimate * hscale

        # calculate min,max base ranges
        t <- range(base[[i]])
        baserange[1] <- min(baserange[1],t[1])
        baserange[2] <- max(baserange[2],t[2])

    }

    # pass 2
    #
    # - plot graphics

    # setup parameters for plot
    if(!add){
      xlim <- if(n==1)
               at + c(-.5, .5)
              else
               range(at) + min(diff(at))/2 * c(-1,1)

      if (is.null(ylim)) {
         ylim <- baserange
      }
    }
    if (is.null(names)) {
        label <- 1:n
    } else {
        label <- names
    }

    boxwidth <- 0.05 * wex

    # setup plot
    if(!add)
      plot.new()
    if(!horizontal) {
      if(!add){
        plot.window(xlim = xlim, ylim = ylim)
        axis(2)
        axis(1,at = at, label=label )
      }

      box()
      for(i in 1:n) {
          # plot left/right density curve
          polygon( c(at[i]-height[[i]], rev(at[i]+height[[i]])),
                   c(base[[i]], rev(base[[i]])),
                   col = col[i], border=border, lty=lty, lwd=lwd)

          if(drawRect){
            # plot IQR
            lines( at[c( i, i)], c(lower[i], upper[i]) ,lwd=lwd, lty=lty)

            # plot 50% KI box
            rect( at[i]-boxwidth/2, q1[i], at[i]+boxwidth/2, q3[i], col=rectCol)

            # plot median point
            points( at[i], med[i], pch=pchMed, col=colMed )
         }
      }

    }
    else {
      if(!add){
        plot.window(xlim = ylim, ylim = xlim)
        axis(1)
        axis(2,at = at, label=label )
      }

      box()
      for(i in 1:n) {
          # plot left/right density curve
          polygon( c(base[[i]], rev(base[[i]])),
                   c(at[i]-height[[i]], rev(at[i]+height[[i]])),
                   col = col[i], border=border, lty=lty, lwd=lwd)

          if(drawRect){
            # plot IQR
            lines( c(lower[i], upper[i]), at[c(i,i)] ,lwd=lwd, lty=lty)

            # plot 50% KI box
            rect( q1[i], at[i]-boxwidth/2, q3[i], at[i]+boxwidth/2,  col=rectCol)

            # plot median point
            points( med[i], at[i], pch=pchMed, col=colMed )
          }
      }
    }
    invisible (list( upper=upper, lower=lower, median=med, q1=q1, q3=q3))
}

library(ineq)
BinaryHeterogeneity <- function(x){
    cutoff_1 <- 1/8
	cutoff_2 <- 3/8
	cutoff_3 <- 5/8
	cutoff_4 <- 7/8
	a1 <- length(which(x<=cutoff_1))
	a2 <- length(which(x>cutoff_1 & x<=cutoff_2))
	a3 <- length(which(x>cutoff_2 & x<=cutoff_3))
	a4 <- length(which(x>cutoff_3 & x<=cutoff_4))
	a5 <- length(which(x>cutoff_4))
	if (a1+a2+a3+a4+a5 <= 2){
		return(NA)
	}else{
		total <- length(x)
		return(1-Gini(c(a1/total,a2/total,a3/total,a4/total,a5/total),corr=T))	
	}
}

MethylationRatio2Class <- function(x){
	cutoff_1 <- 1/8
	cutoff_2 <- 3/8
	cutoff_3 <- 5/8
	cutoff_4 <- 7/8
	x[x<=cutoff_1] <- 0
	x[x>cutoff_1 & x<=cutoff_2] <- 1/4
	x[x>cutoff_2 & x<=cutoff_3] <- 1/2
	x[x>cutoff_3 & x<=cutoff_4] <- 3/4
	x[x>cutoff_4] <- 1
	return(x)
}

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################

alldata <- read.table("../Data/hg_MethylationRatio.txt",header=T,row.names=1)
Methylation <- alldata[,c("hg_8C_9_1","hg_8C_9_2","hg_8C_9_3","hg_8C_9_4","hg_8C_9_5","hg_8C_9_6","hg_8C_9_7","hg_8C_9_8")]
out_name <- "SFig2"
Methylation_class <- apply(Methylation,2,MethylationRatio2Class)
Methylation_ratio <- apply(Methylation,1,mean,na.rm=T)
Methylation_heter <- apply(Methylation,1,BinaryHeterogeneity)

###################################################################################################
################################                PLOT               ################################
###################################################################################################

pdf(paste(out_name,"B.MethylationRatioHeterogeneityCorrelation.pdf",sep=""),width=3.6,height=4)
plot(Methylation_ratio,Methylation_heter,type="p",pch=".",ylim=c(0,1),xlim=c(0,1),col=cccol[7],main="",xlab="Average methylation level",ylab="Methylation heterogeneity score");box(lwd=2)
# legend("topright",c(paste("cor =",round(cor(Methylation_ratio,Methylation_heter,use="pair"),2))),bty="n")
# dev.off()

h1_gene <- which(Methylation_heter<0.2 & Methylation_heter>=0)
h2_gene <- which(Methylation_heter<0.4 & Methylation_heter>=0.2)
h3_gene <- which(Methylation_heter<0.6 & Methylation_heter>=0.4)
h4_gene <- which(Methylation_heter<0.8 & Methylation_heter>=0.6)
h5_gene <- which(Methylation_heter<=1 & Methylation_heter>=0.8)
vioplot(Methylation_ratio[h1_gene],Methylation_ratio[h2_gene],Methylation_ratio[h3_gene],Methylation_ratio[h4_gene],Methylation_ratio[h5_gene],col=cccol80[1:5],border=cccol[8],lwd=2, rectCol="black", colMed="white",names=c(0.2,0.4,0.6,0.8,1.0),horizontal=T);box(lwd=2)
dev.off()

pdf(paste(out_name,"C.MeanDistribution.pdf",sep=""),width=3.6,height=4)
hist(Methylation_ratio,breaks=10,col=cccol[8],border="white",main="Average methylation ratio",xlab="Methylation Ratio");box(lwd=2)
dev.off()

pdf(paste(out_name,"D.BinaryHeterogeneityDistribution.pdf",sep=""),width=3.6,height=4)
hist(Methylation_heter,breaks=10,col=cccol[9],border="white",main="Methylation heterogeneity",xlab="Methylation Heterogeneity");box(lwd=2)
dev.off()
