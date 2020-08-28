###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
library(numDeriv)

sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#64C0AB","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#64C0AB50","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")

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

newton <- function(func = objfun, x0, tol = 1e-5, n.max = 100,...){
    x <- x0
    g <- grad(func, x, ...)
    h <- hessian(func, x, ...)

    n <- 0
    while( max(abs(g))>tol && n<n.max ){
        x <- x-solve(h,g)
        g <- grad(func, x, ...)
        h <- hessian(func, x, ...)
        n <- n+1
    }
    if(n == n.max){
        cat('newton failed to converge\n')
        return(x)
    }
    return(x)
}

fun2 <- function(x,observe){
	# x[1] : u
	# x[2] : d
	# x[3] : p
	predict <- c((-x[1] + 1)**4,
			0.5*(-x[1] + 1)**2*(0.5*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.5*(-x[1] + 1)**2) + 0.25*(-x[1] + 1)**2*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2),
			(1/3)*x[2]*(-x[3] + 1)*(-x[1] + 1)**3 + (2/3)*(0.25*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.25*(-x[1] + 1)**2)*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2),
			1.0*x[2]*(-x[3] + 1)*(-x[1] + 1)*(0.25*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.25*(-x[1] + 1)**2) + 0.25*x[2]*(-x[3] + 1)*(-x[1] + 1)*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2),
			1.0*x[2]**2*(-x[3] + 1)**2*(-x[1] + 1)**2,
			4*x[1]*(-x[1] + 1)**3,
			1.0*x[1]*(-x[1] + 1)*(0.5*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.5*(-x[1] + 1)**2) + 0.5*x[1]*(-x[1] + 1)*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2) + 0.5*(-x[1] + 1)**2*(0.5*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.5*x[1]*(-x[1] + 1)) + 0.25*(-x[1] + 1)**2*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1)) + 0.5*(-x[1] + 1)**2*(0.5*x[1]*(-x[1] + 1) + 0.5*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 0.25*(-x[1] + 1)**2*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)),
			(2/3)*x[2]*x[1]*(-x[3] + 1)*(-x[1] + 1)**2 + (1/3)*(-x[1] + 1)**2*(0.5*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.5*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (1/6)*(-x[1] + 1)**2*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (2/3)*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*x[1]*(-x[1] + 1))*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2) + (2/3)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1))*(0.25*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.25*(-x[1] + 1)**2) + (2/3)*(0.25*x[1]*(-x[1] + 1) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2) + (2/3)*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(0.25*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.25*(-x[1] + 1)**2),
			1.0*x[2]*(-x[3] + 1)*(-x[1] + 1)*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*x[1]*(-x[1] + 1)) + 0.25*x[2]*(-x[3] + 1)*(-x[1] + 1)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1)) + 1.0*x[2]*(-x[3] + 1)*(-x[1] + 1)*(0.25*x[1]*(-x[1] + 1) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 0.25*x[2]*(-x[3] + 1)*(-x[1] + 1)*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 1.0*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2) + 1.0*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(0.25*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.25*(-x[1] + 1)**2),
			4*x[2]*(-x[3] + 1)*(-x[1] + 1)*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 1.0*x[2]*(-x[3] + 1)*(-x[1] + 1)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)),
			6*x[1]**2*(-x[1] + 1)**2,
			0.5*x[1]**2*(0.5*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.5*(-x[1] + 1)**2) + 0.25*x[1]**2*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2) + 1.0*x[1]*(-x[1] + 1)*(0.5*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.5*x[1]*(-x[1] + 1)) + 0.5*x[1]*(-x[1] + 1)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1)) + 1.0*x[1]*(-x[1] + 1)*(0.5*x[1]*(-x[1] + 1) + 0.5*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 0.5*x[1]*(-x[1] + 1)*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 0.5*(-x[1] + 1)**2*(0.5*x[1]**2 + 0.5*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])) + 0.25*(-x[1] + 1)**2*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])),
			(1/3)*x[2]*x[1]**2*(-x[3] + 1)*(-x[1] + 1) + (2/3)*x[1]*(-x[1] + 1)*(0.5*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.5*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (1/3)*x[1]*(-x[1] + 1)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (1/3)*(-x[2] + 1)*(-x[1] + 1)**2*(-x[3]*x[1] + x[3] + x[1]) + (2/3)*(0.25*x[1]**2 + 0.25*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2) + (2/3)*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(0.25*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.25*(-x[1] + 1)**2) + (2/3)*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*x[1]*(-x[1] + 1))*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1)) + (2/3)*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*x[1]*(-x[1] + 1))*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (2/3)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1))*(0.25*x[1]*(-x[1] + 1) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (2/3)*(0.25*x[1]*(-x[1] + 1) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)),
			1.0*x[2]*(-x[3] + 1)*(-x[1] + 1)*(0.25*x[1]**2 + 0.25*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])) + 0.25*x[2]*(-x[3] + 1)*(-x[1] + 1)*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])) + 1.0*(-x[2] + 1)*(0.25*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.25*(-x[1] + 1)**2)*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2)*(-x[3]*x[1] + x[3] + x[1]) + 1.0*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*x[1]*(-x[1] + 1))*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 1.0*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1)) + 1.0*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 1.0*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(0.25*x[1]*(-x[1] + 1) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)),
			2.0*x[2]*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)*(-x[3]*x[1] + x[3] + x[1]) + 4*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)),
			4*x[1]**3*(-x[1] + 1),
			0.5*x[1]**2*(0.5*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.5*x[1]*(-x[1] + 1)) + 0.25*x[1]**2*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1)) + 0.5*x[1]**2*(0.5*x[1]*(-x[1] + 1) + 0.5*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 0.25*x[1]**2*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 1.0*x[1]*(-x[1] + 1)*(0.5*x[1]**2 + 0.5*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])) + 0.5*x[1]*(-x[1] + 1)*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])),
			(1/3)*x[1]**2*(0.5*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.5*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (1/6)*x[1]**2*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (2/3)*x[1]*(-x[2] + 1)*(-x[1] + 1)*(-x[3]*x[1] + x[3] + x[1]) + (2/3)*(0.25*x[1]**2 + 0.25*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1)) + (2/3)*(0.25*x[1]**2 + 0.25*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (2/3)*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*x[1]*(-x[1] + 1)) + (2/3)*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(0.25*x[1]*(-x[1] + 1) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)),
			1.0*(-x[2] + 1)*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*x[1]*(-x[1] + 1))*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1))*(-x[3]*x[1] + x[3] + x[1]) + 1.0*(-x[2] + 1)*(0.25*x[1]*(-x[1] + 1) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(-x[3]*x[1] + x[3] + x[1]) + 1.0*(0.25*x[1]**2 + 0.25*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 1.0*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)),
			4*(-x[2] + 1)*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(-x[3]*x[1] + x[3] + x[1]) + 1.0*(-x[2] + 1)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(-x[3]*x[1] + x[3] + x[1]),
			x[1]**4,
			0.5*x[1]**2*(0.5*x[1]**2 + 0.5*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])) + 0.25*x[1]**2*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])),
			(1/3)*x[1]**2*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]) + (2/3)*(0.25*x[1]**2 + 0.25*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])),
			1.0*(-x[2] + 1)*(0.25*x[1]**2 + 0.25*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(-x[3]*x[1] + x[3] + x[1]),
			1.0*(-x[2] + 1)**2*(-x[3]*x[1] + x[3] + x[1])**2)
	y <- 0
	for (each in seq(length(observe))){
		y <- y + (predict[each]-observe[each])^2

	}
	return (y)
}

ParameterEstimation <- function(observe,out_name,file_name){
	para_maxtrix <- c()
	i <- 1
	while(i<=50){
		x0 <- c(runif(1),runif(1),runif(1))
		tryCatch({y <- newton(func = fun2, x0=x0,observe=observe);
			para_maxtrix <- rbind(para_maxtrix,y);
			i <- i+1},
			error=function(e){i <- i+1})
	}
	selected_parameter_matrix <- para_maxtrix[which(para_maxtrix[,1]<=1 & para_maxtrix[,1]>=0 & para_maxtrix[,2]<=1 & para_maxtrix[,2]>=0 & para_maxtrix[,3]<=1 & para_maxtrix[,3]>=0),]
	# print(selected_parameter_matrix)
	# print(nrow(selected_parameter_matrix))
	if (!is.null(nrow(selected_parameter_matrix))){
		cost_value <- round(apply(selected_parameter_matrix,1,fun2,observe),6)
		optimal_para <- selected_parameter_matrix[which(cost_value==min(cost_value)),]
		if (!is.null(nrow(optimal_para))){
			selected_parameter <- round(apply(optimal_para,2,mean,na.rm=T),6)
		}else{
			selected_parameter <- round(optimal_para,6)
		}	
	}else{
		selected_parameter <- round(selected_parameter_matrix,6)
	}
	write(paste(out_name,selected_parameter[1],selected_parameter[2],selected_parameter[3],sep="\t"),file=file_name,append=TRUE)
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
######## methylation class
hg_methyclass <- read.table("../Data/hg_MethylationClass.txt",header=T,row.names=1)

hg_genes <- row.names(hg_methyclass)
hg_all_cells <- colnames(hg_methyclass)

hg_GV_id <- grep("GV", hg_all_cells)
hg_MII_id <- grep("MII", hg_all_cells)
hg_Sp_id <- grep("Sp", hg_all_cells)
hg_LPN_id <- grep("LPN", hg_all_cells)
hg_2C_id <- grep("2C_6", hg_all_cells)
hg_4C_id <- c(grep("4C_7", hg_all_cells),grep("4C_8", hg_all_cells))
hg_8C_id <- c(grep("8C_6", hg_all_cells),grep("8C_9", hg_all_cells))

runid <- as.numeric(sys_argv[1])
###################################################################################################
################################                PLOT               ################################
###################################################################################################
each_LPN <- "hg_LPN_8"
each_2C <- "hg_2C_6_1"
observe <- MethylationClassClass2ClassSummary(hg_methyclass[,each_LPN],hg_methyclass[,each_2C])
out_name <- paste(each_LPN,each_2C,"ALL",sep="_")
file_name <- paste("NewtonParametersDistribution_LPN2C_",each_LPN,"_",each_2C,"_run",runid,".txt",sep="")
ParameterEstimation(observe,out_name,file_name)
counts <- 1
while(counts <= 20){
	for (sample in seq(0.1,0.9,0.1)){	
		gene_part <- sample(hg_genes,length(hg_genes)*sample)
		observe <- MethylationClassClass2ClassSummary(hg_methyclass[gene_part,each_LPN],hg_methyclass[gene_part,each_2C])
		out_name <- paste(each_LPN,each_2C,counts,sep="_")
		ParameterEstimation(observe,out_name,file_name)
	}
	counts <- counts+1
}

# for (each_2C in hg_2C_id){
# 	for (each_4C in hg_4C_id){
# 		observe <- MethylationClassClass2ClassSummary(hg_methyclass[,each_2C],hg_methyclass[,each_4C])
# 		out_name <- paste(hg_all_cells[each_2C],hg_all_cells[each_4C],"ALL",sep="_")
# 		file_name <- paste("NewtonParametersDistribution_5_2C4C_",each_2C,"_",each_4C,".txt",sep="")
# 		write(paste(out_name,paste(observe,collapse="\t"),sep="\t"),file="Observation_2C4C.txt",append=TRUE)
# 		ParameterEstimation(observe,out_name,file_name)
# 		counts <- 1
# 		while(counts <=200){
# 			gene_part <- sample(hg_genes,length(hg_genes)*0.05)
# 			observe <- MethylationClassClass2ClassSummary(hg_methyclass[gene_part,each_2C],hg_methyclass[gene_part,each_4C])
# 			out_name <- paste(hg_all_cells[each_2C],hg_all_cells[each_4C],counts,sep="_")
# 			write(paste(out_name,paste(observe,collapse="\t"),sep="\t"),file="Observation_2C4C.txt",append=TRUE)
# 			ParameterEstimation(observe,out_name,file_name)
# 			counts <- counts+1
# 		}
# 	}
# } 

# for (each_4C in hg_4C_id){
# 	for (each_8C in hg_8C_id){
# 		observe <- MethylationClassClass2ClassSummary(hg_methyclass[,each_4C],hg_methyclass[,each_8C])
# 		out_name <- paste(hg_all_cells[each_4C],hg_all_cells[each_8C],"ALL",sep="_")
# 		file_name <- paste("NewtonParametersDistribution_5_4C8C_",each_4C,"_",each_8C,".txt",sep="")
# 		write(paste(out_name,paste(observe,collapse="\t"),sep="\t"),file="Observation_4C8C.txt",append=TRUE)
# 		ParameterEstimation(observe,out_name,file_name)
# 		counts <- 1
# 		while(counts <=200){
# 			gene_part <- sample(hg_genes,length(hg_genes)*0.05)
# 			observe <- MethylationClassClass2ClassSummary(hg_methyclass[gene_part,each_4C],hg_methyclass[gene_part,each_8C])
# 			out_name <- paste(hg_all_cells[each_4C],hg_all_cells[each_8C],counts,sep="_")
# 			write(paste(out_name,paste(observe,collapse="\t"),sep="\t"),file="Observation_4C8C.txt",append=TRUE)
# 			ParameterEstimation(observe,out_name,file_name)
# 			counts <- counts+1
# 		}
# 	}
# } 


# pdf("LPN_PCA.pdf",width=6,height=6)
# pca <- prcomp(t(na.omit(hg_methyclass[,hg_LPN_id])))
# plot(as.numeric(pca$x[,1]),as.numeric(pca$x[,2]),pch=16,col=cccol,main='PCA',
# 	xlab=paste('PC1','(',summary(pca)[[6]]['Proportion of Variance',1]*100,'%)',sep=''),
# 	ylab=paste('PC2','(',summary(pca)[[6]]['Proportion of Variance',2]*100,'%)',sep=''))
# box(lwd=2)
# text(as.numeric(pca$x[,1]),as.numeric(pca$x[,2]),hg_all_cells[hg_LPN_id],cex=0.5)
# dev.off()
