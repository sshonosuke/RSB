###-----------------------------------------------------------###
###        R code to replicate Figure 2 in main text          ###
###-----------------------------------------------------------###
##   This code requires "sim-20.RData" and "sim-50.RData" 
##        obtained by running 'Sim-regression.R' 


## Preparation 
p <- 15   # number of covariates 
R <- 300   # number of Monte Carlo replications 
S <- 8    # number of scenarios 

meth <- c("RSB", "SB", "PO", "ZIP", "NB", "ZINB")
M <- length(meth)

Res1 <- array(NA, c(S, M, 3))
Error1 <- array(NA, c(S, M, 3))
Res2 <- array(NA, c(S, M, 3))
Error2 <- array(NA, c(S, M, 3))
dimnames(Res1)[[3]] <- dimnames(Res2)[[3]] <- c("MSE1", "MSE2", "IS")
dimnames(Res1)[[2]] <- dimnames(Res2)[[2]] <- meth



# a=20
load("sim-20.RData")
Res1[,,1] <- 1000*apply(MSE_beta, c(2,3), mean)
Res1[,,2] <- 100*apply(MSE_mu, c(2,3), mean)
Res1[,,3] <- apply(IS, c(2,3), mean)
Error1[,,1] <- 1000*apply(MSE_beta, c(2,3), sd)/sqrt(R*p) 
Error1[,,2] <- 100*apply(MSE_mu, c(2,3), sd)/sqrt(R*p) 
Error1[,,3] <- apply(IS, c(2,3), sd)/sqrt(R*p) 


# a=50
load("sim-50.RData")
Res2[,,1] <- 1000*apply(MSE_beta, c(2,3), mean)
Res2[,,2] <- 100*apply(MSE_mu, c(2,3), mean)
Res2[,,3] <- apply(IS, c(2,3), mean)
Error2[,,1] <- 1000*apply(MSE_beta, c(2,3), sd)/sqrt(R*p) 
Error2[,,2] <- 100*apply(MSE_mu, c(2,3), sd)/sqrt(R*p) 
Error2[,,3] <- apply(IS, c(2,3), sd)/sqrt(R*p) 



## Figure 
meth <- c("RSB", "SB", "PO", "ZIP", "NB", "ZINB")
M <- length(meth)

bb <- seq(-0.1, 0.1, length=M)
col <- c(1, 1, 2, 2, 4, 4)
lty <- c(1, 2, 1, 2, 1, 2)
ww <- qnorm(0.995)     #  Multiplier for Monte Carlo error

Ylab <- c("MSE", "MSE", "interval score")
pos <- c("topleft", "topleft", "topleft")


# Range 
ran <- matrix(NA, 3, 2)
for(j in 1:3){
  ran[j,] <- range(Res1[,,j]+ww*Error1[,,j], Res1[,,j]-ww*Error1[,,j], 
                   Res2[,,j]+ww*Error2[,,j], Res2[,,j]-ww*Error2[,,j])
}
ran[1,2] <- 400
ran[2,2] <- 3000
ran[3,2] <- 200


# make pdf file 
pdf("Figure2.pdf", height=15, width=10, pointsize=15)
par(mfcol=c(3,2))
Main <- c("MSE (b=20)", "SMSE (b=20)", "IS (b=20)")
for(j in 1:3){
  if(j==2){
    plot(cbind(1:S, NA), ylim=ran[j,], xaxt="n", xlab="scenario", ylab=Ylab[j], main=Main[j], log="y")
  }else{
    plot(cbind(1:S, NA), ylim=ran[j,], xaxt="n", xlab="scenario", ylab=Ylab[j], main=Main[j], log="y")
  }
  axis(1, at=1:S, 1:S)
  for(k in 1:M){
    points((1:S)+bb[k], Res1[,k,j], type="l", lty=lty[k], col=col[k])
    for(s in 1:S){
      lines(x=c(s,s)+bb[k], y=Res1[s,k,j]+c(-ww,ww)*Error1[s,k,j], col=col[k], lwd=2)
    }
  }
  legend(pos[j], meth, lty=lty, col=col, ncol=2)
}

Main <- c("MSE (b=50)", "SMSE (b=50)", "IS (b=50)")
for(j in 1:3){
  if(j==2){
    plot(cbind(1:S, NA), ylim=ran[j,], xaxt="n", xlab="scenario", ylab=Ylab[j], main=Main[j], log="y")
  }else{
    plot(cbind(1:S, NA), ylim=ran[j,], xaxt="n", xlab="scenario", ylab=Ylab[j], main=Main[j], log="y")
  }
  axis(1, at=1:S, 1:S)
  for(k in 1:M){
    points((1:S)+bb[k], Res2[,k,j], type="l", lty=lty[k], col=col[k])
    for(s in 1:S){
      lines(x=c(s,s)+bb[k], y=Res2[s,k,j]+c(-ww,ww)*Error2[s,k,j], col=col[k], lwd=2)
    }
  }
  legend(pos[j], meth, lty=lty, col=col, ncol=2)
}
dev.off()

