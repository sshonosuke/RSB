rm(list=ls())
set.seed(1)

## load functions
library(MASS)
source("RSB-TF-function.R")


## simulation settings
scenario <- 1     #  1-4
n <- 100
om1 <- rep(0.05, n)   # zero-inflation ratio
om2 <- rep(0.05, n)   # outlier ratio


## true signal 
if(scenario==1){
  lam <- rep(20, n)
}
if(scenario==2){
  lam <- 25 - 10*ifelse((1:n)>20, 1, 0) + 25*ifelse((1:n)>40, 1, 0) - 20*ifelse((1:n)>60, 1, 0)
}
if(scenario==3){
  mat <- matrix(NA, n, n)
  for(i in 1:n){
    mat[i,] <- 100*exp(-(i-(1:n))^2/(2*10^2))
  }
  lam <- mvrnorm(1, rep(20, n), mat)
}
if(scenario==4){
  gg <- sin(4*(1:n)/n-2) + 2*exp(-30*(4*(1:n)/n-2)^2)
  lam <- 20 + 10*gg
}

  
  
## data generation 
Y <- rpois(n, lam)
ch <- apply(cbind(1-om1-om2, om1, om2), 1, sample, x=c(0, 1, 2), size=1, replace=T)
ch[1:3] <- 0
ch[(n-2):n] <- 0
Y[ch==1] <- 0
Y[ch==2] <- Y[ch==2] + 40
  

## fitting
rfit <- RSB.TF(Y, K=2)     # RSB
fit <- RSB.TF(Y, K=2, robust=F)    # Po


## Estimates and credible intervals  
est1 <- apply(rfit$Lam, 2, median)
CI1 <- apply(rfit$Lam, 2, quantile, prob=c(0.025, 0.975))
est2 <- apply(fit$Lam, 2, median)
CI2 <- apply(fit$Lam, 2, quantile, prob=c(0.025, 0.975))



## Plot (results)
plot(Y, ylab="", main=paste0("scenario ", scenario), ylim=range(CI2, Y), pch=16)
points(lam, type="l", lwd=2)
points(est1, type="l", col="blue")
polygon(c(1:n,rev(1:n)),c(CI1[1,],rev(CI1[2,])),col="#0000ff20",border=NA)
points(est2, type="l", col="red")
polygon(c(1:n,rev(1:n)),c(CI2[1,],rev(CI2[2,])),col="#ff000020",border=NA)
lab <- c("true", "Po", "RSB")
col <- c("black", "red", "blue")
legend("topright", legend=lab, col=col, lty=1) 

