rm(list=ls())

## load dataset
load("CrimeDataset.RData")


## load function
source("RSB-spatial-function.R")


## data 
sum(Y==0)/length(Y)    #  ratio of zero count 
plot(Y/Area, ylim=c(0, 100))

par(mfcol=c(1, 2))
plot(Y, ylab="observed count", xlab="sample index")
plot(Y, ylim=c(0, 20), ylab="observed count", xlab="sample index")



## fitting RSB and Po models
mc <- 5000
bn <- 1000
M <- 50    # number of knots for predictive process (takes a long time to compute with large M) 
set.seed(1)
rfit <- RSB.sp(Y, X, Sp, offset=log(Area), M=M, band=0.05, mc=mc, burn=bn, delta=5000)    # RSB
fit <- RSB.sp(Y, X, Sp, offset=log(Area), M=M, robust=F, band=0.05, mc=mc, burn=bn, delta=5000)    # Po




## DIC (RSB)
Lam.pos <- exp( cbind(1,X)%*%t(rfit$Beta) + t(rfit$Xi) + log(Area) )
Eta.pos <- t(rfit$Eta)
D1 <- sum( apply(-2*dpois(Y, Lam.pos*Eta.pos, log=T), 1, mean) )
hBeta <- apply(rfit$Beta, 2, mean)
hXi <- apply(rfit$Xi, 2, mean)
hEta <- apply(rfit$Eta, 2, mean)
hLam <- exp( as.vector(cbind(1,X)%*%hBeta) + hXi + log(Area) )
D2 <- sum(-2*dpois(Y, hLam*hEta, log=T))
2*D1 - D2


## DIC (Poisson)
Lam.pos <- exp( cbind(1,X)%*%t(fit$Beta) + t(fit$Xi) + log(Area) )
D1 <- sum( apply(-2*dpois(Y, Lam.pos, log=T), 1, mean) )
hBeta <- apply(fit$Beta, 2, mean)
hXi <- apply(fit$Xi, 2, mean)
hLam <- exp( as.vector(cbind(1,X)%*%hBeta) + hXi + log(Area) )
D2 <- sum(-2*dpois(Y, hLam, log=T))
2*D1 - D2




## plot (Eta and S)
Eta <- apply(rfit$Eta, 2, median)
n <- length(Eta)
col <- rep(1, n)


par(mfcol=c(1,2))
plot(Eta, xlab="index", ylab="Eta", cex=0.3, pch=4, col=col)
plot(rfit$S, xlab="MCMC iteration", ylab="s", type="l")



##  Plot (spatial effects)
out <- (1:n)[Eta>10]
ZI <- (1:n)[Eta<0.1]

Xi1 <- apply(rfit$Xi, 2, mean)
Xi2 <- apply(fit$Xi, 2, mean)

par(mfcol=c(1,1))
matplot(cbind(Xi1, Xi2), type="l", lty=1, xlab="area index")
legend("bottomleft", legend=c("RSB", "Po"), lty=1, col=1:2)




