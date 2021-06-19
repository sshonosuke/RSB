rm(list=ls())
set.seed(1)

## load functions
library(pscl)
library(MASS)
source("RSB-reg-function.R")

## function for credible intervals
quant <- function(x){ quantile(x, prob=c(0.025,0.975)) }


## simulation settings
n <- 300   # the number of observations
p <- 5  # the number of predictor variables
om1 <- 0.1   #  zero-inflation ratio
om2 <- 0.1   #  large outlier ratio

aa <- 20    # outlier location (10 or 20)
st <- T     # structured zero-inflation (T) or not (F)


## true parameters
Beta <- rep(0,p)  
Beta[c(1,3,5)] <- c(-0.2, 0.3, 0.6)   # non-zero coefficients
int <- 2   # intercept
Para <- c(int, Beta)


## generation of covariates
rho <- 0.2
mat <- matrix(NA, p, p)
for(k in 1:p){
  for(j in 1:p){ mat[k,j]=rho^(abs(k-j)) }
}

X <- mvrnorm(n, rep(0,p), mat) 
if(st){ 
  PP <- 2/( 1+exp(-(0.5*X[,2]+0.5*X[,3])) )
  om1 <- om1*PP 
}


##  data generation  
Mu <- as.vector( exp(int + X%*%Beta) )
Y <- rpois(n, Mu)
ch <- apply(cbind(1-om1-om2, om1, om2), 1, sample, x=c(0, 1, 2), size=1, replace=T)
Y[ch==1] <- 0
Y[ch==2] <- Y[ch==2] + aa


## fitting existing methods 
#  Poisson 
fit1 <- glm(Y~X, family="poisson")
hbeta1 <- coef(fit1)[-1]
sd1 <- summary(fit1)$coefficients[-1,2]
CI1 <- cbind(hbeta1-1.96*sd1, hbeta1+1.96*sd1)

# quasi-Poisson
fit2 <-  glm(Y~X, family="quasipoisson")
hbeta2 <- coef(fit2)[-1]
sd2 <- summary(fit2)$coefficients[-1,2]
CI2 <- cbind(hbeta2-1.96*sd2, hbeta2+1.96*sd2)

# zero-inflated Poisson
if(min(Y)==0){
  fit3 <-  zeroinfl(Y~X | 1)
  hbeta3 <- coef(fit3)[2:6]
  sd3 <- (summary(fit3)$coefficients)$count[2:6, 2]
  CI3 <- cbind(hbeta3-1.96*sd3, hbeta3+1.96*sd3)
  
  try(fit4 <-  zeroinfl(Y~X | X[,2:3]))
  hbeta4 <- coef(fit4)[2:6]
  sd4 <- (summary(fit4)$coefficients)$count[2:6, 2]
  CI4 <- cbind(hbeta4-1.96*sd4, hbeta4+1.96*sd4)
}

# Poisson-gamma
PG.fit <- PG.reg(Y, X, mc=3000, burn=1000)
hbeta5 <- apply(PG.fit$Beta, 2, mean)[-1]
CI5 <- t( apply(PG.fit$Beta, 2, quant)[,-1] )

# (proposed) RSB
rfit <- RSB.reg(Y, X, mc=3000, burn=1000)
hbeta6 <- apply(rfit$Beta, 2, mean)[-1]
CI6 <- t( apply(rfit$Beta, 2, quant)[,-1] )



## squared error of point estimation
Est <- cbind(hbeta1, hbeta2, hbeta3, hbeta4, hbeta5, hbeta6)
100*sqrt( apply((Est-Beta)^2, 2, mean) )


## confidence/credible intervals 
CI <- list(CI1, CI2, CI3, CI4, CI5, CI6)
CI[[3]][,1] < Beta & CI[[3]][,2] > Beta   # coverage of zero-inflated Poisson
CI[[6]][,1] < Beta & CI[[6]][,2] > Beta   # coverage of RSB


