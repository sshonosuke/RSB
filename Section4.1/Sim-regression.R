###-----------------------------------------------------------###
###      R code to replicate the results in Section 4.1       ###
###        (Simulation study under Poisson regression)        ###                
###-----------------------------------------------------------###

## Preparation 
library(MASS)
source("RSB-function.R")
set.seed(1)

# simulation scenario 
aa <- 20    # outlier location (20 or 50)

# function for CI
quant <- function(x){ quantile(x, prob=c(0.025,0.975)) }

# function for interval score
Interval_score <- function(Interval, Para){
  sum( mean(Interval[2,]-Interval[1,]) + 40*ifelse(Para>Interval[2,], 1, 0)*(Para-Interval[2,]) + 40*ifelse(Para<Interval[1,], 1, 0)*(Interval[1,]-Para) )
}


## Settings
R <- 500  # number of Monte Carlo replications
n <- 200   # number of observations
p <- 15    # number of predictor variables
om1_set <- c(0, rep(0.05, 4), rep(0.1, 3))  #  zero-inflation ratio
om2_set <- c(0, 0, rep(c(0.05, 0.1, 0.15), 2))  #  large outlier ratio
S <- length(om1_set)

mc <- 2000
burn <- 1000

# prior
prior_beta_m <- rep(0, p+1)
prior_beta_v <- diag(rep(100, p+1))


## Replications 
meth <- c("RSB", "SB", "PO", "ZIP", "PG", "ZIPG")
M <- length(meth)
MSE_beta <- array(NA, c(R, S, M))
MSE_mu <- array(NA, c(R, S, M))
CP <- array(NA, c(R, S, M))
AL <- array(NA, c(R, S, M))
IS <- array(NA, c(R, S, M))


for(s in 1:S){
  om1 <- rep(om1_set[s], n)
  om2 <- rep(om2_set[s], n)
  for(r in 1:R){
    # True parameters 
    Beta <- c(0.2, 0.5, 0.3, rep(0, p-3))
    Beta[1] <- Beta[1] + runif(1, -0.2, 0.2)
    Beta[2] <- Beta[2] + runif(1, -0.2, 0.2)
    Beta[3] <- Beta[3] + runif(1, -0.2, 0.2)
    Int <- (0.5)   # intercept
    Coef <- c(Int, Beta)
    
    # covariates
    rho <- 0.2
    mat <- matrix(NA, p, p)
    for(k in 1:p){
      for(j in 1:p){ mat[k,j]=rho^(abs(k-j)) }
    }
    X <- mvrnorm(n, rep(1, p), mat) 
    
    # response values
    tMu <- as.vector( exp(Int + X%*%Beta) )
    Y <- rpois(n, tMu)
    
    # outlier and zero inflation
    ch <- apply(cbind(1-om1-om2, om1, om2), 1, sample, x=c(0, 1, 2), size=1, replace=T)
    Y[ch==1] <- 0
    Y[ch==2] <- Y[ch==2] + aa
    
    
    ## Fitting models
    # RSB
    fit.RSB <- RSB.reg(Y, X, a=1/2, b=1/2, dist="RSB", mc=mc, burn=burn)
    
    # SB 
    fit.SB <- RSB.reg(Y, X, a=1/2, b=1/10, dist="SB", mc=mc, burn=burn)
    
    # Poison distribution 
    fit.Po <- ZIPG.reg(Y, X, ZI=F, mix=F, mc=mc, burn=burn)
    
    # ZIP distribution 
    fit.ZIP <- ZIPG.reg(Y, X, ZI=T, mix=F, mc=mc, burn=burn)
    
    # Poison-gamma distribution (negative binomial)
    fit.PG <- ZIPG.reg(Y, X, ZI=F, mix=T, mc=mc, burn=burn)
    
    # ZIPG distribution 
    fit.ZIPG <- ZIPG.reg(Y, X, ZI=T, mix=T, mc=mc, burn=burn)
    
    ## Evaluation 
    # MSE (beta)
    Est <- NULL
    Est <- cbind(Est, apply(fit.RSB$Beta, 2, mean))
    Est <- cbind(Est, apply(fit.SB$Beta, 2, mean))
    Est <- cbind(Est, apply(fit.Po$Beta, 2, mean))
    Est <- cbind(Est, apply(fit.ZIP$Beta, 2, mean))
    Est <- cbind(Est, apply(fit.PG$Beta, 2, mean))
    Est <- cbind(Est, apply(fit.ZIPG$Beta, 2, mean))
    
    MSE_beta[r,s,] <- apply((Est-Coef)^2, 2, mean) 
    
    # MSE (mu)
    hMu <- NULL
    hMu <- cbind(hMu, apply(exp(cbind(1,X)%*%t(fit.RSB$Beta)), 1, mean))
    hMu <- cbind(hMu, apply(exp(cbind(1,X)%*%t(fit.SB$Beta)), 1, mean))
    hMu <- cbind(hMu, apply(exp(cbind(1,X)%*%t(fit.Po$Beta)), 1, mean))
    hMu <- cbind(hMu, apply(exp(cbind(1,X)%*%t(fit.ZIP$Beta)), 1, mean))
    hMu <- cbind(hMu, apply(exp(cbind(1,X)%*%t(fit.PG$Beta)), 1, mean))
    hMu <- cbind(hMu, apply(exp(cbind(1,X)%*%t(fit.ZIPG$Beta)), 1, mean))
    
    MSE_mu[r,s,] <- apply((hMu-tMu)^2/tMu^2, 2, mean) 
    
    # credible intervals
    CI <- list()
    CI[[1]] <- apply(fit.RSB$Beta, 2, quantile, prob=c(0.025, 0.975))
    CI[[2]] <- apply(fit.SB$Beta, 2, quantile, prob=c(0.025, 0.975))
    CI[[3]] <- apply(fit.Po$Beta, 2, quantile, prob=c(0.025, 0.975))
    CI[[4]] <- apply(fit.ZIP$Beta, 2, quantile, prob=c(0.025, 0.975))
    CI[[5]] <- apply(fit.PG$Beta, 2, quantile, prob=c(0.025, 0.975))
    CI[[6]] <- apply(fit.ZIPG$Beta, 2, quantile, prob=c(0.025, 0.975))
    
    sCP <- NULL
    sAL <- NULL
    sIS <- NULL
    for(l in 1:M){
      sCP <- c(sCP, 100*mean( CI[[l]][1,]<Coef & CI[[l]][2,]>Coef ))
      sAL <- c(sAL, mean( CI[[l]][2,]-CI[[l]][1,] ))
      sIS <- c(sIS, Interval_score(CI[[l]], Coef))
    }
    CP[r,s,] <- sCP
    AL[r,s,] <- sAL
    IS[r,s,] <- sIS
    
    # print
    if(round(r/100)==(r/100)){ 
      print(r)
      print( sqrt(apply(MSE_beta[1:r,s,], 2, mean)) )
      print( sqrt(apply(MSE_mu[1:r,s,], 2, mean)) )
      print( apply(IS[1:r,s,], 2, mean) )
    }
  }
  
  # print
  print(s)
}


# Save results
save(list=ls(), file=paste0("Sim-",aa,".RData")) 

