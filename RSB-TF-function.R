library(pgdraw)
library(MCMCpack)

###-------------------------------------------------------------###
###          Robust Bayesian trend filtering                    ###
###-------------------------------------------------------------###
# y: observations
# K: order of difference
RSB.TF <- function(Y, K=2, a=0.5, gam=0.5, mc=2000, bn=1000, robust=T){
  ## settings
  n <- length(Y)
  Ep <- 10^(-15)    # lower bound of local scale parameter (to avoid numerical error)
  D1 <- matrix(0, n-1, n)   # D matrix
  for(i in 1:(n-1)){
    for(j in 1:n){
      if(i==j){ D1[i, j] <- -1 }
      if((j-i)==1){ D1[i, j] <- 1 }
    }
  }
  D.set <- list()
  D.set[[1]] <- D1
  if(K>1){
    for(k in 2:K){
      D.set[[k]] <- D1[1:(n-k), 1:(n-k+1)]%*%D.set[[k-1]]
    }
  }
  DD <- D.set[[K]]
  
  ##  useful matrices and values
  delta <- 5000
  Kap <- 0.5*(Y-delta)   # key values in PG-augmentation
  M <- t(DD)%*%DD
  cc <- 1    # shape and scale parameters of gamma prior for alpha
  
  ##  MCMC box
  Lam.pos <- matrix(NA, mc, n)
  Eta.pos <- matrix(NA, mc, n)
  S.pos <- rep(0, mc)
  Z.pos <- matrix(NA, mc, n)
  alpha.pos <- c()
  upper <- max(Y)*100    # maximum values for Eta2 (to avoid numerical error)
  
  ##  initial values of theta (using moving median)
  len <- 5
  EY <- c(rep(Y[1], len), Y, rep(Y[n], len))
  theta <- c()
  for(t in 1:n){
    theta[t] <- log( median( EY[t:(t+2*len)] )+0.01 )
  }
  lam <- exp(theta)
  
  ##  other initial values
  resid <- (Y - lam) / sqrt(lam)
  Z <- ifelse(abs(resid)>3, 1, 0)
  Eta <- Eta1 <- Eta2 <- rep(1, n)
  U <- V <- W <- rep(1, n)    # latent variables for Eta2
  S <- 0.05     # mixing probability 
  uu <- xi <- rep(1, n-K)   # local parameters in HS
  lam2 <- alpha <- 1    # global parameters in HS
  
  ##  MCMC replications
  for(k in 1:mc){
    Lam <- exp(theta)
    if(robust){
      # Eta2 (outlier)  
      V <- rgamma(n, 1-a, log(1+Eta2))
      W <- rgamma(n, a+gam, 1+log(1+Eta2))
      U <- rgamma(n, V+W+1, 1+Eta2)
      p.Eta2 <- rgamma(n, 1, U)
      p.Eta2[p.Eta2 > upper] <- upper
      Eta2 <- Z*rgamma(n,Y+1, Lam+U) + (1-Z)*p.Eta2
      # Eta 
      Eta <- (1-Z)*Eta1 + Z*Eta2
      Eta.pos[k,] <- Eta
      # Z 
      lp1 <- log(1-S) + dpois(Y, Eta1*Lam, log=T) 
      lp2 <- log(S) + dpois(Y, Eta2*Lam, log=T)
      AP <- 1/( 1 + exp(lp1-lp2) )
      Z <- rbinom(n, 1, AP)
      Z.pos[k,] <- Z
      # S
      S <- rbeta(1, 1+sum(Z), 1+n-sum(Z))
      S.pos[k] <- S
    }
    
    # omega (normal approximation)
    Psi <- theta + log(Eta) - log(delta)
    if(delta>10^3){
      mv <- pgdraw.moments(Y+delta, Psi)
      Om <- rnorm(n, mv$mu, sqrt(mv$var))
    }else{
      Om <- pgdraw(Y+delta, Psi)
    }
    
    # theta
    Sc <- uu*lam2
    Sc[Sc<Ep] <- Ep
    A <- solve( diag(Om) + t(DD)%*%diag(1/Sc)%*%DD ) 
    B <- Kap - Om*( log(Eta) - log(delta) ) 
    theta <- mvrnorm(1, A%*%B, A)
    Lam.pos[k,] <- exp(theta)
    # uu
    et <- as.vector(DD%*%theta)
    bb <- 0.5*et^2/lam2 + 1/xi
    uu <- rinvgamma(n-K, 1, bb)
    # xi
    xi <- rinvgamma(n-K, 1, 1+1/uu)
    # Lam2
    bb <- 0.5*sum(et^2/uu) + 1/alpha
    lam2 <- rinvgamma(1, (n-K+1)/2, bb)
    # alpha
    alpha <- rinvgamma(1, 1, 1/lam2)
  }
  
  #  Summary
  om <- 1:bn
  Lam.pos <- Lam.pos[-om,]
  Eta.pos <- Eta.pos[-om,]
  Z.pos <- Z.pos[-om,]
  S.pos <- S.pos[-om]
  Res <- list(Lam=Lam.pos, Eta=Eta.pos, Z=Z.pos, S=S.pos, alpha=alpha.pos)
  return(Res)
}













###-------------------------------------------------------------###
###     Robust Bayesian trend filtering (Gaussian prior)        ###
###-------------------------------------------------------------###
# y: observations
# K: order of difference
RSB.TF.Gauss <- function(Y, K=2, a=0.5, gam=0.5, mc=2000, bn=1000, robust=T){
  ## settings
  n <- length(Y)
  D1 <- matrix(0, n-1, n)   # D matrix
  for(i in 1:(n-1)){
    for(j in 1:n){
      if(i==j){ D1[i, j] <- -1 }
      if((j-i)==1){ D1[i, j] <- 1 }
    }
  }
  D.set <- list()
  D.set[[1]] <- D1
  if(K>1){
    for(k in 2:K){
      D.set[[k]] <- D1[1:(n-k), 1:(n-k+1)]%*%D.set[[k-1]]
    }
  }
  DD <- D.set[[K]]
  
  ##  useful matrices and values
  delta <- 5000
  Kap <- 0.5*(Y-delta)   # key values in PG-augmentation
  M <- t(DD)%*%DD
  cc <- 1    # shape and scale parameters of gamma prior for alpha
  
  ##  MCMC box
  Lam.pos <- matrix(NA, mc, n)
  Eta.pos <- matrix(NA, mc, n)
  S.pos <- rep(0, mc)
  Z.pos <- matrix(NA, mc, n)
  alpha.pos <- c()
  upper <- max(Y)*100    # maximum values for Eta2 (to avoid numerical error)
  
  ##  initial values of theta (using moving median)
  len <- 5
  EY <- c(rep(Y[1], len), Y, rep(Y[n], len))
  theta <- c()
  for(t in 1:n){
    theta[t] <- log( median( EY[t:(t+2*len)] )+0.01 )
  }
  lam <- exp(theta)
  
  ##  other initial values
  resid <- (Y - lam) / sqrt(lam)
  Z <- ifelse(abs(resid)>3, 1, 0)
  Eta <- Eta1 <- Eta2 <- rep(1, n)
  U <- V <- W <- rep(1, n)    # latent variables for Eta2
  S <- 0.05     # mixing probability 
  alpha <- 1   # global scale parameter
  
  ##  MCMC replications
  for(k in 1:mc){
    Lam <- exp(theta)
    if(robust){
      # Eta2 (outlier)  
      V <- rgamma(n, 1-a, log(1+Eta2))
      W <- rgamma(n, a+gam, 1+log(1+Eta2))
      U <- rgamma(n, V+W+1, 1+Eta2)
      p.Eta2 <- rgamma(n, 1, U)
      p.Eta2[p.Eta2 > upper] <- upper
      Eta2 <- Z*rgamma(n,Y+1, Lam+U) + (1-Z)*p.Eta2
      # Eta 
      Eta <- (1-Z)*Eta1 + Z*Eta2
      Eta.pos[k,] <- Eta
      # Z 
      lp1 <- log(1-S) + dpois(Y, Eta1*Lam, log=T) 
      lp2 <- log(S) + dpois(Y, Eta2*Lam, log=T)
      AP <- 1/( 1 + exp(lp1-lp2) )
      Z <- rbinom(n, 1, AP)
      Z.pos[k,] <- Z
      # S
      S <- rbeta(1, 1+sum(Z), 1+n-sum(Z))
      S.pos[k] <- S
    }
    
    # omega (normal approximation)
    Psi <- theta + log(Eta) - log(delta)
    if(delta>10^3){
      mv <- pgdraw.moments(Y+delta, Psi)
      Om <- rnorm(n, mv$mu, sqrt(mv$var))
    }else{
      Om <- pgdraw(Y+delta, Psi)
    }
    
    # theta
    diag(1/(uu*lam2))
    A <- solve( diag(Om) + alpha*M )  
    B <- Kap - Om*( log(Eta) - log(delta) ) 
    theta <- mvrnorm(1, A%*%B, A)
    Lam.pos[k,] <- exp(theta)
    # alpha
    ss <- as.vector(DD%*%theta)
    alpha <- rgamma(1, cc+length(ss)/2, cc+sum(ss^2)/2)
    alpha.pos[k] <- alpha
  }
  
  #  Summary
  om <- 1:bn
  Lam.pos <- Lam.pos[-om,]
  Eta.pos <- Eta.pos[-om,]
  Z.pos <- Z.pos[-om,]
  S.pos <- S.pos[-om]
  Res <- list(Lam=Lam.pos, Eta=Eta.pos, Z=Z.pos, S=S.pos, alpha=alpha.pos)
  return(Res)
}




