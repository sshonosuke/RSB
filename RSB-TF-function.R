library(pgdraw)
library(MCMCpack)
library(BayesLogit)


###-----------------------------------------------------------###
###        Function for robust Poisson trend filtering        ###
###        with horseshoe prior and RSB/SB distribution       ### 
###-----------------------------------------------------------###
## INPUT 
# y: observations
# K: order of difference (default is 2)
# robust: If "T", RSB or SB distributions are applied (default is 'T')
# dist: option of robust distribution ("RSB" or "SB")
# a: 1st shape parameter of RSB/SB distribution (default is 0.5) 
# b: 2nd shape parameter of RSB/SB distribution (default is 0.5)
# mc: length of Markov Chain Monte Carlo (default is 2000)
# burn: length of burn-in period (default is 1000)
# Init: initial values of trend values (optional)

## OUTPUT
# Posterior samples of the following parameter and latent variables 
# Lam: trend values 
# Eta: error components for outliers  
# Z: indicator of using RSB/SB distributions 
# S: mixing proportion
# alpha: scale parameter in horseshoe prior 

RSB.TF <- function(Y, K=2, robust=T, dist="RSB", a=0.5, b=0.5, mc=2000, bn=1000, init=NULL){
  ## settings
  n <- length(Y)
  Ep <- 10^(-8)    # lower bound of local scale parameter (to avoid numerical error)
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
  delta <- 3000
  Kap <- 0.5*(Y-delta)   # key values in PG-augmentation
  M <- t(DD)%*%DD
  cc <- 1    # shape and scale parameters of gamma prior for alpha
  a_S <- 0.05*sqrt(n)
  b_S <- 0.95*sqrt(n)
  
  ##  MCMC box
  Lam.pos <- matrix(NA, mc, n)
  Eta.pos <- matrix(NA, mc, n)   
  S.pos <- rep(0, mc)     # mixing proportion 
  Z.pos <- matrix(NA, mc, n)      # mixing indicator
  alpha.pos <- c()
  upper <- 10^8    # maximum values for Eta2 (to avoid numerical error)
  
  ##  initial values of theta (using moving median)
  len <- 5
  EY <- c(rep(Y[1], len), Y, rep(Y[n], len))
  theta <- c()
  for(t in 1:n){
    theta[t] <- log( median( EY[t:(t+2*len)] )+0.01 )
  }
  lam <- exp(theta)
  if(is.null(init)==F){ 
    lam <- init 
    theta <- log(init)
  }
  
  ##  other initial values
  resid <- (Y - lam) / sqrt(lam)
  Z <- ifelse(abs(resid)>3, 1, 0)
  Eta <- Eta1 <- Eta2 <- rep(1, n)
  Eta2 <- Y/lam + 0.01
  U <- V <- W <- rep(1, n)    # latent variables for Eta2
  S <- 0.05     # mixing probability 
  uu <- xi <- rep(1, n-K)   # local parameters in HS
  lam2 <- alpha <- 1    # global parameters in HS
  
  ##  MCMC replications
  for(k in 1:mc){
    Lam <- exp(theta)
    if(robust){
      # Eta2 (outlier)
      if(dist=="RSB"){
        V <- rgamma(n, 1-a, log(1+Eta2))
        W <- rgamma(n, a+b, 1+log(1+Eta2))
        U <- rgamma(n, V+W+1, 1+Eta2)
        p.Eta2 <- rgamma(n, 1, U)
        p.Eta2[p.Eta2 > upper] <- upper
        Eta2 <- Z*rgamma(n, Y+1, Lam+U) + (1-Z)*p.Eta2
      }
      if(dist=="SB"){
        V <- rgamma(n, a+b, 1+Eta2)
        p.Eta2 <- rgamma(n, a, V)
        p.Eta2[p.Eta2 > upper] <- upper
        Eta2 <- Z*rgamma(n, Y+a, Lam+V) + (1-Z)*p.Eta2
      }
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
      S <- rbeta(1, a_S+sum(Z), b_S+n-sum(Z))
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







###-----------------------------------------------------------###
###      Function for Poisson-gamma (negative binomial)       ###  
###          trend filtering with horseshoe prior             ###  
###-----------------------------------------------------------###
## INPUT 
# y: observations
# K: order of difference (default is 2)
# mc: length of Markov Chain Monte Carlo (default is 2000)
# burn: length of burn-in period (default is 1000)

## OUTPUT
## Posterior samples of the following parameter and latent variables 
# Lam: trend values 
# Nu: shape paraemter of latent gamma distribution 

NB.TF <- function(Y, K=2, mc=2000, bn=1000){
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
  M <- t(DD)%*%DD
  cc <- 1    # shape and scale parameters of gamma prior for alpha
  
  ##  MCMC box
  Lam.pos <- matrix(NA, mc, n)     # mean parameter
  Nu.pos <- c()      # overdispersion 
  alpha.pos <- c()     # global scale parameter in HS
  upper <- max(Y)*100    # maximum values for Eta2 (to avoid numerical error)
  
  ##  initial values of theta (using moving median)
  len <- 5
  EY <- c(rep(Y[1], len), Y, rep(Y[n], len))
  theta <- c()
  for(t in 1:n){
    theta[t] <- log( median( EY[t:(t+2*len)] )+0.01 )
  }
  nu <- 1
  
  ##  other initial values
  uu <- xi <- rep(1, n-K)   # local parameters in HS
  lam2 <- alpha <- 1    # global parameters in HS
  delta <- Y + 0.1    # auxiliary variable for sampling of nu 
  
  ##  MCMC replications
  for(k in 1:mc){
    Lam <- exp(theta)
    # nu (dispersion parameter)
    nu <- rgamma(1, cc+sum(delta), cc+sum(log(1+Lam)))
    Nu.pos[k] <- nu
    # delta
    for(i in 1:n){
      if(Y[i]==0){ delta[i]=0 }
      if(Y[i]>0){
        pp <- nu/((1:Y[i])-1+nu)
        delta[i] <- sum( rbinom(length(pp), 1, pp) )
      }
    }
    
    # omega (PG augmentation)
    Kap <- 0.5*(Y-nu)
    Om <- rpg(num=n, h=Y+nu, z=theta)
    
    # theta
    Sc <- uu*lam2
    Sc[Sc<Ep] <- Ep
    A <- solve( diag(Om) + t(DD)%*%diag(1/Sc)%*%DD ) 
    B <- Kap 
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
  Nu.pos <- Nu.pos[-om]
  Res <- list(Lam=Lam.pos, Nu=Nu.pos)
  return(Res)
}









