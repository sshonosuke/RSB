library(MCMCpack)
library(statmod)
library(GIGrvg)
library(nleqslv)



###-----------------------------------------------------------###
###        Function for robust poisson regression             ###  
###               with RSB/SB distribution                    ### 
###-----------------------------------------------------------###
## INPUT 
# Y: vector of response variable 
# X: (n,p)-matrix of covariates (n: sample size, p: number of covariates)
# offset: n-dimensional vector of offset term (optional)
# robust: If "T", RSB or SB distributions are applied (default is 'T')
# dist: option of robust distribution ("RSB" or "SB")
# a: 1st shape parameter of RSB/SB distribution (default is 0.5) 
# b: 2nd shape parameter of RSB/SB distribution (default is 0.5)
# mc: length of Markov Chain Monte Carlo (default is 2000)
# burn: length of burn-in period (default is 1000)
# Init: initial values of regression coefficients (optional)
# S_prior: 2-dimensional vector of shape parameters in beta prior for mixing proportion (default is c(1,1))

## OUTPUT
## Posterior samples of the following parameter and latent variables 
# Beta: regression coefficients 
# Eta: error components for outliers  
# Z: indicator of using RSB/SB distributions 
# S: mixing proportion

RSB.reg <- function(Y, X, offset=NULL, robust=T, dist="RSB", a=0.5, b=0.5, 
                    mc=2000, burn=1000, Init=NULL, S_prior=c(1,1)){
  ## settings
  XX <- cbind(1, X)
  n <- dim(XX)[[1]]
  p <- dim(XX)[[2]]
  if(is.null(offset)){ offset <- rep(0, n) }
  add <- exp(offset)
  
  # hyperparameters for Beta and S
  beta_v <- rep(0, p)      
  beta_Om <- 0.01*diag(p)
  a_S <- S_prior[1]
  b_S <- S_prior[2]
  
  ## Objects to store posterior samples 
  Beta.pos <- matrix(NA, mc, p)
  Eta.pos <- matrix(NA, mc, n)
  S.pos <- rep(0, mc)
  Z.pos <- matrix(NA, mc, n)
  upper <- 10^5   # maximum values for Eta2 (to avoid numerical error)
  
  ##  Initial values
  Beta <- Init
  if(is.null(Beta)){
    Beta <- coef( glm(Y~X, offset=offset, family="poisson") )   
    mu <- as.vector(exp(XX%*%Beta))
    resid <- (Y - mu) / sqrt(mu)
    Z <- ifelse(abs(resid)>3, 1, 0)
    Beta <- coef( glm(Y[Z==0]~X[Z==0,], offset=offset[Z==0], family="poisson") )   
  }
  mu <- as.vector(exp(XX%*%Beta))
  resid <- (Y - mu) / sqrt(mu)
  Z <- ifelse(abs(resid)>3, 1, 0)
  Eta <- Eta1 <- Eta2 <- rep(1, n)
  U <- V <- W <- rep(1, n)    # latent variables for Eta2
  S <- 0.1     # mixing probability 
  
  ##  MCMC iteration 
  for(k in 1:mc){
    Lam <- exp( as.vector(XX%*%Beta) + offset )
    
    if(robust){
      # Eta2 (outlier) 
      if(dist=="RSB"){
        V <- rgamma(n, 1-a, log(1+Eta2))
        W <- rgamma(n, a+b, 1+log(1+Eta2))
        U <- rgamma(n, V+W+1, 1+Eta2)
        p.Eta2 <- rgamma(n, 1, U)   # from prior
        p.Eta2[p.Eta2 > upper] <- upper
        Eta2 <- Z*rgamma(n, Y+1, Lam+U) + (1-Z)*p.Eta2
      }
      if(dist=="SB"){
        V <- rgamma(n, a+b, 1+Eta2)
        p.Eta2 <- rgamma(n, a, V)   # from prior
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
    
    # Beta (independent MH)
    Pos.beta <- function(bb) {
      mu <- as.vector( exp(XX%*%bb) )
      return( t(XX)%*%(Y-add*Eta*mu) )
    }
    vb0 <- Beta    # current value
    vbh <- nleqslv(vb0, Pos.beta)$x     # mode
    vc <- add*Eta*as.vector(exp(XX%*%vbh))
    mS <- t(XX*vc)%*%XX
    mC <- solve( mS+beta_Om )
    vm <- mC%*%( mS%*%vbh + beta_Om%*%beta_v )
    vb1 <- mvrnorm(1, vm, mC)    # proposal
    log.ratio <- t(Y)%*%XX%*%(vb1-vb0)-(add*Eta)%*%(exp(XX%*%vb1)-exp(XX%*%vb0))+0.5*( t(vb1-vbh)%*%mS%*%(vb1-vbh)-t(vb0-vbh)%*%mS%*%(vb0-vbh) )
    prob <- min(exp(log.ratio), 1)
    ch <- rbinom(1, 1, prob)
    Beta <- as.vector( vb0+ch*(vb1-vb0) )
    Beta.pos[k,] <- Beta
  }
  
  ## Summary
  om <- 1:burn
  Beta.pos <- Beta.pos[-om,]
  Eta.pos <- Eta.pos[-om,]
  S.pos <- S.pos[-om]
  Z.pos <- Z.pos[-om,]
  Res <- list(Beta=Beta.pos, Eta=Eta.pos, Z=Z.pos, S=S.pos)
  return(Res)
}









###-----------------------------------------------------------###
###        Function for zero-inflated Poisson-gamma           ###
###            (negative binomial) regression                 ###  
###-----------------------------------------------------------###
## INPUT 
# Y: vector of response variable 
# X: (n,p)-matrix of covariates (n: sample size, p: number of covariates)
# ZI: If 'T', zero inflation component is added (default is 'F')
# mix: If 'T', Poisson-gamma model is applied (default is 'T')  
# mc: length of Markov Chain Monte Carlo (default is 2000)
# burn: length of burn-in period (default is 1000)

## OUTPUT
# Posterior samples of the following parameter and latent variables 
# Beta: regression coefficients 
# Eta: error components for outliers  
# gam: shape parameter for latent gamma distribution  
# S: mixing proportion

ZIPG.reg <- function(Y, X, ZI=F, mix=T, mc=2000, burn=1000){
  ## settings
  XX <- cbind(1, X)
  n <- dim(XX)[[1]]
  p <- dim(XX)[[2]]
  
  # hyperparameter
  vm0 <- rep(0, p)      
  Om <- 0.001*diag(p)
  
  ## MCMC box
  Beta.pos <- matrix(NA, mc, p)
  Eta.pos <- matrix(NA, mc, n)
  gam.pos <- c()
  S.pos <- c()  # zero-inflation probability
  
  ## Initial values
  Beta <- coef( glm(Y~X, family="poisson") )   
  mu <- as.vector(exp(XX%*%Beta))
  Eta <- rep(1, n)
  gam <- 1
  Z <- rep(0, n)
  S <- 0.1
  zero_count <- (Y==0)   # indicator of zero count
  
  ## MCMC iteration
  for(k in 1:mc){
    Mu <- exp( as.vector(XX%*%Beta) )
    
    if(mix){
      # Eta
      Eta_pos <- rgamma(n, gam+Y, gam+Mu)
      Eta_prior <- rgamma(n, gam, gam)
      Eta <- (1-Z)*Eta_pos + Z*Eta_prior
      Eta[Eta<0.0001] <- 0.0001
      Eta.pos[k,] <- Eta
      # gam (random-walk MH)
      gam.new <- gam + 0.3*rnorm(1)
      if( gam.new<0.001 ){ gam.new <- 0.001 }
      LL <- function(x){
        x*log(x) - log(gamma(x)) + (x-1)*mean(log(Eta)) - x*mean(Eta)
      }
      prob <- min(1, exp(n*LL(gam.new) - n*LL(gam)))
      ch <- rbinom(1, 1, prob)
      gam <- (1-ch)*gam + ch*gam.new
      gam.pos[k] <- gam
    }
    
    # Beta (independent MH)
    sXX <- XX[Z==0,]
    sY <- Y[Z==0]
    sEta <- Eta[Z==0]
    func <- function(bb) {
      mu <- exp(sXX%*%bb)
      return( t(sXX)%*%(sY-sEta*mu) )
    }
    vb0 <- Beta
    vbh <- nleqslv(vb0, func)$x
    vc <- sEta*exp(sXX%*%vbh)
    mS <- t(sXX*as.vector(vc))%*%sXX
    mC <- solve( mS+Om )
    vm <- mC%*%( mS%*%vbh + Om%*%vm0 )
    vb1 <- mvrnorm(1, vm, mC)    # proporsal
    log.ratio <- t(sY)%*%sXX%*%(vb1-vb0)-sEta%*%(exp(sXX%*%vb1)-exp(sXX%*%vb0))+0.5*( t(vb1-vbh)%*%mS%*%(vb1-vbh)-t(vb0-vbh)%*%mS%*%(vb0-vbh) )
    prob <- min(exp(log.ratio), 1)
    ch <- rbinom(1, 1, prob)
    Beta <- as.vector( vb0+ch*(vb1-vb0) )
    Beta.pos[k,] <- Beta
    
    # zero-inflation
    if(ZI){
      # Z
      lam <- as.vector( exp(XX%*%Beta) )
      dP <- dpois(Y, lam*Eta)
      Z[zero_count] <- rbinom(sum(zero_count), 1, 1/(1+dP[zero_count]))
      # S
      S <- rbeta(1, 1+sum(Z), 1+n-sum(Z))
      S.pos[k] <- S
    }
  }
  
  ## Summary
  om <- 1:burn
  Beta.pos <- Beta.pos[-om,]
  Eta.pos <- Eta.pos[-om,]
  gam.pos <- gam.pos[-om]
  S.pos <- S.pos[-om]
  Res <- list(Beta=Beta.pos, Eta=Eta.pos, gam=gam.pos, S=S.pos)
  return(Res)
}



