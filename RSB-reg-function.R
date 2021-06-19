library(MCMCpack)
library(statmod)
library(GIGrvg)
library(nleqslv)



###   Robust poisson regression with RSB distribution  ### 
RSB.reg <- function(Y, X, offset=NULL, robust=T, a=0.5, gam=0.5, mc=2000, burn=1000){
  ## settings
  XX <- cbind(1, X)
  n <- dim(XX)[[1]]
  p <- dim(XX)[[2]]
  if(is.null(offset)){ offset <- rep(0, n) }
  add <- exp(offset)
  
  # hyperparameters for Beta
  beta_v <- rep(0, p)      
  beta_Om <- 0.001*diag(p)
  
  ## MCMC box
  Beta.pos <- matrix(NA, mc, p)
  Eta.pos <- matrix(NA, mc, n)
  S.pos <- rep(0, mc)
  Z.pos <- matrix(NA, mc, n)
  upper <- max(Y)*100    # maximum values for Eta2 (to avoid numerical error)
  
  ##  Initial values
  Beta <- coef( glm(Y~X, offset=offset, family="poisson") )   
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
      
      # S
      S <- rbeta(1, 1+sum(Z), 1+n-sum(Z))
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








###  Poisson-gamma regression 
PG.reg <- function(Y, X, mc=2000, burn=1000){
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
  
  ## Initial values
  Beta <- coef( glm(Y~X, family="poisson") )   
  mu <- as.vector(exp(XX%*%Beta))
  Eta <- rep(1, n)
  gam <- 1
  
  ## MCMC iteration
  for(k in 1:mc){
    Mu <- exp( as.vector(XX%*%Beta) )
    # Eta 
    Eta <- rgamma(n, gam+Y, gam+Mu)
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
    
    # Beta (independent MH)
    func <- function(bb) {
      mu <- exp(XX%*%bb)
      return( t(XX)%*%(Y-Eta*mu) )
    }
    vb0 <- Beta
    vbh <- nleqslv(vb0, func)$x
    vc <- Eta*exp(XX%*%vbh)
    mS <- t(XX*as.vector(vc))%*%XX
    mC <- solve( mS+Om )
    vm <- mC%*%( mS%*%vbh + Om%*%vm0 )
    vb1 <- mvrnorm(1, vm, mC)    # proporsal
    log.ratio <- t(Y)%*%XX%*%(vb1-vb0)-Eta%*%(exp(XX%*%vb1)-exp(XX%*%vb0))+0.5*( t(vb1-vbh)%*%mS%*%(vb1-vbh)-t(vb0-vbh)%*%mS%*%(vb0-vbh) )
    prob <- min(exp(log.ratio), 1)
    ch <- rbinom(1, 1, prob)
    Beta <- as.vector( vb0+ch*(vb1-vb0) )
    Beta.pos[k,] <- Beta
  }
  
  ## Summary
  om <- 1:burn
  Beta.pos <- Beta.pos[-om,]
  Eta.pos <- Eta.pos[-om,]
  gam.pos <- gam.pos[-om]
  Res <- list(Beta=Beta.pos, Eta=Eta.pos, gam=gam.pos)
  return(Res)
}








