library(MCMCpack)
library(statmod)
library(GIGrvg)
library(pgdraw)





###   Robust spatial poisson regression     ###
# band: bandwidth in RW-MH for spatial range parameter 
RSB.sp <- function(Y, X, Sp, offset=NULL, robust=T, a=0.5, gam=0.5, band=0.1, M=100, mc=3000, burn=1000, delta=1000){
  # settings
  XX <- cbind(1, X)
  n <- dim(XX)[[1]]
  p <- dim(XX)[[2]]
  knot <- kmeans(Sp, M)$centers
  if(is.null(offset)){ offset <- rep(0, n) }
  Kap <- 0.5*(Y-delta)   # key values in PG-augmentation 
  
  P.mat <- matrix(NA, n, M)
  for(m in 1:M){
    P.mat[,m] <- apply((t(Sp)-knot[m,])^2, 2, sum)
  }
  kDD <- as.matrix(dist(knot)^2)
  
  # hyper-parameters
  C <- 0.01*diag(p)      # precision matrix for beta
  cc <- 1       # shape and scale parameters of gamma prior for Tau 
  
  # MCMC box
  Beta.pos <- matrix(NA, mc, p)   # regression coefficients
  Eta.pos <- matrix(NA, mc, n)    # adjustment term
  S.pos <- rep(NA, mc)      # mixing proportion 
  Z.pos <- matrix(NA, mc, n)    # outlier indicator
  Xi.pos <- matrix(NA, mc, n)   # spatial effect 
  upper <- max(Y)*100     # thresholding values for Eta2 (to avoid numerical error)
  h.pos <- rep(NA, mc)     # bandwidth for spatial effect
  Tau.pos <- rep(NA, mc)      # precision parameter for spatial effect
  
  # Initial values
  Beta <- coef( glm(Y~X, offset=offset, family="poisson") )   
  mu <- exp( as.vector(XX%*%Beta) + offset )
  resid <- (Y - mu) / sqrt(mu)
  Z <- ifelse(abs(resid)>3, 1, 0)
  Eta <- Eta1 <- Eta2 <- rep(1, n)
  U <- V <- W <- rep(1, n)    # latent variables for Eta2
  S <- 0.1     # mixing probability 
  h <- 0.2     # spatial range 
  Xi <- rep(0, n)    # spatial effect
  Tau <- 1      # precision of spatial effect
  
  # MCMC
  for(k in 1:mc){
    Lam <- exp(as.vector(XX%*%Beta) + offset + Xi)
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
    Psi <- as.vector(XX%*%Beta) + offset + Xi + log(Eta) - log(delta)
    if(delta>10^3){
      mv <- pgdraw.moments(Y+delta, Psi)
      Om <- rnorm(n, mv$mu, sqrt(mv$var))
    }else{
      Om <- pgdraw(Y+delta, Psi)
    }
      
    # Beta 
    mat <- solve( t(XX*Om)%*%XX + C )
    bb <- Kap - Om*(offset + Xi + log(Eta) - log(delta)) 
    mm <- as.vector( mat%*%t(XX)%*%bb )
    Beta <- mvrnorm(1, mm, mat)
    Beta.pos[k,] <- Beta
    reg <- as.vector(XX%*%Beta)
    
    # Xi (predictive process)
    H <- exp(-kDD/(2*h^2) )
    dec <- eigen(H)
    IH <- (dec$vectors)%*%diag(1/dec$values)%*%t(dec$vectors)
    P <- exp(-P.mat/(2*h^2))
    D <- P%*%IH
    mat <- t(D*Om)%*%D
    A <- solve( mat + Tau*IH )
    bb <- Kap - Om*( reg + offset + log(Eta) - log(delta) ) 
    mm <- as.vector( A%*%t(D)%*%bb )
    Ep.Xi <- mvrnorm(1, mm, A)
    Xi <- as.vector( D%*%Ep.Xi )  
    Xi.pos[k,] <- Xi
    # Tau
    Tau <- rgamma(1, cc+M/2, cc+t(Ep.Xi)%*%IH%*%Ep.Xi/2)
    Tau.pos[k] <- Tau
    
    # bandwidth 
    bb <- -sum(log(eigen(H)$values[eigen(H)$values>0]))
    new.h <- h + band*rnorm(1)
    new.h[new.h <- 10^(-8)] <- 10^(-8)
    new.H <- exp(-kDD/(2*new.h^2))
    dec <- eigen(new.H)
    new.IH <- (dec$vectors)%*%diag(1/dec$values)%*%t(dec$vectors)
    new.bb <- -sum(log(eigen(new.H)$values[eigen(new.H)$values>0]))
    val1 <- 0.5*bb - 0.5*t(Ep.Xi)%*%IH%*%Ep.Xi/Tau^2
    val2 <- 0.5*new.bb - 0.5*t(Ep.Xi)%*%new.IH%*%Ep.Xi/Tau^2
    prob <- min(1, exp(val2-val1))
    ch <- rbinom(1, 1, prob)
    h <- h + ch*(new.h-h)
    h.pos[k] <- h
    
    if(round(k/500)==(k/500)){ print(k); plot(h.pos)}
  }
  
  # Summary 
  om <- 1:burn
  Beta.pos <- Beta.pos[-om,]
  Eta.pos <- Eta.pos[-om,]
  Xi.pos <- Xi.pos[-om,]
  Z.pos <- Z.pos[-om,]
  h.pos <- h.pos[-om]
  S.pos <- S.pos[-om]
  Tau.pos <- Tau.pos[-om]
  Res <- list(Beta=Beta.pos, Eta=Eta.pos, Xi=Xi.pos, Z=Z.pos, h=h.pos, Tau=Tau.pos, Z=Z.pos, S=S.pos)
  return(Res)
}
