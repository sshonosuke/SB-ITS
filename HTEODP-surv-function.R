library(survival)
library(SparseM)

######  Survival response ######
# Y: censored survival response variable 
# Event: indicator variable for event (1: event, 0: censor)
# X: (n,p)-matrix of biomarkers 
# Z: (n,s)-matrix of adjutment covariates
# PS: propensity score
# Tr: treatment indicator (1 or -1)
# L: number of knots

HTE.ODP.surv <- function(Y, Event, X, Z=NULL, Tr, PS, L=100, alpha.set=NULL, approx=F, maxitr=100){
  p <- dim(X)[2]   # number of biomarkers
  n <- dim(X)[1]   # number of samples
  s <- dim(Z)[2]   # number of adjutment covariates
  if( is.null(Z) ){ s <- 0 }
  if( is.null(alpha.set) ){ alpha.set <- c(0.05,0.1,0.15,0.2) }
  if( is.null(PS) ){ PS <- rep(1,n) }
  ww <- 2*( Tr*PS + (1-Tr)/2 )
  ww <- ( (1/ww) / mean(1/ww) )^(-1)   # normalization
  
  # Summary stat
  Est <- matrix(NA, p, 2)
  NP <- matrix(NA, p, s+1)   # estimtes of nuisance parameters
  for(k in 1:p){
    u1 <- Tr*X[,k]
    u2 <- Tr*Z
    if( is.null(Z)==F ){ fit <- coxph(Surv(Y, Event)~-1+Tr+u1+u2, weights=1/ww) }
    if( is.null(Z)==T ){ fit <- coxph(Surv(Y, Event)~-1+Tr+u1, weights=1/ww) }
    Est[k,] <- summary(fit)$coefficients[2,c(1,3)]
    NP[k,] <- summary(fit)$coefficients[-2,1]
  }
  
  # knots
  ran <- range(Est[,1])
  a.set <- seq(ran[1], ran[2], length=L) 
  a.set <- subset(a.set, a.set!=0)
  L <- length(a.set)
  
  # Initial values of unknown parameters
  Pi <- 0.5
  PP <- rep(1/L,L)
  
  # Estimation (approximation method)
  if(approx==T){
    hbeta <- Est[,1]
    hsig <- Est[,2]
    for(itr in 1:maxitr){
      par <- c(Pi,PP)
      d0 <- dnorm(hbeta, 0, hsig)
      d1 <- matrix(NA, p, L)
      for(k in 1:L){
        d1[,k] <- dnorm(hbeta, a.set[k], hsig)
      }
      md1 <- apply(t(d1)*PP, 2, sum)
      Pw <- t(t(d1)*PP)/md1
      Pz <- (Pi*d0) / (Pi*d0+(1-Pi)*md1)
      
      Pi <- mean(Pz)
      PPw <- (1-Pz)*Pw
      PP <- apply(PPw, 2, sum) / sum(PPw)
      dd <- sqrt( sum(( c(Pi,PP)-par )^2) ) / sqrt(sum(par^2))
      if(dd < 10^(-5) ){ break }
    }
    ODP <- md1/d0
  }
  
  # Estimation (plug-in method)
  if(approx==F){
    tt <- apply(Event/ww*Tr*X, 2, sum)
    if( is.null(Z)==F ){ mat <- Tr*cbind(1, Z) }
    if( is.null(Z)==T ){ mat <- Tr*rep(1,n) }
    ad <- mat%*%t(NP)
    XT <- Tr*X
    CC <- matrix(NA,n,n)
    for(i in 1:n){  CC[i,] <- ifelse(Y>=Y[i],1,0)  }
    CC <- as(CC,"sparseMatrix")
    
    PL <- function(beta){
      val <- apply(Event/ww*log( CC%*%exp(ad+XT*beta)+0.0001 ), 2, sum)
      return( tt*beta-val )
    }
    
    for(itr in 1:maxitr){
      par <- c(Pi,PP)
      L0 <- PL(0)
      L1 <- matrix(NA, p, L)
      for(k in 1:L){ L1[,k] <- PL(a.set[k]) }
      PL1 <- t(log(PP)+t(L1))
      sPL1 <- PL1 - apply(PL1, 1, max)
      Pw <- exp(sPL1) / apply(exp(sPL1), 1, sum)
      ratio <- exp(L1-L0)
      c <- (1-Pi)/Pi*apply(PP*t(ratio), 2, sum)
      Pz <- 1/(1+c)
      
      Pi <- mean(Pz)
      PPw <- (1-Pz)*Pw
      PP <- apply(PPw,2,sum) / sum(PPw)
      dd <- sqrt( sum(( c(Pi,PP)-par )^2) ) / sqrt(sum(par^2))
      if(dd < 10^(-3) ){ break }
    }
    ODP <- apply(PP*t(ratio), 2, sum)
  }
  
  ## Optimal discovery procedure ##
  Lam <- sort(ODP,decreasing=T)
  FDR <- function(lam){ mean(Pz[ODP>=lam]) }
  fdr=c()
  for(r in 1:p){ fdr[r] <- FDR(Lam[r]) }
  
  M <- length(alpha.set)
  SV <- list()
  for(k in 1:M){
    op.lam <- tail(Lam[fdr<alpha.set[k]], 1)
    SV[[k]] <- (1:p)[ODP>=op.lam]
  }
  names(SV) <- paste0("FDR=",alpha.set)
  
  ## Summary ##
  Res <- list(SV=SV, Pi=Pi, PP=PP, Est=Est, knots=a.set)
  return(Res)
}






