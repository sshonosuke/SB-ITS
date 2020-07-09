# SB-ITS
This package implements effieicnt screening methods of predictive biomarkers for individual treatment selection via optimal discovery procedure, as proposed by the following papers.

Sugasawa, S. and Noma, H. (2020). Efficient Screening of Predictive Biomarkers for Individual Treatment Selection. *Biometrics*, to appear (https://arxiv.org/abs/1905.01582)

Functions are implemented in HTEODP-bin-function.R and HTEODP-surv-function.R available in the repository.


# Binary response
Install R function and demo data.
```{r}
load("Example-Data-bin.RData") 
source("HTEODP-bin-function.R")
```

Also install `qvalue' package used for comparison.
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("qvalue")
library(qvalue)
```

Apply the proposed method.
(n: sample size;  p: the number of candidate biomarkers) 

Input of `HTE.ODP.bin`
- `Y`: n-dimensional vector of binary response variables (1 or 0) 
- `X`: (n,p)-matrix of biomarkers 
- `Z`: (n,q)-matrix of adjutment covariates, where q is the number of covariates
- `PS`: n-dimensional vector of propensity scores
- `Tr`: n-dimensional vector of treatment indicators (1 or -1)
- `L`: number of knots (default is 100) 
- `alpha.set`: vector of nominal FDR values
- `approx`: ODP-N is applied if `T`, and ODP-P is applied if `F`
- `maxitr`: maximum number of iterations in EM algorithm (default us 100)

Output of `HTE.ODP.bin`
- `SV`: list of significant biomarkers
- `Pi`: esitmated probability of null
- `PP`: estimated point mass probabilities for non-null distribution
- `Est`: (n,2)-matrix of point estimates and standard errors
- `knots`: set of knots 

```{r}
alpha.set <- c(0.05, 0.1, 0.15, 0.2)
fit.N <- HTE.ODP.bin(Y, X, Z, Tr, PS, L=200, alpha.set=alpha.set, approx=T, maxitr=300)
fit.P <- HTE.ODP.bin(Y, X, Z, Tr, PS, L=200, alpha.set=alpha.set, approx=F, maxitr=300)
```

Significant biomarkers
```{r}
SV.N <- fit.N$SV
SV.P <- fit.P$SV
```
q-value method
```{r}
p <- dim(X)[2]
Stat <- fit.N$Est    # summary statistics
PV1 <- 2*( 1-pnorm( abs(Stat[,1]/Stat[,2]) ) )   
PV2 <- c()
for(k in 1:p){
  f1 <- glm(Y[Tr==1]~X[Tr==1,k]+Z[Tr==1,], family="binomial")
  f2 <- glm(Y[Tr==-1]~X[Tr==-1,k]+Z[Tr==-1,], family="binomial")
  a1 <- summary(f1)$coefficients[2,1] - summary(f2)$coefficients[2,1]
  a2 <- sqrt( summary(f1)$coefficients[2,2]^2 + summary(f2)$coefficients[2,2]^2 )
  PV2[k] <- 2*( 1-pnorm(abs(a1/a2)) )
}

qv1 <- list()
qv2 <- list()
for(k in 1:4){
  qv1[[k]] <- (1:p)[ qvalue(PV1,fdr=alpha.set[k])$significant ]
  qv2[[k]] <- (1:p)[ qvalue(PV2,fdr=alpha.set[k])$significant ]
}
```

Summary of the results 
```{r}
SV <- list(SV.N, SV.P, qv1, qv2)   # list of significant biomarkers from the four methods
NS <- matrix(NA, 4, 4)    # number of significant 
TP <- matrix(NA, 4, 4)    # number of true positives
dimnames(NS)[[1]] <- dimnames(TP)[[1]] <- c("ODP.N", "ODP.P", "qvalue1", "qvalue2")
dimnames(NS)[[2]] <- dimnames(TP)[[2]] <- paste0("FDR=",alpha.set)

for(j in 1:4){
  NS[j,] <- as.numeric(lapply(SV[[j]], length))
  for(k in 1:4){
    TP[j,k] <- length(intersect(SV[[j]][[k]], Ind))
  }
}
FP <- NS - TP   # number of false positives
```




# Survival response 
Install R function and demo data.
```{r}
load("Example-Data-surv.RData") 
source("HTEODP-surv-function.R")
```

Apply the proposed method.
(n: sample size;  p: the number of candidate biomarkers) 

In addition to the input of `HTE.ODP.bin`, `HTE.ODP.surv` requires `Event`. 
- `Event`: n-dimensional vector of event indicator (1: event; 0: censored)

Output of `HTE.ODP.surv` is the same as `HTE.ODP.bin`.
```{r}
alpha.set <- c(0.05, 0.1, 0.15, 0.2)
fit.N <- HTE.ODP.surv(Y, Event, X, Z, Tr, PS, L=200, alpha.set=alpha.set, approx=T, maxitr=300)
fit.P <- HTE.ODP.surv(Y, Event, X, Z, Tr, PS, L=200, alpha.set=alpha.set, approx=F, maxitr=300)
```

q-value method
```{r}
p <- dim(X)[2]
Stat <- fit.N$Est
PV1 <- 2*( 1-pnorm( abs(Stat[,1]/Stat[,2]) ) )   
PV2 <- c()
for(k in 1:p){
  f1 <- coxph( Surv(Y[Tr==1],Event[Tr==1])~X[Tr==1,k]+Z[Tr==1,] )
  f2 <- coxph( Surv(Y[Tr==-1],Event[Tr==-1])~X[Tr==-1,k]+Z[Tr==-1,] )
  a1 <- summary(f1)$coefficients[1] - summary(f2)$coefficients[1]
  a2 <- sqrt( summary(f1)$coefficients[1,3]^2 + summary(f2)$coefficients[1,3]^2 )
  PV2[k] <- 2*( 1-pnorm(abs(a1/a2)) )
}

qv1 <- list()
qv2 <- list()
for(k in 1:4){
  qv1[[k]] <- (1:p)[ qvalue(PV1,fdr=alpha.set[k])$significant ]
  qv2[[k]] <- (1:p)[ qvalue(PV2,fdr=alpha.set[k])$significant ]
}
```

