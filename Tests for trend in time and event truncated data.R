# Tests for trend in time and event truncated data


# Explanation of the different inputs:
## tvec: a vector of the times until events
## cv: coefficient of variation
## a: the point where the integral gets split
## tau: the time truncation point
## sigma: not relevant as we will be using findfailCV() instead of findCV() which can be found in:
### Methods of estimating the coefficient of variation.R

# Test statistics for the tests for trend in event truncated data

# LR
LRfailtestobs <- function(tvec,cv=1){
  n <- length(tvec)
  xvec <- diff(c(0,tvec))
  i <- c(0:(n-1))
  Tn <- tvec[n]
  X1 <- xvec[1]
  LR <- (-1/cv)*(sum(tvec[1:(n-1)]) - ((n-1)/2)*Tn )/(Tn*sqrt(n/12))
  return(LR)
}

# KS
KSfailtestobs <- function(tvec,cv=1){
  n <- length(tvec)
  xvec <- diff(c(0,tvec))
  i <- c(0:(n-1))
  Tn <- tvec[n]
  X1 <- xvec[1]
  KS <- sqrt(n)/cv*max(abs(tvec/Tn-i/n))
  return(KS)
}

# CvM based on continuous V0
CvMfailtestobsnew <- function(tvec,cv=1){
  n <- length(tvec)
  xvec <- diff(c(0,tvec))
  i <- c(0:(n-1))
  Tn <- tvec[n]
  X1 <- xvec[1]
  CvM <- (1/(cv^2))*(sum( -(2*i+1)/(n*Tn)*tvec + 1/(Tn^2)*(tvec^2+tvec*xvec) - (3*i+2)/(3*n*Tn)*xvec + 1/(3*Tn^2)*xvec^2 ) + n/3)
  CvM <- CvM*exp(-8/(cv^1.5*n)) # adjustment
  return(CvM)
}

# CvM based on discrete WO
CvMdiscrete <- function(tvec,cv=1){
  n <- length(tvec)
  i <- c(0:(n-1))
  Tn <- tvec[n]
  CvM <- 1/cv^2*(1/Tn*(sum((tvec^2)/Tn-(2*i+1)/n*tvec)) + n/3)
  return(CvM)
}

# AD based on continuous V0
ADfailtestobs <- function(tvec,cv=1){
  n <- length(tvec)
  xvec <- diff(c(0,tvec))
  mu <- mean(xvec)
  sdest <- sqrt(var(xvec))
  konst <- n/cv^2
  iseq <- 2:(n-1)
  imseq <- 1:(n-2)
  Niseq <- (n-1):2
  Nimseq <- (n-2):1
  qi <- tvec[2:(n-1)]/tvec[n]-iseq*xvec[2:(n-1)]/tvec[n]
  ri <- n*xvec[2:(n-1)]/tvec[n]-1
  qn <- 1-n*xvec[n]/tvec[n]
  r1 <- n*xvec[1]/tvec[n]-1
  rn <- n*xvec[n]/tvec[n]-1
  forsteledd <- r1^2*log(n/(n-1))-r1^2/n 
  sisteledd <- qn^2*log(n/(n-1))-rn^2/n
  sumledd <- sum(qi^2*log(iseq/imseq)+(qi+ri)^2*log(Niseq/Nimseq)-ri^2/n)
  AD <- konst*(forsteledd+sumledd+sisteledd)
  return(AD)
}

# ELR based on continuous V0
ELRafix <- function(tvec,a=0.5,cv=1){
  n <- length(tvec)
  m <- round(n*a)
  xvec <- diff(c(0,tvec))
  i1 <- c(0:(m-1))
  i2 <- c(m:(n-1))
  Tn <- tvec[n]
  X1 <- xvec[1]
  tvec1 <- tvec[i1]
  tvec2 <- tvec[i2]
  i1new <- c(1:m)
  i2new <- c((m+1):n)
  xvec1 <- xvec[i1new]
  xvec2 <- xvec[i2new]
  sum1 <- 1/(n*Tn)*(sum(tvec1)+sum((1/2)*xvec1))
  sum2 <- 1/(n*Tn)*(sum(tvec2)+sum((1/2)*xvec2))
  a <- length(xvec1)/n # adjustment
  ELR <- 1/(sqrt(1/12 - a^2*(1-a)^2))*sqrt(n)/cv*(sum1-sum2 - a^2 + 1/2)
  return(ELR)
}

# ELR based on discrete W0
ELRafixdiscrete <- function(tvec,a=0.5,cv=1){
  n <- length(tvec)
  m <- round(n*a)
  i1 <- c(0:(m-1))
  i2 <- c(m:(n-1))
  Tn <- tvec[n]
  tvec1 <- tvec[i1]
  tvec2 <- tvec[i2]
  i1new <- c(1:m)
  i2new <- c((m+1):n)
  xvec1 <- xvec[i1new]
  sum1 <- 1/(n*Tn)*(sum(tvec1))
  sum2 <- 1/(n*Tn)*(sum(tvec2))
  a <- length(xvec1)/n # adjustment
  ELR <- 1/(sqrt(1/12 - a^2*(1-a)^2))*sqrt(n)/cv*(sum1-sum2 - a^2 + 1/2)
  return(ELR)
}


# Test statistics for the tests for trend in event truncated data
# For our purpose, findCV() = findfailCV() which can be found in:
# Methods of estimating the coefficient of variation.R

# LR
LRtestobs <- function(tvec,tau,sigma,cv){
  CV <- findCV(tvec,tau,sigma,cv)
  n <- length(tvec)
  LR <- (sqrt(12)/(CV*tau*sqrt(n)))*(sum(tvec)-n*tau/2)
  return(LR)
}

# KS
KStestobs <- function(tvec,tau,sigma,cv){
  CV <- findCV(tvec,tau,sigma,cv)
  n <- length(tvec)
  xvec <- diff(c(0,tvec))
  t1vec <- n*c(0,tvec)/tau
  t2vec <- n*c(tvec,tau)/tau
  nvec <- 0:n
  KS <- (1/(CV*sqrt(n)))*max(c(max(abs(nvec-t1vec))),c(max(abs(nvec-t2vec))))
  return(KS)
}  

# CvM
CvMtestobs <- function(tvec,tau,sigma,cv){
  CV <- findCV(tvec,tau,sigma,cv)
  n <- length(tvec)
  xvec <- diff(c(0,tvec))
  konst <- 1/(CV^2*n)
  indseq <- 0:(n-1)
  sumledd <- sum(indseq^2*xvec/tau)-n*sum(indseq*(tvec^2-c(0,tvec[1:(n-1)])^2)/tau^2)
  sisteledd <- n^2/3+n^2*(tvec[n]^2/tau^2-tvec[n]/tau)
  CV <- konst*(sumledd+sisteledd)
  return(CV)
}  

# AD
ADtestobs <- function(tvec,tau,sigma="s",cv=1){
  CV <- findCV(tvec,tau,sigma,cv)
  n <- length(tvec)
  if(n==1)
    return((log(tau^2/((tau-tvec[1])*tvec[1]))-1)/CV^2)
  konst <- 1/(CV^2*n)
  iseq <- 1:(n-1)
  iseq2 <- iseq^2
  Niseq <- (n-1):1
  Niseq2 <- Niseq^2
  tip <- tvec[2:n]
  ti <- tvec[1:(n-1)]
  lntfrac <- log(tip/ti)
  lntautfrac <- log((tau-ti)/(tau-tip))
  sumledd <- sum(Niseq2*lntautfrac+iseq2*lntfrac)
  sisteledd <- n^2*(log(tau/(tau-tvec[1]))+log(tau/tvec[n])-1)
  AD <- konst*(sumledd+sisteledd)
  return(AD)
}

# ELR
ELRtestobs <- function(tvec,tau,sigma,cv,a){
  CV <- findCV(tvec,tau,sigma,cv)
  n <- length(tvec)
  konst <- 1/(CV*tau*sqrt(n)*sqrt((1/12)-a^2*(1-a)^2))
  sumledd <- sum(abs(tvec-a*tau))
  sisteledd <- (0.5-a*(1-a))*n*tau  
  ELR <- konst*(sumledd-sisteledd)
  return(ELR)
}

