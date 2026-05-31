# Improved coefficient of variation (CV) estimates


# Function to calculate estimated CV
findfailCV <- function(tvec){
  n <- length(tvec)
  xvec <- diff(c(0,tvec))
  sdest <- sqrt(var(xvec))
  mu <- mean(xvec)
  CV <- sdest/mu
  return(CV)
}


# CV from Weibull distribution
truecv <- function(beta){
  sqrt((gamma(1+2/beta)-gamma(1+1/beta)^2))/(gamma(1+1/beta))
}


# Method 1
# CV from dividing tvec into three, then taking the mean

cvdiv3 <- function(tvec){
  n <- length(tvec)
  
  tvec1 <- tvec[1:round(n/3)]
  tvec2 <- tvec[(round(n/3+1)):(round((2*n)/3))]
  tvec3 <- tvec[(round((2*n)/3+1)):n]
  
  
  n1 <- length(tvec1)
  xvec1 <- diff(c(0,tvec1))
  sdest1 <- sqrt(var(xvec1))
  mu1 <- mean(xvec1)
  CV1 <- sdest1/mu1
  
  n2 <- length(tvec2)
  xvec2 <- diff(c(max(tvec1),tvec2))
  sdest2 <- sqrt(var(xvec2))
  mu2 <- mean(xvec2)
  CV2 <- sdest2/mu2
  
  n3 <- length(tvec3)
  xvec3 <- diff(c(max(tvec2),tvec3))
  sdest3 <- sqrt(var(xvec3))
  mu3 <- mean(xvec3)
  CV3 <- sdest3/mu3
  
  return(mean(CV1,CV2,CV3))
}


# Method 2
# CV from dividing tvec into three, then scaling

cvstar <- function(tvec){
  n <- length(tvec)
  
  tvec1 <- tvec[1:round(n/3)]
  tvec2 <- tvec[(round(n/3+1)):(round((2*n)/3))]
  tvec3 <- tvec[(round((2*n)/3+1)):n]
  
  xvec1 <- diff(c(0,tvec1))
  xvec2 <- diff(c(max(tvec1),tvec2))
  xvec3 <- diff(c(max(tvec2),tvec3))
  
  m1 <- sum(xvec1)/length(xvec1)
  m2 <- sum(xvec2)/length(xvec2)
  m3 <- sum(xvec3)/length(xvec3)
  
  xvec1star <- xvec1/m1
  xvec2star <- xvec2/m2
  xvec3star <- xvec3/m3
  
  xvecstar <- c(xvec1star,xvec2star,xvec3star)
  sdest <- sqrt(var(xvecstar))
  mu <- mean(xvecstar)
  CV <- sdest/mu
  return(CV)
}


# Method 3
# CV using MLE b

cvb <- function(tvec){
  n <- length(tvec)
  i <- c(0:(n-1))
  Tn <- tvec[n]
  
  b <- -n/(sum(log(tvec[1:(n-1)]/Tn)))
  
  xvec <- diff(c(0,tvec^b))
  sdest <- sqrt(var(xvec))
  mu <- mean(xvec)
  CV <- sdest/mu
  return(CV)
}

# Method 4
# CV from using splines

cvspline <- function(tvec){
  tvec <- predict(lm(1:length(tvec) ~ ns(tvec,df=2)))
  xvec <- diff(c(0,tvec))
  CV <- sd(xvec)/mean(xvec)
  return(CV)
}

# Method 5
# Method 2 by splitting by time

cvstartid <- function(tvec,tau){
  n <- length(tvec)
  n1 <- length(tvec[tvec<(1*tau/4)])
  n2 <- length(tvec[tvec<(3*tau/4)])
  
  if(n1==1)
    n1 <- n1+1
  if(n2==n)
    n2 <- n2-1
  if(n1==n2)
    n2 <- n2+1
  
  tvec1 <- tvec[1:n1]
  tvec2 <- tvec[(n1+1):n2]
  tvec3 <- tvec[(n2+1):n]
  
  xvec1 <- diff(c(0,tvec1))
  xvec2 <- diff(c(max(tvec1),tvec2))
  xvec3 <- diff(c(max(tvec2),tvec3))
  
  m1 <- sum(xvec1)/length(xvec1)
  m2 <- sum(xvec2)/length(xvec2)
  m3 <- sum(xvec3)/length(xvec3)
  
  xvec1star <- xvec1/m1
  xvec2star <- xvec2/m2
  xvec3star <- xvec3/m3
  
  xvecstar <- c(xvec1star,xvec2star,xvec3star)
  sdest <- sqrt(var(xvecstar))
  mu <- mean(xvecstar)
  CV <- sdest/mu
  return(CV)
}

