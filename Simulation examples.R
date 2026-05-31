# Simulation examples for event truncated data


# Level properties
nsim <- 10000
beta <- 1
LRlp <- numeric(0)
KSlp <- numeric(0)
CvMlp <- numeric(0)
ADlp <- numeric(0)
ELRlp <- numeric(0)
set.seed(1)
for(j in 1:100){
  n <- (j)
  testobsLR <- numeric(nsim)
  testobsKS <- numeric(nsim)
  testobsCvM <- numeric(nsim)
  testobsAD <- numeric(nsim)
  testobsELR <- numeric(nsim)
  for(i in 1:nsim){
    tvector <- cumsum(rweibull(n,beta,1))
    testobsLR[i] <- LRfailtestobs(tvector,cv=findfailCV(tvector))
    testobsKS[i] <- KSfailtestobs(tvector,cv=findfailCV(tvector))
    testobsCvM[i] <- CvMfailtestobsnew(tvector,cv=findfailCV(tvector))
    testobsAD[i] <- ADfailtestobs(tvector,cv=findfailCV(tvector))
    testobsELR[i] <- ELRafix(tvector,0.5,cv=findfailCV(tvector))
  }
  LRlp[j] <- mean(abs(testobsLR)>1.96)
  KSlp[j] <- mean(testobsKS>1.358)
  CvMlp[j] <- mean(testobsCvM>0.461)
  ADlp[j] <- mean(testobsAD>2.492)
  ELRlp[j] <- mean(abs(testobsELR)>1.96)
}

plot(x=c(10:100),LRlp[10:100], type="l",ylim=c(0,0.3), col="red", lwd=2,
     main="HPP",
     xlab="number of events per simulation",
     ylab="rejection probability",xaxt="n")
lines(x=c(10:100),KSlp[10:100], col="orange", lwd=2)
lines(x=c(10:100),CvMlp[10:100], col="green", lwd=2)
lines(x=c(10:100),ADlp[10:100], col="blue", lwd=2)
lines(x=c(10:100),ELRlp[10:100], col="purple", lwd=2)
abline(0.05,0)
legend("topright", legend=c("LR","KS","CvM","AD","ELR"),
       col=c("red","orange","green","blue","purple"),lty=1, lwd=2)
axis(1, at=seq(10, 100, by=10), labels = c("10","20","30","40","50","60","70","80","90","100"))





# Power properties, monotonic trend

# Function to generate data event truncated Power Law Weibull TRP 
gendataPLWfeil <- function(alpha=1,beta=1,a=1,b=1,n=20){
  rvec=cumsum(rweibull(n,shape=beta,scale=alpha))
  tvec=(rvec/a)^(1/b)
  return(tvec)
}


nsim <- 10000
beta <- 1
b <- c(seq(0.25,0.75,0.08),seq(0.8,1.2,0.03),seq(1.25,1.5,0.1),seq(1.6,3,0.1),3.5,4,5,6)
LRpp <- numeric(length(b))
KSpp <- numeric(length(b))
CvMpp <- numeric(length(b))
ADpp <- numeric(length(b))
ELRpp <- numeric(length(b))
set.seed(1)
for(j in 1:length(b)){
  bj <- b[j]
  testobsLR <- numeric(nsim)
  testobsKS <- numeric(nsim)
  testobsCvM <- numeric(nsim)
  testobsAD <- numeric(nsim)
  testobsELR <- numeric(nsim)
  for(i in 1:nsim){
    tvector <- gendataPLWfeil(beta=beta,b=bj,n=30)
    tcv <- findfailCV(tvector)
    #tcv <- truecv(beta)
    #tcv <- cvdiv3(tvector)
    #tcv <- cvstar(tvector)
    #tcv <- cvb(tvector)
    #tvc <- cvspline(tvector)
    
    testobsLR[i] <- LRfailtestobs(tvector,cv=tcv)
    testobsKS[i] <- KSfailtestobs(tvector,cv=tcv)
    testobsCvM[i] <- CvMfailtestobsnew(tvector,cv=tcv)
    testobsAD[i] <- ADfailtestobs(tvector,cv=tcv)
    testobsELR[i] <- ELRafix(tvector,0.5,cv=tcv)
  }
  LRpp[j] <- mean(abs(testobsLR)>1.96)
  KSpp[j] <- mean(testobsKS>1.358)
  CvMpp[j] <- mean(testobsCvM>0.461)
  ADpp[j] <- mean(testobsAD>2.492)
  ELRpp[j] <- mean(abs(testobsELR)>1.96)
}

plot(b,LRpp,type="l",col="red",lwd=2,ylim=c(0,1), main="NHPP",
     xlab="b", ylab="rejection probability")
lines(b,KSpp,col="orange",lwd=2)
lines(b,CvMpp,col="green",lwd=2)
lines(b,ADpp,col="blue",lwd=2)
lines(b,ELRpp,col="purple",lwd=2)
segments(0.5, 0.05, 1.5, 0.05,lwd=2)
legend("bottomright", legend=c("LR","KS","CvM","AD","ELR"),
       col=c("red","orange","green","blue","purple"),lty=1, lwd=2)

