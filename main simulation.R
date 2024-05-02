

#Packages.
library(actuar)
library(pracma)

#Number of simulations 250,000.
R<- 2500

#Create vector to store sum of X1+X2+...+X_{n-1}.
vectS <- c(0)

#Number of days holding portfolio.
N <- 30

for (i in 1:R) {
  
  #Generate V1,V2,...Vn
  V <- rexp(N,1)
  
  #Set seed for replication.
  set.seed(i+1)
  
  #Random variable Z >0.
  #Z <- rexp(1,1)
  Zclay <- rgamma(1,1/2.66)
  
  #V divided by Z.
  VdZ <- V/Zclay
  
  #Generator functions for Archimedean copulas.
  psiClay <- 1/((1+VdZ*2.66)^(1/2.66))
  #psiBB7 <- 1 - (1 - (1/(VdZ+1))^(1/1.85)) ^ (1/2.86)
  #psiBB1 <- (1/(VdZ^(1/2.04)+1))^(1/0.83)
  #psiJoe <- 1 - (1-exp(-VdZ))^(1/3.23)
  
  #Get inverse of marginal distribution F by input one of the generator functions from above.
  Finvs <- 1 - as.vector(pt.scaled(psiClay,fitPortReturns$estimate[1],fitPortReturns$estimate[2],fitPortReturns$estimate[3]))
  
  #Vector to store the sum of n-1 values of Finvs.
  vectS[i] <- sum(Finvs[1:length(Finvs)-1])
}

#Check out data.
summary(vectS)
hist(vectS,xlim = c(0,15),breaks = 50)

#Parameters for BB1 copula.
thetaBB1<-.83
betaBB1 <- 2.04

#Parameters for BB7 copula.
thetaBB7<-2.86
betaBB7 <- 1.87

#Parameter for Clayton copula.
thetaClay<-2.66

#Parameter for Joe copula.
thetaJoe <-3.23

#Estimated CDF.
Fncond <- function(x,vectS) {
  s<-c(0) 
  for (i in 1:R) {
    
    #Set seed for replication of Z and Zclay. Used for example 4.11.
    set.seed(i+1)
    #Z <- rexp(1,1)
    Zclay <- rgamma(1,1/2.66)
    
    #F(x-S_{n-1})
    inputF <- pt.scaled(x-vectS[i],fitPortReturns$estimate[1],fitPortReturns$estimate[2],fitPortReturns$estimate[3])
    
    #Inverse generators for each copula.
    #invGeneratorBB1 <- (inputF^(-thetaBB1)-1)^(betaBB1)
    #invGeneratorBB7 <- ((1-(1-inputF)^thetaBB7)^(-betaBB7)) - 1
    invGeneratorClay <- (1/((inputF)^(thetaClay)) - 1)/thetaClay
    #invGeneratorJoe <- (-log(1-(1-inputF)^(thetaJoe)))
    
    s[i] <- exp(-Zclay*invGeneratorClay)
  }
  return(mean(s))
}
#Compute the inverse of Fncond.
Fncondinvs <- function(x) {
  sINV<-c(0) 
  for (i in 1:R) {
    
    #Set seed for replication of Z and Zclay. Used for example 4.11.
    set.seed(i+1)
    #Z <- rexp(1,1)
    Zclay <- rgamma(1,1/2.66)
    
    #F(x-S_{n-1})
    inputF <- pt.scaled(x-vectS[i],fitPortReturns$estimate[1],fitPortReturns$estimate[2],fitPortReturns$estimate[3])
    
    #Inverse generators for each copula.
    #invGeneratorBB1 <- (inputF^(-thetaBB1)-1)^(betaBB1)
    #invGeneratorBB7 <- ((1-(1-inputF)^thetaBB7)^(-betaBB7)) - 1
    invGeneratorClay <- (1/((inputF)^(thetaClay)) - 1)/thetaClay
    #invGeneratorJoe <- (-log(1-(1-inputF)^(thetaJoe)))
    
    sINV[i] <- exp(-Zclay*invGeneratorClay)
  }
  return(1-mean(sINV))
}
#Vectorize Fncondinvs.
FninvscondVect <- Vectorize(Fncondinvs, "x")

#Quantile.
alpha <- 0.95

#Solve Fncond(1/(1-q)-1,vectS) = 0.95 for q, store as qu
qu <- uniroot(f = function(q) alpha - Fncond(1/(1 - q) - 1,vectS), interval = c(0,1),check.conv = TRUE)$root

#This is the VaR estimate
VaR <- 1 /(1 - qu) - 1

#Integrate the inverse distribution function from the VaR estimate to infinity to get the estimate for m.
fun <- function(x) FninvscondVect(x)
m <- integral(fun,VaR,Inf)

#Compute expected shortfall.
ES <- VaR + (m)/(1-alpha)
