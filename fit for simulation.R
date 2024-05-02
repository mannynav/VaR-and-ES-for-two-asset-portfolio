

#Packages.
library(Quandl)
library(quantmod)
library(copula)
library(QRM)
library(PerformanceAnalytics)
library(fitdistrplus)
library(descr)
library(tseries)
library(metRology)
library(mosaic)
library(VineCopula)

#Get yahoo finance. Function from quantmod.
setSymbolLookup(BP='yahoo',AMZ = 'yahoo')

#Get data. Function from quantmod.
data <- getSymbols(c('^DJI','^RUT'))

#Get log returns. Function from performance analytics.
DJI <- as.vector(na.omit(Return.calculate(DJI$DJI.Adjusted,method = c("log"))))
DJI <- DJI[1:4178]
RUT <- as.vector(na.omit(Return.calculate(RUT$RUT.Adjusted,method = c("log"))))
RUT <- RUT[1:4178]

#Check normality of returns.
chart.QQPlot(DJI,distribution = "norm")
jarque.bera.test(DJI)
shapiro.test(DJI)

chart.QQPlot(RUT,distribution = "norm")
jarque.bera.test(RUT)
shapiro.test(RUT)

#Check if returns could follow t-distribution.
summary(DJI)
skewness(DJI)
kurtosis(DJI) #fat tails

summary(RUT)
skewness(RUT)
kurtosis(RUT) #fat tails

#The t-distribution is similar to the Gaussian but with fatter tails, so we'll try to fit the data to the t distribution
fitDJIt<- fitdist(DJI, distr = "t.scaled",start = list(df=100,mean = mean(DJI),sd = sd(DJI)))
plot(fitDJIt)
fitRUTt<- fitdist(RUT, distr = "t.scaled",start = list(df=100,mean = mean(RUT),sd = sd(RUT)))
plot(fitRUTt)

#Plot histograms and their fits.
par(mfrow=c(1,2))
hist(DJI, freq = FALSE , ylim = c(0,60),breaks = 50, xlim = c(-0.08,0.08), xlab = "DJI Returns")
p<- seq(-0.07,0.07, by = 0.0001)
lines(p,dt.scaled(p,fitDJIt$estimate[1],mean = fitDJIt$estimate[2], sd = fitDJIt$estimate[3]), type = "l", col = "red", ylim = c(0,50))
hist(RUT, freq = FALSE , ylim = c(0,60),breaks = 50, xlim = c(-0.08,0.08),xlab = "RUT Returns")
p<- seq(-0.07,0.07, by = 0.0001)
lines(p,dt.scaled(p,fitRUTt$estimate[1],mean = fitRUTt$estimate[2], sd = fitRUTt$estimate[3]), type = "l", col = "red", ylim = c(0,50))
par(mfrow=c(1,1))

#Summary of fits.
summary(fitDJIt)
summary(fitRUTt)

#Kolomogorov-Smifnov tests.
ks.test(DJI,"pnorm", mean = mean(DJI), sd = sd(DJI))
ks.test(DJI, "pt.scaled", df = fitDJIt$estimate[1], mean = fitDJIt$estimate[2], sd = fitDJIt$estimate[3])
ks.test(RUT,"pnorm", mean = mean(RUT), sd = sd(RUT))
ks.test(RUT, "pt.scaled", df = fitRUTt$estimate[1], mean = fitRUTt$estimate[2], sd = fitRUTt$estimate[3])

#Plot returns together.
par(mfrow=c(1,2))
ret <- cbind(DJI,RUT)
plot(ret, xlab = "DJI Daily Returns", ylab = "RUT Daily Returns",main ="DJI and RUT Returns Scatter plot")

#Apply ECDF functions to each column.
ret <- apply(ret,2,edf,adjust = 1)
plot(ret, xlab = "ECDF of DJI Daily Returns", ylab = "ECDF of RUT Daily Returns", main = "DJI and RUT ECDF Scatterplot") #It looks like there could be some more negative dependence based on the bottom left being more dense than the top right
#Compute empirical cumulative distribution function of returns separately.
RUTCDF<- ecdf(RUT)(RUT)
DJICDF<-ecdf(DJI)(DJI)

#Fit copula to the ECDF of returns. From VineCopula package.
CopulaFit<-BiCopSelect(DJICDF,RUTCDF, family = NA)
summary(CopulaFit)

###### Fit for marginal distribution F of X1,X2,...,Xn ######
#Assign weights to returns.
w1<-0.7
w2<- 0.3
weightedReturns <- cbind(w1*DJI,w2*RUT)

#Get weighted portfolio values, these are X1,X,2,...,Xn.
portReturns <- rowSums(weightedReturns)
summary(portReturns)

#Check normality of portfolio.
chart.QQPlot(portReturns,distribution = "norm")
jarque.bera.test(portReturns)
shapiro.test(portReturns)
ks.test(portReturns,"pnorm", mean = mean(portReturns), sd = sd(portReturns))

#Check if portfolio follows t-distribution.
fitPortReturns<- fitdist(portReturns, distr = "t.scaled",start = list(df=100,mean = mean(portReturns),sd = sd(portReturns)))
hist(portReturns)
plot(fitPortReturns)

#P value is greater than 0.05, so we cannot reject the null that portReturns follows a t.scaled distribution.
ks.test(portReturns, "pt.scaled", df = fitPortReturns$estimate[1], mean = fitPortReturns$estimate[2], sd = fitPortReturns$estimate[3])

#Goodness of fit tests.
gofstat(fitPortReturns)

