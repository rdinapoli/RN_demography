library(ReIns)
library(here)
source(here('ABC_routine','growthModels.R'))
load(here('calibratedDates.RData'))
# Observed Data Settings ####
timeRange = c(800,150)

# Model 1 (Logistic Growth) ####
nsim = 1000
K = rep(1,nsim)
nt0 = rtexp(nsim,rate=5,endpoint=1)
r = rexp(nsim,rate=5)
param=data.frame(nt0=nt0,K=K,r=r)
resMatrix = matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nsim)

for (i in 1:nsim)
{
  resMatrix[,i] = logisticGrowthModel(nt0=param$nt0[i],r=param$r[i],K=param$K[i],timeRange=timeRange)$PrDens
}

plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(resMatrix),xlab='Cal BP',ylab='Probability')
apply(resMatrix,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.1))

# Model 2 (Logistic Growth + Palm) ####
set.seed(123)
nsim = 1000
a = rep(1,nsim)
b = rnorm(nsim,mean=0,sd=0.2)
nt0 = rtexp(nsim,rate=5,endpoint=1)
r = rexp(nsim,rate=5)
param=data.frame(nt0=nt0,a=a,b=b,r=r)
resMatrix = matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nsim)

for (i in 1:nsim)
{
  resMatrix[,i] = SingleCovariateLogisticGrowthModel(nt0=param$nt0[i],r=param$r[i],x=palmData$Palm,a=param$a[i],b=param$b[i],timeRange=timeRange)$PrDens
}

plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=c(0,0.005),xlab='Cal BP',ylab='Probability')
apply(resMatrix,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.1))
lines(x=timeRange[1]:timeRange[2],resMatrix[,164],col=2,lwd=2)
param[164,]
lines(x=timeRange[1]:timeRange[2],resMatrix[,38],col=4,lwd=2)
param[38,]


# Model 3 (Logistic Growth + SOI) ####
set.seed(123)
nsim = 1000
a = rep(1,nsim)
b = rnorm(nsim,mean=0,sd=0.4)
nt0 = rtexp(nsim,rate=5,endpoint=1)
r = rexp(nsim,rate=5)
param=data.frame(nt0=nt0,a=a,b=b,r=r)
resMatrix = matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nsim)

for (i in 1:nsim)
{
  resMatrix[,i] = SingleCovariateLogisticGrowthModel(nt0=param$nt0[i],r=param$r[i],x=soi$SOIpr,a=param$a[i],b=param$b[i],timeRange=timeRange)$PrDens
}

plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=c(0,0.005),xlab='Cal BP',ylab='Probability')
apply(resMatrix,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.1))
lines(x=timeRange[1]:timeRange[2],resMatrix[,997],col=2,lwd=2)
param[997,]
lines(x=timeRange[1]:timeRange[2],resMatrix[,998],col=4,lwd=2)
param[998,]

# Model 4 (Logistic Growth + Palm + SOI) ####
set.seed(123)
nsim = 1000
a = rep(1,nsim)
b1 = rnorm(nsim,mean=0,sd=0.2)
b2 = rnorm(nsim,mean=0,sd=0.4)
nt0 = rtexp(nsim,rate=5,endpoint=1)
r = rexp(nsim,rate=5)
param=data.frame(nt0=nt0,a=a,b1=b1,b2=b2,r=r)
resMatrix = matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nsim)

for (i in 1:nsim)
{
  resMatrix[,i] = DoubleCovariateLogisticGrowthModel(nt0=param$nt0[i],r=param$r[i],x1=palmData$Palm,x2=soi$SOIpr,a=param$a[i],b1=param$b1[i],b2=param$b2[i],timeRange=timeRange)$PrDens
}

plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=c(0,0.005),xlab='Cal BP',ylab='Probability')
apply(resMatrix,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.1))
lines(x=timeRange[1]:timeRange[2],resMatrix[,995],col=2,lwd=2)
param[995,]
lines(x=timeRange[1]:timeRange[2],resMatrix[,285],col=4,lwd=2)
param[285,]
