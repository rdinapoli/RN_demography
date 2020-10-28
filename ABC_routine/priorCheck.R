library(reIns)
library(here)
source(here('ABC_routine','growthModels.R'))

# Observed Data Settings
timeRange = c(800,150)

# Generate Priors
nsim = 1000
k = rtexp(nsim,rate=5,endpoint=1)
r = rexp(nsim,rate=10)
param=data.frame(k=k,r=r)

resMatrix = matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nsim)

for (i in 1:nsim)
{
  resMatrix[,i] = logisticModel(k=param$k[i],r=param$r[i],timeRange=timeRange)$PrDens
}

plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(resMatrix),xlab='Cal BP',ylab='Probability')
apply(resMatrix,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.1))

