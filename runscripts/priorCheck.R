library(ReIns)
library(truncnorm)
library(here)
library(lattice)
source(here('src','growthModel.R'))
load(here('R_imagefiles/','variables.RData'))

# Observed Data Settings ####
timeRange = c(800,150)
# Extract Covariates
x1 = palm$PollenPerc
x2 = soi$SOIpr

# Model 1 Prior Check
nsim=1000
nt0 = rtexp(nsim,rate=10,endpoint=0.5)
a = rnorm(nsim,mean=1,sd=0.25)
r = rtexp(nsim,rate=20,endpoint=0.1)
param=data.frame(nt0=nt0,a=a,r=r,b1=0,b2=0)
ppcheck.model1 = matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nsim)

for (i in 1:nsim)
{
  ppcheck.model1[,i] = GrowthModel(nt0=param$nt0[i],r=param$r[i],timeRange=timeRange,b1=param$b1[i],b2=param$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
}

# Model 2 Prior Check
nsim=1000
nt0 = rtexp(nsim,rate=10,endpoint=0.5)
a = rnorm(nsim,mean=1,sd=0.25)
r = rtexp(nsim,rate=20,endpoint=0.1)
b1 = rtruncnorm(nsim,a=-0.015,b=0.015,mean=0,sd=0.005)
b2 = 0
param=data.frame(nt0=nt0,a=a,r=r,b1=b1,b2=b2)

# Check Routine
checkFun = function(x,x1,x2){all((1+x[4]*x1+x[5]*x2) > 0)}
param$check=apply(param,1,checkFun,x1=x1,x2=x2)

while(any(param$check==FALSE))
{
  n = sum(param$check==FALSE)
  i = which((param$check==FALSE))
  param$nt0[i] = rtexp(n,rate=10,endpoint=0.5)
  param$a[i] = rnorm(n,mean=1,sd=0.25)
  param$r[i] = rtexp(n,rate=20,endpoint=0.1)
  param$b1[i] = rtruncnorm(nsim,a=-0.015,b=0.015,mean=0,sd=0.005)
  param$b2[i] = 0
  param$check=apply(param,1,checkFun,x1=x1,x2=x2)
}


ppcheck.model2 = matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nsim)

for (i in 1:nsim)
{
  ppcheck.model2[,i] = GrowthModel(nt0=param$nt0[i],r=param$r[i],timeRange=timeRange,b1=param$b1[i],b2=param$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
}

# Model 3 Prior Check
nsim=1000
nt0 = rtexp(nsim,rate=10,endpoint=0.5)
r = rtexp(nsim,rate=20,endpoint=0.1)
b1 = 0
b2 = rtruncnorm(nsim,mean=0,a=-0.52,b=0.52,sd=0.2)
param=data.frame(nt0=nt0,a=a,r=r,b1=b1,b2=b2)

# Check Routine
checkFun = function(x,x1,x2){all((1+x[4]*x1+x[5]*x2) > 0)}
param$check=apply(param,1,checkFun,x1=x1,x2=x2)

while(any(param$check==FALSE))
{
  n = sum(param$check==FALSE)
  i = which((param$check==FALSE))
  param$nt0[i] = rtexp(n,rate=10,endpoint=0.5)
  param$a[i] = rnorm(n,mean=0,sd=0.25)
  param$r[i] = rtexp(n,rate=20,endpoint=0.1)
  param$b1[i] = 0
  param$b2[i] = rtruncnorm(nsim,mean=0,a=-0.52,b=0.52,sd=0.2)
  param$check=apply(param,1,checkFun,x1=x1,x2=x2)
}


ppcheck.model3 = matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nsim)

for (i in 1:nsim)
{
  ppcheck.model3[,i] = GrowthModel(nt0=param$nt0[i],r=param$r[i],timeRange=timeRange,b1=param$b1[i],b2=param$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
}

# Model 4 Prior Check
nsim=1000
nt0 = rtexp(nsim,rate=10,endpoint=0.5)
a = rnorm(nsim,mean=1,sd=0.25)
r = rtexp(nsim,rate=20,endpoint=0.1)
b1 = rnorm(nsim,mean=0,sd=0.01)
b2 = rnorm(nsim,mean=0,sd=0.2)
param=data.frame(nt0=nt0,a=1,r=r,b1=b1,b2=b2)

# Check Routine
checkFun = function(x,x1,x2){all((1+x[4]*x1+x[5]*x2) > 0)}
param$check=apply(param,1,checkFun,x1=x1,x2=x2)

while(any(param$check==FALSE))
{
  n = sum(param$check==FALSE)
  i = which((param$check==FALSE))
  param$nt0[i] = rtexp(n,rate=10,endpoint=0.5)
  param$a[i] = rnorm(n,mean=1,sd=0.25)
  param$r[i] = rtexp(n,rate=20,endpoint=0.1)
  param$b1[i] = rnorm(n,mean=0,sd=0.01)
  param$b2[i] = rnorm(n,mean=0,sd=0.2)
  param$check=apply(param,1,checkFun,x1=x1,x2=x2)
}

ppcheck.model4 = matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nsim)

for (i in 1:nsim)
{
  ppcheck.model4[,i] = GrowthModel(nt0=param$nt0[i],r=param$r[i],timeRange=timeRange,b1=param$b1[i],b2=param$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
}


jpeg(file=here('figures','temporary','prior_predictive_check.jpeg'),width = 7, height = 7,units='in',res=300)
par(mfrow=c(2,2),mar=c(5,4,2,1))
plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(ppcheck.model1),xlab='Cal BP',ylab='N',main='Model 1 (Logistic, K=1)')
apply(ppcheck.model1,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))

plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(ppcheck.model2),xlab='Cal BP',ylab='N',main='Model 2 (Logistic, K ~ Palm)')
apply(ppcheck.model2,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))

plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(ppcheck.model3),xlab='Cal BP',ylab='N',main='Model 3 (Logistic, K ~ SOI)')
apply(ppcheck.model3,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))

plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(ppcheck.model4),xlab='Cal BP',ylab='N',main='Model 4 (Logistic, K ~ Palm,SOI)')
apply(ppcheck.model4,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))
dev.off()
