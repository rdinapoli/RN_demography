library(ReIns)
library(here)
library(lattice)
source(here('src','growthModel.R'))
load(here('R_imagefiles/','variables.RData'))

# Observed Data Settings ####
timeRange = c(800,150)
# Extract Covariates
x1 = palm$PollenPerc
x2 = soi$SOIpr

# Check Function
checkKFun = function(x,x1,x2){all((x[3]+x[4]*x1+x[5]*x2) > 0)}



# Model 1 Prior Check
nsim=1000
nt0 = rtexp(nsim,rate=10,endpoint=1)
r = rexp(nsim,rate=50) 
a = 1
b1 = 0
b2 = 0
param=data.frame(nt0=nt0,r=r,a=a,b1=b1,b2=b2)
ppcheck.model1_N = ppcheck.model1_PrDens = matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nsim)

param$check=apply(param,1,checkKFun,x1=x1,x2=x2)

while(any(param$check==FALSE))
{
  n = sum(param$check==FALSE)
  i = which((param$check==FALSE))
  param$nt0[i] = rtexp(n,rate=10,endpoint=1)
  param$r[i] = rexp(n,rate=50)
  param$a[i] = 1
  param$b1[i] = 0
  param$b2[i] = 0
  param$check=apply(param,1,checkFun,x1=x1,x2=x2)
}


for (i in 1:nsim)
{
  ppcheck.model1_N[,i] = GrowthModel(nt0=param$nt0[i],r=param$r[i],timeRange=timeRange,b1=param$b1[i],b2=param$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
  ppcheck.model1_PrDens[,i] = GrowthModel(nt0=param$nt0[i],r=param$r[i],timeRange=timeRange,b1=param$b1[i],b2=param$b2[i],x1=x1,x2=x2,a=1,raw = FALSE)$PrDens
}


# Model 2 Prior Check
nsim=1000
param$nt0[i] = rtexp(n,rate=10,endpoint=1)
param$r[i] = rexp(n,rate=50)
a = 1
b1 = rnorm(nsim,mean=0,sd=0.01)
b2 = 0
param=data.frame(nt0=nt0,r=r,a=a,b1=b1,b2=b2)
ppcheck.model2_N = ppcheck.model2_PrDens = matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nsim)

param$check=apply(param,1,checkKFun,x1=x1,x2=x2)

while(any(param$check==FALSE))
{
  n = sum(param$check==FALSE)
  i = which((param$check==FALSE))
  param$nt0[i] = rtexp(n,rate=10,endpoint=1)
  param$r[i] = rexp(n,rate=50)
  param$a[i] = 1
  param$b1[i] = rnorm(n,mean=0,sd=0.01)
  param$b2[i] = 0
  param$check=apply(param,1,checkKFun,x1=x1,x2=x2)
}

for (i in 1:nsim)
{
  ppcheck.model2_N[,i] = GrowthModel(nt0=param$nt0[i],r=param$r[i],timeRange=timeRange,b1=param$b1[i],b2=param$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
  ppcheck.model2_PrDens[,i] = GrowthModel(nt0=param$nt0[i],r=param$r[i],timeRange=timeRange,b1=param$b1[i],b2=param$b2[i],x1=x1,x2=x2,a=1,raw = FALSE)$PrDens
}

# Model 3 Prior Check
nsim=1000
param$nt0[i] = rtexp(n,rate=10,endpoint=1)
param$r[i] = rexp(n,rate=50)
a = 1
b1 = 0
b2 = rnorm(nsim,mean=0,sd=0.2)
param=data.frame(nt0=nt0,r=r,a=a,b1=b1,b2=b2)
ppcheck.model3_N = ppcheck.model3_PrDens = matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nsim)

param$check=apply(param,1,checkKFun,x1=x1,x2=x2)

while(any(param$check==FALSE))
{
  n = sum(param$check==FALSE)
  i = which((param$check==FALSE))
  param$nt0[i] = rtexp(n,rate=10,endpoint=1)
  param$r[i] = rexp(n,rate=50)
  param$a[i] = 1
  param$b1[i] = 0
  param$b2[i] = rnorm(n,mean=0,sd=0.2)
  param$check=apply(param,1,checkKFun,x1=x1,x2=x2)
}


for (i in 1:nsim)
{
  ppcheck.model3_N[,i] = GrowthModel(nt0=param$nt0[i],r=param$r[i],timeRange=timeRange,b1=param$b1[i],b2=param$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
  ppcheck.model3_PrDens[,i] = GrowthModel(nt0=param$nt0[i],r=param$r[i],timeRange=timeRange,b1=param$b1[i],b2=param$b2[i],x1=x1,x2=x2,a=1,raw = FALSE)$PrDens
}

# Model 4 Prior Check
nsim=1000
param$nt0[i] = rtexp(n,rate=10,endpoint=1)
param$r[i] = rexp(n,rate=50)
a = 1
b1 = rnorm(nsim,mean=0,sd=0.01)
b2 = rnorm(nsim,mean=0,sd=0.2)
param=data.frame(nt0=nt0,r=r,a=a,b1=b1,b2=b2)
ppcheck.model4_N = ppcheck.model4_PrDens = matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nsim)

param$check=apply(param,1,checkKFun,x1=x1,x2=x2)

while(any(param$check==FALSE))
{
  n = sum(param$check==FALSE)
  i = which((param$check==FALSE))
  param$nt0[i] = rtexp(n,rate=10,endpoint=1)
  param$r[i] = rexp(n,rate=50)
  param$a[i] = 1
  param$b1[i] = rnorm(n,mean=0,sd=0.01)
  param$b2[i] = rnorm(n,mean=0,sd=0.2)
  param$check=apply(param,1,checkKFun,x1=x1,x2=x2)
}

for (i in 1:nsim)
{
  ppcheck.model4_N[,i] = GrowthModel(nt0=param$nt0[i],r=param$r[i],timeRange=timeRange,b1=param$b1[i],b2=param$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
  ppcheck.model4_PrDens[,i] = GrowthModel(nt0=param$nt0[i],r=param$r[i],timeRange=timeRange,b1=param$b1[i],b2=param$b2[i],x1=x1,x2=x2,a=1,raw = FALSE)$PrDens
  
}


# Plot Prior Check ####


jpeg(file=here('figures','temporary','prior_predictive_check_N.jpeg'),width = 7, height = 7,units='in',res=300)
par(mfrow=c(2,2),mar=c(5,4,2,1))
plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(ppcheck.model1_N),xlab='Cal BP',ylab='N',main='Model 1 (Logistic, K=1)')
apply(ppcheck.model1_N,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))

plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(ppcheck.model2_N),xlab='Cal BP',ylab='N',main='Model 2 (Logistic, K ~ Palm)')
apply(ppcheck.model2_N,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))

plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(ppcheck.model3_N),xlab='Cal BP',ylab='N',main='Model 3 (Logistic, K ~ SOI)')
apply(ppcheck.model3_N,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))

plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(ppcheck.model4_N),xlab='Cal BP',ylab='N',main='Model 4 (Logistic, K ~ Palm,SOI)')
apply(ppcheck.model4_N,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))
dev.off()

jpeg(file=here('figures','temporary','prior_predictive_check_PrDens.jpeg'),width = 7, height = 7,units='in',res=300)
par(mfrow=c(2,2),mar=c(5,4,2,1))
plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(ppcheck.model1_PrDens),xlab='Cal BP',ylab='PrDens',main='Model 1 (Logistic, K=1)')
apply(ppcheck.model1_PrDens,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))

plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(ppcheck.model2_PrDens),xlab='Cal BP',ylab='PrDens',main='Model 2 (Logistic, K ~ Palm)')
apply(ppcheck.model2_PrDens,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))

plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(ppcheck.model3_PrDens),xlab='Cal BP',ylab='PrDens',main='Model 3 (Logistic, K ~ SOI)')
apply(ppcheck.model3_PrDens,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))

plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(ppcheck.model4_PrDens),xlab='Cal BP',ylab='PrDens',main='Model 4 (Logistic, K ~ Palm,SOI)')
apply(ppcheck.model4_PrDens,2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))
dev.off()

