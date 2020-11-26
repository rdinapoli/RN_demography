library(rcarbon)
library(ReIns)
library(foreach)
library(doParallel)
library(doSNOW)
library(here)
load(here('R_imagefiles','variables.RData')) # Load Calibrated Dates and Custom Curves
source(here('src','growthModel.R'))
source(here('src','modelRunner.R'))
source(here('src','randomThin.R'))
source(here('src','calibrate2.R'))

# Observed Data Settings
x.norm=cal.combined.norm
x.nnorm=cal.combined.nnorm
bins = bins
ccurves = list(custom=customCurve)
timeRange = c(800,150)

# Extract Covariates
x1 = palm$PollenPerc
x2 = soi$SOIpr

# ABC Settings
nsim=250000

# Parallelisation Settings
ncores = 10
cl <- makeCluster(ncores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Generate Priors
set.seed(123)
nt0 = rtexp(nsim,rate=10,endpoint=1)
r = rexp(nsim,rate=50)
a = 1
b1 = 0
b2 = rnorm(nsim,mean=0,sd=0.2)
param=data.frame(nt0=nt0,r=r,a=a,b1=b1,b2=b2)

# Check Whether Prior Combination generates negative K
checkKFun = function(x,x1,x2){all((x[3]+x[4]*x1+x[5]*x2) > 0)}
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
param = param[,-ncol(param)]

# Main for Loop Starts Here

reslist <- foreach (i=1:nsim,.packages=c('rcarbon'),.options.snow = opts) %dopar%
  {
    set.seed(i)
    candidateModel = GrowthModel(nt0=param$nt0[i],r=param$r[i],timeRange=timeRange,b1=param$b1[i],b2=param$b2[i],x1=x1,x2=x2,a=param$a[i])
    res=modelRunner(model=candidateModel,x.norm=cal.combined.norm,x.nnorm=cal.combined.nnorm,bins=bins,customCurves=ccurves,timeRange=timeRange,simonly=FALSE)
    return(res)
  }

epsilon=do.call('rbind.data.frame',reslist)
model3 = cbind.data.frame(param,epsilon)
save(model3,file=here('R_imagefiles','abc_model3.RData'))


