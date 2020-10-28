library(rcarbon)
library(ReIns)
library(foreach)
library(doParallel)
library(doSNOW)
library(here)
load('calibratedDates.RData') # Load Calibrated Dates and Custom Curves
source(here('ABC_routine','growthModels.R'))
source(here('ABC_routine','modelRunner.R'))
source(here('ABC_routine','calibrate2.R'))

# Observed Data Settings
x.norm=cal.combined.norm
x.nnorm=cal.combined.nnorm
bins = bins
ccurves = list(custom=customCurve)
timeRange = c(800,150)

# ABC Settings
nsim=100000
tol=0.1

# Parallelisation Settings
ncores = 
cl <- makeCluster(ncores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Generate Priors
k = rtexp(nsim,rate=5,endpoint=1)
r = rexp(nsim,rate=10)
param=data.frame(k=k,r=r)
# Main for Loop Starts Here

reslist <- foreach (i=1:nsim,.packages=c('rcarbon'),.options.snow = opts) %dopar%
  {
    candidateModel = logisticModel(k=param$k[i],r=param$r[i],timeRange=timeRange)
    res=modelRunner(model=candidateModel,x.norm=cal.combined.norm,x.nnorm=cal.combined.nnorm,bins=bins,customCurves=ccurves,timeRange=timeRange,simonly=FALSE)
    return(res)
  }

epsilon=do.call('rbind.data.frame',reslist)
result = cbind.data.frame(param,epsilon)
write.csv(result,file='abc_res.csv')
save(result,file='abc_res.RData')


