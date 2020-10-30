library(rcarbon)
library(ReIns)
library(foreach)
library(doParallel)
library(doSNOW)
library(here)
library(coda)
library(here)
load(here('calibratedDates.RData')) # Load Calibrated Dates and Custom Curves
load(here('ABC_routine','abc_res.RData'))
source(here('ABC_routine','plotPosterior.R'))
source(here('ABC_routine','calibrate2.R'))
source(here('ABC_routine','growthModels.R'))
source(here('ABC_routine','modelRunner.R'))

# Observed Data Setting
x.norm=cal.combined.norm
x.nnorm=cal.combined.nnorm
bins = bins
ccurves = list(custom=customCurve)
timeRange = c(800,150)

# ABC Tolerance level ####
tol=0.001 # 1%, i.e. 1,000 best fit models
nsim = nrow(result)

# Retrieve Posterior Samples ####
post.cal.norm.euc = result[order(result$euc.cal.norm)[1:(nsim*tol)],]
post.cal.nnorm.euc = result[order(result$euc.cal.nnorm)[1:(nsim*tol)],]
post.uncal.norm.euc = result[order(result$euc.uncal.norm)[1:(nsim*tol)],]
post.uncal.nnorm.euc = result[order(result$euc.uncal.nnorm)[1:(nsim*tol)],]
post.cal.norm.nrmse = result[order(result$nrmse.cal.norm)[1:(nsim*tol)],]
post.cal.nnorm.nrmse = result[order(result$nrmse.cal.nnorm)[1:(nsim*tol)],]
post.uncal.norm.nrmse = result[order(result$nrmse.uncal.norm)[1:(nsim*tol)],]
post.uncal.nnorm.nrmse = result[order(result$nrmse.uncal.nnorm)[1:(nsim*tol)],]

posteriors.euc=list(post.uncal.norm.euc=post.uncal.norm.euc,post.uncal.nnorm.euc=post.uncal.nnorm.euc,post.cal.norm.euc=post.cal.norm.euc,post.cal.nnorm.euc=post.cal.nnorm.euc)

posteriors.nrmse=list(post.uncal.norm.nrmse=post.uncal.norm.nrmse,post.uncal.nnorm.nrmse=post.uncal.nnorm.nrmse,post.cal.norm.nrmse=post.cal.norm.nrmse,post.cal.nnorm.nrmse=post.cal.nnorm.nrmse)


ppchecks.euc = vector('list',length=4)
ppchecks.euc[1:4] = list(matrix(NA,nrow=abs(diff(timeRange))+1,ncol=tol*nsim))
ppchecks.nrmse = ppchecks.euc

print('Generating Posterior Predictive Check with EUC posteriors')
pb <- txtProgressBar(max = nsim*tol*4, style = 3)
for (p in 1:4)
{
  for (i in 1:100)
    {
      setTxtProgressBar(pb, (p-1)*100 + i)
      candidateModel = logisticModel(k=posteriors.euc[[p]]$k[i],r=posteriors.euc[[p]]$r[i],timeRange=timeRange)
      res=modelRunner(model=candidateModel,x.norm=cal.combined.norm,x.nnorm=cal.combined.nnorm,bins=bins,customCurves=ccurves,timeRange=timeRange,simonly=TRUE)
      ppchecks.euc[[p]][,i]=res[,p+1]
    }
}
close(pb)

print('Generating Posterior Predictive Check with NRSME posteriors')
pb <- txtProgressBar(max = nsim*tol*4, style = 3)

for (p in 1:4)
{
  for (i in 1:100)
  {
    setTxtProgressBar(pb, (p-1)*100 + i)
    candidateModel = logisticModel(k=posteriors.nrmse[[p]]$k[i],r=posteriors.nrmse[[p]]$r[i],timeRange=timeRange)
    res=modelRunner(model=candidateModel,x.norm=cal.combined.norm,x.nnorm=cal.combined.nnorm,bins=bins,customCurves=ccurves,timeRange=timeRange,simonly=TRUE)
    ppchecks.nrmse[[p]][,i]=res[,p+1]
  }
}

save.image(here('ABC_routine','ppcheckSPD.RData'))




