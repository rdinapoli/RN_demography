library(coda)
library(here)
load(here('ABC_routine','abc_res.RData'))
load(here('ABC_routine','ppcheckSPD.RData'))
source(here('ABC_routine','plotPosterior.R'))
source(here('ABC_routine','growthModels.R'))

# ABC Tolerance level ####
tol=0.01 # 1%, i.e. 1,000 best fit models
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

# Compute 90% HPDI Intervals
options(scipen = 9999)
pdf(file = here('ABC_routine','posterior_logistic_euc.pdf'),width = 5,height = 8)
par(mfrow=c(4,2),mar=c(5,4,2,1))
plotPosterior(post.cal.norm.euc$k,main='Cal Norm Euc (k)',xlim=c(0,0.5))
plotPosterior(post.cal.norm.euc$r,main='Cal Norm Euc (r)',xlim=c(0,0.1))
plotPosterior(post.cal.nnorm.euc$k,main='Cal Non-Norm Euc (k)',xlim=c(0,0.5))
plotPosterior(post.cal.nnorm.euc$r,main='Cal Non-Norm Euc (r)',xlim=c(0,0.1))
plotPosterior(post.uncal.norm.euc$k,main='Uncal Norm Euc (k)',xlim=c(0,0.5))
plotPosterior(post.uncal.norm.euc$r,main='Uncal Norm Euc (r)',xlim=c(0,0.1))
plotPosterior(post.uncal.nnorm.euc$k,main='Uncal Non-Norm Euc (k)',xlim=c(0,0.5))
plotPosterior(post.uncal.nnorm.euc$r,main='Uncal Non-Norm Euc (r)',xlim=c(0,0.1))
dev.off()

pdf(file = here('ABC_routine','posterior_logistic_nrmse.pdf'),width = 5,height = 8)
par(mfrow=c(4,2),mar=c(5,4,2,1))
plotPosterior(post.cal.norm.nrmse$k,main='Cal Norm NRMSE (k)',xlim=c(0,0.5))
plotPosterior(post.cal.norm.nrmse$r,main='Cal Norm NRMSE (r)',xlim=c(0,0.1))
plotPosterior(post.cal.nnorm.nrmse$k,main='Cal Non-Norm NRMSE (k)',xlim=c(0,0.5))
plotPosterior(post.cal.nnorm.nrmse$r,main='Cal Non-Norm NRMSE (r)',xlim=c(0,0.1))
plotPosterior(post.uncal.norm.nrmse$k,main='Uncal Norm NRMSE (k)',xlim=c(0,0.5))
plotPosterior(post.uncal.norm.nrmse$r,main='Uncal Norm NRMSE (r)',xlim=c(0,0.1))
plotPosterior(post.uncal.nnorm.nrmse$k,main='Uncal Non-Norm NRMSE (k)',xlim=c(0,0.5))
plotPosterior(post.uncal.nnorm.nrmse$r,main='Uncal Non-Norm NRMSE (r)',xlim=c(0,0.1))
dev.off()


# Posterior Predictive Checks (raw non-spd based) ####
png(file = here('ABC_routine','ppcheck_raw_euc.png'),width = 550,height = 600)
par(mfrow=c(2,2),mar=c(5,4,2,1))
postPredictiveCheck(model=logisticModel,params=post.cal.norm.euc[,1:2],timeRange = c(800,150),alpha=0.2,type='spaghetti',main='Cal-Norm-Euc',spaghettiSize=500)
postPredictiveCheck(model=logisticModel,params=post.cal.nnorm.euc[,1:2],timeRange = c(800,150),alpha=0.2,type='spaghetti',main='Cal-Non-Norm-Euc',spaghettiSize=500)
postPredictiveCheck(model=logisticModel,params=post.uncal.norm.euc[,1:2],timeRange = c(800,150),alpha=0.2,type='spaghetti',main='Uncal-Norm-Euc',spaghettiSize=500)
postPredictiveCheck(model=logisticModel,params=post.uncal.nnorm.euc[,1:2],timeRange = c(800,150),alpha=0.2,type='spaghetti',main='Uncal-Non-Norm-Euc',spaghettiSize=500)
dev.off()

png(file = here('ABC_routine','ppcheck_raw_nrmse.png'),width = 550,height = 600)
par(mfrow=c(2,2),mar=c(5,4,2,1))
postPredictiveCheck(model=logisticModel,params=post.cal.norm.nrmse[,1:2],timeRange = c(800,150),alpha=0.2,type='spaghetti',main='Cal-Norm-nrmse',spaghettiSize=500)
postPredictiveCheck(model=logisticModel,params=post.cal.nnorm.nrmse[,1:2],timeRange = c(800,150),alpha=0.2,type='spaghetti',main='Cal-Non-Norm-nrmse',spaghettiSize=500)
postPredictiveCheck(model=logisticModel,params=post.uncal.norm.nrmse[,1:2],timeRange = c(800,150),alpha=0.2,type='spaghetti',main='Uncal-Norm-nrmse',spaghettiSize=500)
postPredictiveCheck(model=logisticModel,params=post.uncal.nnorm.nrmse[,1:2],timeRange = c(800,150),alpha=0.2,type='spaghetti',main='Uncal-Non-Norm-nrmse',spaghettiSize=500)
dev.off()

# Posterior Predictive Checks (spd based) ####
index=randomThin(bins)
obs.norm =spd(cal.combined.norm[index],datenormalised=FALSE, spdnormalised = TRUE,timeRange=c(800,150),verbose=FALSE)
obs.nnorm =spd(cal.combined.nnorm[index],datenormalised=FALSE, spdnormalised = TRUE,timeRange=c(800,150),verbose=FALSE)


png(file = here('ABC_routine','ppcheck_spd_euc.png'),width = 550,height = 600)
par(mfrow=c(2,2),mar=c(5,4,2,1))
plotPPCheckSPD(obs.norm,ppchecks.euc[[1]],main='Uncalsample Normalised (EUC)')
plotPPCheckSPD(obs.nnorm,ppchecks.euc[[2]],main='Uncalsample Non-Normalised (EUC)')
plotPPCheckSPD(obs.norm,ppchecks.euc[[3]],main='Calsample Normalised (EUC)')
plotPPCheckSPD(obs.nnorm,ppchecks.euc[[4]],main='Calsample Non-Normalised (EUC)')
dev.off()

png(file = here('ABC_routine','ppcheck_spd_euc.png'),width = 550,height = 600)
par(mfrow=c(2,2),mar=c(5,4,2,1))
plotPPCheckSPD(obs.norm,ppchecks.nrmse[[1]],main='Uncalsample Normalised (NRMSE)')
plotPPCheckSPD(obs.nnorm,ppchecks.nrmse[[2]],main='Uncalsample Non-Normalised (NRMSE)')
plotPPCheckSPD(obs.norm,ppchecks.nrmse[[3]],main='Calsample Normalised (NRMSE)')
plotPPCheckSPD(obs.nnorm,ppchecks.nrmse[[4]],main='Calsample Non-Normalised (NRMSE)')
dev.off()



