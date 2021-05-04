library(here)
library(ggplot2)
library(latex2exp)
library(tidyverse)
library(RColorBrewer)
library(rcarbon)

### Figure 3 Fitted Model ####
load(here('R_imagefiles/','variables.RData'))
source(here('src','growthModel.R'))
load(here('R_imagefiles','ABC_posteriors.RData'))
source(here('src','plotFunctions.R'))

# General Settings 
timeRange = c(800,150)
# Extract Covariates
x1 = palm$PollenPerc
x2 = soi$SOIpr

### PrDens Version 
# Prepare Matrices
model1.matrix = model2.matrix = model3.matrix = model4.matrix = matrix(NA,ncol=nrow(model1.posterior$euc.uncal.nnorm),nrow=length(timeRange[1]:timeRange[2]))

for (i in 1:nrow(model1.posterior$euc.uncal.nnorm))
{
  model1.matrix[,i]=GrowthModel(nt0=model1.posterior$euc.uncal.nnorm$nt0[i],r=model1.posterior$euc.uncal.nnorm$r[i],timeRange=timeRange,b1=model1.posterior$euc.uncal.nnorm$b1[i],b2=model1.posterior$euc.uncal.nnorm$b2[i],x1=x1,x2=x2,a=1,raw = FALSE)$PrDens
  model2.matrix[,i]=GrowthModel(nt0=model2.posterior$euc.uncal.nnorm$nt0[i],r=model2.posterior$euc.uncal.nnorm$r[i],timeRange=timeRange,b1=model2.posterior$euc.uncal.nnorm$b1[i],b2=model2.posterior$euc.uncal.nnorm$b2[i],x1=x1,x2=x2,a=1,raw = FALSE)$PrDens
  model3.matrix[,i]=GrowthModel(nt0=model3.posterior$euc.uncal.nnorm$nt0[i],r=model3.posterior$euc.uncal.nnorm$r[i],timeRange=timeRange,b1=model3.posterior$euc.uncal.nnorm$b1[i],b2=model3.posterior$euc.uncal.nnorm$b2[i],x1=x1,x2=x2,a=1,raw = FALSE)$PrDens
  model4.matrix[,i]=GrowthModel(nt0=model4.posterior$euc.uncal.nnorm$nt0[i],r=model4.posterior$euc.uncal.nnorm$r[i],timeRange=timeRange,b1=model4.posterior$euc.uncal.nnorm$b1[i],b2=model4.posterior$euc.uncal.nnorm$b2[i],x1=x1,x2=x2,a=1,raw = FALSE)$PrDens
}

model1.sum= t(apply(model1.matrix,1,quantile,c(0.05,0.5,0.95)))
model2.sum= t(apply(model2.matrix,1,quantile,c(0.05,0.5,0.95)))
model3.sum= t(apply(model3.matrix,1,quantile,c(0.05,0.5,0.95)))
model4.sum= t(apply(model4.matrix,1,quantile,c(0.05,0.5,0.95)))

pdf(file=here('figures','main_paper','figure3_prdens.pdf'),width = 7, height = 7)
par(mfrow=c(2,2),mar=c(4,4,2,1))
plot(0,0,type='n',xlim=timeRange,ylim=c(0,max(c(model1.sum,model2.sum,model3.sum,model4.sum))),xlab='Cal BP',ylab='Probabilty Mass',main='Model 1')
polygon(c(timeRange[1]:timeRange[2],rev(timeRange[1]:timeRange[2])),c(model1.sum[,1],rev(model1.sum[,3])),col=adjustcolor( brewer.pal(4, 'Set2')[1], alpha.f = 0.2),border=NA)
lines(timeRange[1]:timeRange[2],model1.sum[,2],col=brewer.pal(4, 'Set2')[1])
abline(v=228,lty=3,lwd=1.5)

plot(0,0,type='n',xlim=timeRange,ylim=c(0,max(c(model1.sum,model2.sum,model3.sum,model4.sum))),xlab='Cal BP',ylab='Probabilty Mass',main='Model 2')
polygon(c(timeRange[1]:timeRange[2],rev(timeRange[1]:timeRange[2])),c(model2.sum[,1],rev(model2.sum[,3])),col=adjustcolor( brewer.pal(4, 'Set2')[2], alpha.f = 0.2),border=NA)
lines(timeRange[1]:timeRange[2],model1.sum[,2],col=brewer.pal(4, 'Set2')[2])
abline(v=228,lty=3,lwd=1.5)

plot(0,0,type='n',xlim=timeRange,ylim=c(0,max(c(model1.sum,model2.sum,model3.sum,model4.sum))),xlab='Cal BP',ylab='Probabilty Mass',main='Model 3')
polygon(c(timeRange[1]:timeRange[2],rev(timeRange[1]:timeRange[2])),c(model3.sum[,1],rev(model3.sum[,3])),col=adjustcolor( brewer.pal(4, 'Set2')[3], alpha.f = 0.2),border=NA)
lines(timeRange[1]:timeRange[2],model1.sum[,2],col=brewer.pal(4, 'Set2')[3])
abline(v=228,lty=3,lwd=1.5)

plot(0,0,type='n',xlim=timeRange,ylim=c(0,max(c(model1.sum,model2.sum,model3.sum,model4.sum))),xlab='Cal BP',ylab='Probabilty Mass',main='Model 4')
polygon(c(timeRange[1]:timeRange[2],rev(timeRange[1]:timeRange[2])),c(model4.sum[,1],rev(model4.sum[,3])),col=adjustcolor( brewer.pal(4, 'Set2')[4], alpha.f = 0.2),border=NA)
lines(timeRange[1]:timeRange[2],model1.sum[,2],col=brewer.pal(4, 'Set2')[4])
abline(v=228,lty=3,lwd=1.5)
dev.off()


### Figure 4 Posterior Predictive Check ####
# Load Posterior Predictive Check & Observed SPDs
load(here('R_imagefiles','post_pred_spd.RData'))
# Source R Functions
source(here('src','plotFunctions.R'))
source(here('src','randomThin.R'))
# Observed (thinned) SPD
set.seed(123)
index=randomThin(bins)
obs.norm =spd(cal.combined.norm[index],datenormalised=FALSE, spdnormalised = TRUE,timeRange=c(800,150),verbose=FALSE)
obs.nnorm =spd(cal.combined.nnorm[index],datenormalised=FALSE, spdnormalised = TRUE,timeRange=c(800,150),verbose=FALSE)

pdf(file = here('figures','main_paper','figure4.pdf'),width = 7,height = 7)
par(mfrow=c(2,2),mar=c(4,4,2,1))
plotPPCheckSPD(obs.nnorm,postpcheck.model1.spd[[4]],main='Model 1',obs.col = 'black',sim.fillcol = adjustcolor( brewer.pal(4, 'Set2')[1], alpha.f = 0.2),sim.col = brewer.pal(4, 'Set2')[1],sim.lty = 1,obs.lty = 2,obs.lwd = 1)
abline(v=228,lty=3,lwd=1.5)
plotPPCheckSPD(obs.nnorm,postpcheck.model2.spd[[4]],main='Model 2',obs.col = 'black',sim.fillcol = adjustcolor( brewer.pal(4, 'Set2')[2], alpha.f = 0.2),sim.col = brewer.pal(4, 'Set2')[2],sim.lty = 1,obs.lty = 2,obs.lwd = 1)
abline(v=228,lty=3,lwd=1.5)
plotPPCheckSPD(obs.nnorm,postpcheck.model3.spd[[4]],main='Model 3',obs.col = 'black',sim.fillcol = adjustcolor( brewer.pal(4, 'Set2')[3], alpha.f = 0.2),sim.col = brewer.pal(4, 'Set2')[3],sim.lty = 1,obs.lty = 2,obs.lwd = 1)
abline(v=228,lty=3,lwd=1.5)
plotPPCheckSPD(obs.nnorm,postpcheck.model4.spd[[4]],main='Model 4',obs.col = 'black',sim.fillcol = adjustcolor( brewer.pal(4, 'Set2')[4], alpha.f = 0.2),sim.col = brewer.pal(4, 'Set2')[4],sim.lty = 1,obs.lty = 2,obs.lwd = 1)
abline(v=228,lty=3,lwd=1.5)
dev.off()
