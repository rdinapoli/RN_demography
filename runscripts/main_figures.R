library(here)
library(ggplot2)
library(latex2exp)
library(tidyverse)
library(RColorBrewer)
library(rcarbon)
### Figure 1 Correlation Plots ####


### Figure 2 Marginal Posterior Distributions of the Four Models ####
load(here('R_imagefiles','ABC_posteriors.RData'))
source(here('src','plotFunctions.R'))

# Using uncalsample, non-normalised ggplot way
posterior = rbind.data.frame(model1.posterior$euc.uncal.nnorm,model2.posterior$euc.uncal.nnorm,model3.posterior$euc.uncal.nnorm,model4.posterior$euc.uncal.nnorm)
posterior$model = c(rep('model1',nrow(posterior)/4),rep('model2',nrow(posterior)/4),rep('model3',nrow(posterior)/4),rep('model4',nrow(posterior)/4))
posterior$b1[which(posterior$model%in%c('model1','model3'))]=NA 
posterior$b2[which(posterior$model%in%c('model1','model2'))]=NA                    
posterior=posterior %>% select(nt0,r,b1,b2,model) %>%gather(key='parameter',value='value',-model)

param_names <- list(
  'nt0'=TeX('$N_{t=0}$'),
  'r'=TeX('$r$'),
  'b1'=TeX('$\\beta_{palm}$'),
  'b2'=TeX('$\\beta_{soi}$')
)

param_labeller <- function(variable,value){
  return(param_names[value])
}
dummy2 <- data.frame(X = c("b1", "b2"), Z = c(0, 0))


pdf(file=here('figures','main_paper','figure2_ggplot.pdf'),width = 7.2, height = 6.7)
ggplot(posterior, aes(value,color=model)) +
  geom_density(alpha = 0.3) + facet_wrap(~ parameter, scales = "free",labeller = param_labeller) + xlab('Posterior Estimates') + scale_color_brewer(palette="Set2") +theme(strip.text.x = element_text(size = 10))
dev.off()




# Using uncalsample, non-normalised base-r way
models = unique(posterior$model)

pdf(file=here('figures','main_paper','figure2_baser.pdf'),width = 7, height = 7)
par(mfrow=c(2,2),mar=c(4,4,2,1))
# Nt0
dens1 = density(posterior$value[which(posterior$model=='model1'&posterior$parameter=='nt0')])
dens2 = density(posterior$value[which(posterior$model=='model2'&posterior$parameter=='nt0')],bw=dens1$bw)
dens3 = density(posterior$value[which(posterior$model=='model3'&posterior$parameter=='nt0')],bw=dens1$bw)
dens4 = density(posterior$value[which(posterior$model=='model4'&posterior$parameter=='nt0')],bw=dens1$bw)
plot(dens1,xlab=TeX('$N_{t=0}$'),ylab='Density',main='a',lwd=2,col=brewer.pal(4, 'Set2')[1],ylim=range(c(dens1$y,dens2$y,dens3$y,dens4$y)))
lines(dens2,lwd=2,col=brewer.pal(4, 'Set2')[2])
lines(dens3,lwd=2,col=brewer.pal(4, 'Set2')[3])
lines(dens4,lwd=2,col=brewer.pal(4, 'Set2')[4])
legend('topright',legend=c('Model 1','Model 2','Model 3','Model 4'),lwd=2,col=brewer.pal(4, 'Set2'),bty='n')
# r
dens1 = density(posterior$value[which(posterior$model=='model1'&posterior$parameter=='r')])
dens2 = density(posterior$value[which(posterior$model=='model2'&posterior$parameter=='r')],bw=dens1$bw)
dens3 = density(posterior$value[which(posterior$model=='model3'&posterior$parameter=='r')],bw=dens1$bw)
dens4 = density(posterior$value[which(posterior$model=='model4'&posterior$parameter=='r')],bw=dens1$bw)
plot(dens1,xlab=TeX('$r$'),ylab='Density',main='b',lwd=2,col=brewer.pal(4, 'Set2')[1],ylim=range(c(dens1$y,dens2$y,dens3$y,dens4$y)))
lines(dens2,lwd=2,col=brewer.pal(4, 'Set2')[2])
lines(dens3,lwd=2,col=brewer.pal(4, 'Set2')[3])
lines(dens4,lwd=2,col=brewer.pal(4, 'Set2')[4])

# b1
dens2 = density(posterior$value[which(posterior$model=='model2'&posterior$parameter=='b1')])
dens4 = density(posterior$value[which(posterior$model=='model4'&posterior$parameter=='b1')],bw=dens2$bw)
plot(dens2,xlab=TeX('$\\beta_{palm}$'),ylab='Density',main='c',lwd=2,col=brewer.pal(4, 'Set2')[2],ylim=range(c(dens2$y,dens4$y)))
lines(dens4,lwd=2,col=brewer.pal(4, 'Set2')[4])
abline(v=0,lty=3)

# b2
dens3 = density(posterior$value[which(posterior$model=='model3'&posterior$parameter=='b2')])
dens4 = density(posterior$value[which(posterior$model=='model4'&posterior$parameter=='b2')],bw=dens3$bw)
plot(dens3,xlab=TeX('$\\beta_{soi}$'),ylab='Density',main='d',lwd=2,col=brewer.pal(4, 'Set2')[3],ylim=range(c(dens3$y,dens4$y)))
lines(dens4,lwd=2,col=brewer.pal(4, 'Set2')[4])
abline(v=0,lty=3)
dev.off()

### Figure 3 Fitted Model ####
load(here('R_imagefiles/','variables.RData'))
source(here('src','growthModel.R'))

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

plot(0,0,type='n',xlim=timeRange,ylim=c(0,max(c(model1.sum,model2.sum,model3.sum,model4.sum))),xlab='Cal BP',ylab='Probabilty Mass',main='Model 2')
polygon(c(timeRange[1]:timeRange[2],rev(timeRange[1]:timeRange[2])),c(model2.sum[,1],rev(model2.sum[,3])),col=adjustcolor( brewer.pal(4, 'Set2')[2], alpha.f = 0.2),border=NA)
lines(timeRange[1]:timeRange[2],model1.sum[,2],col=brewer.pal(4, 'Set2')[2])

plot(0,0,type='n',xlim=timeRange,ylim=c(0,max(c(model1.sum,model2.sum,model3.sum,model4.sum))),xlab='Cal BP',ylab='Probabilty Mass',main='Model 3')
polygon(c(timeRange[1]:timeRange[2],rev(timeRange[1]:timeRange[2])),c(model3.sum[,1],rev(model3.sum[,3])),col=adjustcolor( brewer.pal(4, 'Set2')[3], alpha.f = 0.2),border=NA)
lines(timeRange[1]:timeRange[2],model1.sum[,2],col=brewer.pal(4, 'Set2')[3])

plot(0,0,type='n',xlim=timeRange,ylim=c(0,max(c(model1.sum,model2.sum,model3.sum,model4.sum))),xlab='Cal BP',ylab='Probabilty Mass',main='Model 4')
polygon(c(timeRange[1]:timeRange[2],rev(timeRange[1]:timeRange[2])),c(model4.sum[,1],rev(model4.sum[,3])),col=adjustcolor( brewer.pal(4, 'Set2')[4], alpha.f = 0.2),border=NA)
lines(timeRange[1]:timeRange[2],model1.sum[,2],col=brewer.pal(4, 'Set2')[4])
dev.off()


### N Version 
# Prepare Matrices
model1.matrix = model2.matrix = model3.matrix = model4.matrix = matrix(NA,ncol=nrow(model1.posterior$euc.uncal.nnorm),nrow=length(timeRange[1]:timeRange[2]))

for (i in 1:nrow(model1.posterior$euc.uncal.nnorm))
{
  model1.matrix[,i]=GrowthModel(nt0=model1.posterior$euc.uncal.nnorm$nt0[i],r=model1.posterior$euc.uncal.nnorm$r[i],timeRange=timeRange,b1=model1.posterior$euc.uncal.nnorm$b1[i],b2=model1.posterior$euc.uncal.nnorm$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
  model2.matrix[,i]=GrowthModel(nt0=model2.posterior$euc.uncal.nnorm$nt0[i],r=model2.posterior$euc.uncal.nnorm$r[i],timeRange=timeRange,b1=model2.posterior$euc.uncal.nnorm$b1[i],b2=model2.posterior$euc.uncal.nnorm$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
  model3.matrix[,i]=GrowthModel(nt0=model3.posterior$euc.uncal.nnorm$nt0[i],r=model3.posterior$euc.uncal.nnorm$r[i],timeRange=timeRange,b1=model3.posterior$euc.uncal.nnorm$b1[i],b2=model3.posterior$euc.uncal.nnorm$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
  model4.matrix[,i]=GrowthModel(nt0=model4.posterior$euc.uncal.nnorm$nt0[i],r=model4.posterior$euc.uncal.nnorm$r[i],timeRange=timeRange,b1=model4.posterior$euc.uncal.nnorm$b1[i],b2=model4.posterior$euc.uncal.nnorm$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
}

model1.sum= t(apply(model1.matrix,1,quantile,c(0.05,0.5,0.95)))
model2.sum= t(apply(model2.matrix,1,quantile,c(0.05,0.5,0.95)))
model3.sum= t(apply(model3.matrix,1,quantile,c(0.05,0.5,0.95)))
model4.sum= t(apply(model4.matrix,1,quantile,c(0.05,0.5,0.95)))

pdf(file=here('figures','main_paper','figure3_N.pdf'),width = 7, height = 7)
par(mfrow=c(2,2),mar=c(4,4,2,1))
plot(0,0,type='n',xlim=timeRange,ylim=c(0,max(c(model1.sum,model2.sum,model3.sum,model4.sum))),xlab='Cal BP',ylab='N',main='Model 1')
polygon(c(timeRange[1]:timeRange[2],rev(timeRange[1]:timeRange[2])),c(model1.sum[,1],rev(model1.sum[,3])),col=adjustcolor( brewer.pal(4, 'Set2')[1], alpha.f = 0.2),border=NA)
lines(timeRange[1]:timeRange[2],model1.sum[,2],col=brewer.pal(4, 'Set2')[1])

plot(0,0,type='n',xlim=timeRange,ylim=c(0,max(c(model1.sum,model2.sum,model3.sum,model4.sum))),xlab='Cal BP',ylab='N',main='Model 2')
polygon(c(timeRange[1]:timeRange[2],rev(timeRange[1]:timeRange[2])),c(model2.sum[,1],rev(model2.sum[,3])),col=adjustcolor( brewer.pal(4, 'Set2')[2], alpha.f = 0.2),border=NA)
lines(timeRange[1]:timeRange[2],model1.sum[,2],col=brewer.pal(4, 'Set2')[2])

plot(0,0,type='n',xlim=timeRange,ylim=c(0,max(c(model1.sum,model2.sum,model3.sum,model4.sum))),xlab='Cal BP',ylab='N',main='Model 3')
polygon(c(timeRange[1]:timeRange[2],rev(timeRange[1]:timeRange[2])),c(model3.sum[,1],rev(model3.sum[,3])),col=adjustcolor( brewer.pal(4, 'Set2')[3], alpha.f = 0.2),border=NA)
lines(timeRange[1]:timeRange[2],model1.sum[,2],col=brewer.pal(4, 'Set2')[3])

plot(0,0,type='n',xlim=timeRange,ylim=c(0,max(c(model1.sum,model2.sum,model3.sum,model4.sum))),xlab='Cal BP',ylab='N',main='Model 4')
polygon(c(timeRange[1]:timeRange[2],rev(timeRange[1]:timeRange[2])),c(model4.sum[,1],rev(model4.sum[,3])),col=adjustcolor( brewer.pal(4, 'Set2')[4], alpha.f = 0.2),border=NA)
lines(timeRange[1]:timeRange[2],model1.sum[,2],col=brewer.pal(4, 'Set2')[4])
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
plotPPCheckSPD(obs.nnorm,postpcheck.model2.spd[[4]],main='Model 2',obs.col = 'black',sim.fillcol = adjustcolor( brewer.pal(4, 'Set2')[2], alpha.f = 0.2),sim.col = brewer.pal(4, 'Set2')[2],sim.lty = 1,obs.lty = 2,obs.lwd = 1)
plotPPCheckSPD(obs.nnorm,postpcheck.model3.spd[[4]],main='Model 3',obs.col = 'black',sim.fillcol = adjustcolor( brewer.pal(4, 'Set2')[3], alpha.f = 0.2),sim.col = brewer.pal(4, 'Set2')[3],sim.lty = 1,obs.lty = 2,obs.lwd = 1)
plotPPCheckSPD(obs.nnorm,postpcheck.model4.spd[[4]],main='Model 4',obs.col = 'black',sim.fillcol = adjustcolor( brewer.pal(4, 'Set2')[4], alpha.f = 0.2),sim.col = brewer.pal(4, 'Set2')[4],sim.lty = 1,obs.lty = 2,obs.lwd = 1)
dev.off()








