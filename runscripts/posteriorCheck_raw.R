# Load Libraries
library(coda)
library(here)

# Source Scripts
source(here('src','growthModel.R'))

# Load Posteriors and Covariates ####
load(here('R_imagefiles','ABC_posteriors.RData'))
load(here('R_imagefiles/','variables.RData'))

# General Settings ####
timeRange = c(800,150)
# Extract Covariates
x1 = palm$PollenPerc
x2 = soi$SOIpr
nrep= 1000 # Plot top 100 Models


# Run raw models #####


# Model 1 ####
postpcheck.model1= vector('list',length=8)
postpcheck.model1[1:8] = rep(list(matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nrep)),8)

for (k in 1:8)
{
  post= model1.posterior[[k]][1:nrep,]
  for (i in 1:nrep)
  {
    postpcheck.model1[[k]][,i] = GrowthModel(nt0=post$nt0[i],r=post$r[i],timeRange=timeRange,b1=post$b1[i],b2=post$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
  }
}


# Plot results
titles = c('Calsample - Normalised','Calsample - Non-normalised','Uncalsample - Normalised','Uncalsample - Non-normalised')

jpeg(file=here('figures','temporary','posterior_predictive_check_model1_euc.jpeg'),width = 7, height = 7,units='in',res=300)
par(mfrow=c(2,2),mar=c(5,4,2,1))
for(i in 1:4)
{
  plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(postpcheck.model1[[i]]),xlab='Cal BP',ylab='N',main=titles[i])
  apply(postpcheck.model1[[i]],2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))
}
dev.off()

jpeg(file=here('figures','temporary','posterior_predictive_check_model1_nrmse.jpeg'),width = 7, height = 7,units='in',res=300)
par(mfrow=c(2,2),mar=c(5,4,2,1))
for(i in 5:8)
{
  plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(postpcheck.model1[[i]]),xlab='Cal BP',ylab='N',main=titles[i-4])
  apply(postpcheck.model1[[i]],2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))
}
dev.off()


# Model 2 ####
postpcheck.model2= vector('list',length=8)
postpcheck.model2[1:8] = rep(list(matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nrep)),8)

for (k in 1:8)
{
  post= model2.posterior[[k]][1:nrep,]
  for (i in 1:nrep)
  {
    postpcheck.model2[[k]][,i] = GrowthModel(nt0=post$nt0[i],r=post$r[i],timeRange=timeRange,b1=post$b1[i],b2=post$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
  }
}


# Plot results
titles = c('Calsample - Normalised','Calsample - Non-normalised','Uncalsample - Normalised','Uncalsample - Non-normalised')

jpeg(file=here('figures','temporary','posterior_predictive_check_model2_euc.jpeg'),width = 7, height = 7,units='in',res=300)
par(mfrow=c(2,2),mar=c(5,4,2,1))
for(i in 1:4)
{
  plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(postpcheck.model2[[i]]),xlab='Cal BP',ylab='N',main=titles[i])
  apply(postpcheck.model2[[i]],2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))
}
dev.off()

jpeg(file=here('figures','temporary','posterior_predictive_check_model2_nrmse.jpeg'),width = 7, height = 7,units='in',res=300)
par(mfrow=c(2,2),mar=c(5,4,2,1))
for(i in 5:8)
{
  plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(postpcheck.model2[[i]]),xlab='Cal BP',ylab='N',main=titles[i-4])
  apply(postpcheck.model2[[i]],2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))
}
dev.off()






# Model 4 ####
postpcheck.model4 = vector('list',length=8)
postpcheck.model4[1:8] = rep(list(matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nrep)),8)

for (k in 1:8)
{
  post= model4.posterior[[k]][1:nrep,]
  for (i in 1:nrep)
  {
    postpcheck.model4[[k]][,i] = GrowthModel(nt0=post$nt0[i],r=post$r[i],timeRange=timeRange,b1=post$b1[i],b2=post$b2[i],x1=x1,x2=x2,a=1,raw = TRUE)$N
  }
}


# Plot results
titles = c('Calsample - Normalised','Calsample - Non-normalised','Uncalsample - Normalised','Uncalsample - Non-normalised')

jpeg(file=here('figures','temporary','posterior_predictive_check_model4_euc.jpeg'),width = 7, height = 7,units='in',res=300)
par(mfrow=c(2,2),mar=c(5,4,2,1))
for(i in 1:4)
{
  plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(postpcheck.model4[[i]]),xlab='Cal BP',ylab='N',main=titles[i])
  apply(postpcheck.model4[[i]],2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))
}
dev.off()

jpeg(file=here('figures','temporary','posterior_predictive_check_model4_nrmse.jpeg'),width = 7, height = 7,units='in',res=300)
par(mfrow=c(2,2),mar=c(5,4,2,1))
for(i in 5:8)
{
  plot(timeRange[1]:timeRange[2],xlim=timeRange,ylim=range(postpcheck.model4[[i]]),xlab='Cal BP',ylab='N',main=titles[i-4])
  apply(postpcheck.model4[[i]],2,lines,x=timeRange[1]:timeRange[2],col=rgb(0,0,0,0.05))
}
dev.off()

