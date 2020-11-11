# Load Libraries
library(coda)
library(here)
library(MASS)
library(ReIns)
library(latex2exp)
# Load Posteriors ####
load(here('R_imagefiles','ABC_posteriors.RData'))
source(here('src','plotFunctions.R'))

# Model 1 ####
for (i in 1:8)
{
  name = names(model1.posterior)[i]
  fn = paste0('m1_post_',name,'.pdf')
  post = model1.posterior[[i]]
  pdf(file=here('figures','temporary',fn),width = 7, height = 3.5)
  par(mfrow=c(1,2))
  plotMarginalPosterior(post$nt0,xlab=TeX('$N_0$'),main=TeX('$N_0$'))
  lines(density(post$nt0)$x,dtexp(density(post$nt0)$x,rate=10,endpoint=0.5),lty=3,col='darkgrey',lwd=2)
  plotMarginalPosterior(post$r,xlab=TeX('$r$'),main=TeX('$r$'))
  lines(density(post$r)$x,dtexp(density(post$r)$x,rate=20,endpoint=0.5),lty=3,col='darkgrey',lwd=2)
  dev.off()
}

for (i in 1:8)
{
  name = names(model1.posterior)[i]
  fn = paste0('m1_jointpost_',name,'.jpeg')
  post = model1.posterior[[i]]
  jpeg(file=here('figures','temporary',fn),width = 4, height = 4,res=300,units = 'in',pointsize = 12)
  par(mfrow=c(2,2),mar=c(4,5,1,0))
  plotMarginalPosterior(post$nt0,simple=TRUE,xlab=TeX('$N_0$'),main='')
  replicate(1,plot(0,0,axes=F,xlab='',ylab='',type='n'))
  image(kde2d(post$nt0,post$r,n=500),xlab=TeX('$N_0$'),ylab=TeX('$r$'),col=hcl.colors(12, "Blues 3", rev = TRUE))
  plotMarginalPosterior(post$r,simple=TRUE,xlab=TeX('$r$'),main='')
  dev.off()
}



# Model 4 ####

for (i in 1:8)
{
name = names(model4.posterior)[i]
fn = paste0('m4_post_',name,'.pdf')
post = model4.posterior[[i]]
pdf(file=here('figures','temporary',fn),width = 7, height = 7)
par(mfrow=c(2,2))
plotMarginalPosterior(post$nt0,xlab=TeX('$N_0$'),main=TeX('$N_0$'))
lines(density(post$nt0)$x,dtexp(density(post$nt0)$x,rate=10,endpoint=0.5),lty=3,col='darkgrey',lwd=2)
plotMarginalPosterior(post$r,xlab=TeX('$r$'),main=TeX('$r$'))
lines(density(post$r)$x,dtexp(density(post$r)$x,rate=20,endpoint=0.5),lty=3,col='darkgrey',lwd=2)
plotMarginalPosterior(post$b1,xlab=TeX('$\\beta_{palm}$'),main=TeX('$\\beta_{palm}$'))
lines(density(post$b1)$x,dnorm(density(post$b1)$x,mean=0,sd=0.01),lty=3,col='darkgrey',lwd=2)
plotMarginalPosterior(post$b2,xlab=TeX('$\\beta_{soi}$'),main=TeX('$\\beta_{soi}$'))
lines(density(post$b2)$x,dnorm(density(post$b2)$x,mean=0,sd=0.2),lty=3,col='darkgrey',lwd=2)
dev.off()
}


for (i in 1:8)
{
  name = names(model4.posterior)[i]
  fn = paste0('m4_jointpost_',name,'.jpeg')
  post = model4.posterior[[i]]
  jpeg(file=here('figures','temporary',fn),width = 8, height = 8,res=300,units = 'in',pointsize = 12)
  par(mfrow=c(4,4),mar=c(4,5,1,0))
  plotMarginalPosterior(post$nt0,simple=TRUE,xlab=TeX('$N_0$'),main='')
  replicate(3,plot(0,0,axes=F,xlab='',ylab='',type='n'))
  image(kde2d(post$nt0,post$r,n=500),xlab=TeX('$N_0$'),ylab=TeX('$r$'),col=hcl.colors(12, "Blues 3", rev = TRUE))
  plotMarginalPosterior(post$r,simple=TRUE,xlab=TeX('$r$'),main='')
  replicate(2,plot(0,0,axes=F,xlab='',ylab='',type='n'))
  image(kde2d(post$nt0,post$b1,n=500),xlab=TeX('$N_0$'),ylab=TeX('$\\beta_{palm}$'),col=hcl.colors(12, "Blues 3", rev = TRUE))
  image(kde2d(post$r,post$b1,n=500),xlab=TeX('$r$'),ylab=TeX('$\\beta_{palm}$'),col=hcl.colors(12, "Blues 3", rev = TRUE))
  plotMarginalPosterior(post$b1,simple=TRUE,xlab=TeX('$\\beta_{palm}$'),main='')
  replicate(1,plot(0,0,axes=F,xlab='',ylab='',type='n'))
  image(kde2d(post$nt0,post$b2,n=500),xlab=TeX('$N_0$'),ylab=TeX('$\\beta_{soi}$'),col=hcl.colors(12, "Blues 3", rev = TRUE))
  image(kde2d(post$r,post$b2,n=500),xlab=TeX('$r$'),ylab=TeX('$\\beta_{soi}$'),col=hcl.colors(12, "Blues 3", rev = TRUE))
  image(kde2d(post$b1,post$b2,n=500),xlab=TeX('$\\beta_{palm}$'),ylab=TeX('$\\beta_{soi}$'),col=hcl.colors(12, "Blues 3", rev = TRUE))
  plotMarginalPosterior(post$b2,simple=TRUE,xlab=TeX('$\\beta_{soi}$'),main='')
  dev.off()
}



