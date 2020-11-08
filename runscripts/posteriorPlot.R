# Load Libraries
library(coda)
library(here)
library(ReIns)
library(latex2exp)
# Load Posteriors ####
load(here('R_imagefiles','ABC_posteriors.RData'))



# Model 4

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
