plotMarginalPosterior = function(x,hpdi=0.9,simple=FALSE,...)
{
  interval = HPDinterval(mcmc(x),prob = hpdi)
  dens=density(x)
  plot(dens,type='n',...)
  if (simple)
  {
    polygon(x=c(dens$x,rev(dens$x)),y=c(dens$y,rep(0,length(dens$y))),border=NA,col='#0077BC')
  } else {
  hpdi.x = dens$x[which(dens$x>=interval[1]&dens$x<=interval[2])]
  hpdi.y = dens$y[which(dens$x>=interval[1]&dens$x<=interval[2])]
  polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col='lightblue')
  polygon(x=c(dens$x,rev(dens$x)),y=c(dens$y,rep(0,length(dens$y))),border='darkgrey')
  abline(v=median(x),lty=2)
  legend('topright',legend=c(paste0(hpdi*100,'% HPDI(Lo): ',round(interval[1],5)),paste0(hpdi*100,'% HPDI(Hi): ',round(interval[2],5)),paste0('Median: ',round(median(x),5))),bty='n')
  }
}

plotPPCheckSPD = function(obs.spd,ppmat,...)
{
 plot(obs.spd$grid$calBP,obs.spd$grid$PrDens,xlim=rev(range(obs.spd$grid$calBP)),ylim=c(0,max(c(as.numeric(obs.spd$grid$PrDens),as.numeric(ppmat)))),type='n',xlab='Cal BP',ylab='SPD',...)
 lo = apply(ppmat,1,quantile,0.025)
 hi = apply(ppmat,1,quantile,0.975)
 median = apply(ppmat,1,median)
 polygon(c(obs.spd$grid$calBP,rev(obs.spd$grid$calBP)),c(lo,rev(hi)),border=NA,col='lightgrey')
 lines(obs.spd$grid$calBP,median,lwd=2,lty=2)
 lines(obs.spd$grid$calBP,obs.spd$grid$PrDens,lwd=2,col='red')
 legend('topleft',legend=c('Observed SPD','Median PPC','95% PPC Interval'),lwd=c(2,2,10),col=c('red','black','lightgrey'),lty=c(1,2,2))
 }

