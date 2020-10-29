plotPosterior = function(x,hpdi=0.9,...)
{
  interval = HPDinterval(mcmc(x),prob = hpdi)
  dens=density(x)
  plot(dens,type='n',...)
  hpdi.x = dens$x[which(dens$x>=interval[1]&dens$x<=interval[2])]
  hpdi.y = dens$y[which(dens$x>=interval[1]&dens$x<=interval[2])]
  polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col='lightblue')
  polygon(x=c(dens$x,rev(dens$x)),y=c(dens$y,rep(0,length(dens$y))),border='darkgrey')
  abline(v=median(x),lty=2)
  legend('topright',legend=c(paste0(hpdi*100,'% HPDI(Lo): ',round(interval[1],5)),paste0(hpdi*100,'% HPDI(Hi): ',round(interval[2],5)),paste0('Median: ',round(median(x),5))),bty='n')
}

postPredictiveCheck = function(model=logisticModel,params,type,timeRange=c(800,150),spaghettiSize=100,alpha=0.2,...)
{
  nsim=nrow(x)
  ppMat = matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nsim)
  
  for (i in 1:nsim)
  {
    args = as.list(params[i,])
    args[[length(args)+1]] = timeRange
    names(args)[length(args)]='timeRange'
    ppMat[,i]=do.call(model,args)$PrDens
  }
  plot(timeRange[1]:timeRange[2],xlim=timeRange,xlab='Cal BP',ylab='Probability',ylim=c(0,max(ppMat)),...)
  if (type=='envelope')
  {
  lo = apply(ppMat,1,quantile,0.025)
  hi = apply(ppMat,1,quantile,0.975)
  polygon(c(c(timeRange[1]:timeRange[2]),c(timeRange[2]:timeRange[1])),c(lo,rev(hi)),col='lightblue',border=NA)
  lines(timeRange[1]:timeRange[2],apply(ppMat,1,mean),lwd=2,lty=2)
  }
  
  if (type=='spaghetti')
  {
    apply(ppMat[,sample(nsim,size=spaghettiSize)],2,lines,x=c(timeRange[1]:timeRange[2]),col=rgb(0,0,0,alpha))
  }
}