tolDiagPlot = function(x,param,tolcol,minsim=100)
{
  require(RColorBrewer)
  nsim = nrow(x)
  x = x[order(x[tolcol]),]
  p = x[param][,1]
  
  xx = nsim:minsim
  
  if (is.numeric(p))
  {
    lo = hi = med= numeric(length(xx))
    for (i in 1:length(xx))
    {
      lo[i]=quantile(p[1:xx[i]],0.025)
      hi[i]=quantile(p[1:xx[i]],0.975)
      med[i]=median(p[1:xx[i]])
    }
  plot(xx,med,lwd=2,xlim=rev(range(xx)),ylim=range(p),xlab='Number of simulations',ylab=param,type='l')
  lines(xx,lo,lty=2)
  lines(xx,hi,lty=2)
  }
  
  if (is.character(p))
  {
    n.cat = length(unique(p))
    cats = sort(unique(p))
    mat = matrix(0,nrow=n.cat,ncol=length(xx))
    for (i in 1:length(xx))
    {
      mat[,i] = sapply(cats,function(x,y){sum(x==y)},y=p[1:xx[i]])/xx[i]
    }
    plot(xx,mat[1,],lwd=2,xlim=rev(range(xx)),ylim=c(0,1),xlab='Number of simulations',ylab=param,type='n')
    for (i in 1:n.cat)
    {
      lines(xx,mat[i,],lwd=2,col=brewer.pal(n.cat, 'Set2')[i])
    }
    legend('topleft',legend=cats,lwd=2,col=brewer.pal(n.cat, 'Set2'))
  }
}