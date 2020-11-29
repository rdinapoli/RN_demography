GrowthModel = function(nt0,r,x1,x2,a,b1,b2,timeRange,raw=FALSE)
{
  CalBP = timeRange[1]:timeRange[2]
  tt = 0:length(CalBP)
  pop = numeric(length=length(tt))
  pop[1] = nt0
  for (i in 2:length(tt))
  {
    pop[i] = pop[i-1] * exp(r*(1-pop[i-1]/exp(a+b1*x1[i-1]+b2*x2[i-1])))
  }
  pop = pop[-1]
  if (raw)
  {
    d = data.frame(CalBP=CalBP,N=pop)
  } else
  {
    PrDens = pop/sum(pop)
    d = data.frame(CalBP=CalBP,PrDens=PrDens)
    class(d)='CalGrid'
  }
  return(d)
}


