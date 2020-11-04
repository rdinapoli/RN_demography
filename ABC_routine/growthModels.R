logisticModel = function(k,r,timeRange)
{
  # k ... proportion of carrying capacity population at start date
  # r ... growth rate
  CalBP = timeRange[1]:timeRange[2]
  tt = 1:length(CalBP)
  pop = 1/(1+((1-k)/k)*exp(-r*tt))
  PrDens = pop/sum(pop)
  d = data.frame(CalBP=CalBP,PrDens=PrDens)
  class(d)='CalGrid'
  return(d)
}


logisticGrowthModel = function(nt0=0.01,r,K=1,timeRange,raw=FALSE)
{
  CalBP = timeRange[1]:timeRange[2]
  tt = 0:length(CalBP)
  pop = numeric(length=length(tt))
  pop[1] = nt0
  for (i in 2:length(tt))
  {
    pop[i] = pop[i-1] * exp(r*(1-pop[i-1]/K))
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


SingleCovariateLogisticGrowthModel = function(nt0,r,x,a,b,timeRange,raw=FALSE)
{
  CalBP = timeRange[1]:timeRange[2]
  tt = 0:length(CalBP)
  pop = numeric(length=length(tt))
  pop[1] = nt0
  for (i in 2:length(tt))
  {
    pop[i] = pop[i-1] * exp(r*(1-pop[i-1]/(a+b*x[i-1])))
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

DoubleCovariateLogisticGrowthModel = function(nt0,r,x1,x2,a,b1,b2,timeRange,raw=FALSE)
{
  CalBP = timeRange[1]:timeRange[2]
  tt = 0:length(CalBP)
  pop = numeric(length=length(tt))
  pop[1] = nt0
  for (i in 2:length(tt))
  {
    pop[i] = pop[i-1] * exp(r*(1-pop[i-1]/(a+b1*x1[i-1]+b2*x2[i-1])))
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


# Sanity Checks
# m1a=logisticGrowthModel(nt0=0.01,r=0.04,K=1,timeRange=c(800,150))
# s = 3456
# m1b=logisticGrowthModel(nt0=0.01*s,r=0.04,K=1*s,timeRange=c(800,150))
# all.equal(m1a$PrDens,m1b$PrDens)

# onefnoisegen<-function(omega=-1,np=1000,m=100,M=200){
#   mp=1
#   X<-rnorm(np,0,1)
#   np=length(X)
#   Y=fft(X)
#   tau=(1:(np/2))/(np/2)
#   tau=c(tau,sort(tau,decreasing=T))
#   redY=Y*(tau^omega)/(tau[mp]^omega)
#   Xi=as.double(fft(redY,inverse=TRUE))
#   Xi=((Xi-min(Xi))/abs(min(Xi)-max(Xi)))*abs(M-m)+m
#   return(Xi)
# }
# 
# set.seed(12345)
# x = onefnoisegen(omega=-2,np=length(800:150),m=-3,M=3)
# m1a=SingleCovariateLogisticGrowthModel(nt0=0.01,r=0.04,x=x,a=1,b=0.3,timeRange=c(800,150))
# m1b=SingleCovariateLogisticGrowthModel(nt0=0.01*s,r=0.04,x=x,a=1*s,b=0.3*s,timeRange=c(800,150))
# m1c=SingleCovariateLogisticGrowthModel(nt0=0.01*s,r=0.04,x=x*s,a=1*s,b=0.3,timeRange=c(800,150))
# all.equal(m1a$PrDens,m1b$PrDens)
# all.equal(m1a$PrDens,m1c$PrDens)
# 
# 
# set.seed(6789)
# x1 = onefnoisegen(omega=-2,np=length(800:150),m=-3,M=3)
# x2 = onefnoisegen(omega=-0.5,np=length(800:150),m=-10,M=10)
# m1a=DoubleCovariateLogisticGrowthModel(nt0=0.01,r=0.04,x1=x1,x2=x2,a=1,b1=0.3,b2=-0.01,timeRange=c(800,150))
# m1b=DoubleCovariateLogisticGrowthModel(nt0=0.01*s,r=0.04,x1=x1,x2=x2,a=1*s,b1=0.3*s,b2=-0.01*s,timeRange=c(800,150))
# m1c=DoubleCovariateLogisticGrowthModel(nt0=0.01*s,r=0.04,x1=x1*s,x2=x2*s,a=1*s,b1=0.3,b2=-0.01,timeRange=c(800,150))
# all.equal(m1a$PrDens,m1b$PrDens)
# all.equal(m1a$PrDens,m1c$PrDens)


