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