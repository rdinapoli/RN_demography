modelRunner = function(model,x.norm,x.nnorm,bins,customCurves,timeRange,simonly=FALSE)
{
  # check model
  if (any(is.na(model$PrDens)|is.nan(model$PrDens)|model$PrDens==Inf))
  {
    return(data.frame(euc.cal.norm=NA,euc.cal.nnorm=NA,euc.uncal.norm=NA,euc.uncal.nnorm=NA,nrmse.cal.norm=NA,nrmse.cal.nnorm=NA,nrmse.uncal.norm=NA,nrmse.uncal.nnorm=NA))
  }
  
  # generate Target SPD via random thinning
  index=randomThin(bins)
  if (!simonly)
  {
  target.spd.norm =spd(x.norm[index],datenormalised=FALSE, spdnormalised = TRUE,timeRange=timeRange,verbose=FALSE)
  target.spd.nnorm =spd(x.nnorm[index],datenormalised=FALSE, spdnormalised = TRUE,timeRange=timeRange,verbose=FALSE)
  }
  
  #Back calibrate model and sample radiocarbon dates
  curves = unique(x.norm$metadata$CalCurve)
  n.curves = length(curves)
  calsample.norm = uncalsample.norm =  calsample.nnorm = uncalsample.nnorm = vector('list',length=n.curves) #List of Uncalibrated CalGrids from the model and Sampled Errors
  idCounter = 0 
  for (i in 1:n.curves)
  {
    if (curves[[i]]%in%c("intcal13","intcal20","shcal13","shcal20","marine13","marine20","intcal13nhpine16","shcal13shkauri16"))
    {
      d = uncalibrate(model,verbose=F,calCurves = curves[[i]])
      n.samples = sum(x.norm[index]$metadata$CalCurve==curves[[i]])
      errors = sample(x.norm$metadata$Error[which(x.norm$metadata$CalCurve==curves[[i]])],size=n.samples,replace=TRUE)
      uncal.samples = sample(d$CRA,size=n.samples,prob=d$PrDens,replace=TRUE)
      cal.samples = sample(d$CRA,size=n.samples,prob=d$Raw,replace=TRUE)
      ids = idCounter + (1:n.samples)
      uncalsample.norm[[i]] = calibrate2(uncal.samples,errors,calCurves = curves[[i]],ids = ids,normalised=TRUE,verbose=FALSE)
      uncalsample.nnorm[[i]] = calibrate2(uncal.samples,errors,calCurves = curves[[i]],ids = ids,normalised=FALSE,verbose=FALSE)
      calsample.norm[[i]] = calibrate2(cal.samples,errors,calCurves = curves[[i]],ids = ids,normalised=TRUE,verbose=FALSE)
      calsample.nnorm[[i]] = calibrate2(cal.samples,errors,calCurves = curves[[i]],ids = ids,normalised=FALSE,verbose=FALSE)
    } else {
      d = uncalibrate(model,verbose=F,calCurves = customCurves[[curves[i]]])
      n.samples = sum(x.norm[index]$metadata$CalCurve==curves[i])
      errors = sample(x.norm$metadata$Error[which(x.norm$metadata$CalCurve==curves[[i]])],size=n.samples[i],replace=TRUE)
      uncal.samples = sample(d$CRA,size=n.samples,prob=d$PrDens,replace=TRUE)
      cal.samples = sample(d$CRA,size=n.samples,prob=d$Raw,replace=TRUE)
      ids = idCounter + (1:n.samples)
      uncalsample.norm[[i]] = calibrate2(uncal.samples,errors,calCurves =customCurves[[curves[i]]],ids = ids,normalised=TRUE,verbose=FALSE)
      uncalsample.nnorm[[i]] = calibrate2(uncal.samples,errors,calCurves = customCurves[[curves[i]]],ids = ids,normalised=FALSE,verbose=FALSE)
      calsample.norm[[i]] = calibrate2(cal.samples,errors,calCurves =customCurves[[curves[i]]],ids = ids,normalised=TRUE,verbose=FALSE)
      calsample.nnorm[[i]] = calibrate2(cal.samples,errors,calCurves = customCurves[[curves[i]]],ids = ids,normalised=FALSE,verbose=FALSE)
    }
    idCounter = max(ids)
  }
  
  #Create candidate SPDs
  uncalsample.norm <- do.call("combine", uncalsample.norm)
  uncalsample.nnorm <- do.call("combine", uncalsample.nnorm)
  calsample.norm <- do.call("combine", calsample.norm)
  calsample.nnorm <- do.call("combine", calsample.nnorm)
  
  spd.uncalsample.norm = spd(uncalsample.norm,timeRange=timeRange,spdnormalised = TRUE,verbose=FALSE)
  spd.uncalsample.nnorm = spd(uncalsample.nnorm,timeRange=timeRange,spdnormalised = TRUE,verbose=FALSE)
  spd.calsample.norm = spd(calsample.norm,timeRange=timeRange,spdnormalised = TRUE,verbose=FALSE)
  spd.calsample.nnorm = spd(calsample.nnorm,timeRange=timeRange,spdnormalised = TRUE,verbose=FALSE)
  
  if (simonly)
  {
    res = data.frame(calBP=spd.uncalsample.norm$grid$calBP,uncalsample.norm=spd.uncalsample.norm$grid$PrDens,uncalsample.nnorm=spd.uncalsample.nnorm$grid$PrDens,calsample.norm=spd.calsample.norm$grid$PrDens,calsample.nnorm=spd.calsample.nnorm$grid$PrDens)
    return(res)
  } else {
  
  #Compute Epsilon
  euc.cal.norm=sqrt(sum((target.spd.norm$grid$PrDens-spd.calsample.norm$grid$PrDens)^2))
  euc.cal.nnorm=sqrt(sum((target.spd.nnorm$grid$PrDens-spd.calsample.nnorm$grid$PrDens)^2))
  euc.uncal.norm=sqrt(sum((target.spd.norm$grid$PrDens-spd.uncalsample.norm$grid$PrDens)^2))
  euc.uncal.nnorm=sqrt(sum((target.spd.nnorm$grid$PrDens-spd.uncalsample.nnorm$grid$PrDens)^2))
  
  nrmse.cal.norm=sqrt(sum((target.spd.norm$grid$PrDens-spd.calsample.norm$grid$PrDens)^2)/length(model$CalBP))/sd(target.spd.norm$grid$PrDens)
  nrmse.cal.nnorm=sqrt(sum((target.spd.nnorm$grid$PrDens-spd.calsample.nnorm$grid$PrDens)^2)/length(model$CalBP))/sd(target.spd.norm$grid$PrDens)
  nrmse.uncal.norm=sqrt(sum((target.spd.norm$grid$PrDens-spd.uncalsample.norm$grid$PrDens)^2)/length(model$CalBP))/sd(target.spd.norm$grid$PrDens)
  nrmse.uncal.nnorm=sqrt(sum((target.spd.nnorm$grid$PrDens-spd.uncalsample.nnorm$grid$PrDens)^2)/length(model$CalBP))/sd(target.spd.norm$grid$PrDens)

  #return results
  return(data.frame(euc.cal.norm=euc.cal.norm,euc.cal.nnorm=euc.cal.nnorm,euc.uncal.norm=euc.uncal.norm,euc.uncal.nnorm=euc.uncal.nnorm,nrmse.cal.norm=nrmse.cal.norm,nrmse.cal.nnorm=nrmse.cal.nnorm,nrmse.uncal.norm=nrmse.uncal.norm,nrmse.uncal.nnorm=nrmse.uncal.nnorm))
  }
}

