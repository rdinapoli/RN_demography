calibrate2 <- function(x, errors, ids=NA, dateDetails=NA, calCurves='intcal20', resOffsets=0 , resErrors=0, timeRange=c(55000,0), normalised=TRUE, F14C=FALSE, calMatrix=FALSE, eps=1e-5, ncores=1, verbose=TRUE, ...){
  
  if (ncores>1&!requireNamespace("doSNOW", quietly=TRUE)){	
    warning("the doSNOW package is required for multi-core processing; ncores has been set to 1")
    ncores=1
  }	
  
  # age and error checks
  if (length(x) != length(errors)){
    stop("Ages and errors (and ids/date details/offsets if provided) must be the same length.")
  }
  if (!is.na(ids[1]) & (length(x) != length(ids))){
    stop("Ages and errors (and ids/details/offsets if provided) must be the same length.")
  }
  if (any(is.na(x))|any(is.na(errors))){
    stop("Ages or errors contain NAs")
  }
  if (is.null(calCurves)|anyNA(calCurves)){
    stop("calCurves is NULL or contain NAs")
  }
  
  if (!any(class(calCurves) %in% c("matrix","data.frame"))& any(calCurves %in% c("intcal13","shcal13","marine13","intcal13nhpine16","shcal13shkauri16")) & (max(timeRange) > 50000 | min(timeRange)<0))
  {
    timeRange=c(50000,0)
    warning("timeRange value not supported with the selected curve(s); calibrating using timeRange=c(50000,0)")
  }
  
  if (any(class(calCurves) %in% c("matrix","data.frame")))
  {
    if (max(calCurves[,1])<timeRange[1] | min(calCurves[,1]) < timeRange[2] )
    {
      timeRange=c(max(calCurves[,1]),min(calCurves[,1]))
      warning(paste0("timeRange value not supported with the selected curve(s); calibrating using timeRange=c(",timeRange[1],",",timeRange[1],")"))
    }
  }
  
  if (F14C==TRUE&normalised==FALSE)
  {
    normalised=TRUE
    warning("normalised cannot be FALSE when F14C is set to TRUE, calibrating with normalised=TRUE")
  }
  # calCurve checks and set-up
  if (any(class(calCurves) %in% c("matrix","data.frame"))){
    cctmp <- as.matrix(calCurves)
    if (ncol(cctmp)!=3 | !all(sapply(cctmp,is.numeric))){
     # stop("The custom calibration curve must have just three numeric columns.")
    } else {
      colnames(cctmp) <- c("CALBP","C14BP","Error")
      if (max(cctmp[,2]) < max(x) | min(cctmp[,2]) > min(x)){
       # stop("The custom calibration curve does not cover the input age range.")
      }
      #cclist <- vector(mode="list", length=1)
      #cclist[[1]] <- cctmp
      #names(cclist) <- "custom"
      calCurves <- rep("custom",length(x))
      cclist2 <- vector(mode="list", length=1)
      names(cclist2) <- "custom"
      calBPrange = seq(max(cctmp[,1]),min(cctmp[,1]),-1)
      if (F14C)
      {
        F14 <- exp(cctmp[,2]/-8033) 
        F14Error <-  F14*cctmp[,3]/8033 
        calf14 <- approx(cctmp[,1], F14, xout=calBPrange)$y 
        calf14error <-  approx(cctmp[,1], F14Error, xout=calBPrange)$y 
        cclist2[[1]] <- list(calf14=calf14,calf14error=calf14error,calBPrange=calBPrange,calBPindex=which(calBPrange<=timeRange[1]&calBPrange>=timeRange[2]))
      } else {
        cclist2[[1]] = list(mu=stats::approx(cctmp[,1], cctmp[,2], xout = calBPrange)$y,tau2 = stats::approx(cctmp[,1], cctmp[,3], xout = calBPrange)$y^2,calBPrange=calBPrange,calBPindex=which(calBPrange<=timeRange[1]&calBPrange>=timeRange[2]))
      } 
    } 
  } else if (!all(calCurves %in% c("intcal13","intcal20","shcal13","shcal20","marine13","marine20","intcal13nhpine16","shcal13shkauri16","normal"))){
    stop("calCurves must be a character vector specifying one or more known curves or a custom three-column matrix/data.frame (see ?calibrate.default).")
  } else {
    tmp <- unique(calCurves)
    if (length(calCurves)==1){ calCurves <- rep(calCurves,length(x)) }
    #cclist <- vector(mode="list", length=length(tmp))
    #names(cclist) <- tmp
    cclist2 <- vector(mode="list", length=length(tmp))
    names(cclist2) <- tmp
    for (a in 1:length(tmp)){
      cctmp <- read_cal_curve_from_file(tmp[a])
      #cclist[[tmp[a]]] <- cctmp
      calBPrange = seq(max(cctmp[,1]),min(cctmp[,1]),-1)
      if (F14C)
      {
        F14 <- exp(cctmp[,2]/-8033) 
        F14Error <-  F14*cctmp[,3]/8033 
        calf14 <- approx(cctmp[,1], F14, xout=calBPrange)$y 
        calf14error <-  approx(cctmp[,1], F14Error, xout=calBPrange)$y 
        cclist2[[tmp[a]]] <- list(calf14=calf14,calf14error=calf14error,calBPrange=calBPrange,calBPindex=which(calBPrange<=timeRange[1]&calBPrange>=timeRange[2]))
      } else {
        cclist2[[tmp[a]]] = list(mu=stats::approx(cctmp[,1], cctmp[,2], xout = calBPrange)$y,tau2 = stats::approx(cctmp[,1], cctmp[,3], xout = calBPrange)$y^2,calBPrange=calBPrange,calBPindex=which(calBPrange<=timeRange[1]&calBPrange>=timeRange[2]))}
    }
  }
  ## container and reporting set-up
  reslist <- vector(mode="list", length=2)
  sublist <- vector(mode="list", length=length(x))
  if (calMatrix){
    calmBP <- seq(timeRange[1],timeRange[2],-1)
    calmat <- matrix(ncol=length(x), nrow=length(calmBP))
    rownames(calmat) <- calmBP
    calmat[] <- 0
  }
  if (is.na(ids[1])){
    ids <- as.character(1:length(x))
  } else {
    ids <- as.character(ids)
    if (any(duplicated(ids))){ stop("The values in the ids argument must be unique or left as defaults.") }
  }
  if (length(resOffsets)==1){ resOffsets <- rep(resOffsets,length(x)) }
  if (length(resErrors)==1){ resErrors <- rep(resErrors,length(x)) }
  names(sublist) <- ids
  names(reslist) <- c("metadata","grids")
  opts = NULL
  if (verbose){
    print("Calibrating radiocarbon ages...")
    flush.console()
    if (ncores>1){ print(paste("Running in parallel (standard calibration only) on ",getDoParWorkers()," workers...",sep=""))}
    pb <- txtProgressBar(min=0, max=length(x), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }
  # calibration
  `%dofun%` <- `%do%`
  if (ncores>1){
    # parallellised
    cl <- snow::makeCluster(ncores)
    registerDoSNOW(cl)
    `%dofun%` <- `%dopar%`
  }
  b <- NULL # Added to solve No Visible Binding for Global Variable NOTE
  age = x - resOffsets
  error = sqrt(errors^2 + resErrors^2)
  
  sublist <- foreach (b=1:length(x),.options.snow = opts) %dofun% {
    if (verbose & ncores==1) {setTxtProgressBar(pb, b)}
    #age <- x[b] - resOffsets[b]
    #error <- sqrt(errors[b]^2 + resErrors[b]^2)
    if (F14C==FALSE)
    {  
      dens <- BP14C_calibration(age[b], error[b], cclist2[[calCurves[b]]]$mu, cclist2[[calCurves[b]]]$tau2, eps)
      print(sum(dens))
    }
    if (F14C==TRUE)
    {
      dens <- F14C_calibration(age[b], error[b], cclist2[[calCurves[b]]]$calf14, cclist2[[calCurves[b]]]$calf14error, eps)
    }
    if (normalised){
      dens <- normalise_densities(dens, eps)
    }
    calBPrange = cclist2[[calCurves[b]]]$calBPrange
    calBP = calBPrange[cclist2[[calCurves[b]]]$calBPindex]
    PrDens = dens[cclist2[[calCurves[b]]]$calBPindex]
    #if (anyNA(PrDens)){stop("One or more dates are outside the calibration range")}
    res = list(calBP = calBP[PrDens > 0],PrDens = PrDens[PrDens > 0])
    return(res)
  }
  if (ncores>1){
    snow::stopCluster(cl)
  }
  names(sublist) <- ids
  if (calMatrix){
    for (a in 1:length(sublist)){
      calmat[as.character(sublist[[a]]$calBP),a] <- sublist[[a]]$PrDens
    }
  }
  ## clean-up and results
  if (length(x)>1 & verbose){ close(pb) }
  df <- data.frame(DateID=ids, CRA=x, Error=errors, Details=dateDetails, CalCurve=calCurves,ResOffsets=resOffsets, ResErrors=resErrors, StartBP=timeRange[1], EndBP=timeRange[2], Normalised=normalised, F14C=F14C, CalEPS=eps, stringsAsFactors=FALSE)
  reslist[["metadata"]] <- df
  if (calMatrix){
    reslist[["grids"]] <- NA
    reslist[["calmatrix"]] <- calmat
  } else {
    reslist[["grids"]] <- lapply(sublist,function(x){tmp=data.frame(calBP=x[[1]],PrDens=x[[2]]);class(tmp)=append("calGrid",class(tmp));return(tmp)})
    reslist[["calmatrix"]] <- NA
  }
  class(reslist) <- c("CalDates",class(reslist))
  if (verbose){ print("Done.") }
  return(reslist)
}


# normalise densities to 1
normalise_densities <- function(dens,eps) {
  dens <- dens/sum(dens)
  dens[dens < eps] <- 0
  dens <- dens/sum(dens)
  return(dens)
}

# calibrates in F14C space
F14C_calibration <- function(age, error, calf14,calf14error, eps) {
  f14age <- exp(age/-8033) 
  f14err <- f14age*error/8033 
  p1 <- (f14age - calf14)^2 
  p2 <- 2 * (f14err^2 + calf14error^2) 
  p3 <- sqrt(f14err^2 + calf14error^2) 
  dens <- exp(-p1/p2)/p3 
  dens[dens < eps] <- 0
  return(dens)
}

# calibrates in 14C BP space
BP14C_calibration <- function(age, error, mu, tau2, eps) {
  tau <- error^2 + tau2
  dens <- dnorm(age, mean=mu, sd=sqrt(tau))
  dens[dens < eps] <- 0
  return(dens)
}

# reads a cal curve file from extdata
read_cal_curve_from_file <- function(calCurves) {
  calCurveFile <- paste(system.file("extdata", package="rcarbon"), "/", calCurves,".14c", sep="")
  options(warn=-1)
  calcurve <- readLines(calCurveFile, encoding="UTF-8")
  calcurve <- calcurve[!grepl("[#]",calcurve)]
  calcurve.con <- textConnection(calcurve)
  calcurve <- as.matrix(read.csv(calcurve.con, header=FALSE, stringsAsFactors=FALSE))[,1:3]
  close(calcurve.con)
  options(warn=0)
  colnames(calcurve) <- c("CALBP","C14BP","Error")
  return(calcurve)
}

