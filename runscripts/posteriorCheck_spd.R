# Load Libraries ####
library(coda)
library(here)
library(rcarbon)

# Load Posteriors ####
load(here('R_imagefiles','ABC_posteriors.RData'))
load(here('R_imagefiles/','variables.RData'))

# Source Functions ####
source(here('src','growthModel.R'))
source(here('src','calibrate2.R'))
source(here('src','modelRunner.R'))
source(here('src','randomThin.R'))

# General Settings ####
timeRange = c(800,150)
x.norm=cal.combined.norm
x.nnorm=cal.combined.nnorm
bins = bins
ccurves = list(custom=customCurve)
# Extract Covariates
x1 = palm$PollenPerc
x2 = soi$SOIpr
nrep= 100 # Plot top 100 Models

# Run Models #### 

# Model 1 #
postpcheck.model1.spd = vector('list',length=8)
postpcheck.model1.spd[1:8] = rep(list(matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nrep)),8)
index = c(4,5,2,3,4,5,2,3) # Matching Index so that relevant columns of the object res is extracted

pb <- txtProgressBar(max = nrep*8, style = 3)

for (k in 1:8)
{
  post= model1.posterior[[k]][1:nrep,]
  for (i in 1:nrep)
  {
    setTxtProgressBar(pb, i+(nrep*(k-1)))
    set.seed(as.numeric(row.names(post))[i]) #set the same seed as in the original run
    model = GrowthModel(nt0=post$nt0[i],r=post$r[i],timeRange=timeRange,b1=post$b1[i],b2=post$b2[i],x1=x1,x2=x2,a=1,raw = FALSE)
    postpcheck.model1.spd[[k]][,i]=modelRunner(model=model,x.norm=cal.combined.norm,x.nnorm=cal.combined.nnorm,bins=bins,customCurves=ccurves,timeRange=timeRange,simonly=TRUE)[,index[k]]
  }
}



# Model 2 #
postpcheck.model2.spd = vector('list',length=8)
postpcheck.model2.spd[1:8] = rep(list(matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nrep)),8)
index = c(4,5,2,3,4,5,2,3) # Matching Index so that relevant columns of the object res is extracted

pb <- txtProgressBar(max = nrep*8, style = 3)

for (k in 1:8)
{
  post= model2.posterior[[k]][1:nrep,]
  for (i in 1:nrep)
  {
    setTxtProgressBar(pb, i+(nrep*(k-1)))
    set.seed(as.numeric(row.names(post))[i]) #set the same seed as in the original run
    model = GrowthModel(nt0=post$nt0[i],r=post$r[i],timeRange=timeRange,b1=post$b1[i],b2=post$b2[i],x1=x1,x2=x2,a=1,raw = FALSE)
    postpcheck.model2.spd[[k]][,i]=modelRunner(model=model,x.norm=cal.combined.norm,x.nnorm=cal.combined.nnorm,bins=bins,customCurves=ccurves,timeRange=timeRange,simonly=TRUE)[,index[k]]
  }
}

# Model 3 #
postpcheck.model3.spd = vector('list',length=8)
postpcheck.model3.spd[1:8] = rep(list(matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nrep)),8)
index = c(4,5,2,3,4,5,2,3) # Matching Index so that relevant columns of the object res is extracted

pb <- txtProgressBar(max = nrep*8, style = 3)

for (k in 1:8)
{
  post= model3.posterior[[k]][1:nrep,]
  for (i in 1:nrep)
  {
    setTxtProgressBar(pb, i+(nrep*(k-1)))
    set.seed(as.numeric(row.names(post))[i]) #set the same seed as in the original run
    model = GrowthModel(nt0=post$nt0[i],r=post$r[i],timeRange=timeRange,b1=post$b1[i],b2=post$b2[i],x1=x1,x2=x2,a=1,raw = FALSE)
    postpcheck.model3.spd[[k]][,i]=modelRunner(model=model,x.norm=cal.combined.norm,x.nnorm=cal.combined.nnorm,bins=bins,customCurves=ccurves,timeRange=timeRange,simonly=TRUE)[,index[k]]
  }
}

# Model 4 #
postpcheck.model4.spd = vector('list',length=8)
postpcheck.model4.spd[1:8] = rep(list(matrix(NA,nrow=abs(diff(timeRange))+1,ncol=nrep)),8)
index = c(4,5,2,3,4,5,2,3) # Matching Index so that relevant columns of the object res is extracted

pb <- txtProgressBar(max = nrep*8, style = 3)

for (k in 1:8)
{
  post= model4.posterior[[k]][1:nrep,]
  for (i in 1:nrep)
  {
    setTxtProgressBar(pb, i+(nrep*(k-1)))
    set.seed(as.numeric(row.names(post))[i]) #set the same seed as in the original run
    model = GrowthModel(nt0=post$nt0[i],r=post$r[i],timeRange=timeRange,b1=post$b1[i],b2=post$b2[i],x1=x1,x2=x2,a=1,raw = FALSE)
    postpcheck.model4.spd[[k]][,i]=modelRunner(model=model,x.norm=cal.combined.norm,x.nnorm=cal.combined.nnorm,bins=bins,customCurves=ccurves,timeRange=timeRange,simonly=TRUE)[,index[k]]
  }
}

save(postpcheck.model1.spd,postpcheck.model2.spd,postpcheck.model4.spd,file=here('R_imagefiles','post_pred_spd.RData'))
# save(postpcheck.model1.spd,postpcheck.model2.spd,postpcheck.model3.spd,postpcheck.model4.spd,file=here('R_imagefiles','post_pred_spd.RData'))




