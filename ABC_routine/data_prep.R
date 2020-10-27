library(rcarbon)
library(readxl)
library(dplyr)
library(here)
# Source modified calibration script
source('./ABC_routine/calibrate2.R')

# Read Data
RNdates=read.csv('rapanui_DiNapoli_etal.csv')

# Sub-setting
RNdates = subset(RNdates, Context %in% c("Residential","Agricultural","Burial"))
# Create SiteID field for making life easier
RNdates$SiteID = as.numeric(as.factor(RNdates$CorrectSiteDesignation))

# Calibrate separately
RNdates.terrestrial = subset(RNdates,CalCurve=='shcal20')
RNdates.mixed = subset(RNdates,CalCurve=='Bone')

# Calibration Settings
normalised=TRUE
DeltDaR = -214
DeltaRErr = 16
customCurve = mixCurves('shcal20','marine20',p=0.5,resOffsets=DeltDaR,resErrors=DeltaRErr)

# Calibrate Terrestrial and Marine first
cal.terrestrial = calibrate2(RNdates.terrestrial$Age,RNdates.terrestrial$SD,ids=RNdates.terrestrial$Laboratory.ID,calCurves='shcal20',normalised=normalised)
cal.mixed = calibrate2(RNdates.mixed$Age,RNdates.mixed$SD,ids=RNdates.mixed$Laboratory.ID,calCurves=customCurve,normalised=normalised)

#Combine Calibrations
cal.combined = rcarbon::combine(cal.terrestrial,cal.mixed)

# Match dates order in cal.combined with original SiteID for binning
index=match(RNdates$Laboratory.ID,cal.combined$metadata$DateID)
cal.combined=cal.combined[index] #reorder to match SiteID sequence

# Binning, with h=50 as in the original
bins = binPrep(sites = RNdates$SiteID,ages = cal.combined,h=50)

# Create a curve sampler for pre-ABC random thinning, in case a given bin has dates from different different curves
binCurveSelector=table(bins,cal.combined$metadata$CalCurve)
#which(apply(binCurveSelector,1,function(x){return(sum(x>0))})>1) #damn, just one bin having this issue...
binCurveSelector = prop.table(binCurveSelector,1)

# Create SPD
rn.spd =spd(cal.combined,bins=bins,timeRange=c(800,0))
plot(rn.spd)

# Check sensitivity to mixed dates
index=match(RNdates.terrestrial$Laboratory.ID,cal.terrestrial$metadata$DateID)
cal.terrestrial=cal.terrestrial[index] #reorder to match SiteID sequence
bins.terrestrial = binPrep(sites = RNdates.terrestrial$SiteID,ages = cal.terrestrial,h=50)
rn.spd.terrestrial = spd(cal.terrestrial,bins=bins.terrestrial,timeRange=c(800,0))
lines(rn.spd.terrestrial$grid,col=2,lty=2)

# Stack SPD to examine the general impact of different contexts
contextSPD=stackspd(cal.combined,bins = bins,group = RNdates$Context,timeRange=c(800,0))
plot(contextSPD,legend.arg=list(cex=0.8))


# Save Dates, Bins, curves, and binCurveselector 
save(cal.combined,bins,binCurveSelector,customCurves,file='calibratedDates.RData')





