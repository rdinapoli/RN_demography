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
DeltDaR = -214
DeltaRErr = 16
customCurve = mixCurves('shcal20','marine20',p=0.5,resOffsets=DeltDaR,resErrors=DeltaRErr)

# Calibrate Terrestrial and Marine first
cal.terrestrial.norm = calibrate2(RNdates.terrestrial$Age,RNdates.terrestrial$SD,ids=RNdates.terrestrial$Laboratory.ID,calCurves='shcal20',normalised=TRUE)
cal.mixed.norm = calibrate2(RNdates.mixed$Age,RNdates.mixed$SD,ids=RNdates.mixed$Laboratory.ID,calCurves=customCurve,normalised=TRUE)
cal.terrestrial.nnorm = calibrate2(RNdates.terrestrial$Age,RNdates.terrestrial$SD,ids=RNdates.terrestrial$Laboratory.ID,calCurves='shcal20',normalised=FALSE)
cal.mixed.nnorm = calibrate2(RNdates.mixed$Age,RNdates.mixed$SD,ids=RNdates.mixed$Laboratory.ID,calCurves=customCurve,normalised=FALSE)


#Combine Calibrations
cal.combined.norm = rcarbon::combine(cal.terrestrial.norm,cal.mixed.norm)
cal.combined.nnorm = rcarbon::combine(cal.terrestrial.nnorm,cal.mixed.nnorm)

# Match dates order in cal.combined with original SiteID for binning
index=match(RNdates$Laboratory.ID,cal.combined.norm$metadata$DateID)
cal.combined.norm=cal.combined.norm[index] #reorder to match SiteID sequence
cal.combined.nnorm=cal.combined.nnorm[index] #reorder to match SiteID sequence

# Binning, with h=50 as in the original
bins = binPrep(sites = RNdates$SiteID,ages = cal.combined.norm,h=50)

# Create a curve sampler for pre-ABC random thinning, in case a given bin has dates from different different curves
binCurveSelector=table(bins,cal.combined.norm$metadata$CalCurve)
binCurveSelector = prop.table(binCurveSelector,1)

# Create SPD
rn.spd.norm =spd(cal.combined.norm,bins=bins,timeRange=c(800,150))
rn.spd.nnorm =spd(cal.combined.nnorm,bins=bins,timeRange=c(800,150))
# par(mfrow=c(2,1))
# plot(rn.spd.norm)
# plot(rn.spd.nnorm)


# Stack SPD to examine the general impact of different contexts
# contextSPD.norm=stackspd(cal.combined.norm,bins = bins,group = RNdates$Context,timeRange=c(800,150))
# contextSPD.nnorm=stackspd(cal.combined.nnorm,bins = bins,group = RNdates$Context,timeRange=c(800,150))
# plot(contextSPD.norm,legend.arg=list(cex=0.8))
# plot(contextSPD.nnorm,legend.arg=list(cex=0.8))


# Save Dates, Bins, curves, and binCurveselector 
save(cal.combined.norm,cal.combined.nnorm,bins,binCurveSelector,customCurve,file='calibratedDates.RData')





