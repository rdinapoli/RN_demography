# Load required libraries ####
library(rcarbon)
library(here)

# Source modified calibration script ####
source(here('src','calibrate2.R'))

# Read 14C Dates ####
RNdates=read.csv(here('data','rapanui_DiNapoli_etal.csv'))

# Sub-setting to Anthrophic Contexts
RNdates = subset(RNdates, Context %in% c("Residential","Agricultural","Burial"))

# Assign SiteID 
RNdates$SiteID = as.numeric(as.factor(RNdates$CorrectSiteDesignation))

# Calibrate separately
RNdates.terrestrial = subset(RNdates,CalCurve=='shcal20')
RNdates.mixed = subset(RNdates,CalCurve=='Bone')

# Calibration Settings for mixed terrestial/marine samples
DeltDaR = -214
DeltaRErr = 16
customCurve = mixCurves('shcal20','marine20',p=0.5,resOffsets=DeltDaR,resErrors=DeltaRErr)

# Calibrate (both normalised and non-normalised)
cal.terrestrial.norm = calibrate2(RNdates.terrestrial$Age,RNdates.terrestrial$SD,ids=RNdates.terrestrial$Laboratory.ID,calCurves='shcal20',normalised=TRUE)
cal.mixed.norm = calibrate2(RNdates.mixed$Age,RNdates.mixed$SD,ids=RNdates.mixed$Laboratory.ID,calCurves=customCurve,normalised=TRUE)
cal.terrestrial.nnorm = calibrate2(RNdates.terrestrial$Age,RNdates.terrestrial$SD,ids=RNdates.terrestrial$Laboratory.ID,calCurves='shcal20',normalised=FALSE)
cal.mixed.nnorm = calibrate2(RNdates.mixed$Age,RNdates.mixed$SD,ids=RNdates.mixed$Laboratory.ID,calCurves=customCurve,normalised=FALSE)

#Combine Calibrations
cal.combined.norm = rcarbon::combine(cal.terrestrial.norm,cal.mixed.norm)
cal.combined.nnorm = rcarbon::combine(cal.terrestrial.nnorm,cal.mixed.nnorm)

# Match dates order with original SiteID for binning
index=match(RNdates$Laboratory.ID,cal.combined.norm$metadata$DateID)
cal.combined.norm=cal.combined.norm[index] #reorder to match SiteID sequence
cal.combined.nnorm=cal.combined.nnorm[index] #reorder to match SiteID sequence

# Binning, with h=50 as in the original
bins = binPrep(sites = RNdates$SiteID,ages = cal.combined.norm,h=50)

# Read Environmental ####
# Linear Interpolation of Palm% from Table S5 of Lima et al 2020
tableS5=read.csv(here('data','tableS5.csv'))
tableS5$CalBP = BCADtoBP(tableS5$BCCE)
palm = data.frame(CalBP=800:150)
palm$PollenPerc = approx(x=tableS5$CalBP[which(!is.na(tableS5$Palm))],y=tableS5$Palm[which(!is.na(tableS5$Palm))],xout=palm$CalBP)$y

# SOI from Yan et al 2011
soi = read.csv(here('yan2011soipr.txt'),skip=109,sep='')
soi = subset(soi,Age.AD.<1900)
soi$CalBP = BCADtoBP(soi$Age.AD.)
soi = subset(soi,CalBP<=800&CalBP>=150)
soi = soi[,c(3,2)]



# Save Dates, Bins, curves, and binCurveselector 
save(cal.combined.norm,cal.combined.nnorm,customCurve,bins,soi,palm,file=here('R_imagefiles','variables.RData'))





