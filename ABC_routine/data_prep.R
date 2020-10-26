library(rcarbon)
library(readxl)
library(dplyr)
library(here)

# Read and Combine Data
RNdates_lima = read.csv(here('data','date_si_lima.csv'))
RNdates_diNapoli=read.csv('rapanui_DiNapoli_etal.csv')
RNdates = merge(RNdates_diNapoli,RNdates_lima,by.x='Laboratory.ID',by.y='Laboratory.ID',x.all=TRUE)

#Optional Sub-setting
# RNdates = subset(RNdates, Context == "Residential" | Context == "Agricultural")
# Create SiteID field for making life easier
RNdates$SiteID = as.numeric(as.factor(RNdates$SiteName))
## Should this be: RNdates$SiteID = as.numeric(as.factor(RNdates$CorrectSiteDesignation)) ?

# Calibrate separately
RNdates.terrestrial = subset(RNdates,pmarine==0)
RNdates.marine = subset(RNdates,pmarine==1)
RNdates.mixed = subset(RNdates,pmarine>0&pmarine<1)

# Normalise?
normalised=TRUE

# Calibrate Terrestrial and Marine first
cal.terrestrial = calibrate(RNdates.terrestrial$Age,RNdates.terrestrial$SD,ids=RNdates.terrestrial$Laboratory.ID,calCurves='shcal20',normalised=normalised)
cal.marine = calibrate(RNdates.marine$Age,RNdates.marine$SD,ids=RNdates.marine$Laboratory.ID,calCurves='marine20',resOffsets=-83,resErrors=34,normalised=normalised)

# Calibrate Mixed... Each unique terrestrial/marine combo has its own curve, and some dates cannot be calibrated
customCurves = vector('list',length=nrow(RNdates.mixed))
check = vector(length=nrow(RNdates.mixed)) #store whether the date was successfully calibrated or not via tryCatch

for (i in 1:nrow(RNdates.mixed))
{
	if (i==1)
	{
		customCurves[[i]] = mixCurves('shcal20','marine20',p=1-RNdates.mixed$pmarine[i],resOffsets=-87,resErrors=34)
		cal.mixed = tryCatch(calibrate(RNdates.mixed$Age[i],RNdates.mixed$SD[i],ids=RNdates.mixed$Laboratory.ID[i],calCurves=customCurves[[i]],normalised=normalised),error=function(e){return(NULL)})
		check[i] = ifelse(is.null(cal.mixed),FALSE,TRUE)
		cal.mixed$metadata$CalCurve = paste0('custom',i)
	}

	if (i>1)
	{

 		customCurves[[i]] = mixCurves('shcal20','marine20',p=1-RNdates.mixed$pmarine[i],resOffsets=-87,resErrors=34)
		cal.mixed.tmp = tryCatch(calibrate(RNdates.mixed$Age[i],RNdates.mixed$SD[i],ids=RNdates.mixed$Laboratory.ID[i],calCurves=customCurves[[i]],normalised=normalised),error=function(e){return(NULL)})
		check[i] = ifelse(is.null(cal.mixed.tmp),FALSE,TRUE)
		if (check[i]==TRUE){
		  cal.mixed.tmp$metadata$CalCurve = paste0('custom',i)
		  cal.mixed = rcarbon::combine(cal.mixed,cal.mixed.tmp)
		}
	}
}

# which dates were excluded? 
print(RNdates.mixed[which(!check),])

#Combine Calibrations
cal.combined = rcarbon::combine(cal.terrestrial,cal.marine,cal.mixed)

# Match dates order in cal.combined with original SiteID for binning
RNdatesWithoutExcludedDates = subset(RNdates,!Laboratory.ID%in%RNdates.mixed[which(!check),1])
index=match(RNdatesWithoutExcludedDates$Laboratory.ID,cal.combined$metadata$DateID)
cal.combined=cal.combined[index] #reorder to match SiteID sequence

# Binning, with h=50 as in the original
bins = binPrep(sites = RNdatesWithoutExcludedDates$SiteID,ages = cal.combined,h=50)

# Create a curve sampler for pre-ABC random thinning, in case a given bin has dates from different different curves
binCurveSelector=table(bins,cal.combined$metadata$CalCurve)
#which(apply(binCurveSelector,1,function(x){return(sum(x>0))})>1) #damn, just one bin having this issue...
binCurveSelector = prop.table(binCurveSelector,1)

# Create SPD
rn.spd =spd(cal.combined,bins=bins,timeRange=c(800,0),runm=100)
plot(rn.spd)

# Save Dates, Bins, curves, and binCurveselector 
save(cal.combined,bins,binCurveSelector,customCurves,file='calibratedDates.RData')





