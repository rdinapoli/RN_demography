library(rcarbon)
library(here)
load(here('R_imagefiles','variables.RData')) # Load Calibrated Dates and Custom Curves
source(here('src','randomThin.R'))

# Load data
x.norm=cal.combined.norm
x.nnorm=cal.combined.nnorm
bins = bins
ccurves = list(custom=customCurve)
timeRange = c(800,150)

# Create SPDs
set.seed(123)
index=randomThin(bins)
#Observed SPDs
spd.norm =spd(x.norm[index],datenormalised=FALSE, spdnormalised = TRUE,timeRange=timeRange, verbose=FALSE)
spd.nnorm =spd(x.nnorm[index],datenormalised=FALSE, spdnormalised = TRUE,timeRange=timeRange,verbose=FALSE)
#Observed SPDs with taphonomic correction
spd.norm.trans=transformSPD(spd.norm)
spd.nnorm.trans=transformSPD(spd.nnorm)

pdf(file = here('figures','supplemental','taphonomic_correction_test.pdf'),width = 8,height = 5)
par(mfrow=c(2,2))
plot(spd.nnorm,main="non-normalized")
plot(spd.nnorm.trans,main="non-normalized taphonomic correction")
plot(spd.norm, main="normalized")
plot(spd.norm.trans, main="normalized taphonomic correction")
par(mfrow=c(1,1))
dev.off()