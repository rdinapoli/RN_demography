library(rcarbon)
library(ggplot2)
library(gridExtra)
library(here)
setwd(here())

##################################
##### LOAD DATA AND CALIBRATE ####
##################################

#Load corrected dataset
RN_c14 <- read.csv("rapanui_DiNapoli_etal.csv")

# having trouble calibrating bone dates with mixed curve, removing them and marine dates for now
RN_c14 <- subset(RN_c14, CalCurve == "shcal20")

# restrict to just dates from Residential and Agricultural contexts
RN_c14 <- subset(RN_c14, Context == "Residential" | Context == "Agricultural")

#Define correct time range
RN_timerange <- c(800,0)

#calibrate dates without normalization
RN_caldates <- calibrate(x=RN_c14$Age, RN_c14$SD, calCurves='shcal20', normalised=F,  calMatrix = F, verbose = F)

#calibrate dates with normalization
RN_caldates_norm <- calibrate(x=RN_c14$Age, RN_c14$SD, calCurves='shcal20', normalised=T,  calMatrix = F, verbose = F)

#use corrected site bins
RN_bins <- binPrep(sites=RN_c14$CorrectSiteDesignation, ages=RN_c14$Age, h=50)

#create SPD with non-normalized dates
RN.spd <- spd(RN_caldates, timeRange = RN_timerange, bins=RN_bins, spdnormalised = T, runm=100)

#create SPD with normalized dates
RN.spd_norm <- spd(RN_caldates_norm, timeRange = RN_timerange, bins=RN_bins, spdnormalised = T, runm=100)

#compare normalized and non-normalized spd
par(mfrow=c(1,2))
plot(RN.spd, calendar="BCAD", main="SPD w/o normalization")
plot(RN.spd_norm, calendar="BCAD", main="SPD w/normalization")
par(mfrow=c(1,1))


#########################
### LINEAR MODELS #######
#########################

set.seed(12345)
linear_model1 <- modelTest(x=RN_caldates, bins=RN_bins, datenormalised = F, 
                        spdnormalised = TRUE, errors=RN_c14$SD, 
                        runm=100,timeRange=RN_timerange, model="linear",nsim=5000, method='uncalsample', 
                        ncores=3, raw=TRUE)
set.seed(12345)
linear_model2 <- modelTest(x=RN_caldates, bins=RN_bins, datenormalised = F, 
                           spdnormalised = TRUE, errors=RN_c14$SD, 
                           runm=100,timeRange=RN_timerange, model="linear",nsim=5000, method='calsample', 
                           ncores=3, raw=TRUE)
set.seed(12345)
linear_model3 <- modelTest(x=RN_caldates_norm, bins=RN_bins, datenormalised = T, 
                           spdnormalised = TRUE, errors=RN_c14$SD, 
                           runm=100,timeRange=RN_timerange, model="linear",nsim=5000, method='uncalsample', 
                           ncores=3, raw=TRUE)
set.seed(12345)
linear_model4 <- modelTest(x=RN_caldates_norm, bins=RN_bins, datenormalised = T, 
                           spdnormalised = TRUE, errors=RN_c14$SD, 
                           runm=100,timeRange=RN_timerange, model="linear",nsim=5000, method='calsample', 
                           ncores=3, raw=TRUE)

#### plot linear models

par(mfrow=c(2,2))
plot(linear_model1, calendar="BCAD", main="non-normalized, uncalsample", xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
plot(linear_model2, calendar="BCAD", main="non-normalized, calsample", xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
plot(linear_model3, calendar="BCAD", main="normalized, uncalsample", xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
plot(linear_model4, calendar="BCAD", main="normalized, calsample", xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
par(mfrow=c(1,1))
########################
### UNIFORM MODELS ####
########################

#test hypothesis that there was uniform population between 1400 and 1722
set.seed(12345)
uni_model1 <- modelTest(x=RN_caldates, bins=RN_bins, datenormalised = FALSE, 
                        spdnormalised = TRUE, errors=RN_c14$SD,  
                        runm=100,timeRange=c(550, 100), model="uniform",nsim=5000, method='uncalsample', 
                        ncores=3, raw=TRUE)

set.seed(12345)
uni_model2 <- modelTest(x=RN_caldates, bins=RN_bins, datenormalised = FALSE, 
                        spdnormalised = TRUE, errors=RN_c14$SD,  
                        runm=100,timeRange=c(550, 100), model="uniform",nsim=5000, method='calsample', 
                        ncores=3, raw=TRUE)
set.seed(12345)
uni_model3 <- modelTest(x=RN_caldates_norm, bins=RN_bins, datenormalised = TRUE, 
                        spdnormalised = TRUE, errors=RN_c14$SD, 
                        runm=100,timeRange=c(550, 100), model="uniform",nsim=5000, method='uncalsample', 
                        ncores=3, raw=TRUE)

set.seed(12345)
uni_model4 <- modelTest(x=RN_caldates_norm, bins=RN_bins, datenormalised = TRUE, 
                        spdnormalised = TRUE, errors=RN_c14$SD, 
                        runm=100,timeRange=c(550, 100), model="uniform",nsim=5000, method='calsample', 
                        ncores=3, raw=TRUE)

#### plot uniform models

par(mfrow=c(2,2))
plot(uni_model1, calendar="BCAD", main="non-normalized, uncalsample", xlim=c(1400, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
plot(uni_model2, calendar="BCAD", main="non-normalized, calsample", xlim=c(1400, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
plot(uni_model3, calendar="BCAD", main="normalized, uncalsample", xlim=c(1400, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
plot(uni_model4, calendar="BCAD", main="normalized, calsample", xlim=c(1400, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
par(mfrow=c(1,1))


#########################################
# LIMA ET AL'S LOGISTIC GROWTH MODELS ###
#########################################
# note that we're using the SPD not rate, but the form of the equations matches 
# the equations they present in their paper

# in order to do Lima's models, we need finer-scale interpolations of palm cover

# Load Lima et al.'s supplementary data with Palm cover and SOI variables
dat1<-read.csv('RAPANUINW.csv')
# create linear interpolation of palms for all years
yr <- dat1$YEAR
palms <- dat1$PalmRR
palm_interp <- data.frame(yr, palms)

#linear interpolation
palm_interp_lin <- data.frame(
  with(palm_interp,
       approx(yr, palms, xout=seq(1150, 1950, by = 1), method='linear')),
  method='approx()'
)
palms <- palm_interp_lin$y
#SOI index
# do linear interpolation for now but should update

SOI <- dat1$SOI1
SOI_interp <- data.frame(yr, SOI)

SOI_interp_lin <- data.frame(
  with(SOI_interp,
       approx(yr, SOI, xout=seq(1150, 1950, by = 1), method='linear')),
  method='approx()'
)
SOI <- SOI_interp_lin$y

#######################
# Create dataframe's to fit models to

#create dataframe of values for non-normalized dates
grd_rn <- RN.spd$grid
x <- grd_rn$calBP
y <- grd_rn$PrDens
#compute Nt-1 with a 30 yr backsight
Nt_1 <- function(year,backsight)
{
  obs=rep(NA,length(year))
  
  for (i in 1:30)
  {
    obs[i] = 0
  }
  for (i in 31:c(length(obs)))
  {
    t_1 = year[i-backsight]
    obs[i] = t_1
  }
  return(obs)
}

minus1 <- Nt_1(y, 30)  

new_dat <- data.frame(x, y, palms, minus1, SOI )
#change NA to 0, probably need to think about impacts of this
new_dat[is.na(new_dat)] <- 0

#create dataframe of values for normalized dates
grd_rn_norm <- RN.spd_norm$grid
x_norm <- grd_rn_norm$calBP
y_norm <- grd_rn_norm$PrDens
#compute Nt-1 with a 30 yr backsite

minus1_norm <- Nt_1(y_norm, 30)  

new_dat_norm <- data.frame(x_norm, y_norm, palms, minus1_norm, SOI )
#change NA to 0
new_dat_norm[is.na(new_dat_norm)] <- 0

##################################################
# Build logistic models for non-normalized dates #
##################################################

#basic logistic model
log1 <- nls(y ~ minus1 * exp(b * (1 - (minus1/k))), data=new_dat, start=list(b=0.47, k=0.002), trace=TRUE)
summary(log1)
log1_fit <- data.frame(calBP=x,PrDens=predict(log1))

#their model with palm cover, had to increase d to get a good fit
log2 <- nls(y ~ minus1 * exp(b * (1 - (minus1/(k+d*palms)))), data=new_dat, start=list(b=0.47, k=0.002, d=0.0001), trace=TRUE)
summary(log2)
log2_fit <- data.frame(calBP=x,PrDens=predict(log2))

#their model with SOI index
log3 <- nls(y ~ minus1 * exp(b * (1 - (minus1/(k+d*SOI)))), data=new_dat, start=list(b=0.47, k=0.002, d=0.001), trace=TRUE)
summary(log3)
log3_fit <- data.frame(calBP=x,PrDens=predict(log3))

#their model with palms and SOI index
#can't get this model to fit
log4 <- nls(y ~ minus1 * exp(b * (1 - (minus1/(k+d*palms+e*SOI)))), data=new_dat, start=list(b=0.47, k=0.002, d=0.001, e=-0.001), trace=TRUE)
#summary(log4)
#log4_fit <- data.frame(calBP=x,PrDens=predict(log4))

##############################################
# Build logistic models for normalized dates #
##############################################

#basic logistic model
log1_norm <- nls(y_norm ~ minus1_norm * exp(b * (1 - (minus1_norm/k))), data=new_dat_norm, start=list(b=0.47, k=0.002), trace=TRUE)
summary(log1_norm)
log1_fit_norm <- data.frame(calBP=x_norm,PrDens=predict(log1_norm))

#their model with palm cover 
#had to decrease d to get a good fit
log2_norm <- nls(y_norm ~ minus1_norm * exp(b * (1 - (minus1_norm/(k+d*palms)))), data=new_dat_norm, start=list(b=0.47, k=0.002, d=0.0001), trace=TRUE)
summary(log2_norm)
log2_fit_norm <- data.frame(calBP=x_norm,PrDens=predict(log2_norm))

#their model with SOI index
# wont fit with a normalized spd, had to decrease d to get good fit
log3_norm <- nls(y_norm ~ minus1_norm * exp(b * (1 - (minus1_norm/(k+d*SOI)))), data=new_dat_norm, start=list(b=0.47, k=0.002, d=0.0001), trace=TRUE)
summary(log3_norm)
log3_fit_norm <- data.frame(calBP=x_norm,PrDens=predict(log3_norm))

#their model with palms and SOI index
#can't get this model to fit
log4_norm <- nls(y_norm ~ minus1_norm * exp(b * (1 - (minus1_norm/(k+d*palms+e*SOI)))), data=new_dat_norm, start=list(b=0.47, k=0.002, d=0.001, e=-0.001), trace=TRUE)
#summary(log4_norm)
#log4_fit_norm <- data.frame(calBP=x_norm,PrDens=predict(log4_norm))


#####################################
# Run logistic models in rcarbon ####
#####################################

#####################################
# General logistic model

#logistic growth model with non-normalized dates and with method set to 'uncalsample'
set.seed(12345)
log1_model1 <- modelTest(x=RN_caldates, bins=RN_bins, datenormalised = FALSE, 
                        spdnormalised = TRUE, errors=RN_c14$SD, predgrid=log1_fit, 
                        runm=100,timeRange=RN_timerange, model="custom",nsim=5000, method='uncalsample', 
                        ncores=3, raw=TRUE)
#logistic growth model with non-normalized dates and with method set to 'calsample'
set.seed(12345)
log1_model2 <- modelTest(x=RN_caldates, bins=RN_bins, datenormalised = FALSE, 
                        spdnormalised = TRUE, errors=RN_c14$SD, predgrid=log1_fit, 
                        runm=100,timeRange=RN_timerange, model="custom",nsim=5000, method='calsample', 
                        ncores=3, raw=TRUE)
#logistic growth model normalized dates and method set to uncalsample
set.seed(12345)
log1_model3 <- modelTest(x=RN_caldates_norm, bins=RN_bins, datenormalised = TRUE, 
                        spdnormalised = TRUE, errors=RN_c14$SD, predgrid=log1_fit_norm, 
                        runm=100,timeRange=RN_timerange, model="custom",nsim=5000, method='uncalsample', 
                        ncores=3, raw=TRUE)

#logistic growth model normalized dates and method set to calsample
set.seed(12345)
log1_model4 <- modelTest(x=RN_caldates_norm, bins=RN_bins, datenormalised = TRUE, 
                        spdnormalised = TRUE, errors=RN_c14$SD, predgrid=log1_fit_norm, 
                        runm=100,timeRange=RN_timerange, model="custom",nsim=5000, method='calsample', 
                        ncores=3, raw=TRUE)

#print summaries of models
summary(log1_model1)
summary(log1_model2)
summary(log1_model3)
summary(log1_model4)

#compare results

par(mfrow=c(2,2))
plot(log1_model1, calendar='BCAD', main='non-normalized, uncalsample', ylim=c(0,0.0035), xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5) #Lima et al.'s proposed 'collapses'
plot(log1_model2, calendar='BCAD', main='non-normalized, calsample', ylim=c(0,0.0035), xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
plot(log1_model3, calendar='BCAD', main='normalized, uncalsample', ylim=c(0,0.0035), xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
plot(log1_model4, calendar='BCAD', main='normalized, calsample', ylim=c(0,0.0035), xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
par(mfrow=c(1,1))


# #rates of change
# par(mfrow=c(2,2))
# plot(log1_model1, type='roc', calendar='BCAD', main='non-normalized, uncalsample')#, ylim=c(0,0.0035))
# abline(v=1722, col='red')
# abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
# plot(log1_model3, type='roc', calendar='BCAD', main='normalized, uncalsample')#, ylim=c(0,0.0035))
# abline(v=1722, col='red')
# abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
# plot(log1_model2, type='roc', calendar='BCAD', main='non-normalized, calsample')#, ylim=c(0,0.0035))
# abline(v=1722, col='red')
# abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
# plot(log1_model4, type='roc', calendar='BCAD', main='normalized, calsample')#, ylim=c(0,0.0035))
# abline(v=1722, col='red')
# abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
# par(mfrow=c(1,1))

#####################################
# Logistic model with palms as lateral perturbation effect

#logistic growth model with palms with non-normalized dates and with method set to 'uncalsample'
set.seed(12345)
log2_model1 <- modelTest(x=RN_caldates, bins=RN_bins, datenormalised = FALSE, 
                         spdnormalised = TRUE, errors=RN_c14$SD, predgrid=log2_fit, 
                         runm=100,timeRange=RN_timerange, model="custom",nsim=5000, method='uncalsample', 
                         ncores=3, raw=TRUE)
#logistic growth model with palms with non-normalized dates and with method set to 'calsample'
set.seed(12345)
log2_model2 <- modelTest(x=RN_caldates, bins=RN_bins, datenormalised = FALSE, 
                         spdnormalised = TRUE, errors=RN_c14$SD, predgrid=log2_fit, 
                         runm=100,timeRange=RN_timerange, model="custom",nsim=5000, method='calsample', 
                         ncores=3, raw=TRUE)
#logistic growth model with palms normalized dates and method set to uncalsample
set.seed(12345)
log2_model3 <- modelTest(x=RN_caldates_norm, bins=RN_bins, datenormalised = TRUE, 
                         spdnormalised = TRUE, errors=RN_c14$SD, predgrid=log2_fit_norm, 
                         runm=100,timeRange=RN_timerange, model="custom",nsim=5000, method='uncalsample', 
                         ncores=3, raw=TRUE)
#logistic growth model with palms normalized dates and method set to calsample
set.seed(12345)
log2_model4 <- modelTest(x=RN_caldates_norm, bins=RN_bins, datenormalised = TRUE, 
                         spdnormalised = TRUE, errors=RN_c14$SD, predgrid=log2_fit_norm, 
                         runm=100,timeRange=RN_timerange, model="custom",nsim=5000, method='calsample', 
                         ncores=3, raw=TRUE)

#print summaries of models
summary(log2_model1)
summary(log2_model2)
summary(log2_model3)
summary(log2_model4)

#compare results
par(mfrow=c(2,2))
plot(log2_model1, calendar='BCAD', main='non-normalized, uncalsample', ylim=c(0,0.0035), xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
plot(log2_model2, calendar='BCAD', main='non-normalized, calsample', ylim=c(0,0.0035), xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
plot(log2_model3, calendar='BCAD', main='normalized, uncalsample', ylim=c(0,0.0035), xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
plot(log2_model4, calendar='BCAD', main='normalized, calsample', ylim=c(0,0.0035), xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5, xlim=c(1150, 1750))
par(mfrow=c(1,1))


#####################################
# Logistic model with SOI as lateral perturbation effect

#logistic growth model with SOI with non-normalized dates and with method set to 'uncalsample'
set.seed(12345)
log3_model1 <- modelTest(x=RN_caldates, bins=RN_bins, datenormalised = FALSE, 
                         spdnormalised = TRUE, errors=RN_c14$SD, predgrid=log3_fit, 
                         runm=100,timeRange=RN_timerange, model="custom",nsim=5000, method='uncalsample', 
                         ncores=3, raw=TRUE)
#logistic growth model with SOI with non-normalized dates and with method set to 'calsample'
set.seed(12345)
log3_model2 <- modelTest(x=RN_caldates, bins=RN_bins, datenormalised = FALSE, 
                         spdnormalised = TRUE, errors=RN_c14$SD, predgrid=log3_fit, 
                         runm=100,timeRange=RN_timerange, model="custom",nsim=5000, method='calsample', 
                         ncores=3, raw=TRUE)
#logistic growth model with SOI normalized dates and method set to uncalsample
set.seed(12345)
log3_model3 <- modelTest(x=RN_caldates_norm, bins=RN_bins, datenormalised = TRUE, 
                         spdnormalised = TRUE, errors=RN_c14$SD, predgrid=log3_fit_norm, 
                         runm=100,timeRange=RN_timerange, model="custom",nsim=5000, method='uncalsample', 
                         ncores=3, raw=TRUE)

#logistic growth model with SOI normalized dates and method set to calsample
set.seed(12345)
log3_model4 <- modelTest(x=RN_caldates_norm, bins=RN_bins, datenormalised = TRUE, 
                         spdnormalised = TRUE, errors=RN_c14$SD, predgrid=log3_fit_norm, 
                         runm=100,timeRange=RN_timerange, model="custom",nsim=5000, method='calsample', 
                         ncores=3, raw=TRUE)

#print summaries of models
summary(log3_model1)
summary(log3_model2)
summary(log3_model3)
summary(log3_model4)

#compare results
par(mfrow=c(2,2))
plot(log3_model1, calendar='BCAD', main='non-normalized, uncalsample', ylim=c(0,0.0035), xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
plot(log3_model2, calendar='BCAD', main='non-normalized, calsample', ylim=c(0,0.0035), xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
plot(log3_model3, calendar='BCAD', main='normalized, uncalsample', ylim=c(0,0.0035), xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
plot(log3_model4, calendar='BCAD', main='normalized, calsample', ylim=c(0,0.0035), xlim=c(1150, 1750))
abline(v=1722, col='red')
abline(v=c(1430, 1550, 1640, 1700), col=gray(level=0.1, alpha=1), lty=2, lwd=0.5)
par(mfrow=c(1,1))



###################################################
# Plot relationship between SPDs and variables ####
###################################################


#plot SPD over time as a function of palm cover and SOI index (similar to Lima et al. figure 2)
#non-normalized SPD as a function of palms and SOI
p1 <- ggplot(data=new_dat, aes(x=x, y=y)) + 
  geom_rect(aes(xmin=520, xmax=400, ymin=0, ymax=Inf), color=NA, fill='gray60', alpha=0.009) +
  geom_rect(aes(xmin=310, xmax=250, ymin=0, ymax=Inf), color=NA, fill='gray60', alpha=0.009) +
  geom_line(aes(color=palms), size=3) + 
  scale_color_viridis_c() +
  labs(subtitle="a) non-normalized SPD") +
  xlab("cal BP")+
  ylab("SPD") +
  xlim(800,200) +
  geom_vline(xintercept=228, linetype='dashed', color='black', size=1)
p2 <- ggplot(data=new_dat, aes(x=x, y=y)) + 
  geom_rect(aes(xmin=520, xmax=400, ymin=0, ymax=Inf), color=NA, fill='gray60', alpha=0.009) +
  geom_rect(aes(xmin=310, xmax=250, ymin=0, ymax=Inf), color=NA, fill='gray60', alpha=0.009) +
  geom_line(aes(color=SOI), size=3) + 
  scale_color_viridis_c(option='magma')+
  labs(subtitle="b)") +
  xlab("cal BP")+
  ylab("SPD") +
  xlim(800,200) +
  geom_vline(xintercept=228, linetype='dashed', color='black', size=1)
#non-normalized SPD as a function of palms and SOI
p3 <- ggplot(data=new_dat_norm, aes(x=x_norm, y=y_norm)) + 
  geom_rect(aes(xmin=520, xmax=400, ymin=0, ymax=Inf), color=NA, fill='gray60', alpha=0.009) +
  geom_rect(aes(xmin=310, xmax=250, ymin=0, ymax=Inf), color=NA, fill='gray60', alpha=0.009) +
  geom_line(aes(color=palms), size=3) + 
  scale_color_viridis_c() +
  labs(subtitle="c) normalized SPD") +
  xlab("cal BP")+
  ylab("SPD") +
  xlim(800,200) +
  geom_vline(xintercept=228, linetype='dashed', color='black', size=1)
p4 <- ggplot(data=new_dat_norm, aes(x=x_norm, y=y_norm)) + 
  geom_rect(aes(xmin=520, xmax=400, ymin=0, ymax=Inf), color=NA, fill='gray60', alpha=0.009) +
  geom_rect(aes(xmin=310, xmax=250, ymin=0, ymax=Inf), color=NA, fill='gray60', alpha=0.009) +
  geom_line(aes(color=SOI), size=3) + 
  scale_color_viridis_c(option='magma')+
  labs(subtitle="d)") +
  xlab("cal BP")+
  ylab("SPD") +
  xlim(800,200) +
  geom_vline(xintercept=228, linetype='dashed', color='black', size=1)
grid.arrange(p1,p2,p3,p4, nrow=2)


# tiff('spds_palm_SOI.tiff', width=12, height=8, compression='lzw', units='in', res=600)
# grid.arrange(p1,p2,p3,p4, nrow=2)
# dev.off()

#####################################
###### Plot the SHCal20 curve #######
#####################################

# The SHCal20 radiocarbon calibration curve
# For reader reference we also include a plot of the SHCal20 calibration curve here (Hogg et al. 2010).
# In the graph of the SHCal20 curve below, we have shaded the two critical regions of interest for Lima et al.'s analysis. 
# The light red region shows the relatively steep portion of the curve from ca. 1475-1500 AD, which is followed by a plateau 
# and several 'wiggles' (blue shading).

#This code is based on the rcarbon documentation at https://rdrr.io/github/ahb108/rcarbon/src/R/calibration-helpers.R
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

# Plot the curve with critical areas shaded
sh20 <- read_cal_curve_from_file('shcal20')
plot(sh20[,2]~sh20[,1], xlim=c(1050, 0), ylim=c(100,1200), type='l', xaxt='n', col='blue', xlab="cal. AD", ylab="C14 BP")
segments(sh20[,1], sh20[,2] - sh20[,3], sh20[,1], sh20[,2] + sh20[,3], col='blue')
lablist<-as.vector(seq(50, 1050, by=100))
axis(1, at=seq(50, 1050, by=100), labels=BPtoBCAD(lablist))
rect(575,50,450,1300, col=rgb(1,0,0,0.1), border=F)
rect(450,50,0,1300, col=rgb(0,0,1,0.1), border=F)
