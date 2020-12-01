library(rcarbon)
library(here)

# Setting for simulations
time_range <- c(800, 200)
n_cores <- 8
n_sim <- 1000
set.seed(123)

# Generate uniformly distributed calendar dates and uncalibrate them
Cals <- data.frame(calBP=800:200,PrDens=runif(length(800:200)))
Cals <- as.CalGrid(Cals)
Cals <- uncalibrate(Cals, calCurves = 'shcal20')

# Convert to dataframe with constant error 
test <- data.frame(CRAs=Cals$CRA, Ers=rep(30, times=length(Cals$CRA))) 

# Calibrate with and without normalization
test_norm <- calibrate(test$CRAs, test$Ers, calCurves = 'shcal20', normalised = T)
test_nnorm <- calibrate(test$CRAs, test$Ers, calCurves = 'shcal20', normalised = F)

# Test calsample and uncalsample with normalized and non-normalized dates compared to uniform null model
uni_norm_cal <- modelTest(x=test_norm, errors=test$Ers, normalised = T, spdnormalised = T, nsim=n_sim, timeRange = time_range, 
                           model='uniform', method='calsample', ncores=n_cores)
uni_nnorm_cal <- modelTest(x=test_nnorm, errors=test$Ers, normalised = F, spdnormalised = T, nsim=n_sim, timeRange = time_range, 
                            model='uniform', method='calsample', ncores=n_cores)
uni_norm_uncal <- modelTest(x=test_norm, errors=test$Ers, normalised = T, spdnormalised = T, nsim=n_sim, timeRange = time_range, 
                             model='uniform', method='uncalsample', ncores=n_cores)
uni_nnorm_uncal <- modelTest(x=test_nnorm, errors=test$Ers, normalised = F, spdnormalised = T, nsim=n_sim, timeRange = time_range, 
                              model='uniform', method='uncalsample', ncores=n_cores)

jpeg(file=here('figures','temporary','back-calibration_simulation.jpeg'), width=8, height=8, units='in', res=300)
par(mfrow=c(2,2))
plot(uni_norm_cal, lwd=2, main='normalized, calsample')
plot(uni_nnorm_cal, lwd=2,  main='non-normalized, calsample')
plot(uni_norm_uncal, lwd=2, main='normalized, uncalsample')
plot(uni_nnorm_uncal, lwd=2, main='non-normalized, uncalsample')
par(mfrow=c(1,1))
dev.off()

tiff(file=here('figures','temporary','back-calibration_simulation.tiff'),compression = 'lzw', width=8, height=8, units='in', res=600)
par(mfrow=c(2,2))
plot(uni_norm_cal, lwd=2, main='normalized, calsample')
plot(uni_nnorm_cal, lwd=2,  main='non-normalized, calsample')
plot(uni_norm_uncal, lwd=2, main='normalized, uncalsample')
plot(uni_nnorm_uncal, lwd=2, main='non-normalized, uncalsample')
par(mfrow=c(1,1))
dev.off()