# Load Posterior Predictive Check & Observed SPDs ####
load(here('R_imagefiles','post_pred_spd.RData'))
load(here('R_imagefiles/','variables.RData'))

# Source R Functions ####
source(here('src','plotFunctions.R'))
source(here('src','randomThin.R'))

# Observed (thinned) SPD ####
index=randomThin(bins)
obs.norm =spd(cal.combined.norm[index],datenormalised=FALSE, spdnormalised = TRUE,timeRange=c(800,150),verbose=FALSE)
obs.nnorm =spd(cal.combined.nnorm[index],datenormalised=FALSE, spdnormalised = TRUE,timeRange=c(800,150),verbose=FALSE)


# Plot Posterior Checks


# Model 1
pdf(file = here('figures','temporary','posterior_predictive_check_model1_spd_euc.pdf'),width = 7,height = 8)
par(mfrow=c(2,2),mar=c(5,4,2,1))
plotPPCheckSPD(obs.norm,postpcheck.model1.spd[[1]],main='Calsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model1.spd[[2]],main='Calsample - non-normalised')
plotPPCheckSPD(obs.norm,postpcheck.model1.spd[[3]],main='Uncalsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model1.spd[[4]],main='Uncalsample - non-normalised')
dev.off()

pdf(file = here('figures','temporary','posterior_predictive_check_model1_spd_nrmse.pdf'),width = 7,height = 8)
par(mfrow=c(2,2),mar=c(5,4,2,1))
plotPPCheckSPD(obs.norm,postpcheck.model1.spd[[5]],main='Calsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model1.spd[[6]],main='Calsample - non-normalised')
plotPPCheckSPD(obs.norm,postpcheck.model1.spd[[7]],main='Uncalsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model1.spd[[8]],main='Uncalsample - non-normalised')
dev.off()


# Model 2
pdf(file = here('figures','temporary','posterior_predictive_check_model2_spd_euc.pdf'),width = 7,height = 8)
par(mfrow=c(2,2),mar=c(5,4,2,1))
plotPPCheckSPD(obs.norm,postpcheck.model2.spd[[1]],main='Calsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model2.spd[[2]],main='Calsample - non-normalised')
plotPPCheckSPD(obs.norm,postpcheck.model2.spd[[3]],main='Uncalsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model2.spd[[4]],main='Uncalsample - non-normalised')
dev.off()

pdf(file = here('figures','temporary','posterior_predictive_check_model2_spd_nrmse.pdf'),width = 7,height = 8)
par(mfrow=c(2,2),mar=c(5,4,2,1))
plotPPCheckSPD(obs.norm,postpcheck.model2.spd[[5]],main='Calsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model2.spd[[6]],main='Calsample - non-normalised')
plotPPCheckSPD(obs.norm,postpcheck.model2.spd[[7]],main='Uncalsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model2.spd[[8]],main='Uncalsample - non-normalised')
dev.off()


# Model 3
pdf(file = here('figures','temporary','posterior_predictive_check_model3_spd_euc.pdf'),width = 7,height = 8)
par(mfrow=c(2,2),mar=c(5,4,2,1))
plotPPCheckSPD(obs.norm,postpcheck.model3.spd[[1]],main='Calsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model3.spd[[2]],main='Calsample - non-normalised')
plotPPCheckSPD(obs.norm,postpcheck.model3.spd[[3]],main='Uncalsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model3.spd[[4]],main='Uncalsample - non-normalised')
dev.off()

pdf(file = here('figures','temporary','posterior_predictive_check_model3_spd_nrmse.pdf'),width = 7,height = 8)
par(mfrow=c(2,2),mar=c(5,4,2,1))
plotPPCheckSPD(obs.norm,postpcheck.model3.spd[[5]],main='Calsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model3.spd[[6]],main='Calsample - non-normalised')
plotPPCheckSPD(obs.norm,postpcheck.model3.spd[[7]],main='Uncalsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model3.spd[[8]],main='Uncalsample - non-normalised')
dev.off()


# Model 4
pdf(file = here('figures','temporary','posterior_predictive_check_model4_spd_euc.pdf'),width = 7,height = 8)
par(mfrow=c(2,2),mar=c(5,4,2,1))
plotPPCheckSPD(obs.norm,postpcheck.model4.spd[[1]],main='Calsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model4.spd[[2]],main='Calsample - non-normalised')
plotPPCheckSPD(obs.norm,postpcheck.model4.spd[[3]],main='Uncalsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model4.spd[[4]],main='Uncalsample - non-normalised')
dev.off()

pdf(file = here('figures','temporary','posterior_predictive_check_model4_spd_nrmse.pdf'),width = 7,height = 8)
par(mfrow=c(2,2),mar=c(5,4,2,1))
plotPPCheckSPD(obs.norm,postpcheck.model4.spd[[5]],main='Calsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model4.spd[[6]],main='Calsample - non-normalised')
plotPPCheckSPD(obs.norm,postpcheck.model4.spd[[7]],main='Uncalsample - normalised')
plotPPCheckSPD(obs.nnorm,postpcheck.model4.spd[[8]],main='Uncalsample - non-normalised')
dev.off()