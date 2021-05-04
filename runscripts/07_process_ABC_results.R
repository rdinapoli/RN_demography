# Load Libraries ####
library(here)

# Load Raw Results ####
load(here('R_imagefiles','abc_model1.RData'))
load(here('R_imagefiles','abc_model2.RData'))
load(here('R_imagefiles','abc_model3.RData'))
load(here('R_imagefiles','abc_model4.RData'))

# ABC settings ####
tol=0.004 #1000 best models
nsim = 250000

# Extract Posteriors for each model ####

# model 1
model1.posterior = vector('list',length=8)
names(model1.posterior) = c('euc.cal.norm','euc.cal.nnorm','euc.uncal.norm','euc.uncal.nnorm','nrmse.cal.norm','nrmse.cal.nnorm','nrmse.uncal.norm','nrmse.uncal.nnorm')
model1.posterior[[1]] = model1[order(model1$euc.cal.norm)[1:(nsim*tol)],]
model1.posterior[[2]] = model1[order(model1$euc.cal.nnorm)[1:(nsim*tol)],]
model1.posterior[[3]] = model1[order(model1$euc.uncal.norm)[1:(nsim*tol)],]
model1.posterior[[4]] = model1[order(model1$euc.uncal.nnorm)[1:(nsim*tol)],]
model1.posterior[[5]] = model1[order(model1$nrmse.cal.norm)[1:(nsim*tol)],]
model1.posterior[[6]] = model1[order(model1$nrmse.cal.nnorm)[1:(nsim*tol)],]
model1.posterior[[7]] = model1[order(model1$nrmse.uncal.norm)[1:(nsim*tol)],]
model1.posterior[[8]] = model1[order(model1$nrmse.uncal.nnorm)[1:(nsim*tol)],]

# # model 2
model2.posterior = vector('list',length=8)
names(model2.posterior) = c('euc.cal.norm','euc.cal.nnorm','euc.uncal.norm','euc.uncal.nnorm','nrmse.cal.norm','nrmse.cal.nnorm','nrmse.uncal.norm','nrmse.uncal.nnorm')
model2.posterior[[1]] = model2[order(model2$euc.cal.norm)[1:(nsim*tol)],]
model2.posterior[[2]] = model2[order(model2$euc.cal.nnorm)[1:(nsim*tol)],]
model2.posterior[[3]] = model2[order(model2$euc.uncal.norm)[1:(nsim*tol)],]
model2.posterior[[4]] = model2[order(model2$euc.uncal.nnorm)[1:(nsim*tol)],]
model2.posterior[[5]] = model2[order(model2$nrmse.cal.norm)[1:(nsim*tol)],]
model2.posterior[[6]] = model2[order(model2$nrmse.cal.nnorm)[1:(nsim*tol)],]
model2.posterior[[7]] = model2[order(model2$nrmse.uncal.norm)[1:(nsim*tol)],]
model2.posterior[[8]] = model2[order(model2$nrmse.uncal.nnorm)[1:(nsim*tol)],]

# model 3
model3.posterior = vector('list',length=8)
names(model3.posterior) = c('euc.cal.norm','euc.cal.nnorm','euc.uncal.norm','euc.uncal.nnorm','nrmse.cal.norm','nrmse.cal.nnorm','nrmse.uncal.norm','nrmse.uncal.nnorm')
model3.posterior[[1]] = model3[order(model3$euc.cal.norm)[1:(nsim*tol)],]
model3.posterior[[2]] = model3[order(model3$euc.cal.nnorm)[1:(nsim*tol)],]
model3.posterior[[3]] = model3[order(model3$euc.uncal.norm)[1:(nsim*tol)],]
model3.posterior[[4]] = model3[order(model3$euc.uncal.nnorm)[1:(nsim*tol)],]
model3.posterior[[5]] = model3[order(model3$nrmse.cal.norm)[1:(nsim*tol)],]
model3.posterior[[6]] = model3[order(model3$nrmse.cal.nnorm)[1:(nsim*tol)],]
model3.posterior[[7]] = model3[order(model3$nrmse.uncal.norm)[1:(nsim*tol)],]
model3.posterior[[8]] = model3[order(model3$nrmse.uncal.nnorm)[1:(nsim*tol)],]

# model 4
model4.posterior = vector('list',length=8)
names(model4.posterior) = c('euc.cal.norm','euc.cal.nnorm','euc.uncal.norm','euc.uncal.nnorm','nrmse.cal.norm','nrmse.cal.nnorm','nrmse.uncal.norm','nrmse.uncal.nnorm')
model4.posterior[[1]] = model4[order(model4$euc.cal.norm)[1:(nsim*tol)],]
model4.posterior[[2]] = model4[order(model4$euc.cal.nnorm)[1:(nsim*tol)],]
model4.posterior[[3]] = model4[order(model4$euc.uncal.norm)[1:(nsim*tol)],]
model4.posterior[[4]] = model4[order(model4$euc.uncal.nnorm)[1:(nsim*tol)],]
model4.posterior[[5]] = model4[order(model4$nrmse.cal.norm)[1:(nsim*tol)],]
model4.posterior[[6]] = model4[order(model4$nrmse.cal.nnorm)[1:(nsim*tol)],]
model4.posterior[[7]] = model4[order(model4$nrmse.uncal.norm)[1:(nsim*tol)],]
model4.posterior[[8]] = model4[order(model4$nrmse.uncal.nnorm)[1:(nsim*tol)],]




# Model Selection ####
model1$model='m1'
model2$model='m2'
model3$model='m3'
model4$model='m4'
 
## Combine models
allmodels = rbind.data.frame(model1,model2,model3,model4)
 
# Extract Posterior with tol=0.004, the 1000 best models
tol=0.004
allmodels.posterior = vector('list',length=8)
nims=nrow(allmodels)
names(allmodels.posterior) = c('euc.cal.norm','euc.cal.nnorm','euc.uncal.norm','euc.uncal.nnorm','nrmse.cal.norm','nrmse.cal.nnorm','nrmse.uncal.norm','nrmse.uncal.nnorm')
allmodels.posterior[[1]] = allmodels[order(allmodels$euc.cal.norm)[1:(nsim*tol)],]
allmodels.posterior[[2]] = allmodels[order(allmodels$euc.cal.nnorm)[1:(nsim*tol)],]
allmodels.posterior[[3]] = allmodels[order(allmodels$euc.uncal.norm)[1:(nsim*tol)],]
allmodels.posterior[[4]] = allmodels[order(allmodels$euc.uncal.nnorm)[1:(nsim*tol)],]
allmodels.posterior[[5]] = allmodels[order(allmodels$nrmse.cal.norm)[1:(nsim*tol)],]
allmodels.posterior[[6]] = allmodels[order(allmodels$nrmse.cal.nnorm)[1:(nsim*tol)],]
allmodels.posterior[[7]] = allmodels[order(allmodels$nrmse.uncal.norm)[1:(nsim*tol)],]
allmodels.posterior[[8]] = allmodels[order(allmodels$nrmse.uncal.nnorm)[1:(nsim*tol)],]

save(model1.posterior,model2.posterior,model3.posterior,model4.posterior,allmodels.posterior,file=here('R_imagefiles','ABC_posteriors.RData'))

