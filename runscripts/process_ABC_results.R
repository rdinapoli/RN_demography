# Load Libraries ####
library(here)

# Load Raw Results ####
load(here('R_imagefiles','abc_model4.RData'))
model4 = result
# load(here('R_imagefiles','abc_model3.RData'))
# model3 = result
# load(here('R_imagefiles','abc_model2.RData'))
# model2 = result
# load(here('R_imagefiles','abc_model1.RData'))
# model1 = result

# ABC settings ####
tol=0.004
nsim = 250000

# Extract Posteriors ####
model4.posterior = vector('list',length=8)
names(model4.posterior) = c('euc.cal.norm','euc.cal.nnorm','euc.uncal.norm','euc.uncal.nnorm','nrmse.cal.norm','nrmse.cal.nnorm','nrmse.uncal.norm','nrmse.uncal.nnorm')
model4.posterior[[1]] = result[order(model4$euc.cal.norm)[1:(nsim*tol)],]
model4.posterior[[2]] = result[order(model4$euc.cal.nnorm)[1:(nsim*tol)],]
model4.posterior[[3]] = result[order(model4$euc.uncal.norm)[1:(nsim*tol)],]
model4.posterior[[4]] = result[order(model4$euc.uncal.nnorm)[1:(nsim*tol)],]
model4.posterior[[5]] = result[order(model4$nrmse.cal.norm)[1:(nsim*tol)],]
model4.posterior[[6]] = result[order(model4$nrmse.cal.nnorm)[1:(nsim*tol)],]
model4.posterior[[7]] = result[order(model4$nrmse.uncal.norm)[1:(nsim*tol)],]
model4.posterior[[8]] = result[order(model4$nrmse.uncal.nnorm)[1:(nsim*tol)],]

# Save to R image file
save(model4.posterior,file=here('R_imagefiles','ABC_posteriors.RData'))
