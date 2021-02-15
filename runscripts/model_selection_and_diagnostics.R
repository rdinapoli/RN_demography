library(here)
load(here('R_imagefiles','ABC_posteriors.RData'))
source(here('src','modelSelection.R'))
source(here('src','tolDiagPlot.R'))

# Compute Bayes Factor Matrix
modelSelection(allmodels.posterior[[1]]) #euc.cal.norm
modelSelection(allmodels.posterior[[2]]) #euc.cal.nnorm
modelSelection(allmodels.posterior[[3]]) #euc.uncal.norm
modelSelection(allmodels.posterior[[4]]) #euc.uncal.nnorm, results presented in main text
modelSelection(allmodels.posterior[[5]]) #nrmse.cal.norm
modelSelection(allmodels.posterior[[6]]) #nrmse.cal.nnorm
modelSelection(allmodels.posterior[[7]]) #nrmse.uncal.norm
modelSelection(allmodels.posterior[[8]]) #nrmse.uncal.nnorm

# Sample Diagnostic Plots ####

# Model Selection
tolDiagPlot(x=allmodels.posterior[[1]],param = 'model',tolcol='euc.cal.norm')

# Posterior
tolDiagPlot(model2.posterior[[1]],param='b1',tolcol='euc.cal.norm')
tolDiagPlot(model4.posterior[[1]],param='b1',tolcol='euc.cal.norm')
tolDiagPlot(model4.posterior[[1]],param='b2',tolcol='euc.cal.norm')
tolDiagPlot(model1.posterior[[1]],param='r',tolcol='euc.cal.norm')
