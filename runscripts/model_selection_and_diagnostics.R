library(here)
load(here('R_imagefiles','ABC_posteriors.RData'))
source(here('src','modelSelection.R'))
source(here('src','tolDiagPlot.R'))

# Compute Bayes Factor Matrix
modelSelection(allmodels.posterior[[1]])
modelSelection(allmodels.posterior[[2]])
modelSelection(allmodels.posterior[[3]])
modelSelection(allmodels.posterior[[4]]) #euc.uncal.nnorm, results presented in main text
modelSelection(allmodels.posterior[[5]])
modelSelection(allmodels.posterior[[6]])
modelSelection(allmodels.posterior[[7]])
modelSelection(allmodels.posterior[[8]])

# Sample Diagnostic Plots ####

# Model Selection
tolDiagPlot(x=allmodels.posterior[[1]],param = 'model',tolcol='euc.cal.norm')

# Posterior
tolDiagPlot(model2.posterior[[1]],param='b1',tolcol='euc.cal.norm')
tolDiagPlot(model4.posterior[[1]],param='b1',tolcol='euc.cal.norm')
tolDiagPlot(model4.posterior[[1]],param='b2',tolcol='euc.cal.norm')
tolDiagPlot(model1.posterior[[1]],param='r',tolcol='euc.cal.norm')
