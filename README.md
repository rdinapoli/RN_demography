# RN_demography
This repository contains data and code to reproduce analyses in "Evaluating models of climate-induced population dynamics using radiocarbon-based paleodemography and Approximate Bayesian Computation: an example from Rapa Nui (Easter Island)"

# Notes on files/folders

## Data

All data files are contained in the data folder, which contains: 

rapanui_DiNapoli_etal.csv is the radiocarbon data with corrected site codes and classification of dates by archaeological context

RAPANUINW.csv is supplementary data from Lima et al. (2020) used to fit their nls models, incluiding time series data on palm cover and SOI index

rapnveget.csv is supplementary data from Lima et al. (2020) used to fit their linear regression models

date_si_lima.csv is the supplementary radiocarbon data from Lima et al. (2020)

yan2011soipr.txt is the Southern Oscillation Index (SOI) reconstruction from Yan et al. (2011).

## src

This folder contains custom R functions for running the ABC analysis:

calibrate2.R is a modification to the rcarbon calibrate() function that allows for...

growthModel.R is a function for running logistic growth models 1-4.

modelRunner.R is a function for executing the approximate Bayesian computation (ABC) simulations and calculation of distance measures.

modelSelection.R is a function for comparing ABC results using Bayes factors.

plotFunctions.R is used to reproduce figured in the main text and supplementary information.

randomThin.R is a function for generating the target SPDs in the ABC routine using a random thinning method.

tolDiagPlot.R is a function for comparing how the Bayes factors may change with different tolerance values.

## runscripts

This folder conaints the code to reproduce the results in the main text and supplementary information:

priorCheck.R runs the prior predictive checks for the 4 growth models.

data_prep.R Prepares the data used in the analyses, including calibrating and binning the radiocarbon dates, extracting the SOI values from Yan et al. (2011), and performing a linear interpolation of palm values from Lima et al. (2020).

submit_ABC_model1.R through model4.R execute the ABC routine for each of the growth models. Note that each of these take a very long time to run.

## Shiny app

An interactive Shiny app is available for users to explore the demographic model results under different parameterizations. This app can be executed within Rstudio or accessed online here: https://dinapoli.shinyapps.io/RN_Demographic_Growth_Model/
