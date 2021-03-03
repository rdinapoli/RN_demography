## A robust approach for modeling human demography and paleoenvironmental interactions using Approximate Bayesian Computation: an example from Rapa Nui (Easter Island): source code, data, scripts, and Shiny app
 
This repository contains data and scripts to reproduce analyses in "A robust approach for evaluating models of climate-induced population dynamics using radiocarbon-based paleodemography and Approximate Bayesian Computation: an example from Rapa Nui (Easter Island)"

The main workflow is described below and outputs are stored as R image files in the R_images(add link) directory and figures in the figures(add link) directory.

## Data Sets and Data Preparation

All data files used in the paper can be found in the [data](./data) directory, which contains: 
rapanui_DiNapoli_etal.csv is the current radiocarbon data for Rapa Nui with corrected site codes and classification of dates by archaeological context.
RAPANUINW.csv is supplementary data from Lima et al. (2020) used to fit their nls models, incluiding time series data on palm cover and SOI index.
rapnveget.csv is supplementary data from Lima et al. (2020) used to fit their linear regression models.
date_si_lima.csv is the supplementary radiocarbon data used by Lima et al. (2020).
tableS5.csv is supplementary information from Lima et al. (2020).
yan2011soipr.txt is the Southern Oscillation Index (SOI) reconstruction from Yan et al. (2011).

The file data_prep.R in the [runscripts](./runscripts) directory contains script for calibrating and binning the radiocarbon dates, generating the observed SPDs, and paleoenvironmental time series (extracting the SOI values from Yan et al. (2011), and performing a linear interpolation of palm values from Lima et al. (2020)). The output of this script is saved in the variables.Rdata file in the [R_imagefiles](./R_imagefiles) directory. Note that data preparation uses a custom calibration function, calibrate2.R in the [src](./src) directory that allows the use of custom mixed curves.

## Paleodemographic modeling using Approximate Bayesian Computation (ABC)

The scripts for defining the functions for the demographic models, running the models, and performing multi-model comparison can be found in the [src](./src) directory.
The file growthModel.R is a function defining logistic growth models 1-4. The file randomThin.R is a function for generating the target SPDs in the ABC routine using a random thinning method. The file modelRunner.R is a function for executing the ABC model-fitting simulations and returns a data.frame of distance measures quantifying the fit between each realization of the model and the observed SPD. The file modelSelection.R contains a function for comparing ABC results using Bayes factors. The file plotFunctions.R contains functions to plot the marginal posteriors of the model parameters and SPD posterior predictive checks.

The [runscripts](./runscripts) folder contains the scripts to run the above functions using the Rapa Nui data and reproduce the results in the main text and supplementary information. The file priorCheck.R runs the prior predictive checks for the 4 growth models and reproduces Figure S2. The files submit_ABC_model1.R through submit_ABC_model4.R execute the ABC routine for each of the growth models. Note that **each of these take a long time to complete**. The current code is written for parallel runs on 10 cores, which took approximately 24 hours for each model. Users should be aware of these settings in attempting to reproduce or expand on the code provided here. The results of each model run are saved in the [R_imagefiles](./R_imagefiles) directory as abc_model1.Rdata, etc. The file process_ABC_results.R contains code to extract the 1000 best fitting runs of each model, and the results are saved as ABC_posteriors.RData in the R_imagefiles directory. The file model_selection_and_diagnostics.R loads these ABC posteriors and calculates the Bayes factor matrices for model comparison, and produces the results shown in the supplementary tables and table 2 in the main text. The file posteriorCheck_raw.R generates plots of the posterior predictive checks for each model iteration shown in the supplementary information. The file posteriorPlot.R uses the ABC_posteriors.RData file and the plotFunctions.R script to generate the marginal posterior plots for each model's parameters in the supplementary information and main text. The file posteriorCheck_spd.R uses these ABC results and generates posterior predictive checks for the fit between the observed and simulated SPDs, and plot_posteriorCheck_spd.R generated plots of these as shown in the supplementary information and main text. The file main_figures.R generates Figures 3-4 in the main text. The file CorrelTest.R explores the direct correlation between the observed SPD and paleoenvironmental proxies and generates supplementary Figure S1 and Figure 2 in the main text.

In summary, the workflow to reproduce these analyses requires first running data_prep.R, then running priorCheck.R and submit_ABC_model1.R through submit_ABC_model4.R. Using the output, the user then runs process_ABC_results.R, model_selection_and_diagnostics.R, posteriorCheck_raw.R, posteriorPlot.R, posteriorCheck_spd.R, plot_posteriorCheck_spd.R, CorrelTest.R, and main_figures.R 

## Shiny app

An interactive Shiny app is available for users to explore the demographic model results under different parameterizations. This app can be executed within Rstudio using the file in the [growthmodel_shiny](./growthmodel_shiny) directory or accessed online here: https://dinapoli.shinyapps.io/RN_Demographic_Growth_Model/

## File Structure

## R Settings

## License
CC-BY 3.0
