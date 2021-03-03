## Approximate Bayesian Computation of radiocarbon and paleoenvironmental record suggests resilience of Rapa Nui (Easter Island) communities to ecological change: source code, data, scripts, and Shiny app
 
This repository contains data and scripts to reproduce analyses in "Approximate Bayesian Computation of radiocarbon and paleoenvironmental record suggests resilience of Rapa Nui (Easter Island) communities to ecological change"

The main workflow is described below and outputs are stored as R image files in the [R_imagefiles](./R_imagefiles) directory and figures in the [figures](./figures) directory.

## Data Sets and Data Preparation

All data files used in the paper can be found in the [data](./data) directory, which contains: 
 - `rapanui_DiNapoli_etal.csv` is the current radiocarbon data for Rapa Nui with corrected site codes and classification of dates by archaeological context.
 - `RAPANUINW.csv` is supplementary data from Lima et al. (2020) used to fit their nls models, incluiding time series data on palm cover and SOI index.
 - `rapnveget.csv` is supplementary data from Lima et al. (2020) used to fit their linear regression models.
 - `date_si_lima.csv` is the supplementary radiocarbon data used by Lima et al. (2020).
 - `tableS5.csv` is supplementary information from Lima et al. (2020).
 - `yan2011soipr.txt` is the Southern Oscillation Index (SOI) reconstruction from Yan et al. (2011).

The file `data_prep.R` in the [runscripts](./runscripts) directory contains script for calibrating and binning the radiocarbon dates, generating the observed SPDs, and paleoenvironmental time series (extracting the SOI values from Yan et al. (2011), and performing a linear interpolation of palm values from Lima et al. (2020)). The output of this script is saved in the `variables.Rdata` file in the [R_imagefiles](./R_imagefiles) directory. Note that data preparation uses a custom fast calibration function, `calibrate2.R` in the [src](./src) directory that allows the use of custom mixed curves.

## Paleodemographic modeling using Approximate Bayesian Computation (ABC)

The scripts for defining the functions for the demographic models, running the models, and performing multi-model comparison can be found in the [src](./src) directory.
The file `growthModel.R` is a function defining logistic growth models 1-4. The file `randomThin.R` is a function for generating the target SPDs in the ABC routine using a random thinning method. The file `modelRunner.R` is a function for executing the ABC model-fitting simulations and returns a data.frame of distance measures quantifying the fit between each realization of the model and the observed SPD. The file `modelSelection.R` contains a function for comparing ABC results using approximate Bayes factors. The file `plotFunctions.R` contains functions to plot the marginal posteriors of the model parameters and SPD posterior predictive checks.

The [runscripts](./runscripts) folder contains the scripts to run the above functions using the Rapa Nui data and reproduce the results in the main text and supplementary information. The file `priorCheck.R` runs the prior predictive checks for the four growth models and reproduces Figure S2. The files `submit_ABC_model1.R` through `submit_ABC_model4.R` execute the ABC routine for each of the growth models. Note that **each of these take a long time to complete**. The current code is written for parallel runs on 10 cores, which took approximately 24 hours for each model. Users should be aware of these settings in attempting to reproduce or expand on the code provided here. The results of each model run are saved in the [R_imagefiles](./R_imagefiles) directory as `abc_model1.Rdata`, etc. The file `process_ABC_results.R` contains code to extract the 1000 best fitting runs of each model, and the results are saved as `ABC_posteriors.RData` in the [R_imagefiles](./R_imagefiles) directory. The file `model_selection_and_diagnostics.R` loads these ABC posteriors and calculates the approximate Bayes factor matrices for model comparison, and produces the results shown in the supplementary tables and table 2 in the main text. The file `posteriorCheck_raw.R` generates plots of the posterior predictive checks for each model iteration shown in the supplementary information. The file `posteriorPlot.R` uses the `ABC_posteriors.RData` file and the `plotFunctions.R` script to generate the marginal posterior plots for each model's parameters in the supplementary information and main text. The file `posteriorCheck_spd.R` uses these ABC results and generates posterior predictive checks for the fit between the observed and simulated SPDs, and `plot_posteriorCheck_spd.R` generated plots of these as shown in the supplementary information and main text. The file `main_figures.R` generates Figures 3-4 in the main text. The file `CorrelTest.R` explores the direct correlation between the observed SPD and paleoenvironmental proxies and generates supplementary Figure S1 and Figure 2 in the main text.

In summary, the workflow to reproduce these analyses requires first running `data_prep.R`, then running `priorCheck.R` and `submit_ABC_model1.R` through `submit_ABC_model4.R`. Using the output, the user then runs `process_ABC_results.R`, `model_selection_and_diagnostics.R`, `posteriorCheck_raw.R`, `posteriorPlot.R`, `posteriorCheck_spd.R`, `plot_posteriorCheck_spd.R`, `CorrelTest.R`, and `main_figures.R`. 

## Shiny app

An interactive Shiny app is available for users to explore the demographic model results under different parameterizations. This app can be executed within Rstudio using the file in the [growthmodel_shiny](./growthmodel_shiny) directory or accessed online here: https://dinapoli.shinyapps.io/RN_Demographic_Growth_Model/

## File Structure
```
├── data
│   ├── date_si_lima.csv
│   ├── rapanui_DiNapoli_etal.csv
│   ├── RAPANUINW.xlsx
│   ├── rapnveget.xlsx
│   ├── tableS5.csv
│   └── yan2011soipr.txt
├── figures
│   ├── main_paper
│   │   ├── Figure2_obs_SPD.jpeg
│   │   ├── Figure2_obs_SPD.tiff
│   │   ├── figure3_prdens.pdf
│   │   └── figure4.pdf
│   ├── supplemental
│   │   ├── correlation_results.jpeg
│   │   ├── correlation_results.tiff
│   │   ├── m1_jointpost_euc.cal.nnorm.jpeg
│   │   ├── m1_jointpost_euc.cal.norm.jpeg
│   │   ├── m1_jointpost_euc.uncal.nnorm.jpeg
│   │   ├── m1_jointpost_euc.uncal.norm.jpeg
│   │   ├── m1_jointpost_nrmse.cal.nnorm.jpeg
│   │   ├── m1_jointpost_nrmse.cal.norm.jpeg
│   │   ├── m1_jointpost_nrmse.uncal.nnorm.jpeg
│   │   ├── m1_jointpost_nrmse.uncal.norm.jpeg
│   │   ├── m1_post_euc.cal.nnorm.pdf
│   │   ├── m1_post_euc.cal.norm.pdf
│   │   ├── m1_post_euc.uncal.nnorm.pdf
│   │   ├── m1_post_euc.uncal.norm.pdf
│   │   ├── m1_post_nrmse.cal.nnorm.pdf
│   │   ├── m1_post_nrmse.cal.norm.pdf
│   │   ├── m1_post_nrmse.uncal.nnorm.pdf
│   │   ├── m1_post_nrmse.uncal.norm.pdf
│   │   ├── m2_jointpost_euc.cal.nnorm.jpeg
│   │   ├── m2_jointpost_euc.cal.norm.jpeg
│   │   ├── m2_jointpost_euc.uncal.nnorm.jpeg
│   │   ├── m2_jointpost_euc.uncal.norm.jpeg
│   │   ├── m2_jointpost_nrmse.cal.nnorm.jpeg
│   │   ├── m2_jointpost_nrmse.cal.norm.jpeg
│   │   ├── m2_jointpost_nrmse.uncal.nnorm.jpeg
│   │   ├── m2_jointpost_nrmse.uncal.norm.jpeg
│   │   ├── m2_post_euc.cal.nnorm.pdf
│   │   ├── m2_post_euc.cal.norm.pdf
│   │   ├── m2_post_euc.uncal.nnorm.pdf
│   │   ├── m2_post_euc.uncal.norm.pdf
│   │   ├── m2_post_nrmse.cal.nnorm.pdf
│   │   ├── m2_post_nrmse.cal.norm.pdf
│   │   ├── m2_post_nrmse.uncal.nnorm.pdf
│   │   ├── m2_post_nrmse.uncal.norm.pdf
│   │   ├── m3_jointpost_euc.cal.nnorm.jpeg
│   │   ├── m3_jointpost_euc.cal.norm.jpeg
│   │   ├── m3_jointpost_euc.uncal.nnorm.jpeg
│   │   ├── m3_jointpost_euc.uncal.norm.jpeg
│   │   ├── m3_jointpost_nrmse.cal.nnorm.jpeg
│   │   ├── m3_jointpost_nrmse.cal.norm.jpeg
│   │   ├── m3_jointpost_nrmse.uncal.nnorm.jpeg
│   │   ├── m3_jointpost_nrmse.uncal.norm.jpeg
│   │   ├── m3_post_euc.cal.nnorm.pdf
│   │   ├── m3_post_euc.cal.norm.pdf
│   │   ├── m3_post_euc.uncal.nnorm.pdf
│   │   ├── m3_post_euc.uncal.norm.pdf
│   │   ├── m3_post_nrmse.cal.nnorm.pdf
│   │   ├── m3_post_nrmse.cal.norm.pdf
│   │   ├── m3_post_nrmse.uncal.nnorm.pdf
│   │   ├── m3_post_nrmse.uncal.norm.pdf
│   │   ├── m4_jointpost_euc.cal.nnorm.jpeg
│   │   ├── m4_jointpost_euc.cal.norm.jpeg
│   │   ├── m4_jointpost_euc.uncal.nnorm.jpeg
│   │   ├── m4_jointpost_euc.uncal.norm.jpeg
│   │   ├── m4_jointpost_nrmse.cal.nnorm.jpeg
│   │   ├── m4_jointpost_nrmse.cal.norm.jpeg
│   │   ├── m4_jointpost_nrmse.uncal.nnorm.jpeg
│   │   ├── m4_jointpost_nrmse.uncal.norm.jpeg
│   │   ├── m4_post_euc.cal.nnorm.pdf
│   │   ├── m4_post_euc.cal.norm.pdf
│   │   ├── m4_post_euc.uncal.nnorm.pdf
│   │   ├── m4_post_euc.uncal.norm.pdf
│   │   ├── m4_post_nrmse.cal.nnorm.pdf
│   │   ├── m4_post_nrmse.cal.norm.pdf
│   │   ├── m4_post_nrmse.uncal.nnorm.pdf
│   │   ├── m4_post_nrmse.uncal.norm.pdf
│   │   ├── posterior_predictive_check_model1_euc.jpeg
│   │   ├── posterior_predictive_check_model1_nrmse.jpeg
│   │   ├── posterior_predictive_check_model1_spd_euc.pdf
│   │   ├── posterior_predictive_check_model1_spd_nrmse.pdf
│   │   ├── posterior_predictive_check_model2_euc.jpeg
│   │   ├── posterior_predictive_check_model2_nrmse.jpeg
│   │   ├── posterior_predictive_check_model2_spd_euc.pdf
│   │   ├── posterior_predictive_check_model2_spd_nrmse.pdf
│   │   ├── posterior_predictive_check_model3_euc.jpeg
│   │   ├── posterior_predictive_check_model3_nrmse.jpeg
│   │   ├── posterior_predictive_check_model3_spd_euc.pdf
│   │   ├── posterior_predictive_check_model3_spd_nrmse.pdf
│   │   ├── posterior_predictive_check_model4_euc.jpeg
│   │   ├── posterior_predictive_check_model4_nrmse.jpeg
│   │   ├── posterior_predictive_check_model4_spd_euc.pdf
│   │   ├── posterior_predictive_check_model4_spd_nrmse.pdf
│   │   ├── prior_predictive_check_N.jpeg
│   │   └── prior_predictive_check_PrDens.jpeg
│   └── supplementary
├── growthmodel_shiny
│   └── growthModelShiny.R
├── README.md
├── R_imagefiles
│   ├── abc_model1.RData
│   ├── abc_model2.RData
│   ├── abc_model3.RData
│   ├── abc_model4.RData
│   ├── ABC_posteriors.RData
│   ├── post_pred_spd.RData
│   └── variables.RData
├── RN_demography.Rproj
├── runscripts
│   ├── CorrelTest.R
│   ├── data_prep.R
│   ├── main_figures.R
│   ├── model_selection_and_diagnostics.R
│   ├── plot_posteriorCheck_spd.R
│   ├── posteriorCheck_raw.R
│   ├── posteriorCheck_spd.R
│   ├── posteriorPlot.R
│   ├── priorCheck.R
│   ├── process_ABC_results.R
│   ├── submit_ABC_model1.R
│   ├── submit_ABC_model2.R
│   ├── submit_ABC_model3.R
│   └── submit_ABC_model4.R
└── src
    ├── calibrate2.R
    ├── growthModel.R
    ├── modelRunner.R
    ├── modelSelection.R
    ├── plotFunctions.R
    ├── randomThin.R
    └── tolDiagPlot.R
```
## R Settings

```
attached base packages:
[1] parallel  stats     graphics  grDevices utils    
[6] methods   base     

other attached packages:
 [1] shiny_1.5.0        coda_0.19-4        doSNOW_1.0.19     
 [4] snow_0.4-3         doParallel_1.0.15  iterators_1.0.13  
 [7] foreach_1.5.1      ReIns_1.0.10       RColorBrewer_1.1-2
[10] forcats_0.5.0      stringr_1.4.0      dplyr_1.0.2       
[13] purrr_0.3.4        readr_1.4.0        tidyr_1.1.2       
[16] tibble_3.0.6       tidyverse_1.3.0    latex2exp_0.4.0   
[19] here_0.1           gridExtra_2.3      ggplot2_3.3.2     
[22] rcarbon_1.4.1     

loaded via a namespace (and not attached):
 [1] nlme_3.1-148         fs_1.4.2            
 [3] lubridate_1.7.9.2    httr_1.4.2          
 [5] rprojroot_2.0.2      tools_4.0.3         
 [7] backports_1.2.1      R6_2.5.0            
 [9] rpart_4.1-15         DBI_1.1.0           
[11] mgcv_1.8-31          colorspace_2.0-0    
[13] withr_2.3.0          sp_1.4-5            
[15] tidyselect_1.1.0     compiler_4.0.3      
[17] cli_2.3.0            rvest_0.3.6         
[19] xml2_1.3.2           scales_1.1.1        
[21] spatstat.data_2.0-0  spatstat_1.64-1     
[23] goftest_1.2-2        digest_0.6.27       
[25] spatstat.utils_2.0-0 rmarkdown_2.6       
[27] pkgconfig_2.0.3      htmltools_0.5.0     
[29] fastmap_1.0.1        dbplyr_1.4.4        
[31] rlang_0.4.10         readxl_1.3.1        
[33] rstudioapi_0.13      generics_0.1.0      
[35] startup_0.14.1       jsonlite_1.7.2      
[37] magrittr_2.0.1       Matrix_1.2-18       
[39] Rcpp_1.0.5           munsell_0.5.0       
[41] abind_1.4-5          lifecycle_1.0.0     
[43] stringi_1.5.3        yaml_2.2.1          
[45] grid_4.0.3           blob_1.2.1          
[47] promises_1.1.1       crayon_1.4.1        
[49] deldir_0.2-10        lattice_0.20-41     
[51] haven_2.3.1          splines_4.0.3       
[53] tensor_1.5           hms_0.5.3           
[55] knitr_1.31           pillar_1.4.7        
[57] codetools_0.2-16     reprex_0.3.0        
[59] glue_1.4.2           evaluate_0.14       
[61] modelr_0.1.8         httpuv_1.5.4        
[63] vctrs_0.3.6          cellranger_1.1.0    
[65] gtable_0.3.0         polyclip_1.10-0     
[67] assertthat_0.2.1     xfun_0.21           
[69] mime_0.10            xtable_1.8-4        
[71] broom_0.7.0          later_1.1.0.1       
[73] survival_3.2-3       ellipsis_0.3.1 
```
## License
CC-BY 3.0
