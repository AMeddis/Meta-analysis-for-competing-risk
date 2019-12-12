Code and example data for manuscript "Meta-analysis of clinical trials with competing time-to-event endpoints" by Alessandra Meddis, Aurélien Latouche, Bingqing Zhou, Stefan Michiels, Jason Fine.  


For questions, comments or remarks about the code please contact A. Meddis (alessandra.meddis@curie.fr / alessandra.meddis@gmail.com).


The code was written using :
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)

attached base packages:
grid, stats, graphics, grDevices, utils, datasets, methods, base     

packages useful for the code:
etm_1.0.5  ggplot2_3.2.1   dplyr_0.8.3  crrSC_1.1   meta_4.9-7      

cmprsk_2.2-9    survival_2.41-3  prodlim_2019.11.13

(*** If you don't succeed to install this last prodlim version, use: 
 library(devtools)
 install_github("tagteam/prodlim")
 packageVersion("prodlim")   ***)


Data is confidential, thus we provide a recreated data set (compR.Rdata) similar to the one used for the application, with same number of trials for each treatment group. Informations on trial id (trial), group of treatment (grpCT), time of event (failure_time), type of event (failure_type), treatment information (treated/control) and age (continuous variable) are introduced in the data set.


To reproduce the analysis described in the manuscript, just run the Meta-analysis_for_compR.R script. More in details, the code is used to: 

   -draw forest plots
   
   -fit a stratified Fine-Gray model for subdistribution hazards
   
   -perform landmark approaches for subdistribution hazards
   
   -test the proportionality for subdistribution hazard ratios (PSH.test, Schoenfeld residual, cumulative subdistribution hazards)
   
   -evaluate interaction between treatment effect and individual-level covariates

We also provide Meta-analysis_ex.pdf as an example of results of the analysis using the compR data (on local relapse).



