# Meta-analysis-for-competing-risk
We propose a formal guidelines to conduct Individual Patient Data  meta-analysis with competing risk endpoints. 
We propose formal guidelines to conduct an Individual Patient Data meta-analysis with competing risk endpoints. We consider a one-step approach stratified by study.

We provide:

- compR: a recreated data set with trial id (trial), group of treatment (grpCT), time of event (failure_time), type of event (failure_type), treatment information (treated/control) and age (continuous variable)
- Meta-analysis_for_compR: R script with functions for
 
         -drawing forest plots
         
         -fit a stratified Fine-Gray model for subdistribution hazards
         
         -perform landmark approaches for subdistribution hazards
         
         -test the proportionality for subdistribution hazard ratios (PSH.test, Schoenfeld residual, cumulative subdistribution hazards)
         
         -evaluate interaction between treatment effect and individual-level covariates
 
- Meta-analysis_ex: results of the proposed analysis using the compR data (analysis on local relapse).
 

The code was written using :
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)

attached base packages:
grid, stats, graphics, grDevices, utils, datasets, methods, base     

packages useful for the code:
etm_1.0.5  ggplot2_3.2.1   dplyr_0.8.3  crrSC_1.1   meta_4.9-7      

cmprsk_2.2-9    survival_2.41-3  prodlim_2019.11.13


Data is confidential, thus we provide a recreated data set (compR.Rdata) similar to the one used for the application, with same number of trials for each treatment group. Informations on trial id (trial), group of treatment (grpCT), time of event (failure_time), type of event (failure_type), treatment information (treated/control) and age (continuous variable) are introduced in the data set.

For questions, comments or remarks about the code please contact A. Meddis (alessandra.meddis@curie.fr / alessandra.meddis@gmail.com).

 
###### Supplementary material and code for the paper "Meta-analysis of clinical trials with competing time-to-event endpoints" by A. Meddis, A. Latouche, B. Zhou, S. Michiels and J. Fine .
