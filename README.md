# Meta-analysis-for-competing-risk
We propose a formal guidelines to conduct Individual Patient Data  meta-analysis with competing risk endpoints. 
We propose formal guidelines to conduct an Individual Patient Data meta-analysis with competing risk endpoints. We consider a one-step approach stratified by study.

We provide:

- compR: a data set with trial id (trial), group of treatment (grpCT), time of event (failure_time), type of event (failure_type), treatment information (treated/control) and age (continuous variable)
- Meta-analysis for compR: R script with functions for
 
         -drawing forest plots
         
         -fit a stratified Fine-Gray model for subdistribution hazards
         
         -perform landmark approaches for subdistribution hazards
         
         -test the proportionality for subdistribution hazard ratios (PSH.test, Schoenfeld residual, cumulative subdistribution hazards)
         
         -evaluate interaction between treatment effect and individual-level covariates
 
- Meta-analysis_example: results of the proposed analysis using the compR data (analysis on local relapse).
 


Package used for the code:
-survival, comprisk, crrSC
-meta,etm
-grid,dplyr,ggplot2

 
