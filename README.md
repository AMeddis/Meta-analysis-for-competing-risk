# Meta-analysis-for-competing-risk
We propose a formal guidelines to conduct Individual Patient Data  meta-analysis with competing risk endpoints. 
We consider a one-steg approach with fixed effect model (mild heterofenity handled with stratified approach)

We provide:

- compR: a data set with trial id (trial), group of treatment (grpCT), time of event (failure_time), type of event (failure_type), treatment information (treated/control) and age (continuous variable)
- Meta-analysis for competing risks: script with all the function in R we implemented for:
 
         -forest plot (treatment effect and quantification of heterogeneity)
         
         -stratified Fine-Gray model for subdistribution hazards 
         
         -landmark approaches for subdistribution hazards
         
         -test of proportionality for subdistribution hazards ratios (PSH.test, Schoenfeld residual, cumulative subdistribution hazards)
         
         -interaction between treatment effect and inidividual level variable
 
- Meta-analysis_results: results of the proposed analysis for a typical analysis using the compR data (anlaysis on local relapse).



 
