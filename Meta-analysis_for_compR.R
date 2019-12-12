#######################################################################################
#
#                                                                                     
#   Filename    :	Meta-analysis_for_compR.R    												  
#                                                                                     
#   Project     :       BiomJ article "Meta-analysis of clinical trials with competing time-to-event endpoints"                                                            
#   Authors     :       Alessandra Meddis, Aurélien Latouche, Bingqing Zhou, Stefan Michiels, Jason Fine                                                               
#   Date        :       19/11/2019
#   Purpose     :       As described in BiomJ article
#																				  
#   R Version   :       R-3.5.0 (2018-04-23)                                                               
#  
#
#   Input data files  :    compR.Rdata                                                           
#   Output data files :    ---
#
#   Required R packages :  etm, ggplot2, dplyr, crrSC, meta, cmprsk, survival, grid , prodlim
#
#
#######################################################################################



##Analysis of an Individual patient data meta-analysis in a context of competing events:
#1.Forest plot for subdistribution hazards ratios and cause specific hazards ratios.
#2.Stratifed Fine-Gray model for fixed treatment effect
#3.Landmark of subdistribution hazards 
#4.Proportionality assumption: PSH.test, Schoenfeld's residuals plot and cumulative subdistribution hazards plot
#5.Intercation between treatment effect and inividual level variable for a startifed Fine-Gray model

#####packages
library(dplyr) 
library(grid)
library(ggplot2)

library(prodlim) #stacked plot cum inc
library(cmprsk)  #SHR (section3)
library(survival) #CSHR (section 3)
library(meta)  #Forest plot (Figures 3 & 4)
library(crrSC) #stratified SHR (section 3.1)
library(etm) #Schoenfeld residuals plot and cumulative SH (Figures in supplementary material)

rm(list = ls())
#Example data
load("compR.Rdata")
#trial: trial identifiant (1,2,3...,23)
#trt:treatment information (treated/control)
#grpCT: group of treatment (Adjuvant/Concomitant/Concomitant+adjuvant/Induction )
#patid: patients id
#failure_time: time of event 
#failure_type: type of event (Death without failure/Distant failure/ Loco-regional failure/No event )
#age (continuous variable)
#failure_type.n: coded failure_type (0: censored-no event, 1:local relapse-Loco-regional failure
#                                    2: distant relapse-Distant failure, 3:death without failure) 
#trt_temp: coded treatment variable (1:treated,0: control)


#---------------------------------------------------------------------------
# Section 3: Inidividual patient data meta-analysis for competing endpoints #
#---------------------------------------------------------------------------#

### Section 3.1: Regression model for competing risks 

#-----stratified Fine-Gray model for fixed effect (distant relapse)----------
#library(crrSC)
mod_str=crrs(ftime=compR$failure_time, fstatus=compR$failure_type.n,
             cov1=compR$trt_temp, strata=compR$trial, failcode=2,
             cencode=0)

SHR_str=c(exp(mod_str$coef),
          round(exp(mod_str$coef-1.96*sqrt(mod_str$var[1,1])),3),
          round(exp(mod_str$coef+1.96*sqrt(mod_str$var[1,1])),3))

#--------Test of proportionality for SHR---------

##PSH.test

#library(crrSC)
compR$trt_temp<-ifelse(compR$trt=="Control", 0,1)
#goodness of fit for each stratum (trial) id
gof<-function(id){
  #considero il trial id
  temp<-compR[compR$trial==id,]
  #f(t)=t^2
  gofx2<-psh.test(time=temp$failure_time, fstatus=temp$failure_type.n, z=temp$trt_temp,
                  D=c(1),tf=function(x) c(x^2))
  #f(t)=t
  gofx<-psh.test(time=temp$failure_time, fstatus=temp$failure_type.n, z=temp$trt_temp,
                 D=c(1),tf=function(x) c(x))
  #f(t)=log(t)
  goflog<-psh.test(time=temp$failure_time, fstatus=temp$failure_type.n, z=temp$trt_temp,
                   D=c(1),tf=function(x) c(log(x)))
  
  #pvalue for each test
  c(id,gofx2[1,5],gofx[1,5],goflog[1,5])
}

trial<-unique(compR$trial)

#data with trial id and pvalue for each test
gof_trial<-do.call(rbind,lapply(trial,FUN = gof))

#which trial reject (pvalue<0.05) the null hypothesis for t^2 (non proportionality)
np_indx2<-list(trial=which(gof_trial[,2]<0.05),pvalue=gof_trial[which(gof_trial[,2]<0.05),2])

#which trial reject (pvalue<0.05) the null hypothesis for t (non proportionality)
np_indx<-list(trial=which(gof_trial[,3]<0.05),pvalue=gof_trial[which(gof_trial[,3]<0.05),3])

#which trial reject (pvalue<0.05) the null hypothesis for log(t) (non proportionality)
np_indlogx<-list(trial=which(gof_trial[,4]<0.05),pvalue=gof_trial[which(gof_trial[,4]<0.05),4])



##Schoenfeld residuals plot and cumulative SH

##For the Cumulative subdistribution hazards we estimate the cumulative incidence function
#considering the specific case of a multi-state model (A-J model)
#where the only possible transitions are initial state - event with event=local,distant or death

#library(etm)

## create the etm database
compR$from<-rep(8,nrow(compR)) #we define 8 as initial state
compR$id<-compR$patid
compR$time<-compR$failure_time
compR$to<-compR$failure_type.n

#define the transition matrix
tra<-matrix(FALSE,4,4)
#put TRUE to 1-->2,3,4 because we are in competing risk settings 
#and in 1 is the initial state and 2,3,4 the possible events
tra[1,2:4]<-TRUE

#for each trial id we estimate the cumulative SHR (distant relapse)
cumSH<-function(id){
  #data for each trial
  temp<-compR[compR$trial==id,]
  #A-J model for treat 1 or 0
  #the possible state are 8,1,2,3 (initial state,local, distant, death)
  AJ.z1<-etm(temp[temp$trt_temp == 1, ], c("8","1","2","3"), tra, "0", 0)
  AJ.z0<-etm(temp[temp$trt_temp == 0, ], c("8","1","2","3"), tra, "0", 0)
  ###time of events (for both treated and control)
  times <- sort(c(AJ.z0$time, AJ.z1$time))
  #estimation of probability of transition (form initial state to distant relapse)
  cif.z0 <- trprob(AJ.z0, tr.choice = "8 2", timepoints = times)
  cif.z1 <- trprob(AJ.z1, tr.choice = "8 2", timepoints = times)
  #cum sub hazards from the cumulative incidence function
  sub.haz.z0 <- cumsum(1 -((1 - cif.z0) /
                             (1 - c(0, cif.z0[-length(cif.z0)]))))
  sub.haz.z1 <- cumsum(1 -((1 - cif.z1) /
                             (1 - c(0, cif.z1[-length(cif.z1)]))))
  
  data.frame(rep(id,length(times)),sub.haz.z0,sub.haz.z1,times)
}

cumSH_data<-do.call(rbind,lapply(trial,FUN = cumSH))
colnames(cumSH_data)<-c("trial","cumSH0","cumSH1","times")

# plot cumSH and Schoenfeld residuals (Figure in supplementary material)
#trial 16

par(mfrow=c(1,2))
##cumSH
tmp<-cumSH_data[cumSH_data$trial==16,] 
temp<-compR[compR$trial==16,] #for the SCh plot
plot(tmp$cumSH0, tmp$cumSH1, lwd = 2, type = "s",
     xlab = expression(hat(Lambda)(t, "Z = 0")),
     ylab = expression(hat(Lambda)(t, "Z = 1")), 
     main=paste("cumSHplot-trial=",16,"nev=",nrow(filter(temp,to==2)),""))

shfit.lr <- with(compR[compR$trial==16,], crr(time, failure_type.n ,trt_temp ,
                                              failcode = 1, cencode = 0))
abline(a = 0, b = exp(shfit.lr$coef), col = "darkgray", lwd = 2)

#Sch residuals plot
cox.rel.mv <- coxph(Surv(failure_time, failure_type.n == 2) ~ trt_temp, temp)
plot(cox.zph(cox.rel.mv, transform='identity'), ylab="beta(t)",
     main=paste("schoenfeld residual-trial=",16,""))

#dev.off()


### Section 3.3: Effect of follow-up

#--------Landmark of treatment effect (distant relapse)------

# Landmark for the entire data set considering as landmark times the FUP of the studies
#(difference in prediction if the follow-up is longer/shorter)

#library(dplyr)
#library(crrSC)
#library(ggplot2)

#calculate the FUP for each study
FUP<-tapply(compR$failure_time,compR$trial,median)
#create the landmark times sequence
#before min(FUP) we define more 4 ftimes where to estimate the treatment effect(SHR)
ldmt<-c(seq(1,min(FUP)-1,length.out = 4),sort(FUP))

#estimate of the SHR for each landmark time
SHR_l<-sapply(1:(length(ldmt)-1), FUN=function(k){
  #select subjects at risk for the specific landmark time (ldmt)
  admc<-filter(compR, failure_time>ldmt[k])
  #once we defined the database we use a stratified Fine-Gray model
  mod_str=crrs(ftime=admc$failure_time, fstatus=admc$failure_type.n,
               cov1=admc$trt_temp, failcode=2, strata=admc$trial, #distant relapse(filcode=2)
               cencode=0)
  #SHR and c.i.
  SHR_s=exp(mod_str$coef)
  inf=round(exp(mod_str$coef-1.96*sqrt(mod_str$var[1,1])),3)
  sup=round(exp(mod_str$coef+1.96*sqrt(mod_str$var[1,1])),3)
  c(SHR_s,inf,sup)})

#we create the database for the ggplot
SHR_l<-t(SHR_l)
SHR_l<-data.frame(time=ldmt[1:(length(ldmt)-1)], SHR=SHR_l[,1], 
                  inf=SHR_l[,2],sup=SHR_l[,3])

ggplot(SHR_l, 
       aes(x=time,ymin=inf, 
           middle=SHR,  ymax=sup )) +
  geom_errorbar(stat="identity")+
  xlim(0,15)+
  geom_point(aes(x=time,y=SHR), lwd=1, col="red")+
  geom_abline(slope=0,intercept = 1,lty=3)


### Section 3.4: Assesstreatment-by-covariate interaction

#-----interactions age:treatment (distant relapse)----------
##design matrix
M<-model.matrix( ~ 0+ compR$trt +compR$age + compR$trt:compR$age )

#stratified Fine-Gray model
mod_str=crrs(ftime=compR$failure_time, fstatus=compR$failure_type.n,
             cov1=M[,2:4], strata=compR$trial, failcode=2, # distant failure
             cencode=0)                   
#beta_trt<-c(coef=round(exp(mod_str$coef[1]),3),
#            inf=round(exp(mod_str$coef[1]-1.96*sqrt(mod_str$var[1,1])),3),
#            sup=round(exp(mod_str$coef[1]+1.96*sqrt(mod_str$var[1,1])),3))

#interaction term
beta_int<-c(coef=round(exp(mod_str$coef[3]),3),
            inf=round(exp(mod_str$coef[3]-1.96*sqrt(mod_str$var[3,3])),3),
            sup=round(exp(mod_str$coef[3]+1.96*sqrt(mod_str$var[3,3])),3))


#------interaction age:treatment modalities (distant relapse)-------------

#we want to understant if the age is significant in the treatment efficacy for a specific treatment type
#thus we first stratify on the treatment type (e.g. adjuvant) and 
#we employ the stratified (on trials) Fine-Gray model in this subgroup
#with the treatment variable and the interaction term.

#first stratification (treatment modality = adjuvant)
adj<-filter(compR,grpCT=="Adjuvant")
#design Matrix
M2<-model.matrix( ~ 0+ adj$trt +adj$age + adj$trt:adj$age )
#stratified Fine-Gray model
mod_adj<-crrs(ftime=adj$failure_time, fstatus=adj$failure_type.n,
              cov1=M2[,2:4],
              strata=adj$trial, failcode=2, #distant relapse
              cencode=0)
#interaction trt:age for adjuvant chemotherapy       
beta_adj_int<-c(coef=round(exp(mod_adj$coef[3]),3),
                inf=round(exp(mod_adj$coef[3]-1.96*sqrt(mod_adj$var[3,3])),3),
                sup=round(exp(mod_adj$coef[3]+1.96*sqrt(mod_adj$var[3,3])),3))

#beta_adj_int




#----------------------------------------------------------------
#     Section4: Statistical analysis and results                 #
#----------------------------------------------------------------#


###############
#Figure 1 and 2: stacked plot of the cumulative incidence function
###############

#Figure 1 (all individuals)
ajfit <- prodlim(Hist(failure_time,failure_type.n) ~ 1, data = compR)
#postscript(file="results/Figure_1.eps")
plot(ajfit,cause ="stacked",type="cuminc",legend = FALSE,
     lty=c(1,1,1),col=c("black","azure4","azure3"),atrisk=FALSE)
legend("topleft", c("local relapse","distant relapse","death without failure"), 
       lty=c(1,1,1), lwd=c(2,2,2),col=c("black","azure4","azure3"))
#dev.off()

#Figure 2 (two different arms: control and treatment)

#2a (on the left) people treated by chemo
ajfit2 <- prodlim(Hist(failure_time,failure_type.n) ~ 1, data = filter(compR, trt_temp==0))
#postscript(file="results/Figure_2a.eps")
plot(ajfit2,cause ="stacked",type="cuminc",legend = FALSE, 
     lty=c(1,1,1),col=c("black","azure4","azure3"),atrisk=FALSE)
legend("topleft", c("local relapse","distant relapse","death without failure"), 
       lty=c(1,1,1), lwd=c(2,2,2),col=c("black","azure4","azure3"))
title("Cumulative incidence for people treated by chemo")
#dev.off()

#2b (on the right) people treated by chemo + radio
ajfit2t <- prodlim(Hist(failure_time,failure_type.n) ~ 1, data = filter(compR, trt_temp==1))
#postscript(file="results/Figure_2b.eps")
plot(ajfit2t,cause ="stacked",type="cuminc",legend = FALSE,
     lty=c(1,1,1),col=c("black","azure4","azure3"),atrisk=FALSE)
legend("topleft", c("local relapse","distant relapse","death without failure"), 
       lty=c(1,1,1), lwd=c(2,2,2),col=c("black","azure4","azure3"))
title("Cumulative incidence for people treated by chemo+radio")
#dev.off()

###############
#Figure 3 and 4: forest plot for SHR and CSHR
###############

##We create a function that create a forest plot for the desired quantity:

#@importFrom cmprsk survival 
#@importFrom meta 


#@param ftime: vector of failure/censoring times
#@param fstatus: vector of type of events (0 for censoring)
#@param trt: vector of the treatment (treated/control)
#@param failcode: code of fstatus that denotes the failure type of interest
#@param group: variable which define subgroups of subjects
#@param strata: stratification covariate (trial)
#@param quantity: a string which define the quantity we are interested on "SHR" or "CSHR"

#@return object is an object of metagen (meta library):

ForestPlot<-function(ftime,fstatus,trt, failcode, group=NULL, strata, quantity="SHR"){
  library(cmprsk) #SHR
  library(survival) #CSHR
  library(meta)  #Forest plot
  library(grid)
  
  ##create the data frame
  data=data.frame(ftime=ftime,status=fstatus,trt=trt,trial=strata)
  
  if(!is.null(group))
  {data$group<-group}
  
  trial=unique(strata)
  
  ##initialize the object
  res=c()
  
  if(quantity=="SHR"){
    #for each strata we evaluate the SHR (crr function)
    for(i in 1:length(trial)) {
      #select data for each strata
      temp<-data[data$trial==trial[i],]
      #model for SHR
      model=crr(ftime=temp$ftime, fstatus=temp$status,
                cov1=temp$trt, failcode=failcode, 
                cencode=0)
      #SHR and conf int and the standard error(useful for the quantification of heterogeneity)
      res=rbind(res,
                c(round(exp(model$coef),3),
                  round(exp(model$coef-1.96*sqrt(model$var[1,1])),3),
                  round(exp(model$coef+1.96*sqrt(model$var[1,1])),3),round(sqrt(model$var[1,1]),3)))
    }
    
    #add the index of strata
    res=data.frame(trial=trial, res)
    colnames(res)=c("trial", "SHR", "inf", "sup","SE")
    fplot_quantity="Subdistribution HR"
  }
  
  
  else if(quantity=="CSHR"){
    ##for each strata we estimate the CSHR (coxph function)
    for(i in 1:length(trial)) {
      #stratification
      temp<-data[data$trial==trial[i],]
      #model for CSHR
      model<-coxph(formula = Surv(ftime, status == failcode) ~ trt,
                   data= temp) 
      #CSHR and conf int  and SE (for the quantification of heterogeneity)
      res=rbind(res,
                c(round(exp(model$coef),3),
                  round(exp(model$coef-1.96*sqrt(model$var[1,1])),3),
                  round(exp(model$coef+1.96*sqrt(model$var[1,1])),3),round(sqrt(model$var[1,1]),3)))
      
    }
    #add the index of strata
    res=data.frame(trial=trial, res)
    colnames(res)=c("trial", "CSHR", "inf", "sup","SE")
    fplot_quantity="Cause-specific HR"
  }
  else{stop("quantity of interest is not well specified")}
  
  if(!is.null(group)){
    #we want to add the information on the groups
    tabFP=merge(res, unique(data.frame(trial=data$trial, n_g=data$group)),
                by="trial")
    #we need numeric variable for group definition
    tabFP=data.frame(tabFP, grp=as.numeric(tabFP$n_g))
    
    ###we proceed with the Forest plot (meta funct)
    
    treateffect<-eval(parse(text=paste("tabFP$",quantity,"", sep="")))
    dFP=metagen(TE=log(treateffect),  # Treatment effect, i.e. log hazard ratios
                seTE=tabFP$SE,     # Strandard error of treatment estimate
                sm="HR",         # Summary measure
                data=tabFP, studlab =trial, 
                comb.random=F, byvar=n_g)
    
    ##plot it
    forest(dFP, sortvar=tabFP$grp, col.study=tabFP$grp, 
           rightcols=c("effect", "ci"),
           leftcols=c("studlab"),
           leftlabs=c("Study"),smlab=paste("",fplot_quantity," "),rightlabs = c(paste("",quantity," "),"CI-95%"),
           colgap.forest.left = unit(4,"cm"),resid.hetstat=FALSE,
           bylab="",print.I2.ci = TRUE)
    
    return(dFP)
  }
  
  else{
    
    ###we proceed with the Forest plot (meta funct)
    
    treateffect<-eval(parse(text=paste("res$",quantity,"", sep="")))
    FP=metagen(TE=log(treateffect),  # Treatment effect, i.e. log hazard ratios
               seTE=res$SE,     # Strandard error of treatment estimate
               sm="HR",         # Summary measure
               comb.random=F)   #fixed effect model
    
    ##plot it
    #alternative of code in case of error with unit("4cm") (depending on the package version)
    forest(FP, 
           rightcols=c("effect", "ci"),leftcols=c("studlab"),
           leftlabs=c("Study"),smlab=fplot_quantity,rightlabs = c(paste("",quantity," "),"CI-95%"),
           colgap.forest.left = unit(4,"cm"),resid.hetstat=FALSE,
           bylab="",print.I2.ci = TRUE)
    return(FP)
  }
  
  
}

##Figure 3 
#on the left (SHR for local relapse)
#postscript(file="../results/Figure_3a.eps", width=7, height=9)
FP_SHR_lr<-ForestPlot(ftime=compR$failure_time,fstatus=compR$failure_type.n,trt=compR$trt_temp,
           failcode=1, strata=compR$trial, group= compR$grpCT, quantity="SHR")
#dev.off()
#on the right (CSHR for local relapse)
#postscript(file="../results/Figure_3b.eps", width=7, height=9)
FP_CSHR_lr<-ForestPlot(ftime=compR$failure_time,fstatus=compR$failure_type.n,trt=compR$trt_temp,
                      failcode=1, strata=compR$trial, group= compR$grpCT, quantity="CSHR")
#dev.off()

##Figure 4
#on the left (SHR for distant relapse)
#postscript(file="../results/Figure_4a.eps", width=7, height=9)
FP_SHR_dr<-ForestPlot(ftime=compR$failure_time,fstatus=compR$failure_type.n,trt=compR$trt_temp,
           failcode=2, strata=compR$trial, group= compR$grpCT, quantity="SHR")
#dev.off()
#on the rigth (CSHR for distant relapse)
#postscript(file="../results/Figure_4b.eps", width=7, height=9)
FP_SHR_lr<-ForestPlot(ftime=compR$failure_time,fstatus=compR$failure_type.n,trt=compR$trt_temp,
                      failcode=2, strata=compR$trial, group= compR$grpCT, quantity="CSHR")
#dev.off()


###############
#Figure 5 and 6: Landmark of treatment effect 
###############

#Landmark for each study (time-varying trt effect--> non proportionality) 
#we consider the landmark times in order to have enough observations

#library(dplyr)
#library(cmprsk) 
#library(ggplot2)

trial<-unique(compR$trial)

#function for each trial i 
#where failcode defines the endpoint of interest 
SHRt_trial<-function(i,failcode){
  #define the database for the strata i 
  temp<-compR[compR$trial==trial[i],]
  #define the landmark time (sequence of 8 times)
  tstar<-seq(min(temp$failure_time),
             quantile(filter(temp,failure_type.n==2)$failure_time,0.75),length.out = 8)
  dit<-length(tstar)
  #for each time we have SHR in this strata (trial) i 
  a<-sapply(1:dit, FUN=function(k){
    #subjects at risk
    admc<-filter(temp, failure_time>=tstar[k])
    #Fine-Gray model
    mod_str=crr(ftime=admc$failure_time, fstatus=admc$failure_type.n,
                cov1=admc$trt_temp, failcode=failcode, 
                cencode=0)
    #SHR(tstar) and c.i.
    SHR_s=exp(mod_str$coef)
    inf=round(exp(mod_str$coef-1.96*sqrt(mod_str$var[1,1])),3)
    sup=round(exp(mod_str$coef+1.96*sqrt(mod_str$var[1,1])),3)
    c(SHR_s,inf,sup)})
  
  cbind(time=tstar,t(a),trial=rep(i,length(tstar)))
}


#Local relapse (Figure 5)

#call the function for each trial and we create the database
SHRt_t<-do.call(rbind,lapply(trial,1,FUN = SHRt_trial))

#prepare data for ggplot
colnames(SHRt_t)=c('time',"SHR","inf", "sup","trial")
SHRt_t<-data.frame(SHRt_t)

#plot for each trials
#postscript(file="../results/Figure_5.eps", width=9,height=9)
ggplot(SHRt_t, 
       aes(x=time,ymin=inf,
           middle=SHR,  ymax=sup )) +
  geom_errorbar(stat="identity")+
  geom_point(aes(x=time,y=SHR), lwd=1, col="red")+
  geom_abline(slope=0,intercept = 1,lty=3)+
  ylim(0,6)+
  facet_wrap( ~ trial, ncol=4)
#dev.off()

#Distant relapse (Figure 6)

#call the function for each trial and we create the database
SHRt_t<-do.call(rbind,lapply(trial,2,FUN = SHRt_trial))

#prepare data for ggplot
colnames(SHRt_t)=c('time',"SHR","inf", "sup","trial")
SHRt_t<-data.frame(SHRt_t)


#plot for each trials
#postscript(file="../results/Figure_6.eps", width=9,height=9)
ggplot(SHRt_t, 
       aes(x=time,ymin=inf,
           middle=SHR,  ymax=sup )) +
  geom_errorbar(stat="identity")+
  geom_point(aes(x=time,y=SHR), lwd=1, col="red")+
  geom_abline(slope=0,intercept = 1,lty=3)+
  ylim(0,6)+
  facet_wrap( ~ trial, ncol=4)

#dev.off()

#Note: in the paper we plot the landmark of SHR for trials where non-proportionality was detected by the PSH.test
#using the same function, considering the subset : 
#filter(SHRt_t,trial %in% c(8,22)) for Figure 5
#filter(SHRt_t,trial==4) for Figure 6

#Figure 5
# SHR_822<-filter(SHRt_t,trial %in% c(8,22))
# SHR_822$trial<-ifelse(SHR_822$trial==8,"trial 8","trial 22")
# #plot for 8 and 22
# ggplot(SHR_822,
#        aes(x=time,y=SHR, ymin=inf, ymax=sup)) +
#   geom_linerange(lwd=1.5)+
#   geom_point(col="red",lwd=2)+
#   # geom_errorbar()+ 
#   geom_abline(slope=0,intercept = 1,lty=2)+
#   ylim(0,4)+
#   xlim(0,3)+
#   facet_wrap( ~ trial, ncol=2)




                   





