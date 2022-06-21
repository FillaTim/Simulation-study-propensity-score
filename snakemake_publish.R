#This script contains the RCode for simulation of propensity score methods.
#This script was called for different settings by a snakemake file and run 
#on the high performance cluster at HHU
#author: Tim Filla
############################

#use arguments given by the sankemake command
args = commandArgs(trailingOnly=TRUE)
#sample size
sasi<-as.numeric(args[1])
#effect of the covariates on the outcome
beta_cov<-bias<-as.numeric(args[2])
#effect of the treatment on the outcome
beta_tr<-as.numeric(args[3])
#linear or non-linear propensity score model
prop_mod<-as.character(args[4])
#which scenario should be used 
#1 and 2 are good overlap, 3 and 4 are bad overlap
#1 and 3 are high treatment prevalence, 2 and 4 are low treatment prevalence
scenario<-as.numeric(args[5])
k<-as.numeric(args[6])
outname<-args[7]
if(bias<0){
  bias<-NULL
}
#add the libPath, where the needed packages are stored
.libPaths(c("/gpfs/project/tifil100/"))
library(survival)
library(cobalt)
library(WeightIt)
library(optweight)
library(ebal)
library(optimr)

#source the functions for robust variance estimation
setwd("/gpfs/project/tifil100/snakemake_binary/RCode")
source("variance_estimator_function.R")
#load the matrix containing the coefficients for the propensity score model
load("ps_coefficient_matrix.RData")
modelling_approach_functions<-function(dataset,treatment,survey_weights=NULL,grad=1,linkf="logit"){
  p_glob<-length(which(dataset[,treatment]==1))/length(dataset[,treatment])
  #grad states whether higher moment of continuous covariates should
  #be included in the logistic regression model
  #cvar is a binary vector that states whether the i-th
  #variable is continuous (cvar=1) or not (cvar=0)
  cvar<-rep(0,length(names(dataset))-1)
  for(i in 1:length(names(dataset))){
    #we use the cutpoint of more then 10 different values
    #for a continuous covariate
    if(length(unique(dataset[,names(dataset)[i]]))>10){
      cvar[i]<-1
    }
  }
  
  #first we paste all names of the covariates and +, which results in
  #namen_formular=name1+name2+name3+name4+name5+....
  #Hereby we check if the name corresponds to the treatment name
  #and if true we exclude it.
  namen_formular<-c("")
  for(j in 1:length(names(dataset))){
    if(names(dataset)[j]!=treatment){
      namen_formular<-paste(namen_formular,names(dataset)[j],"+",sep="")
    }
  }
  #we exclude the last "+"
  namen_formular<-substr(namen_formular,1,nchar(namen_formular)-1)
  #if there is any continuous covariate and grad >1 we add the term 
  #+I(namexy^p) whereby p ranges from 2 to grad
  if(grad>1 &any(cvar==1) ){
    for(p in 2:grad){
      for(k in which(cvar==1)){
        namen_formular<-paste(namen_formular,"+","I(",names(dataset)[k],"^",p,")",sep="")
      }
    }
  }
  #Finally we need to write the name of the reference category
  #with ~
  #The final formula looks like that
  #treatment~name1+name2+name3+...+nameK+namec1^2+namecK^2+...+namec1^grad+namecK^grad
  namen_formular<-paste(treatment,"~",namen_formular)
  #calculate the parameter of the logistic regression model
  if(is.null(survey_weights)){
    modell<-glm(as.formula(namen_formular),data=dataset,family="binomial"(link=linkf))
  }
  else{
    modell<-glm(as.formula(namen_formular),data=dataset,weights=survey_weights,family="binomial"(link=linkf))
  }
  beta<-modell$coefficients
  
  #get the propensity score
  ps<-predict(modell,newdata=dataset,type="response")
  Z<-dataset[,treatment]
  
  #calculate the weights from the propensity score by using the different 
  #weighting functions
  #matching weights
  wmw<-sapply(1:length(ps),function(x) min(c(1-ps[x],ps[x]))/(Z[x]*ps[x]+(1-Z[x])*(1-ps[x])))
  #stabilized iptw
  wstiptwate<-sapply(1:length(ps),function(x) ifelse(Z[x]==0,(1-p_glob)/(1-ps[x]),(p_glob)/ps[x]))
  #iptw with target population ATE
  wiptwate<-sapply(1:length(ps),function(x) ifelse(Z[x]==0,1/(1-ps[x]),1/ps[x]))
  #iptw with target population ATT
  wiptwatt<-sapply(1:length(ps),function(x) ifelse(Z[x]==0,ps[x]/(1-ps[x]),1))
  #overlap weights
  woverlap<-sapply(1:length(ps),function(x) ifelse(Z[x]==1,1-ps[x],ps[x]))
  return(list("overlap"=woverlap,"mw"=wmw,"iptwate"=wiptwate,"iptwatt"=wiptwatt,"stiptwate" =wstiptwate)) 
}
gmmw<-function(dataset,treatment="treatment"){
  require(MASS)
  bal.loss<-function(beta.curr,sample.weights0=sample.weights){
    sample.weights<-sample.weights0
    ##Generate theta and probabilities.
    theta.curr<-as.vector(X%*%beta.curr)
    probs.curr<-(1+exp(-theta.curr))^-1
    probs.curr<-pmin(1-probs.min,probs.curr)
    probs.curr<-pmax(probs.min,probs.curr)
    ##Generate weights.
    w.curr<-sapply(1:length(probs.curr),function(z) min(c(probs.curr[z],1-probs.curr[z]))*(treat[z]*(1-probs.curr[z])-(1-treat[z])*probs.curr[z])/(probs.curr[z]*(1-probs.curr[z])))
    w.curr[which(treat==1)]<-w.curr[which(treat==1)]/sum(w.curr[which(treat==1)])
    w.curr[which(treat==0)]<-w.curr[which(treat==0)]/sum(abs(w.curr[which(treat==0)]))
    ##Generate mean imbalance.
    Xprimew <- t(sample.weights*X)%*%(w.curr)
    
    #calculate the loss
    loss<-t(Xprimew)%*%Xprimew
    return(loss)
  } 
  cr_weights<-function(beta.curr){
    ##Generate theta and probabilities.
    theta.curr<-as.vector(X%*%beta.curr)
    probs.curr<-(1+exp(-theta.curr))^-1
    probs.curr<-pmin(1-probs.min,probs.curr)
    probs.curr<-pmax(probs.min,probs.curr)
    ##Generate weights.
    
    w.curr<-sapply(1:length(probs.curr),function(z) min(c(probs.curr[z],1-probs.curr[z]))*(treat[z]*(1-probs.curr[z])-(1-treat[z])*probs.curr[z])/(probs.curr[z]*(1-probs.curr[z])))
    w.curr[which(treat==1)]<-w.curr[which(treat==1)]
    w.curr[which(treat==0)]<-(-1)*w.curr[which(treat==0)]
    return(w.curr)
  }
  
  #create treatment vector used in bal.loss
  treat<-dataset[[treatment]]
  dataset[,-which(names(dataset)==treatment)]<-scale(dataset[,-which(names(dataset)==treatment)])
  #calculate logistic regression model with all covariates in dataset_sim[[i]]
  #as explanatory variables
  mod<-glm(tr~.,data=dataset,family="binomial")
  X<-model.matrix(mod)
  thetanew<-as.matrix(X)%*%mod$coefficients
  probs.curr<-(1+exp(-thetanew))^-1
  beta_0<-mod$coefficients
  
  #sample.weights are kept constant 1 for all observations
  sample.weights<-rep(1,dim(X)[1])
  
  #take XprimeX.inv as the inverse of the estimated covariance matrix
  #XprimeX.inv <- ginv(t(sample.weights^0.5 * X) %*% (sample.weights^0.5 *X))
  
  #optimise the parameter beta_0 on a scalar range of (.8,1.1)
  alpha.func<-function(alpha) bal.loss(beta_0*alpha)
  
  #if the probability for one observation to receive treatment is 0
  #this might cause trouble due to weights with value infinity
  #therefor we want the minimum value to be 0.00001 and
  #the same for a maximum value of 1-0.00001
  probs.min<-c(0.00001)
  beta.curr<-beta_0*optimize(alpha.func,interval=c(.8,1.1))$min
  ##find the parameter beta that best fits to the moment conditions
  gmm.init<-beta.curr
  L<-optim(gmm.init,bal.loss,method="BFGS",hessian=TRUE)
  
  #calculate the weights for this parameter.
  gmmw<-cr_weights(L$par)
  theta<-as.vector(X%*%L$par)
  probs_res<-(1+exp(-theta))^-1
  return(list("gmmw"=gmmw,"L"=L,"init"=gmm.init,"fitted_values"=probs_res,"X"=X))
}
ESS_0<-function(w,treatment=x_data$tr){
  return(sum(w[which(treatment==0)])^2/sum(w[which(treatment==0)]^2))
}
ESS_1<-function(w,treatment=x_data$tr){
  return(sum(w[which(treatment==1)])^2/sum(w[which(treatment==1)]^2))
}

#set the value of the coefficients in the propensity score model
#depending on the chosen scenario
p_vec<-Erg_final_Mat[scenario,]
p_vec_null<-p_vec

#create the covariates
#we need the number of events of a binary variable with prevalence 0.08,0.19,0.33
sasi_8<-round(sasi/100*8)
sasi_19<-round(sasi/100*19)
sasi_33<-round(sasi/100*33)

#set the number of iterations
iter<-500

#create empty matrices to store the results
Tr_prev<-matrix(NA,iter,2)
Cox_result_matrix<-as.data.frame(matrix(NA,iter,15))
ESS_result_matrix_0<-as.data.frame(matrix(NA,iter,15))
ESS_result_matrix_1<-as.data.frame(matrix(NA,iter,15))
Bal_bin_raw<-Bal_bin_iptw_ate<-Bal_bin_iptwst<-Bal_bin_iptw_att<-Bal_bin_mw_ml<-Bal_bin_mw_gmm<-matrix(NA,iter,11)
Bal_bin_overlap<-Bal_bin_cbpsj_ate<-Bal_bin_cbpsj_att<-Bal_bin_cbpso_att<-Bal_bin_cbpso_ate<-matrix(NA,iter,11)
Bal_bin_ebal_ate<-Bal_bin_ebal_att<-Bal_bin_vabal_att<-Bal_bin_vabal_ate<-matrix(NA,iter,11)
Bal_con_raw<-Bal_con_iptw_ate<-Bal_con_iptwst<-Bal_con_iptw_att<-Bal_con_mw_ml<-Bal_con_mw_gmm<-matrix(NA,iter,3)
Bal_con_overlap<-Bal_con_cbpsj_ate<-Bal_con_cbpsj_att<-Bal_con_cbpso_att<-Bal_con_cbpso_ate<-matrix(NA,iter,3)
Bal_con_ebal_ate<-Bal_con_ebal_att<-Bal_con_vabal_att<-Bal_con_vabal_ate<-matrix(NA,iter,3)
Corrected_rob_standard_error<-Rob_standard_error<-as.data.frame(matrix(NA,iter,15))
names(Cox_result_matrix)<-names(ESS_result_matrix_0)<-names(ESS_result_matrix_1)<-
  names(Rob_standard_error)<-names(Corrected_rob_standard_error)<-c("raw","iptw_ate_ml","iptwst_ate_ml","iptw_att","cbpsj_ate","cbpso_ate",
                                                                    "cbpsj_att","cbpso_att","ebal_ate","vabal_ate","ebal_att","vabal_att",
                                                                    "overlap","mw_ml","mw_gmm")

#start the iteration
for(i in 1:iter){
  p_vec<-p_vec_null
  set.seed(1902+k+i)
  
  #create the covariates with three continuous independent standard gaussian distributed variables
  x_con<-MASS:::mvrnorm(n=sasi,mu=rep(0,3),Sigma = diag(1,3))
  
  #create 11 binary covariates and store them in x_bin
  x_bin<-matrix(0,nrow=sasi,ncol=11)
  #8%
  x_bin[sample(1:sasi,sasi_8,replace=FALSE),1]<-1
  x_bin[sample(1:sasi,sasi_8,replace=FALSE),2]<-1
  x_bin[sample(1:sasi,sasi_8,replace=FALSE),3]<-1
  #19%
  x_bin[sample(1:sasi,sasi_19,replace=FALSE),4]<-1
  x_bin[sample(1:sasi,sasi_19,replace=FALSE),5]<-1
  x_bin[sample(1:sasi,sasi_19,replace=FALSE),6]<-1
  x_bin[sample(1:sasi,sasi_19,replace=FALSE),7]<-1
  #33%
  x_bin[sample(1:sasi,sasi_33,replace=FALSE),8]<-1
  x_bin[sample(1:sasi,sasi_33,replace=FALSE),9]<-1
  x_bin[sample(1:sasi,sasi_33,replace=FALSE),10]<-1
  x_bin[sample(1:sasi,sasi_33,replace=FALSE),11]<-1
  
  #combine the continuous and binary covariates for propensity score estimation
  x_mat<-cbind(1,x_con,x_bin)
  
  #if the propensity score model is linear, we can remove the beta coefficients 16...21 as those are set to 0
  if(prop_mod=="linear"){
    p_vec<-p_vec[-c(16:21)]
    pr_vec<-1/(1+exp(-x_mat%*%p_vec))
  }
  #if the propensity score model is non-linear, we create the non linear variables and calculate the propensity score
  if(prop_mod=="non_linear"){
    x_mat_non_lin<-x_mat
    x_mat_non_lin<-cbind(x_mat_non_lin,x_mat[,2]^2)
    x_mat_non_lin<-cbind(x_mat_non_lin,x_mat[,3]*x_mat[,7])
    x_mat_non_lin<-cbind(x_mat_non_lin,exp(x_mat[,4]))
    x_mat_non_lin<-cbind(x_mat_non_lin,x_mat[,5]*x_mat[,12])
    x_mat_non_lin<-cbind(x_mat_non_lin,x_mat[,6]*x_mat[,8])
    x_mat_non_lin<-cbind(x_mat_non_lin,x_mat[,9]*x_mat[,13])
    pr_vec<-1/(1+exp(-x_mat_non_lin%*%p_vec))
  }
  
  #sample the treatment variable with P(T=1|X=x_i)= pr_vec_i
  tr<-sapply(1:length(pr_vec),function(z){
    return(sample(0:1,1,replace=TRUE,prob=c(1-pr_vec[z],pr_vec[z])))
  })
  
  #create the final data frame by adding the treatment variable and removing the column with only 1
  x<-cbind(tr,x_mat[,-1])
  x_data<-as.data.frame(x)
  
  #set the effect vector for outcome creation
  beta_ef<-c(beta_tr,rep(beta_cov,ncol(x_mat)-1))
  #create the survival times
  R=as.numeric((-log(runif(length(tr)))/(0.00002*exp(t(beta_ef)%*%t(x))))^(1/2))
  
  #create a survival object for cox Regression
  sobj<-Surv(R,event=rep(1,length(R)))
  
  #calculate the weights
  L<-modelling_approach_functions(x_data,"tr",grad=1)
  
  
  #CBPS
  cbps_ate_over_res<-try(CBPS:::CBPS(tr~.,data=x_data,method="over",ATT=0,standardize = FALSE))
  cbps_ate_just_res<-try(CBPS:::CBPS(tr~.,data=x_data,method="exact",ATT=0,standardize = FALSE))
  cbps_att_over_res<-try(CBPS:::CBPS(tr~.,data=x_data,method="over",ATT=1,standardize = FALSE))
  cbps_att_just_res<-try(CBPS:::CBPS(tr~.,data=x_data,method="exact",ATT=1,standardize = FALSE))
  cbps_ate_over<-cbps_ate_over_res$weights
  cbps_att_over<-cbps_att_over_res$weights
  cbps_ate_just<-cbps_ate_just_res$weights
  cbps_att_just<-cbps_att_just_res$weights
  
  
  #EBAL/VBAL
  ebal_ate<-try(weightit(tr~.,data=x_data,method="ebal",moments = 1,estimand="ATE"))$weights
  ebal_att<-try(weightit(tr~.,data=x_data,method="ebal",moments = 1,estimand="ATT"))$weights
  vabal_ate<-try(weightit(tr~.,data=x_data,method="optweight",moments = 1,estimand="ATE"))$weights
  vabal_att<-try(weightit(tr~.,data=x_data,method="optweight",moments = 1,estimand="ATT"))$weights
  
  #sometimes exact 0 values occur for ebal or vabal, but this causes problem with the cox regression function
  #therefor we replace it with 10^-13
  ebal_ate[which(ebal_ate==0)]<-0.0000000000001
  ebal_att[which(ebal_att==0)]<-0.0000000000001
  vabal_ate[which(vabal_ate==0)]<-0.0000000000001
  vabal_att[which(vabal_att==0)]<-0.0000000000001
  
  #if the convergence failed- this means all weights in one group for ATE or all weights in the control group (ATT)
  #are equal, then we set the class to "try-error". Now we will always check if class is equal to try-error
  #and only store results if that is not true.
  if(length(unique(round(ebal_ate,5)))<=2 |length(unique(round(ebal_ate[which(x_data$tr==1)],5)))<=2){
    class(ebal_ate)<-"try-error"
  }
  if(length(unique(round(ebal_att,5)))<=2 |length(unique(round(ebal_att[which(x_data$tr==0)],5)))<=2){
    class(ebal_att)<-"try-error"
  }
  if(length(unique(round(vabal_ate,5)))<=2 |length(unique(round(vabal_ate[which(x_data$tr==1)],5)))<=2){
    class(vabal_ate)<-"try-error"
  }
  if(length(unique(round(vabal_att,5)))<=5 |length(unique(round(vabal_att[which(x_data$tr==0)],5)))<=2){
    class(vabal_att)<-"try-error"
  }
  
  #gmmw
  gmmw_res<-gmmw(x_data,"tr")
  weights_gmmw<-gmmw_res$gmmw
  
  #calculate the cox regressions and store the results
  cox_raw<-coxph(sobj~tr,data=x_data)
  cox_iptwst<-coxph(sobj~tr,data=x_data,weights=L$stiptwate,robust = TRUE)
  cox_iptw_ate<-coxph(sobj~tr,data=x_data,weights=L$iptwate)
  cox_iptw_att<-coxph(sobj~tr,data=x_data,weights=L$iptwatt)
  cox_overlap<-coxph(sobj~tr,data=x_data,weights=L$overlap)
  cox_mw_ml<-coxph(sobj~tr,data=x_data,weights=L$mw)
  cox_mw_gmm<-coxph(sobj~tr,data=x_data,weights=weights_gmmw)
  cox_cbpsj_ate<-coxph(sobj~tr,data=x_data,weights=cbps_ate_just)
  cox_cbpsj_att<-coxph(sobj~tr,data=x_data,weights=cbps_att_just)
  cox_cbpso_ate<-coxph(sobj~tr,data=x_data,weights=cbps_ate_over)
  cox_cbpso_att<-coxph(sobj~tr,data=x_data,weights=cbps_att_over)
  
  #here we check if ebal or vabal really converged
  if(class(ebal_ate)!="try-error"){
    cox_ebal_ate<-coxph(sobj~tr,data=x_data,weights=ebal_ate)
  }
  if(class(ebal_att)!="try-error"){
    cox_ebal_att<-coxph(sobj~tr,data=x_data,weights=ebal_att)
  }
  if(class(vabal_ate)!="try-error"){
    cox_vabal_ate<-coxph(sobj~tr,data=x_data,weights=vabal_ate)
  }
  if(class(vabal_att)!="try-error"){
    cox_vabal_att<-coxph(sobj~tr,data=x_data,weights=vabal_att)
  }
  
  #store the results for those methods that always converge
  Cox_result_matrix[i,c(1:8,13:15)]<-c(coef(cox_raw),coef(cox_iptw_ate),coef(cox_iptwst),coef(cox_iptw_att),coef(cox_cbpsj_ate),coef(cox_cbpso_ate),
                                       coef(cox_cbpsj_att),coef(cox_cbpso_att),
                                       coef(cox_overlap),coef(cox_mw_ml),coef(cox_mw_gmm))
  #add the results of ebal or vabal if they really converged
  if(class(ebal_ate)!="try-error"){
    Cox_result_matrix[i,9]<-coef(cox_ebal_ate)
  }
  if(class(ebal_att)!="try-error"){
    Cox_result_matrix[i,11]<-coef(cox_ebal_att)
  }
  if(class(vabal_ate)!="try-error"){
    Cox_result_matrix[i,10]<-coef(cox_vabal_ate)
  }
  if(class(vabal_att)!="try-error"){
    Cox_result_matrix[i,12]<-coef(cox_vabal_att)
  }
  
  #store the results of the robust variance estimator
  Rob_standard_error[i,c(1:8,13:15)]<-c(sqrt(diag(cox_raw$var)),sqrt(diag(cox_iptw_ate$var)),sqrt(diag(cox_iptwst$var)),sqrt(diag(cox_iptw_att$var)),sqrt(diag(cox_cbpsj_ate$var)),sqrt(diag(cox_cbpso_ate$var)),
                                        sqrt(diag(cox_cbpsj_att$var)),sqrt(diag(cox_cbpso_att$var)),
                                        sqrt(diag(cox_overlap$var)),sqrt(diag(cox_mw_ml$var)),sqrt(diag(cox_mw_gmm$var)))
  if(class(ebal_ate)!="try-error"){
    Rob_standard_error[i,9]<-sqrt(diag(cox_ebal_ate$var))
  }
  if(class(ebal_att)!="try-error"){
    Rob_standard_error[i,10]<-sqrt(diag(cox_ebal_att$var))
  }
  if(class(vabal_ate)!="try-error"){
    Rob_standard_error[i,11]<-sqrt(diag(cox_vabal_ate$var))
  }
  if(class(vabal_att)!="try-error"){
    Rob_standard_error[i,12]<-sqrt(diag(cox_vabal_att$var))
  }
  
  #calculate the corrected robust variance estimator
  #create data in form so that it can be used by the variance estimator functions
  x_data_cox<-x_data
  x_data_cox$censored<-rep(1,nrow(x_data))
  x_data_cox$R<-R
  gmmw_res$X<-as.data.frame(gmmw_res$X)
  gmmw_res$X$R<-R
  gmmw_res$X$tr<-x_data$tr
  gmmw_res$X$censored<-rep(1,nrow(x_data))
  
  #v_cbpsj_ate<-var_est_cbps(x_data_cox,indA="tr",indX=names(x_data)[-1],indTime="R",indStatus = "censored",psfit=cbps_ate_just_res$fitted.values,coefficients = cbps_ate_just_res$coefficients,weights=cbps_ate_just_res$weights)
  #v_cbpsj_att<-var_est_cbps_att(x_data_cox,indA="tr",indX=names(x_data)[-1],indTime="R",indStatus = "censored",psfit=cbps_att_just_res$fitted.values,coefficients = cbps_att_just_res$coefficients,weights=cbps_att_just_res$weights)
  v_overlap<-var_est_overlap(x_data_cox,indA="tr",indX=names(x_data)[-1],indTime="R",indStatus = "censored")
  v_mw_ml<-var_est_mw(x_data_cox,indA="tr",indX=names(x_data)[-1],indTime="R",indStatus = "censored")
  v_iptw_ate<-var_est(x_data_cox,indA="tr",indX=names(x_data)[-1],indTime="R",indStatus = "censored")
  v_iptw_ate_st<-var_est_iptw_st(x_data_cox,indA="tr",indX=names(x_data)[-1],indTime="R",indStatus = "censored")
  v_iptw_att<-var_est_iptw_att(x_data_cox,indA="tr",indX=names(x_data)[-1],indTime="R",indStatus = "censored")
  v_mw_gmm<-var_est_mw_gmm(gmmw_res$X,indA="tr",indX=names(x_data)[-1],indTime="R",indStatus = "censored",psfit=gmmw_res$fitted_values,coefficients = gmmw_res$L$par,weights=gmmw_res$gmmw)
  
  #store the results of the corrected robust variance estimator
  Corrected_rob_standard_error$iptw_ate_ml[i]<-v_iptw_ate
  Corrected_rob_standard_error$iptwst_ate_ml[i]<-v_iptw_ate_st
  Corrected_rob_standard_error$iptw_att[i]<-v_iptw_att
  Corrected_rob_standard_error$overlap[i]<-v_overlap
  Corrected_rob_standard_error$mw_ml[i]<-v_mw_ml
  Corrected_rob_standard_error$mw_gmm[i]<-v_mw_gmm
  #Corrected_rob_standard_error$cbpsj_att[i]<-v_cbpsj_att
  #Corrected_rob_standard_error$cbpsj_ate[i]<-v_cbpsj_ate
  
  #calculate the effective sample sizes in treatment and control group and store the results
  ESS_result_matrix_0[i,c(1:8,13:15)]<-c(length(which(x_data$tr==0)),ESS_0(L$iptwate),ESS_0(L$stiptwate),ESS_0(L$iptwatt),ESS_0(cbps_ate_just),ESS_0(cbps_ate_over),
                                         ESS_0(cbps_att_just),ESS_0(cbps_att_over),ESS_0(L$overlap),ESS_0(L$mw),ESS_0(weights_gmmw))
  ESS_result_matrix_1[i,c(1:8,13:15)]<-c(length(which(x_data$tr==1)),ESS_1(L$iptwate),ESS_1(L$stiptwate),ESS_1(L$iptwatt),ESS_1(cbps_ate_just),ESS_1(cbps_ate_over),
                                         ESS_1(cbps_att_just),ESS_1(cbps_att_over), ESS_1(L$overlap),ESS_1(L$mw),ESS_1(weights_gmmw))
  if(class(ebal_ate)!="try-error"){
    ESS_result_matrix_1[i,9]<-ESS_1(ebal_ate)
    ESS_result_matrix_0[i,9]<-ESS_0(ebal_ate)
  }
  if(class(ebal_att)!="try-error"){
    ESS_result_matrix_1[i,11]<-ESS_1(ebal_att)
    ESS_result_matrix_0[i,11]<-ESS_0(ebal_att)
  }
  if(class(vabal_ate)!="try-error"){
    ESS_result_matrix_1[i,10]<-ESS_1(vabal_ate)
    ESS_result_matrix_0[i,10]<-ESS_0(vabal_ate)
  }
  if(class(vabal_att)!="try-error"){
    ESS_result_matrix_1[i,12]<-ESS_1(vabal_att)
    ESS_result_matrix_0[i,12]<-ESS_0(vabal_att)
  }
  
  #balance results
  #calculate the weighted standardized difference for binary covariates (function stand_diff)
  #and continuous covariates (function stand_diff_con)
  stand_diff<-function(w){
    pc<-colSums(x_data[which(tr==0),-c(1:4)]*w[which(tr==0)])/sum(w[which(tr==0)])
    pt<-colSums(x_data[which(tr==1),-c(1:4)]*w[which(tr==1)])/sum(w[which(tr==1)])
    return((pc-pt)/(sqrt((pc*(1-pc)+pt*(1-pt))/2)))
  }
  stand_diff_con<-function(w){
    m_d<-colSums(x_data[which(tr==0),c(2:4)]*w[which(tr==0)])/sum(w[which(tr==0)])-
      colSums(x_data[which(tr==1),c(2:4)]*w[which(tr==1)])/sum(w[which(tr==1)])
    s_c<-sum(w[which(tr==0)])/(sum(w[which(tr==0)])^2-sum(w[which(tr==0)]^2))*
      sum(w[which(tr==0)]*(x_data[which(tr==0),c(2:4)]-colSums(x_data[which(tr==0),c(2:4)]*w[which(tr==0)])/sum(w[which(tr==0)]))^2)
    s_t<-sum(w[which(tr==1)])/(sum(w[which(tr==1)])^2-sum(w[which(tr==1)]^2))*
      sum(w[which(tr==1)]*(x_data[which(tr==1),c(2:4)]-colSums(x_data[which(tr==1),c(2:4)]*w[which(tr==1)])/sum(w[which(tr==1)]))^2)
    return(m_d/sqrt((s_c+s_t)/2))
  }
  x_data_bal<-x_data
  x_data_bal[,2:4]<-scale(x_data_bal[,2:4],center=FALSE)
  #for the raw data all weights are set to 1
  md_bin_raw<-stand_diff(rep(1,length(L$overlap)))
  md_bin_iptw_ate<-stand_diff(L$iptwate)
  md_bin_iptw_att<-stand_diff(L$iptwatt)
  md_bin_cbpsj_ate<-stand_diff(cbps_ate_just)
  md_bin_cbpso_ate<-stand_diff(cbps_ate_over)
  md_bin_cbpsj_att<-stand_diff(cbps_att_just)
  md_bin_cbpso_att<-stand_diff(cbps_att_over)
  if(class(ebal_ate)!="try-error"){
    md_bin_ebal_ate<-stand_diff(ebal_ate)
    Bal_bin_ebal_ate[i,]<-md_bin_ebal_ate
  }
  if(class(ebal_att)!="try-error"){
    md_bin_ebal_att<-stand_diff(ebal_att)
    Bal_bin_ebal_att[i,]<-md_bin_ebal_att
  }
  if(class(vabal_ate)!="try-error"){
    md_bin_vabal_ate<-stand_diff(vabal_ate)
    Bal_bin_vabal_ate[i,]<-md_bin_vabal_ate
  }
  if(class(vabal_att)!="try-error"){
    md_bin_vabal_att<-stand_diff(vabal_att)
    Bal_bin_vabal_att[i,]<-md_bin_vabal_att
  }
  md_bin_overlap<-stand_diff(L$overlap)
  md_bin_mw_ml<-stand_diff(L$mw)
  md_bin_mw_gmm<-stand_diff(weights_gmmw)
  Bal_bin_cbpsj_ate[i,]<-md_bin_cbpsj_ate
  Bal_bin_cbpso_att[i,]<-md_bin_cbpso_att
  Bal_bin_cbpsj_att[i,]<-md_bin_cbpsj_ate
  Bal_bin_cbpso_ate[i,]<-md_bin_cbpso_ate
  Bal_bin_iptw_ate[i,]<-md_bin_iptw_ate
  Bal_bin_iptwst[i,]<-md_bin_iptw_ate
  Bal_bin_iptw_att[i,]<-md_bin_iptw_att
  Bal_bin_mw_ml[i,]<-md_bin_mw_ml
  Bal_bin_mw_gmm[i,]<-md_bin_mw_gmm
  Bal_bin_overlap[i,]<-md_bin_overlap
  Bal_bin_raw[i,]<-md_bin_raw
  
  #continuous covariates
  #for the raw data all weights are set to 1
  md_con_raw<-stand_diff_con(rep(1,length(L$overlap)))
  md_con_iptw_ate<-stand_diff_con(L$iptwate)
  md_con_iptw_att<-stand_diff_con(L$iptwatt)
  md_con_cbpsj_ate<-stand_diff_con(cbps_ate_just)
  md_con_cbpso_ate<-stand_diff_con(cbps_ate_over)
  md_con_cbpsj_att<-stand_diff_con(cbps_att_just)
  md_con_cbpso_att<-stand_diff_con(cbps_att_over)
  if(class(ebal_ate)!="try-error"){
    md_con_ebal_ate<-stand_diff_con(ebal_ate)
    Bal_con_ebal_ate[i,]<-md_con_ebal_ate
  }
  if(class(ebal_att)!="try-error"){
    md_con_ebal_att<-stand_diff_con(ebal_att)
    Bal_con_ebal_att[i,]<-md_con_ebal_att
  }
  if(class(vabal_ate)!="try-error"){
    md_con_vabal_ate<-stand_diff_con(vabal_ate)
    Bal_con_vabal_ate[i,]<-md_con_vabal_ate
  }
  if(class(vabal_att)!="try-error"){
    md_con_vabal_att<-stand_diff_con(vabal_att)
    Bal_con_vabal_att[i,]<-md_con_vabal_att
  }
  md_con_overlap<-stand_diff_con(L$overlap)
  md_con_mw_ml<-stand_diff_con(L$mw)
  md_con_mw_gmm<-stand_diff_con(weights_gmmw)
  
  Bal_con_cbpsj_ate[i,]<-md_con_cbpsj_ate
  Bal_con_cbpso_ate[i,]<-md_con_cbpso_ate
  Bal_con_cbpsj_att[i,]<-md_con_cbpsj_att
  Bal_con_cbpso_att[i,]<-md_con_cbpso_att
  Bal_con_iptw_ate[i,]<-md_con_iptw_ate
  Bal_con_iptwst[i,]<-md_con_iptw_ate
  Bal_con_iptw_att[i,]<-md_con_iptw_att
  Bal_con_mw_ml[i,]<-md_con_mw_ml
  Bal_con_mw_gmm[i,]<-md_con_mw_gmm
  Bal_con_overlap[i,]<-md_con_overlap
  Bal_con_raw[i,]<-md_con_raw
  Tr_prev[i,]<-table(tr)
}

#return a large list with all information regarding, balance, ess, cox regression estimators, robust standard error
#and corrected robust standard errors
Result_list<-list("bc_raw"=Bal_con_raw,"bc_overlap"=Bal_con_overlap,"bc_mw_gmm"=Bal_con_mw_gmm,"bc_mw_ml"=Bal_con_mw_ml,"bc_cbpsj_ate"=Bal_con_cbpsj_ate,"bc_cbpsj_att"=Bal_con_cbpsj_att,
                  "bc_cbpso_ate"=Bal_con_cbpso_ate,"bc_cbpso_att"=Bal_con_cbpso_att,"bc_ebal_ate"=Bal_con_ebal_ate,"bc_ebal_att"=Bal_con_ebal_att,"bc_vabal_ate"=Bal_con_vabal_ate,"bc_vabal_att"=Bal_con_vabal_att,
                  "bc_iptw_ate"=Bal_con_iptw_ate,"bc_iptw_att"=Bal_con_iptw_att,"bc_iptwst"=Bal_con_iptwst,
                  "bb_raw"=Bal_bin_raw,"bb_overlap"=Bal_bin_overlap,"bb_mw_gmm"=Bal_bin_mw_gmm,"bb_mw_ml"=Bal_bin_mw_ml,"bb_cbpsj_ate"=Bal_bin_cbpsj_ate,"bb_cbpsj_att"=Bal_bin_cbpsj_att,
                  "bb_cbpso_ate"=Bal_bin_cbpso_ate,"bb_cbpso_att"=Bal_bin_cbpso_att,"bb_ebal_ate"=Bal_bin_ebal_ate,"bb_ebal_att"=Bal_bin_ebal_att,"bb_vabal_ate"=Bal_bin_vabal_ate,"bb_vabal_att"=Bal_bin_vabal_att,
                  "bb_iptw_ate"=Bal_bin_iptw_ate,"bb_iptw_att"=Bal_bin_iptw_att,"bb_iptwst"=Bal_bin_iptwst,
                  "ESS_0"=ESS_result_matrix_0,"ESS_1"=ESS_result_matrix_1,"Cox_result"=Cox_result_matrix,"Rob_standard_error"=Rob_standard_error,"Corrected_Rob_standard_error"=Corrected_rob_standard_error,
                  "prevalence"=Tr_prev)

setwd(outname)
outputname<-paste("first_result_binary","_sample_size_",sasi,"_bias_",bias,"_bias_treatment_",beta_tr,"_prop_mod_",prop_mod,"_scenario_",scenario,"_k_",k,".RData",sep="")
save(Result_list,file=outputname)

