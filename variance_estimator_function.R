var_est<-function(data,indA,indTime,indStatus,indX,ties="breslow"){
  dat = data
  n = nrow(dat)
  dat$id = 1:n
  dat$A = dat[, indA]
  dat$time = dat[, indTime]
  dat$status = dat[, indStatus]
  nX = length(indX) + 1
  covX0 = dat[, indX]
  A = dat$A
  psmd = glm(A ~ ., family = "binomial", data = as.data.frame(cbind(A, 
                                                                    covX0)))
  psfit = predict(psmd, type = "response")
  covX0<-model.matrix(psmd)
  nX<-dim(covX0)[2]
  covX0<-covX0[,-1]
  dat$wt = dat$A/psfit + (1 - dat$A)/(1 - psfit)
  w<-dat$wt
  prevA = mean(dat$A)
  fit <- coxph(Surv(time, status) ~ A + cluster(id), weights = dat$wt, 
               data = dat, ties = ties)
  logHR = fit$coefficients
  covX<-cbind(1,covX0)
  delta<-dat$status
  S0<-S1<-rep(NA,nrow(data))
  #de contains the derivation of the propensity score estimation
  
  #derivation for the coefficient beta:
  #################
  dterm1<-dterm2<-dterm3<-dS1S0<-dS1<-dS0<-dw<-de<-matrix(NA,nrow(data),nX)
  
  for(i in 1:nrow(data)){
    de[i,]<-as.numeric(psfit[i]*(1-psfit[i])*covX[i,])
    dw[i,]<-dat$A[i]*(-1)/psfit[i]^2*de[i,]+(1-dat$A[i])*1/(1-psfit[i])^2*de[i,]
  }
  for(i in 1:nrow(data)){
    indr<-which(dat$time>=dat$time[i]& dat$status[i]==1)
    if(length(indr)==0) next
    S0[i]<-sum(dat$wt[indr]*exp(A[indr]*logHR))
    S1[i]<-sum(dat$wt[indr]*exp(A[indr]*logHR)*A[indr])
    if(length(indr)>1){
      dS1[i,]<-colSums(dw[indr,]*exp(A[indr]*logHR)*A[indr])
      dS0[i,]<-colSums(dw[indr,]*exp(A[indr]*logHR))
    }
    else{
      dS1[i,]<-dw[indr,]*exp(A[indr]*logHR)*A[indr]
      dS0[i,]<-dw[indr,]*exp(A[indr]*logHR)
    }
  }
  #dS1/dS0 derivation of S1/S0
  for(i in 1:nrow(data)){
    dS1S0[i,]<-(dS1[i,]*S0[i]-S1[i]*dS0[i,])/S0[i]^2
    dterm1[i,]<-dat$wt[i]*delta[i]*(-dS1S0[i,])+dw[i,]*delta[i]*(A[i]-S1[i]/S0[i])
  }
  #derivation of the second term product rule and
  #quotient rule
  term21<-rep(NA,nrow(data))
  term22<-rep(NA,nrow(data))
  for(i in 1:nrow(data)){
    indr<-which(data$time<=data$time[i]& dat$status[i]==1)
    term21[i]<-dat$wt[i]*A[i]*exp(A[i]*logHR)
    term22[i]<-sum(delta[indr]*dat$wt[indr]/S0[indr])
  }
  dterm2<-dterm21<-dterm22<-matrix(NA,nrow(data),nX)
  for(i in 1:nrow(data)){
    indr<-which(dat$time<=dat$time[i]& dat$status[i]==1)
    dterm21[i,]<-dw[i,]*A[i]*exp(A[i]*logHR)
    if(length(indr)==0) next
    if(length(indr)>1){
      dterm22[i,]<-colSums((dw[indr,]*delta[indr]*S0[indr]-
                              dS0[indr,]*dat$wt[indr]*delta[indr])/S0[indr]^2,na.rm=T)
    }
    else{
      dterm22[i,]<-(dw[indr,]*delta[indr]*S0[indr]-
                      dS0[indr,]*dat$wt[indr]*delta[indr])/S0[indr]^2
    }
    dterm2[i,]<-term21[i]*dterm22[i,]+dterm21[i,]*term22[i]
  }
  #term 3 we have got a product rule and a quotient rule
  #with a product rule in the quotient rule
  #term31*summe(term32)
  #term32=term321*term322/term333
  term31<-term32<-term33<-rep(NA,nrow(data))
  for(i in 1:nrow(data)){
    indr<-which(dat$time<=dat$time[i]& dat$status[i]==1)
    term31[i]<-dat$wt[i]*exp(A[i]*logHR)
    term32[i]<-sum(delta[indr]*dat$wt[indr]*S1[indr]/S0[indr]^2)
  }
  dterm3<-dterm31<-dterm32<-matrix(NA,nrow(data),nX)
  for(i in 1:nrow(data)){
    indr<-which(dat$time<=dat$time[i]& dat$status[i]==1)
    dterm31[i,]<-dw[i,]*exp(A[i]*logHR)
    if(length(indr)==0) next
    if(length(indr)>1){
      dterm32[i,]<-colSums(((dw[indr,]*delta[indr]*S1[indr]+
                               dat$wt[indr]*delta[indr]*dS1[indr,])*S0[indr]^2-
                              2*dS0[indr,]*dat$wt[indr]*delta[indr]*S1[indr])/S0[indr]^4)
    }
    else{
      dterm32[i,]<-((dw[indr,]*delta[indr]*S1[indr]+
                       dat$wt[indr]*delta[indr]*dS1[indr,])*S0[indr]^2-
                      2*dS0[indr,]*dat$wt[indr]*delta[indr]*S1[indr])/S0[indr]^4
    }
    dterm3[i,]<-term31[i]*dterm32[i,]+dterm31[i,]*term32[i]
  }
  ##################
  #derivation for the coefficient theta
  dS1theta<-S1
  dS0theta<-S1
  dtheta<-rep(NA,nrow(dat))
  for(i in 1:nrow(dat)){
    indr<-which(dat$time<=dat$time[i]& dat$status[i]==1)
    dtheta[i]<-dat$wt[i]*delta[i]*((-S1[i]*S0[i]+S1[i]^2)/S0[i]^2)
  }

  ###################
  #derivation of logistic regression
  covX<-as.matrix(covX)
  dml<-matrix(NA,nX,nX)
  sumsquare <- function(u) {
    u %*% t(u)
  }
  
  A22mat1 = apply(covX, 1, sumsquare) %*% (psfit * (1 - psfit))
  for(i in 1:ncol(covX)){
    if(i==1){
      A22mat<-colSums(dw*covX[,i]*A-(1-A)*covX[,i]*dw)
    }
    else{
      A22mat<-rbind(A22mat,colSums(dw*covX[,i]*A-(1-A)*covX[,i]*dw))
    }
  }
  
  dml = matrix(apply(A22mat1, 1, sum), nX, nX)

  
  
  #################
  #calculate 
  #pi und eta
  pi2<-pi<-matrix(NA,nrow(dat),nX)
  eta<-rep(NA,nrow(dat))
  for(i in 1:nrow(dat)){
    indr<-which(dat$time<=dat$time[i]& dat$status==1)
    pi[i,]<-(A[i]-psmd$fitted.values[i])*covX[i,]
    pi2[i,]<-dat$A[i]*covX[i,]*dat$wt[i]-(1-dat$A[i])*covX[i,]*dat$wt[i]
    eta[i]<-sum(dat$wt[i]*delta[i]*(A[i]-S1[i]/S0[i]),na.rm=T)-
      dat$wt[i]*A[i]*exp(A[i]*logHR)*sum(delta[indr]*dat$wt[indr]/S0[indr],na.rm=T)+
      dat$wt[i]*exp(A[i]*logHR)*sum(delta[indr]*dat$wt[indr]*S1[indr]/S0[indr]^2,na.rm=T)
  }
  ###############
  #variance calculation
  Omega<-t(cbind(eta,pi))
  sumsquare <- function(u) {
    u %*% t(u)
  }
  oot = apply(t(Omega), 1, sumsquare)
  B = matrix(apply(oot, 1, sum), nX + 1, nX + 1)
  AA<-rbind(c(-sum(dtheta,na.rm=T),-colSums(dterm1,na.rm=T)),cbind(0,dml))
  var_est<-solve(AA)%*%B%*%t(solve(AA))
  stderr_est<-sqrt(var_est[1,1])
  return(stderr_est)
}
var_est_overlap<-function(data,indA,indTime,indStatus,indX,ties="breslow"){
  dat = data
  n = nrow(dat)
  dat$id = 1:n
  dat$A = dat[, indA]
  dat$time = dat[, indTime]
  dat$status = dat[, indStatus]
  nX = length(indX) + 1
  covX0 = dat[, indX]
  A = dat$A
  psmd = glm(A ~ ., family = "binomial", data = as.data.frame(cbind(A, 
                                                                    covX0)))
  psfit = predict(psmd, type = "response")
  covX0<-model.matrix(psmd)
  nX<-dim(covX0)[2]
  covX0<-covX0[,-1]
  dat$wt = dat$A *(1-psfit) + (1-dat$A)*psfit
  w<-dat$wt
  prevA = mean(dat$A)
  fit <- coxph(Surv(time, status) ~ A + cluster(id), weights = dat$wt, 
               data = dat, ties = ties)
  logHR = fit$coefficients
  covX<-cbind(1,covX0)
  delta<-dat$status
  S0<-S1<-rep(NA,nrow(data))
  #de contains the derivation of the propensity score estimation
  
  #derivation for the coefficient beta:
  #################
  dterm1<-dterm2<-dterm3<-dS1S0<-dS1<-dS0<-dw<-de<-matrix(NA,nrow(data),nX)
  
  for(i in 1:nrow(data)){
    de[i,]<-as.numeric(psfit[i]*(1-psfit[i])*covX[i,])
    dw[i,]<-(-2)*A[i]*de[i,]+de[i,]
  }
  for(i in 1:nrow(data)){
    indr<-which(dat$time>=dat$time[i]& dat$status[i]==1)
    if(length(indr)==0) next
    S0[i]<-sum(dat$wt[indr]*exp(A[indr]*logHR))
    S1[i]<-sum(dat$wt[indr]*exp(A[indr]*logHR)*A[indr])
    if(length(indr)>1){
      dS1[i,]<-colSums(dw[indr,]*exp(A[indr]*logHR)*A[indr])
      dS0[i,]<-colSums(dw[indr,]*exp(A[indr]*logHR))
    }
    else{
      dS1[i,]<-dw[indr,]*exp(A[indr]*logHR)*A[indr]
      dS0[i,]<-dw[indr,]*exp(A[indr]*logHR)
    }
  }
  #dS1/dS0 derivation of S1/S0
  for(i in 1:nrow(data)){
    dS1S0[i,]<-(dS1[i,]*S0[i]-S1[i]*dS0[i,])/S0[i]^2
    dterm1[i,]<-dat$wt[i]*delta[i]*(-dS1S0[i,])+dw[i,]*delta[i]*(A[i]-S1[i]/S0[i])
  }
  ##################
  #derivation for the coefficient theta
  dS1theta<-S1
  dS0theta<-S1
  dtheta<-rep(NA,nrow(dat))
  for(i in 1:nrow(dat)){
    dtheta[i]<-dat$wt[i]*delta[i]*((-S1[i]*S0[i]+S1[i]^2)/S0[i]^2)
  }

  ###################
  #derivation of logistic regression
  covX<-as.matrix(covX)
  dml<-matrix(NA,nX,nX)
  sumsquare <- function(u) {
    u %*% t(u)
  }
  
  A22mat1 = apply(covX, 1, sumsquare) %*% (psfit * (1 - psfit))
  for(i in 1:ncol(covX)){
    if(i==1){
      A22mat<-colSums(dw*covX[,i]*A-(1-A)*covX[,i]*dw)
    }
    else{
      A22mat<-rbind(A22mat,colSums(dw*covX[,i]*A-(1-A)*covX[,i]*dw))
    }
  }
  
  dml = matrix(apply(A22mat1, 1, sum), nX, nX)
  
  
  #################
  #calculate 
  #pi und eta
  pi2<-pi<-matrix(NA,nrow(dat),nX)
  eta<-rep(NA,nrow(dat))
  for(i in 1:nrow(dat)){
    indr<-which(dat$time<=dat$time[i] &dat$status==1)
    pi[i,]<-(A[i]-psmd$fitted.values[i])*covX[i,]
    pi2[i,]<-dat$A[i]*covX[i,]*dat$wt[i]-(1-dat$A[i])*covX[i,]*dat$wt[i]
    eta[i]<-sum(dat$wt[i]*delta[i]*(A[i]-S1[i]/S0[i]),na.rm=T)-
      dat$wt[i]*A[i]*exp(A[i]*logHR)*sum(delta[indr]*dat$wt[indr]/S0[indr])+
      dat$wt[i]*exp(A[i]*logHR)*sum(delta[indr]*dat$wt[indr]*S1[indr]/S0[indr]^2)
  }
  ###############
  #variance calcualtion
  Omega<-t(cbind(eta,pi))
  sumsquare <- function(u) {
    u %*% t(u)
  }
  oot = apply(t(Omega), 1, sumsquare)
  B = matrix(apply(oot, 1, sum), nX + 1, nX + 1)
  AA<-rbind(c(-sum(dtheta,na.rm=T),-colSums(dterm1,na.rm=T)),cbind(0,dml))
  var_est<-solve(AA)%*%B%*%t(solve(AA))
  stderr_est<-sqrt(var_est[1,1])
  return(stderr_est)
}
var_est_iptw_att<-function(data,indA,indTime,indStatus,indX,ties="breslow"){
  dat = data
  n = nrow(dat)
  dat$id = 1:n
  dat$A = dat[, indA]
  dat$time = dat[, indTime]
  dat$status = dat[, indStatus]
  nX = length(indX) + 1
  covX0 = dat[, indX]
  A = dat$A
  psmd = glm(A ~ ., family = "binomial", data = as.data.frame(cbind(A,covX0)))
  psfit = predict(psmd, type = "response")
  covX0<-model.matrix(psmd)
  nX<-dim(covX0)[2]
  covX0<-covX0[,-1]
  dat$wt = dat$A + (1 - dat$A)*psfit/(1 - psfit)
  w<-dat$wt
  prevA = mean(dat$A)
  fit <- coxph(Surv(time, status) ~ A + cluster(id), weights = dat$wt, 
               data = dat, ties = ties)
  logHR = fit$coefficients
  covX<-cbind(1,covX0)
  delta<-dat$status
  S0<-S1<-rep(NA,nrow(data))
  
  #de contains the derivation of the propensity score estimation
  
  #derivation for the coefficient beta:
  #################
  dterm1<-dS1S0<-dS1<-dS0<-dw<-de<-matrix(NA,nrow(data),nX)
  
  for(i in 1:nrow(data)){
    de[i,]<-as.numeric(psfit[i]*(1-psfit[i])*covX[i,])
    dw[i,]<-(1-dat$A[i])*(psfit[i]/(1-psfit[i])^2*de[i,]+de[i,]/(1-psfit[i]))
  }
  for(i in 1:nrow(data)){
    indr<-which(dat$time>=dat$time[i])
    S0[i]<-sum(dat$wt[indr]*exp(A[indr]*logHR))
    S1[i]<-sum(dat$wt[indr]*exp(A[indr]*logHR)*A[indr])
    if(length(indr)>1){
      dS1[i,]<-colSums(dw[indr,]*exp(A[indr]*logHR)*A[indr])
      dS0[i,]<-colSums(dw[indr,]*exp(A[indr]*logHR))
    }
    else{
      dS1[i,]<-dw[indr,]*exp(A[indr]*logHR)*A[indr]
      dS0[i,]<-dw[indr,]*exp(A[indr]*logHR)
    }
  }
  #dS1/dS0 derivation of S1/S0
  for(i in 1:nrow(data)){
    dS1S0[i,]<-(dS1[i,]*S0[i]-S1[i]*dS0[i,])/S0[i]^2
    dterm1[i,]<-dat$wt[i]*delta[i]*(-dS1S0[i,])+dw[i,]*delta[i]*(A[i]-S1[i]/S0[i])
  }
  ##################
  #derivation for the coefficient theta
  dS1theta<-S1
  dS0theta<-S1
  dtheta<-rep(NA,nrow(dat))
  for(i in 1:nrow(dat)){
    indr<-which(dat$time<=dat$time[i])
    dtheta[i]<-dat$wt[i]*delta[i]*((-S1[i]*S0[i]+S1[i]^2)/S0[i]^2)
  }

  ###################
  #derivation of logistic regression
  covX<-as.matrix(covX)
  dml<-matrix(NA,nX,nX)
  sumsquare <- function(u) {
    u %*% t(u)
  }
  
  A22mat1 = apply(covX, 1, sumsquare) %*% (psfit * (1 - psfit))
  for(i in 1:ncol(covX)){
    if(i==1){
      A22mat<-colSums(dw*covX[,i]*A-(1-A)*covX[,i]*dw)
    }
    else{
      A22mat<-rbind(A22mat,colSums(dw*covX[,i]*A-(1-A)*covX[,i]*dw))
    }
  }
  
  dml = matrix(apply(A22mat1, 1, sum), nX, nX)
  
  
  #################
  #calculate 
  #pi und eta
  pi2<-pi<-matrix(NA,nrow(dat),nX)
  eta<-rep(NA,nrow(dat))
  for(i in 1:nrow(dat)){
    indr<-which(dat$time<=dat$time[i] &dat$status==1)
    pi[i,]<-(A[i]-psmd$fitted.values[i])*covX[i,]
    pi2[i,]<-dat$A[i]*covX[i,]*dat$wt[i]-(1-dat$A[i])*covX[i,]*dat$wt[i]
    eta[i]<-sum(dat$wt[i]*delta[i]*(A[i]-S1[i]/S0[i]),na.rm=T)-
      dat$wt[i]*A[i]*exp(A[i]*logHR)*sum(delta[indr]*dat$wt[indr]/S0[indr])+
      dat$wt[i]*exp(A[i]*logHR)*sum(delta[indr]*dat$wt[indr]*S1[indr]/S0[indr]^2)
  }
  
  ###############
  #variance calcualtion
  Omega<-t(cbind(eta,pi))
  sumsquare <- function(u) {
    u %*% t(u)
  }
  oot = apply(t(Omega), 1, sumsquare)
  B = matrix(apply(oot, 1, sum), nX + 1, nX + 1)
  AA<-rbind(c(-sum(dtheta,na.rm=T),-colSums(dterm1,na.rm=T)),cbind(0,dml))
  var_est<-solve(AA)%*%B%*%t(solve(AA))
  stderr_est<-sqrt(var_est[1,1])
  return(stderr_est)
}
var_est_mw_gmm<-function(data,indA,indTime,indStatus,indX,ties="breslow",psfit,coefficients,weights){
  dat = data
  n = nrow(dat)
  dat$id = 1:n
  dat$A = dat[, indA]
  dat$time = dat[, indTime]
  dat$status = dat[, indStatus]
  nX = length(indX)
  covX0 = dat[, indX]
  A = dat$A
  if(any(psfit)==1){
    psfit[which(psfit==1)]<-0.999999
  }
  if(any(psfit)==0){
    psfit[which(psfit==0)]<-0.0000001
  }
  nX<-dim(covX0)[2]
  dat$wt = weights
  w<-dat$wt
  fit <- coxph(Surv(time, status) ~ A + cluster(id), weights = dat$wt, 
               data = dat, ties = ties)
  logHR = fit$coefficients
  covX<-covX0
  delta<-dat$status
  S0<-S1<-rep(NA,nrow(data))
  
  #de contains the derivation of the propensity score estimation
  k<-1400
  min_approx<-(-1)/k*log(1/2*(exp(-k*(psfit))+exp(-k*(1-psfit))))
  min_real<-sapply(1:length(psfit),function(z){
    return(min(c(psfit[z],1-psfit[z])))
  })
  #derivation for the coefficient beta:
  #################
  d_min_approx<-dterm1<-dS1S0<-dS1<-dS0<-dw<-de<-matrix(NA,nrow(data),nX)
  for(i in 1:nrow(data)){
    de[i,]<-as.numeric(psfit[i]*(1-psfit[i])*covX[i,])
    d_min_approx[i,]<-(-1/k)*1/(1/2*(exp(-k*(psfit[i]))+exp(-k*(1-psfit[i]))))/2*(exp(-k*psfit[i])*(-k*de[i,])+exp(-k*(1-psfit[i]))*(k*de[i,]))
    dw[i,]<-dat$A[i]*(-min_approx[i]*1/psfit[i]^2*de[i,]+d_min_approx[i,]/psfit[i])+(1-dat$A[i])*
      (min_approx[i]/(1-psfit[i])^2*de[i,]+d_min_approx[i,]/(1-psfit[i]))
  }
  for(i in 1:nrow(data)){
    indr<-which(dat$time>=dat$time[i]& dat$status[i]==1)
    if(length(indr)==0) next
    S0[i]<-sum(dat$wt[indr]*exp(A[indr]*logHR))
    S1[i]<-sum(dat$wt[indr]*exp(A[indr]*logHR)*A[indr])
    if(length(indr)>1){
      dS1[i,]<-colSums(dw[indr,]*exp(A[indr]*logHR)*A[indr])
      dS0[i,]<-colSums(dw[indr,]*exp(A[indr]*logHR))
    }
    else{
      dS1[i,]<-dw[indr,]*exp(A[indr]*logHR)*A[indr]
      dS0[i,]<-dw[indr,]*exp(A[indr]*logHR)
    }
  }
  #dS1/dS0 derivation of S1/S0
  for(i in 1:nrow(data)){
    dS1S0[i,]<-(dS1[i,]*S0[i]-S1[i]*dS0[i,])/S0[i]^2
    dterm1[i,]<-dat$wt[i]*delta[i]*(-dS1S0[i,])+dw[i,]*delta[i]*(A[i]-S1[i]/S0[i])
  }
  dS1theta<-S1
  dS0theta<-S1
  dtheta<-rep(NA,nrow(dat))
  for(i in 1:nrow(dat)){
    dtheta[i]<-dat$wt[i]*delta[i]*((-S1[i]*S0[i]+S1[i]^2)/S0[i]^2)
  }

  ###################
  #derivation of logistic regression
  covX<-as.matrix(covX)
  dml<-matrix(NA,nX ,nX )
  sumsquare <- function(u) {
    u %*% t(u)
  }
  A22mat1 = apply(covX0, 1, sumsquare) %*% (psfit * (1 - psfit))
  A22 = matrix(apply(A22mat1, 1, sum), nX, nX)
  A22mat<-NULL
  for(i in 1:ncol(covX0)){
    if(i==1){
      A22mat<-colSums(dw*covX0[,i]*A-dw*covX0[,i]*(1-A))
    }
    else{
      A22mat<-rbind(A22mat,colSums(dw*covX0[,i]*A-(1-A)*covX0[,i]*dw))
    }
  }
  dml<-A22mat
  
  #################
  #calculate 
  #pi und eta
  pi<-matrix(NA,nrow(dat),nX)
  pi2<-matrix(NA,nrow(dat),nX)
  eta<-rep(NA,nrow(dat))
  for(i in 1:nrow(dat)){
    indr<-which(dat$time<=dat$time[i] &dat$status==1)
    pi[i,]<-(A[i]-psfit[i])*covX[i,]
    pi2[i,]<-dat$A[i]*covX[i,]*dat$wt[i]-(1-dat$A[i])*covX[i,]*dat$wt[i]
    eta[i]<-sum(dat$wt[i]*delta[i]*(A[i]-S1[i]/S0[i]),na.rm=T)-
      dat$wt[i]*A[i]*exp(A[i]*logHR)*sum(delta[indr]*dat$wt[indr]/S0[indr])+
      dat$wt[i]*exp(A[i]*logHR)*sum(delta[indr]*dat$wt[indr]*S1[indr]/S0[indr]^2)
  }
  
  ###############
  #variance calcualtion
  Omega<-t(cbind(eta,pi2))
  sumsquare <- function(u) {
    u %*% t(u)
  }
  oot = apply(t(Omega), 1, sumsquare)
  for(j in 1:ncol(Omega)){
    if(j==1){
      B<-matrix(Omega[,j],ncol=1)%*%matrix(Omega[,j],nrow=1)
    }
    else{
      B<-B+matrix(Omega[,j],ncol=1)%*%matrix(Omega[,j],nrow=1)
    }
  }
  AA<-rbind(c(sum(dtheta,na.rm=T),colSums(dterm1,na.rm=T)),cbind(0,dml))
  var_est<-solve(-AA)%*%B%*%solve(t(-AA))
  stderr_est<-sqrt(var_est[1,1])
  return(stderr_est)
}
var_est_mw<-function(data,indA,indTime,indStatus,indX,ties="breslow"){
  dat = data
  n = nrow(dat)
  dat$id = 1:n
  dat$A = dat[, indA]
  dat$time = dat[, indTime]
  dat$status = dat[, indStatus]
  nX = length(indX) + 1
  covX0 = dat[, indX]
  A = dat$A
  psmd = glm(A ~ ., family = "binomial", data = as.data.frame(cbind(A, 
                                                                    covX0)))
  psfit = predict(psmd, type = "response")
  covX0<-model.matrix(psmd)
  nX<-dim(covX0)[2]
  covX0<-covX0[,-1]
  k<-1400
  min_approx<-(-1)/k*log(1/2*(exp(-k*(psfit))+exp(-k*(1-psfit))))
  w_min_approx<-dat$A*min_approx/psfit+(1-dat$A)*min_approx/(1-psfit)
  dat$wt = sapply(1:length(dat$A),function(z){
    return(dat$A[z]*min(c(psfit[z],1-psfit[z]))/(psfit[z]) + (1-dat$A[z])*min(c(psfit[z],1-psfit[z]))/(1-psfit[z]))
  })
  w<-dat$wt
  prevA = mean(dat$A)
  fit <- coxph(Surv(time, status) ~ A + cluster(id), weights = dat$wt, 
               data = dat, ties = ties)
  logHR = fit$coefficients
  covX<-cbind(1,covX0)
  delta<-dat$status
  S0<-S1<-rep(NA,nrow(data))
  #de contains the derivation of the propensity score estimation
  
  #derivation for the coefficient beta:
  #################
  dterm1<-dS1S0<-dS1<-dS0<-dw<-d_min_approx<-de<-matrix(NA,nrow(data),nX)
  
  for(i in 1:nrow(data)){
    de[i,]<-as.numeric(psfit[i]*(1-psfit[i])*covX[i,])
    d_min_approx[i,]<-(-1/k)*1/(1/2*(exp(-k*(psfit[i]))+exp(-k*(1-psfit[i]))))/2*(exp(-k*psfit[i])*(-k*de[i,])+exp(-k*(1-psfit[i]))*(k*de[i,]))
    dw[i,]<-dat$A[i]*(-min_approx[i]*1/psfit[i]^2*de[i,]+d_min_approx[i,]/psfit[i])+(1-dat$A[i])*
      (min_approx[i]/(1-psfit[i])^2*de[i,]+d_min_approx[i,]/(1-psfit[i]))
  }
  for(i in 1:nrow(data)){
    indr<-which(dat$time>=dat$time[i]& dat$status[i]==1)
    if(length(indr)==0) next
    S0[i]<-sum(dat$wt[indr]*exp(A[indr]*logHR))
    S1[i]<-sum(dat$wt[indr]*exp(A[indr]*logHR)*A[indr])
    if(length(indr)>1){
      dS1[i,]<-colSums(dw[indr,]*exp(A[indr]*logHR)*A[indr])
      dS0[i,]<-colSums(dw[indr,]*exp(A[indr]*logHR))
    }
    else{
      dS1[i,]<-dw[indr,]*exp(A[indr]*logHR)*A[indr]
      dS0[i,]<-dw[indr,]*exp(A[indr]*logHR)
    }
  }
  #dS1/dS0 derivation of S1/S0
  for(i in 1:nrow(data)){
    dS1S0[i,]<-(dS1[i,]*S0[i]-S1[i]*dS0[i,])/S0[i]^2
    dterm1[i,]<-dat$wt[i]*delta[i]*(-dS1S0[i,])+dw[i,]*delta[i]*(A[i]-S1[i]/S0[i])
  }
  ##################
  #derivation for the coefficient theta
  dS1theta<-S1
  dS0theta<-S1
  dtheta<-rep(NA,nrow(dat))
  for(i in 1:nrow(dat)){
    indr<-which(dat$time<=dat$time[i])
    dtheta[i]<-dat$wt[i]*delta[i]*((-S1[i]*S0[i]+S1[i]^2)/S0[i]^2)
  }

  ###################
  #derivation of logistic regression
  covX<-as.matrix(covX)
  dml<-matrix(NA,nX,nX)
  sumsquare <- function(u) {
    u %*% t(u)
  }
  
  A22mat1 = apply(covX, 1, sumsquare) %*% (psfit * (1 - psfit))
  for(i in 1:ncol(covX)){
    if(i==1){
      A22mat<-colSums(dw*covX[,i]*A-(1-A)*covX[,i]*dw)
    }
    else{
      A22mat<-rbind(A22mat,colSums(dw*covX[,i]*A-(1-A)*covX[,i]*dw))
    }
  }
  
  dml = matrix(apply(A22mat1, 1, sum), nX, nX)

  
  
  #################
  #calculate 
  #pi und eta
  pi2<-pi<-matrix(NA,nrow(dat),nX)
  eta<-rep(NA,nrow(dat))
  for(i in 1:nrow(dat)){
    indr<-which(dat$time<=dat$time[i] &dat$status==1)
    pi[i,]<-(A[i]-psmd$fitted.values[i])*covX[i,]
    pi2[i,]<-dat$A[i]*covX[i,]*dat$wt[i]-(1-dat$A[i])*covX[i,]*dat$wt[i]
    eta[i]<-sum(dat$wt[i]*delta[i]*(A[i]-S1[i]/S0[i]),na.rm=T)-
      dat$wt[i]*A[i]*exp(A[i]*logHR)*sum(delta[indr]*dat$wt[indr]/S0[indr])+
      dat$wt[i]*exp(A[i]*logHR)*sum(delta[indr]*dat$wt[indr]*S1[indr]/S0[indr]^2)
  }
  
  ###############
  #variance calcualtion
  Omega<-t(cbind(eta,pi))
  sumsquare <- function(u) {
    u %*% t(u)
  }
  oot = apply(t(Omega), 1, sumsquare)
  B = matrix(apply(oot, 1, sum), nX + 1, nX + 1)
  AA<-rbind(c(-sum(dtheta,na.rm=T),-colSums(dterm1,na.rm=T)),cbind(0,dml))
  var_est<-solve(AA)%*%B%*%t(solve(AA))
  stderr_est<-sqrt(var_est[1,1])
  return(stderr_est)
}
var_est_iptw_st<-function(data,indA,indTime,indStatus,indX,ties="breslow"){
  dat = data
  n = nrow(dat)
  dat$id = 1:n
  dat$A = dat[, indA]
  dat$time = dat[, indTime]
  dat$status = dat[, indStatus]
  nX = length(indX) + 1
  covX0 = dat[, indX]
  A = dat$A
  psmd = glm(A ~ ., family = "binomial", data = as.data.frame(cbind(A,covX0)))
  psfit = predict(psmd, type = "response")
  covX0<-model.matrix(psmd)
  nX<-dim(covX0)[2]
  covX0<-covX0[,-1]
  pr = mean(dat$A)
  dat$wt = dat$A*pr/psfit + (1 - dat$A)*(1-pr)/(1 - psfit)
  w<-dat$wt
  w_unst<-ifelse(dat$A==1,w/pr,w/(1-pr))
  fit <- coxph(Surv(time, status) ~ A + cluster(id), weights = dat$wt, 
               data = dat, ties = ties)
  logHR = fit$coefficients
  covX<-cbind(1,covX0)
  delta<-dat$status
  S0<-S1<-rep(NA,nrow(data))
  pr<-sum(dat$A)/length(dat$A)
  #de contains the derivation of the propensity score estimation
  
  #derivation for the coefficient beta:
  #################
  dw<-de<-matrix(NA,nrow(data),nX)
  dterm1<-dS1S0<-dS1<-dS0<-matrix(NA,nrow(data),nX+1)
  dw_p<-rep(NA,nrow(dat))
  for(i in 1:nrow(data)){
    de[i,]<-as.numeric(psfit[i]*(1-psfit[i])*covX[i,])
    dw[i,]<-dat$A[i]*pr*(-1)/psfit[i]^2*de[i,]+(1-dat$A[i])*(1-pr)/(1-psfit[i])^2*de[i,]
    dw_p[i]<-dat$A[i]/psfit[i]-(1-dat$A[i])/(1-psfit[i])
  }
  dw<-cbind(dw,dw_p)
  for(i in 1:nrow(data)){
    indr<-which(dat$time>=dat$time[i]& dat$status[i]==1)
    if(length(indr)==0) next
    S0[i]<-sum(dat$wt[indr]*exp(A[indr]*logHR))
    S1[i]<-sum(dat$wt[indr]*exp(A[indr]*logHR)*A[indr])
    if(length(indr)>1){
      dS1[i,]<-colSums(dw[indr,]*exp(A[indr]*logHR)*A[indr])
      dS0[i,]<-colSums(dw[indr,]*exp(A[indr]*logHR))
    }
    else{
      dS1[i,]<-dw[indr,]*exp(A[indr]*logHR)*A[indr]
      dS0[i,]<-dw[indr,]*exp(A[indr]*logHR)
    }
  }
  #dS1/dS0 derivation of S1/S0
  for(i in 1:nrow(data)){
    dS1S0[i,]<-(dS1[i,]*S0[i]-S1[i]*dS0[i,])/S0[i]^2
    dterm1[i,]<-dat$wt[i]*delta[i]*(-dS1S0[i,])+dw[i,]*delta[i]*(A[i]-S1[i]/S0[i])
  }
  ##################
  #derivation for the coefficient theta
  dS1theta<-S1
  dS0theta<-S1
  dtheta<-rep(NA,nrow(dat))
  for(i in 1:nrow(dat)){
    dtheta[i]<-dat$wt[i]*delta[i]*((-S1[i]*S0[i]+S1[i]^2)/S0[i]^2)
  }

  ###################
  #derivation of logistic regression
  covX<-as.matrix(covX)
  dml<-matrix(NA,nX,nX)
  sumsquare <- function(u) {
    u %*% t(u)
  }
  
  A22mat1 = apply(covX, 1, sumsquare) %*% (psfit * (1 - psfit))
  for(i in 1:ncol(covX)){
    if(i==1){
      A22mat<-colSums(dw*covX[,i]*A-(1-A)*dw*covX[,i])
    }
    else{
      A22mat<-rbind(A22mat,colSums(dw*covX[,i]*A-(1-A)*dw*covX[,i]))
    }
  }
  
  dml = matrix(apply(A22mat1, 1, sum), nX, nX)
  
  
  #################
  #calculate 
  #pi und eta
  pi2<-pi<-matrix(NA,nrow(dat),nX)
  eta<-rep(NA,nrow(dat))
  for(i in 1:nrow(dat)){
    indr<-which(dat$time<=dat$time[i] &dat$status==1)
    pi[i,]<-(A[i]-psmd$fitted.values[i])*covX[i,]
    pi2[i,]<-dat$A[i]*covX[i,]*dat$wt[i]-(1-dat$A[i])*covX[i,]*dat$wt[i]
    eta[i]<-sum(dat$wt[i]*delta[i]*(A[i]-S1[i]/S0[i]),na.rm=T)-
      dat$wt[i]*A[i]*exp(A[i]*logHR)*sum(delta[indr]*dat$wt[indr]/S0[indr])+
      dat$wt[i]*exp(A[i]*logHR)*sum(delta[indr]*dat$wt[indr]*S1[indr]/S0[indr]^2)
  }
  
  diff_pr<-rep(NA,length(dat$A))
  for(i in 1:nrow(dat)){
    diff_pr[i]<-dat$A[i]-pr
  }
  ###############
  #variance calcualtion
  Omega<-t(cbind(eta,pi,diff_pr))
  sumsquare <- function(u) {
    u %*% t(u)
  }
  oot = apply(t(Omega), 1, sumsquare)
  B = matrix(apply(oot, 1, sum), nX + 2, nX + 2)
  AA<-rbind(c(-sum(dtheta,na.rm=T),-colSums(dterm1,na.rm=T)),rbind(cbind(0,dml,0),c(rep(0,nX+1),nrow(data))))
  var_est<-solve(AA)%*%B%*%t(solve(AA))
  var_est_1<-solve(AA[1,1])%*%B[1,1]%*%t(solve(AA[1,1]))
  stderr_est<-sqrt(var_est[1,1])
  return(stderr_est)
}
