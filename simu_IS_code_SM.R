##########################################################################
## Loading packages
## Install these packages if you don't yet have them.
##########################################################################
library(quantreg)
library(survival)
library(emplik)
library(BB)
library(dplyr)
library(prodlim)
library(foreach)
library(doParallel)
library(doRNG)
library(readr)

##########################################################################
## This function is used for generating polynomial basis function
##########################################################################
basis_func <- function(tseq){
  poly=cbind(1,log(tseq),sqrt(tseq),1/sqrt(tseq))
  #poly=cbind(1,log(tseq),sqrt(tseq),1/sqrt(tseq),1/tseq)
  #poly=cbind(1,1/sqrt(tseq),tseq,tseq^2)
  return(poly)}

##########################################################################
## This function is used for calculating survival probability using KM
##########################################################################
KM_func <- function(T,censor,wgt=1){ 
  deathtime <- unique(sort(T[censor[]==1]))
  nrisk <- ndeath <- rep(0,length(deathtime))
  for(i in 1:length(deathtime)){
    nrisk[i] <- sum((deathtime[i]<=T)*wgt)
    ndeath[i] <- sum((T==deathtime[i])*censor*wgt)
  }
  prodobj <- 1-ndeath/nrisk
  survp <- rep(0,length(deathtime))
  for(i in 1:length(deathtime)){survp[i]=prod(prodobj[1:i])}
  return(data.frame(cbind(deathtime,ndeath,nrisk,survp)))}

###########################################################################
## This function is required to generate simulation data
## We denote our method as "IS" and Li et al(2016)'s method as "Li," respectively.
## To compare with Li's method, data was generated the same way
###########################################################################
da_gen_S2 <- function(n, eta0, C0, tvec0, tau0){
  X <- rbinom(n,1,0.5) #time-fixed covariate
  U <- runif(n,0,1)
  lambda <- runif(n,0.5,1.5)
  T <- sqrt(-log(U)/lambda) #failure time
  eta <- rbinom(n,1,0.9)
  C <- runif(n,0,C0)*eta+4*(1-eta) #censoring time
  Y <- pmax(pmin(T,C),0.0001)
  delta <- (T<=C)*1 #failure indicator
  da <- data.frame(cbind(Y,T,delta,X))
  names(da) <- c("Y","T","delta","X")
  
  Umat <- c() 
  nvec <- Yvec <- tvec <- c()
  nvisit <- length(tvec0) #number of visits
  visit <- c()
  for(i in 1:n){
    p_come <- 0.5*(X[i]==0)+0.7*I(X[i]==1)
    come_visit <- rbinom(nvisit,1,p_come) #Whether the participant visited or not
    t0seq <- tvec0
    ni <- length(t0seq)
    Z_t_n <- sqrt(-log(1-tau0)/lambda[i]/t0seq^2+1)-1
    bseq <- sqrt(t0seq)
    Z_t <- log(Z_t_n)/bseq
    
    come_visit[t0seq>Y[i]] <- 0
    Z_t[come_visit==0]=999 #missing data
    base <- basis_func(t0seq)
    Umat_c <- cbind(base,base*X[i],base*Z_t)
    Umat <- rbind(Umat,Umat_c)
    nvec <- c(nvec,ni)
    Yvec <- c(Yvec,log(pmax(Y[i]-t0seq,0.00001)))
    tvec <- c(tvec,t0seq)
    visit <- rbind(visit,come_visit)
  }
  return(list(da=da,nvec=nvec,Umat=Umat,Yvec=Yvec,tvec=tvec,visit=visit))}

##########################################################################
## This function is used for estimation
##########################################################################
est_func_S2 <- function(obj,wt.p,ne,tpseq,tau0){
  da <- obj$da 
  Umat <- obj$Umat  
  tvec <- obj$tvec 
  nvec <- obj$nvec
  Y_reg <- obj$Yvec 
  U_reg <- obj$Umat 
  visit <- obj$visit
  
  Gest <- KM_func(da$Y,1-da$delta,wt.p) #wt.p=1 fitting
  SC <- stepfun(Gest$deathtime,c(1,Gest$survp)) #survival function for censoring time (fitting point)
  wt_reg_Li <- wt_reg_IS <- wt_reg_sd <- wt_reg_add <- come_reg <- c()
  
  for(i in 1:n0){
    ni <- nvec[i]
    index <- sum(nvec[1:(i-1)])
    if(i==1){index=0}
    U_mat_c <- matrix(Umat[(index+1):(index+ni),],ni)
    tseq <- tvec[(index+1):(index+ni)]
    come_c <- visit[i,] #visit indicator
    invw.curr <- pmax(SC(da$Y[i])/SC(tseq),0.0001)
    wt.curr <- (da$delta[i]*(da$Y[i]>tseq)/invw.curr*come_c)*wt.p[i]
    wt.curr2 <- (da$delta[i]/invw.curr)
    pseudo1 <- -apply(U_mat_c*wt.curr,2,sum)
    pseudo2 <- 2*apply(U_mat_c*(da$Y[i]>tseq)*come_c*wt.p[i]*tau0,2,sum)
    
    come_reg <- c(come_reg,come_c)
    Y_reg <- c(Y_reg,M,M)
    U_reg <- rbind(U_reg,rbind(pseudo1,pseudo2))
    wt_reg_Li <- c(wt_reg_Li,wt.curr) 
    wt_reg_add <- c(wt_reg_add,1,1)
    wt_reg_IS <- c(wt_reg_IS,wt.curr2) 
  }
  
  wt_reg_Li <- c(wt_reg_Li,wt_reg_add) #weight for Li's method
  wt_reg_IS <- c(wt_reg_IS) #weight for IS method
  come_reg_idw <- c(come_reg) #visit indicator
  
  #evaluate the quantile model at tpseq
  Li_fit <- rq.wfit(U_reg,Y_reg,weights=wt_reg_Li)$coef
  if(sum(abs(Li_fit))>=5000){Li_fit=Li_fit*NA}
  
  #Point estimators obtained using Li's method
  pbase <- basis_func(tpseq)
  base_est <- as.vector(pbase%*%Li_fit[1:ppoly])
  alpha_est <- as.vector(pbase%*%Li_fit[(ppoly+1):(2*ppoly)])
  beta_est <- as.vector(pbase%*%Li_fit[(2*ppoly+1):(3*ppoly)])
  Li_est <- cbind(base_est,alpha_est,beta_est) 
  
  #Point estimators obtained using IS method
  objectF <- function(beta){
    beta <- as.matrix(beta)
    result <- (1/n0)*t(U*come) %*% (W*(pnorm((U%*%beta-logY)/sqrt(diag(U %*% H %*% t(U)))))-tau0)
    result <- as.vector(result)
  }
  
  U <- obj$Umat 
  come <-come_reg_idw 
  W <- wt_reg_IS #censoring weight
  logY <- obj$Yvec #log of observed failure time Y
  nc <- length(Li_fit)
  H <- diag(1/(12*n0), nc, nc)
  
  betastart <- c(Li_fit)
  IS_fit <- dfsane(betastart, objectF, control = list(tol=5*1e-03))
  idw_fit <- IS_fit$par
  
  tryCatch({
    if (IS_fit$convergence == 0){
      idw_pbase <- basis_func(tpseq)
      idw_base_est <- as.vector(idw_pbase%*%idw_fit[1:ppoly])
      idw_alpha_est <- as.vector(idw_pbase%*%idw_fit[(ppoly+1):(2*ppoly)])
      idw_beta_est <- as.vector(idw_pbase%*%idw_fit[(2*ppoly+1):(3*ppoly)])
      IS_est <- cbind(idw_base_est,idw_alpha_est,idw_beta_est)
    } else {
      IS_est <- c(NA,NA,NA)
    }
  },error=function(e){
    IS_est <- c(NA,NA,NA)
  })
  
  ##########################################################################
  ## This function is used for covariance estimation using IS method
  ##########################################################################
  
  # pre-calculation
  UHU_d <- diag(U %*% H %*% t(U))
  Ubeta_Y <- U %*% idw_fit-logY
  
  result_ismb <- c()
  for (j in 1:ne){
    
    eta_p1 <- rexp(n0,1)
    eta_p <- rep(eta_p1, each=12)
    
    Gest_sd <- KM_func(da$Y,1-da$delta,eta_p1)
    SC_sd <- stepfun(Gest_sd$deathtime,c(1,Gest_sd$survp))
    
    wt_reg_idw_sd <- c()
    
    for(m in 1:n0){
      ni <- nvec[m]
      index <- sum(nvec[1:(m-1)])
      if(m==1){index=0}
      tseq <- tvec[(index+1):(index+ni)]
      
      invw_c_sd <- pmax(SC_sd(da$Y[m])/SC_sd(tseq),0.0001)
      wt_c_sd <- (da$delta[m]/invw_c_sd)
      wt_reg_idw_sd <- c(wt_reg_idw_sd,wt_c_sd)
    }
    result <- (1/n0)*t(eta_p*U*come) %*% 
      {wt_reg_idw_sd*(pnorm(Ubeta_Y/sqrt(UHU_d)))-tau0}
    result_ismb <- cbind(result_ismb,result)
  }
  
  v <- cov(t(result_ismb))
  a_beta <- (1/n0)*t(U*come*W*as.vector(dnorm(Ubeta_Y/sqrt(UHU_d))))%*%(U/sqrt(UHU_d))
  
  inva_beta <- solve(a_beta) 
  sigma <- t(inva_beta) %*% v %*% inva_beta
  base_sigma <- sigma[1:ppoly,1:ppoly]
  X_sigma <- sigma[(ppoly+1):2*ppoly,(ppoly+1):2*ppoly]
  Z_t_sigma <- sigma[(2*ppoly+1):3*ppoly,(2*ppoly+1):3*ppoly]
  
  pbase <- basis_func(tpseq)
  v_base <- pbase%*%base_sigma%*%t(pbase)
  v_X <- pbase%*%X_sigma%*%t(pbase)
  v_Z_t <- pbase%*%Z_t_sigma%*%t(pbase)
  sd_idw <- cbind(sqrt(v_base),sqrt(v_X),sqrt(v_Z_t))
  return(list(est=Li_est,         #point estimator by Li's method
              idw_est=IS_est,     #point estimator by IS method
              sd_idw=sd_idw       #covariance estimator by IS method
  ))
}

##########################################################################
## Fitting
##########################################################################
ppoly <- 4 #number of basis functions
M <- 1E8  
n0 <- 400 #sample size
tvec0 <- c(seq(0.1,0.25,0.05),seq(0.3,1,0.1)) #visit time in the data
tau0 <- 0.25 #quantile level
exp0 <- 1.5; eta0 <- 0.9; C0 <- 4;

##########################################################################
## Dataset for simulation
##########################################################################
obj <- da_gen_S2(n0, eta0, C0, tvec0, tau0)

##########################################################################
## Estimation
# point estimator: est(Li), idw_est(IS)
# variance estimator: sd_idw(IS)
########################################################################## 
tpseq <- c(0.1) #time to estimate the regression coefficient
ne <- 200 #The repeated time for estimating variance using IS method
fit <- est_func_S2(obj, rep(1,n0), ne, tpseq, tau0)

##########################################################################
## Conduct the simulation 1000 times.
##########################################################################
n1 <- 1000 #number of simulation 
n0 <- 800 #sample size

numCores <- detectCores() - 1
myCluster <- makeCluster(numCores)
registerDoParallel(myCluster)

system.time(
  sim_func_S2 <-
    foreach (i=1:n1, .options.RNG = 1234, .packages = c('survival','BB','quantreg','emplik','prodlim')) %dorng% {
      
      obj <- da_gen_S2(n0, eta0, C0, tvec0, tau0)
      fit <- est_func_S2(obj, rep(1,n0), ne, tpseq, tau0)
      
      my_list <- list("li" = fit$est,
                      "idw" = fit$idw_est, 
                      "idw_sd" = fit$sd_idw)
      
      return(my_list)
    })
stopCluster(myCluster)
closeAllConnections()

## initialization
par_est_li <- par_est_is <- sd_est_is <- matrix(rep(0,n1*3), nrow=n1)

## Saves the results in csv format to a specified directory
for(i in 1:n1){
  par_est_li[i,1:3] <- sim_func_S2[[i]]$li
  par_est_is[i,1:3] <- sim_func_S2[[i]]$idw
  sd_est_is[i,1:3] <- sim_func_S2[[i]]$idw_sd
}
res_sim <- cbind(par_est_li,par_est_is,sd_est_is)
write.csv(res_sim, file=paste0("S2_sim re H12_0213_", n1,"_sample_",n0,"_tau",tau0, "_t",tpseq,"_ne", ne,".csv"))

