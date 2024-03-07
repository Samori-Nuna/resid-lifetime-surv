rm(list=ls())
require(msm)
library(JM)
library(mvtnorm)
library(MASS)
library(splines)
library(fpca)
library(quantreg)
library(survival)
library(emplik)
library(BB)
library(dplyr)
library(prodlim)
library(foreach)
library(doParallel)
library(doRNG)
library(diagonals)


##########################################################################
## This function is used for calculating survival probability using KM
##########################################################################
SC.func=function(T,censor,wgt=1){
  deathtime=unique(sort(T[censor[]==1]))
  nrisk=ndeath=rep(0,length(deathtime))
  for(i in 1:length(deathtime)){
    nrisk[i]=sum((deathtime[i]<=T)*wgt)
    ndeath[i]=sum((T==deathtime[i])*censor*wgt)
  }
  prodobj=1-ndeath/nrisk
  survp=rep(0,length(deathtime))
  for(i in 1:length(deathtime)){survp[i]=prod(prodobj[1:i])}
  return(data.frame(cbind(deathtime,ndeath,nrisk,survp)))}

###########################################################################
## This function is used to generate simulation dataset 
###########################################################################

da.genrt=function(N){
  n=N ##number of sample size
  K=21
  # parameters for the f(t)
  betas <- cbind(-2,0.2)
  sigma.y <- 0.3 # measurement error standard deviation
  D=matrix(c(0.1,0.04,0.04,0.1),nrow=2) # random effects bi = (bi1,bi2)^T ~ N(0,D])
  bb <- mvrnorm(n, rep(0, nrow(D)), D)
  randomeffect<-bb
  
  # parameters for the survival model
  betasm=0.4
  sigma.t=1.2
  
  #parameters for censoring
  censoring=0.12 ## small -> censoring rate small
  Xcoef=1
  X.id=runif(n,0,1)
  X=rep(X.id,each=K)
  
  # design matrices for the longitudinal measurement model
  times <- c(replicate(n, c(0,seq(0.5,10,0.5))))
  times[which(times!=0)] <- times[which(times!=0)]+ rnorm(n*(K-1),0,0.1)
  Ai=betas[,1]
  Bi=betas[,2]
  
  id <- rep(1:n, each = K)
  eta.y <-as.vector( Ai+Bi*times+randomeffect[id,1]+randomeffect[id,2]*times)
  y <- rnorm(n * K, eta.y, sigma.y)
  
  # simulate event times
  invS <- function (t, u, i) 
  {
    h <- function (s) 
    {     
      exp(betasm*(Ai+Bi*s+randomeffect[i,1]+randomeffect[i,2]*s)+Xcoef*X.id[i])
    }
    integrate(h, lower = 0, upper = t)$value - u
  }
  u <- exp(rnorm(n,0,sigma.t))
  trueTimes <- numeric(n)
  for (i in 1:n) 
  {
    Up <- 50
    tries <- 5
    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    while(inherits(Root, "try-error") && tries > 0) 
    {
      tries <- tries - 1
      Up <- Up + 50
      Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    }
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 50
  }
  na.ind <- !is.na(trueTimes)
  trueTimes <- trueTimes[na.ind]
  n <- length(trueTimes)
  trueTimes <- rep(trueTimes,each=K)
  
  #Ctimes <- -log(runif(n, 0, 1))/(censoring*exp(X.id[na.ind]*Xcoef))
  Ctimes <- runif(n,0,18)
  Ctimes<-rep(Ctimes,each=K)
  
  Time=pmin(L,trueTimes)
  Time=pmin(Time,Ctimes)
  event=as.numeric(pmin(L,trueTimes) <= Ctimes)
  
  #Time=pmin(trueTimes,Ctimes)
  #Time=pmin(L,Time)
  #event=as.numeric(trueTimes <= Ctimes)
  
  long.na.ind <- rep(na.ind, each = K)
  y <- y[long.na.ind]
  times<-times[long.na.ind] 
  X<-X[long.na.ind]
  
  ind <- times <= Time
  id <- id[long.na.ind][ind]
  y <- y[ind]
  X <- X[ind]
  times=times[ind]
  Time=Time[ind]
  trueTimes=trueTimes[ind]
  Ctimes=Ctimes[ind]
  event=event[ind]
  data=cbind(id,y,X,event,trueTimes,Ctimes,Time,times)
  data.id=data [!duplicated(data[,1]),]
  data.id=as.data.frame(data.id)
  plot(survfit(Surv(Time,event)~1,data=data.id))
  max(data.id[,7])
  sum(data.id[,4])/dim(data.id)[1]
  return(list(data=data,X.id=X.id,Xcoef=Xcoef,Ai=Ai,Bi=Bi,bb=bb,betasm=betasm,sigma.t=sigma.t))
}

###########################################################################
## This function is used to calculate AUC
###########################################################################

AUC<- function(data, qt,tau){
  fz=0
  fm=0
  data=as.data.frame(data)
  data.id=data[!duplicated(data[,1]),]
  wt.p=rep(1,dim(data.id)[1])
  Gest=SC.func(data.id$Time,1-data.id$event,wt.p)
  SC=stepfun(Gest$deathtime,c(1,Gest$survp)) #survival function for censoring time
  for(i in qt[,1]){
    for(j in qt[,1]){
      fz=sum(fz,data.id[data.id$id==i,]$event*SC(data.id[data.id$id==i,]$Time)^(-2)*I(data.id[data.id$id==i,]$Time<data.id[data.id$id==j,]$Time&data.id[data.id$id==i,]$Time<tau)*I( qt[qt[,1]==i,2]<qt[qt[,1]==j,2]))
      fm=sum(fm,data.id[data.id$id==i,]$event*SC(data.id[data.id$id==i,]$Time)^(-2)*I(data.id[data.id$id==i,]$Time<data.id[data.id$id==j,]$Time&data.id[data.id$id==i,]$Time<tau))
    }
  }
  return(fz/fm)
}

###########################################################################
## This function is used to fit our proposed FPCA-based model
###########################################################################
fpca.fit<- function(data){
  fdata=as.data.frame(data)
  fdata=fdata[which(fdata$times<=10 & fdata$times>0),]
  td=fdata
  ID=td[,1]
  measurement=td[,2]
  obstime=td[,8]/10
  tempdata=cbind(ID,measurement,obstime)
  
  ### candidate models for fitting
  M.set<-c(4,5,6,7)
  r.set<-c(3)
  
  ### parameters for fpca.mle
  ini.method="EM"
  basis.method="bs"
  sl.v=rep(0.5,10)
  max.step=50
  grid.l=seq(0,1,0.01)
  grids=seq(0,1,0.001)
  
  ### fit candidate models by fpca.mle
  try4=try(fpca.mle(tempdata, M.set,r.set,ini.method, basis.method,sl.v,max.step,grid.l,grids))
  if (inherits(try4, "try-error")) {next }
  result<-fpca.mle(tempdata, M.set,r.set,ini.method, basis.method,sl.v,max.step,grid.l,grids)
  summary(result)
  
  ### rescaled grid
  grids.new<-result$grid
  
  ### model selection result: the true model M=5, r=3 is selected with the smallest CV score among all converged models
  M<-result$selected_model[1]
  r<-result$selected_model[2]
  
  ## estimated eigenvalues 
  evalest<-result$eigenvalues ## estimated
  
  ## estimated error variance 
  sig2est<-result$error_var ## estimated
  
  ## estimated eigenfunctions
  eigenfest<-result$eigenfunctions
  
  ## estimated mean curve
  muest<-result$fitted_mean
  
  ##look at the CV scores and convergence for each model: note that model (M=5, r=4) does not converge.
  result$cv_scores ##CV
  result$converge ##convergence
  
  ## construct super dataset with extracted FPCA scores
  splinef=cbind(1,ns(t0vec,nbasis-1))
  tempdata=0
  tempfdata=0
  mwfdata=c() 
  for(t in t0vec)
  {
    tempdata=fdata[which(fdata$Time>t & fdata$times<t),]
    if(dim(tempdata)[1]==0 || dim(tempdata)[1]==1) { next }
    tempfdata=cbind(tempdata$id,tempdata$y,tempdata$times/10)
    fpcs<-fpca.score(tempfdata,grids.new,muest,evalest,eigenfest,sig2est,r)
    tempdata.id=tempdata[!duplicated(tempdata$id),]
    mwfdata=rbind(mwfdata,cbind(id=tempdata.id$id,Z1=fpcs[,1],Z2=fpcs[,2],Z3=fpcs[,3],tempdata.id[,3:6],Time=tempdata.id[,7],Shat=rep(0,dim(tempdata.id)[1]),tij=rep(t,dim(tempdata.id)[1])))
    
  }
  mwfdata=mwfdata[order(mwfdata$id),]
  mwfdata.id=mwfdata[!duplicated(mwfdata$id),]
  
  ## arrange the super dataset into the forms used in estimation
  Umat=c()
  Yvec=tvec=c()
  visit.mat=c()
  dataall=data.frame(id=rep(mwfdata.id[,1],each=length(t0vec)),Z1=999,Z2=999,Z3=999,X=999,Time=999,event=999,tij=rep(t0vec,dim(mwfdata.id)[1]))
  for(k in mwfdata.id[,1])
  {
    for(t in t0vec)
    {
      if(sum(mwfdata$id==k & mwfdata$tij==t)){
        dataall[dataall$id==k & dataall$tij==t,c(2,3,4)]= mwfdata[mwfdata$id==k & mwfdata$tij==t,c(2,3,4)]}
    }
    dataall[dataall$id==k,5]= mwfdata[which(mwfdata$id==k)[1] ,]$X
    dataall[dataall$id==k,6]= mwfdata[which(mwfdata$id==k)[1] ,]$Time
    dataall[dataall$id==k,7]= mwfdata[which(mwfdata$id==k)[1] ,]$event
    Umat=rbind(Umat,cbind(splinef,splinef*dataall[dataall$id==k,]$X[1],splinef*dataall[dataall$id==k,]$Z1,splinef*dataall[dataall$id==k,]$Z2,splinef*dataall[dataall$id==k,]$Z3))
    Yvec=c(Yvec,log(pmax(mwfdata.id[mwfdata.id[,1]==k,]$Time-t0vec,0.00001)))
    tvec=c(tvec,t0vec)
    visit.mat=rbind(visit.mat, t0vec%in% mwfdata[mwfdata$id==k,]$tij*1)
  }
  
  da=data.frame(mwfdata.id$Time,mwfdata.id$event,mwfdata.id$X)
  names(da)=c("Y","delta","X")
  ni.vec=rep(length(t0vec),dim(mwfdata.id)[1])
  obj=list(da=da,ni.vec=ni.vec,Umat=Umat,Yvec=Yvec,tvec=tvec,visit.mat=visit.mat,dataall=dataall,mwfdata = mwfdata , mwfdata.id = mwfdata.id)
  return(obj)
}

fpca.analysis <- function(obj)
{
  ########################
  ## estimation
  ########################
  #wt.p is weight used for re-sampling; for fitting only, set wt.p to 1
  #obj is a list generated by da.genrt()
  
  da=obj$da; 
  Umat=obj$Umat;  
  tvec=obj$tvec; 
  ni.vec=obj$ni.vec
  Y.reg=obj$Yvec; 
  U.reg=obj$Umat; 
  visit.mat=obj$visit.mat
  dataall=obj$dataall
  mwfdata = obj$mwfdata
  mwfdata.id = obj$mwfdata.id
  wt.p=rep(1,dim(da)[1])
  
  Gest=SC.func(da$Y,1-da$delta,wt.p)
  SC=stepfun(Gest$deathtime,c(1,Gest$survp)) #survival function for censoring time
  wt.reg.1=wt.reg.2=wt.reg.add=c()
  n0=dim(da)[1]
  for(i in 1:n0){
    ni=ni.vec[i]
    index=sum(ni.vec[1:(i-1)])
    if(i==1){index=0}
    U.mat.curr=matrix(Umat[(index+1):(index+ni),],ni)
    tseq=tvec[(index+1):(index+ni)]
    come.cur=visit.mat[i,] #the visit indicators for the ith subject
  
    invw.curr=pmax(SC(da$Y[i])/pmax(SC(tseq),0.000001),0.0001)
    wt.curr=(da$delta[i]*(da$Y[i]>tseq)/invw.curr*come.cur)*wt.p[i]
    wt.curr2=(da$delta[i]/invw.curr)
    pseudo1=-apply(U.mat.curr*wt.curr,2,sum)
    pseudo2=2*apply(U.mat.curr*(da$Y[i]>tseq)*come.cur*wt.p[i]*tau0,2,sum)
    
    Y.reg=c(Y.reg,M,M)
    U.reg=rbind(U.reg,rbind(pseudo1,pseudo2))
    wt.reg.1=c(wt.reg.1,wt.curr) #lin et al방법의 점추정을 위한 weight
    wt.reg.add=c(wt.reg.add,1,1)
    wt.reg.2=c(wt.reg.2,wt.curr2) #induced smoothing방법의 점추정을 위한 weight
    
  wt.reg = c(wt.reg.1,wt.reg.add)
    
    #evaluate the quantile model at tpseq
  try=try(rq.wfit(U.reg,Y.reg,weights=wt.reg)$coef)
  if (inherits(try, "try-error")) {gamma.fit=rep(NA,20) }else {
    gamma.fit=rq.wfit(U.reg,Y.reg,weights=wt.reg)$coef
    non.conv=(sum(abs(gamma.fit))>=200)|(is.na(sum(gamma.fit)))
    if(sum(abs(gamma.fit))>=5000){gamma.fit=gamma.fit*NA}}

  return(list(gamma.fit = gamma.fit, mwfdata = mwfdata, mwfdata.id =mwfdata.id, come.cur = come.cur , wt.reg = wt.reg, wt.reg.2 = wt.reg.2, da=da))
  }
}

# length(fit.fpca$wt.reg); summary(fit.fpca$wt.reg)
# length(fit.fpca$wt.reg.2); summary(fit.fpca$wt.reg.2)

##simulation

set.seed(2019122036)
NN = 1
nbasis = 4
M = 1E8
L = 6
N = 400
tau0 = 0.5
t0vec = seq(2, L, 1)
vauc.fpca = c()
vmean.fpca = c()

for (nsim in 1:NN) {
  dataobj = da.genrt(N)
  data = dataobj$data
  
  temp.fpca = fpca.fit(data)
  
  auc.fpca = c()
  mean.fpca = matrix(NA, ncol = 5, nrow = 1)
  
  for (tpseq in c(2, 2.5, 3, 3.5, 4)) {
    tryCatch({
      fit.fpca = fpca.analysis(temp.fpca)
      IDall = unique(fit.fpca$data[which(fit.fpca$data$Time >= tpseq), 1])
      qt.fpca = c()
      qt.fpcaauc = c()
      qt.true = c()
      
      data.id = data[!duplicated(data[, 1]), ]
      qt.true = cbind(IDall, pmin(data.id[data.id[, 1] %in% IDall, 5] - tpseq, L - tpseq))
      
      for (k in IDall) {
        ND.fpca = fit.fpca$data[fit.fpca$data[, 1] == k, ]
        ND.fpca = ND.fpca[which(abs(ND.fpca$tij - tpseq) == min(abs(ND.fpca$tij - tpseq)))[1], ]
        qt = exp(fit.fpca$coef[1] + fit.fpca$coef[2] * ND.fpca$X + fit.fpca$coef[3] * ND.fpca$Z1 + fit.fpca$coef[4] * ND.fpca$Z2 + fit.fpca$coef[5] * ND.fpca$Z3)
        qt.fpcaauc = rbind(qt.fpcaauc, c(k, qt))
        qt = exp(fit.fpca$coef[1] + fit.fpca$coef[2] * ND.fpca$X + fit.fpca$coef[3] * ND.fpca$Z1 + fit.fpca$coef[4] * ND.fpca$Z2 + fit.fpca$coef[5] * ND.fpca$Z3)
        qt = ifelse(qt < range(qt.true[, 2])[2] & qt > range(qt.true[, 2])[1], qt, NA)
        qt.fpca = rbind(qt.fpca, c(k, qt))
      }
      
      mean.fpca = cbind(mean.fpca, mean(abs(qt.true[, 2] - qt.fpca[, 2]), na.rm = T))
      auc.fpca = cbind(auc.fpca, AUC(data, qt.fpcaauc, L))
      
    }, error = function(e) {
      cat("Error occurred at nsim:", nsim, "tpseq:", tpseq, "\n")
      error_log[[length(error_log) + 1]] = list(nsim = nsim, tpseq = tpseq, error = e, data = data)
    })
  }
  
  vmean.fpca = rbind(vmean.fpca, mean.fpca)
  vauc.fpca = rbind(vauc.fpca, auc.fpca)
}

# fit.fpca$coef # NULL
# dim(fit.fpca$mwfdata.id)

#idw.temp.est: induced smoothing 방법으로 구한 점추정치

objectF = function(beta){
  beta = as.matrix(beta)
  
  result = (1/N)*t(U*come) %*% (W*(pnorm((U %*% beta-logY)/sqrt(diag(U %*% H %*% t(U)))))-tau0)
  result = as.vector(result)
}

##########################################################################################
### ni.vec 차원때문에 계산 안되는거 해결 ### 
#a = ni.vec*t(U*come)
#b = 5 * t(U*come)

# dim(U) : 610 20
#dim((pnorm((U%*%a-logY)/sqrt(diag(U %*% H %*% t(U))))))
#dim(W)
#dim((1/N)*ni.vec*t(U*come))
#dim((W*(pnorm((U%*%betastart-logY)/sqrt(diag(U %*% H %*% t(U)))))-tau0))
#########################################################################################


obj = temp.fpca
U=obj$Umat #U: covariate W,Z에 fractional basis function을 곱한 형태
ni.vec = obj$ni.vec
t=obj$tvec
come=fit.fpca$come.cur#visit indicator
W=fit.fpca$wt.reg.2 #censoring weight
logY=obj$Yvec #log of Y(observed failure time)
nc=length(fit.fpca$gamma.fit)
H = diag(1/(nc*N), nc, nc) # 12 개수 변경

betastart=as.matrix(fit.fpca$gamma.fit)
is.fit = dfsane(betastart, objectF, control = list(tol=5*1e-03))
print(is.fit)
idw.fit = is.fit$par  ## 최적화된 변수 값들


## lin et al point estimation


gamma.fit = fit.fpca$gamma.fit

splinef=cbind(1,ns(t0vec,nbasis-1))[which(abs(t0vec-tpseq)==min(abs(t0vec-tpseq)))[1],]
alpha0.est=as.vector(splinef%*%gamma.fit[1:nbasis])
alpha1.est=as.vector(splinef%*%gamma.fit[(nbasis+1):(2*nbasis)])
beta1.est=as.vector(splinef%*%gamma.fit[(2*nbasis+1):(3*nbasis)])
beta2.est=as.vector(splinef%*%gamma.fit[(3*nbasis+1):(4*nbasis)])
beta3.est=as.vector(splinef%*%gamma.fit[(4*nbasis+1):(5*nbasis)])
temp.est=cbind(alpha0.est,alpha1.est,beta1.est,beta2.est,beta3.est)


## idw method point estimation

tryCatch({
  if (is.fit$convergence == 0){  ## 수렴했는지
    idw.splinef=cbind(1,ns(t0vec,nbasis-1))[which(abs(t0vec-tpseq)==min(abs(t0vec-tpseq)))[1],]
    idw.alpha0.est=as.vector(splinef%*%idw.fit[1:nbasis])
    idw.alpha1.est=as.vector(splinef%*%idw.fit[(nbasis+1):(2*nbasis)])
    idw.beta1.est=as.vector(splinef%*%idw.fit[(2*nbasis+1):(3*nbasis)])
    idw.beta2.est=as.vector(splinef%*%idw.fit[(3*nbasis+1):(4*nbasis)])
    idw.beta3.est=as.vector(splinef%*%idw.fit[(4*nbasis+1):(5*nbasis)])
    idw.temp.est=cbind(idw.alpha0.est,idw.alpha1.est,idw.beta1.est,idw.beta2.est,idw.beta3.est)
  } else {
    idw.temp.est = c(NA,NA,NA)
  }
},error=function(e){
  idw.temp.est = c(NA,NA,NA)
}) 

temp.est
idw.temp.est

print(fit.fpca)


##########################################################
##########################################################
### data random 하게 바꿔보기 ###
##########################################################
##########################################################

results <- list()
NN = 1
nbasis = 4
M = 1E8
L = 6
N = 400
tau0 = 0.5
t0vec = seq(2, L, 1)
sample_num <-sample(1:100,10)

for (i in 1:10) {  # 100번의 반복을 위한 for문
  set.seed(i)  # i를 이용하여 seed 설정
  vauc.fpca = c()
  vmean.fpca = c()
  temp.est=c()
  idw.temp.est=c()
  for (nsim in 1:NN) {
    dataobj = da.genrt(N)
    data = dataobj$data
    
    temp.fpca = fpca.fit(data)
    
    auc.fpca = c()
    mean.fpca = matrix(NA, ncol = 5, nrow = 1)
    
    for (tpseq in c(2, 2.5, 3, 3.5, 4)) {
      tryCatch({
        fit.fpca = fpca.analysis(temp.fpca)
        IDall = unique(fit.fpca$data[which(fit.fpca$data$Time >= tpseq), 1])
        qt.fpca = c()
        qt.fpcaauc = c()
        qt.true = c()
        
        data.id = data[!duplicated(data[, 1]), ]
        qt.true = cbind(IDall, pmin(data.id[data.id[, 1] %in% IDall, 5] - tpseq, L - tpseq))
        
        for (k in IDall) {
          ND.fpca = fit.fpca$data[fit.fpca$data[, 1] == k, ]
          ND.fpca = ND.fpca[which(abs(ND.fpca$tij - tpseq) == min(abs(ND.fpca$tij - tpseq)))[1], ]
          qt = exp(fit.fpca$coef[1] + fit.fpca$coef[2] * ND.fpca$X + fit.fpca$coef[3] * ND.fpca$Z1 + fit.fpca$coef[4] * ND.fpca$Z2 + fit.fpca$coef[5] * ND.fpca$Z3)
          qt.fpcaauc = rbind(qt.fpcaauc, c(k, qt))
          qt = exp(fit.fpca$coef[1] + fit.fpca$coef[2] * ND.fpca$X + fit.fpca$coef[3] * ND.fpca$Z1 + fit.fpca$coef[4] * ND.fpca$Z2 + fit.fpca$coef[5] * ND.fpca$Z3)
          qt = ifelse(qt < range(qt.true[, 2])[2] & qt > range(qt.true[, 2])[1], qt, NA)
          qt.fpca = rbind(qt.fpca, c(k, qt))
        }
        
        mean.fpca = cbind(mean.fpca, mean(abs(qt.true[, 2] - qt.fpca[, 2]), na.rm = T))
        auc.fpca = cbind(auc.fpca, AUC(data, qt.fpcaauc, L))
        
      }, error = function(e) {
        cat("Error occurred at nsim:", nsim, "tpseq:", tpseq, "\n")
        error_log[[length(error_log) + 1]] = list(nsim = nsim, tpseq = tpseq, error = e, data = data)
      })
    }
    
    vmean.fpca = rbind(vmean.fpca, mean.fpca)
    vauc.fpca = rbind(vauc.fpca, auc.fpca)
  }
  
  obj = temp.fpca
  U=obj$Umat #U: covariate W,Z에 fractional basis function을 곱한 형태
  ni.vec = obj$ni.vec
  t=obj$tvec
  come=fit.fpca$come.cur#visit indicator
  W=fit.fpca$wt.reg.2 #censoring weight
  logY=obj$Yvec #log of Y(observed failure time)
  nc=length(fit.fpca$gamma.fit)
  H = diag(1/(nc*N), nc, nc) # 12 개수 변경
  
  betastart=as.matrix(fit.fpca$gamma.fit)
  is.fit = dfsane(betastart, objectF, control = list(tol=5*1e-03))
  idw.fit = is.fit$par  ## 최적화된 변수 값들
  
  
  ## lin et al point estimation
  
  
  gamma.fit = fit.fpca$gamma.fit
  
  
  splinef=cbind(1,ns(t0vec,nbasis-1))[which(abs(t0vec-tpseq)==min(abs(t0vec-tpseq)))[1],]
  alpha0.est=as.vector(splinef%*%gamma.fit[1:nbasis])
  alpha1.est=as.vector(splinef%*%gamma.fit[(nbasis+1):(2*nbasis)])
  beta1.est=as.vector(splinef%*%gamma.fit[(2*nbasis+1):(3*nbasis)])
  beta2.est=as.vector(splinef%*%gamma.fit[(3*nbasis+1):(4*nbasis)])
  beta3.est=as.vector(splinef%*%gamma.fit[(4*nbasis+1):(5*nbasis)])
  temp.est=cbind(alpha0.est,alpha1.est,beta1.est,beta2.est,beta3.est)
  
  
  ## idw method point estimation
  
  tryCatch({
    if (is.fit$convergence == 0){  ## 수렴했는지
      idw.splinef=cbind(1,ns(t0vec,nbasis-1))[which(abs(t0vec-tpseq)==min(abs(t0vec-tpseq)))[1],]
      idw.alpha0.est=as.vector(splinef%*%idw.fit[1:nbasis])
      idw.alpha1.est=as.vector(splinef%*%idw.fit[(nbasis+1):(2*nbasis)])
      idw.beta1.est=as.vector(splinef%*%idw.fit[(2*nbasis+1):(3*nbasis)])
      idw.beta2.est=as.vector(splinef%*%idw.fit[(3*nbasis+1):(4*nbasis)])
      idw.beta3.est=as.vector(splinef%*%idw.fit[(4*nbasis+1):(5*nbasis)])
      idw.temp.est=cbind(alpha0.est,alpha1.est,beta1.est,beta2.est,beta3.est)
    } else {
      idw.temp.est = c(NA,NA,NA)
    }
  },error=function(e){
    idw.temp.est = c(NA,NA,NA)
  })
  
  # 결과 저장
  results[[i]] <- list(
    i = i,
    convergence = is.fit$convergence == 0,
    temp.est = temp.est,
    idw.temp.est = idw.temp.est
  )
  
}

##########################################################
###
### 초기값 random 하게 설정해보기 
###
##########################################################

# set.seed(10)

betastart=as.matrix(fit.fpca$gamma.fit) ## 이 부분 임의로 설정
is.fit = dfsane(betastart, objectF, control = list(tol=5*1e-03))
print(is.fit)
idw.fit = is.fit$par  ## 최적화된 변수 값들


## lin et al point estimation


gamma.fit = fit.fpca$gamma.fit


splinef=cbind(1,ns(t0vec,nbasis-1))[which(abs(t0vec-tpseq)==min(abs(t0vec-tpseq)))[1],]
alpha0.est=as.vector(splinef%*%gamma.fit[1:nbasis])
alpha1.est=as.vector(splinef%*%gamma.fit[(nbasis+1):(2*nbasis)])
beta1.est=as.vector(splinef%*%gamma.fit[(2*nbasis+1):(3*nbasis)])
beta2.est=as.vector(splinef%*%gamma.fit[(3*nbasis+1):(4*nbasis)])
beta3.est=as.vector(splinef%*%gamma.fit[(4*nbasis+1):(5*nbasis)])
temp.est=cbind(alpha0.est,alpha1.est,beta1.est,beta2.est,beta3.est)


## idw method point estimation

tryCatch({
  if (is.fit$convergence == 0){  ## 수렴했는지
    idw.splinef=cbind(1,ns(t0vec,nbasis-1))[which(abs(t0vec-tpseq)==min(abs(t0vec-tpseq)))[1],]
    idw.alpha0.est=as.vector(splinef%*%idw.fit[1:nbasis])
    idw.alpha1.est=as.vector(splinef%*%idw.fit[(nbasis+1):(2*nbasis)])
    idw.beta1.est=as.vector(splinef%*%idw.fit[(2*nbasis+1):(3*nbasis)])
    idw.beta2.est=as.vector(splinef%*%idw.fit[(3*nbasis+1):(4*nbasis)])
    idw.beta3.est=as.vector(splinef%*%idw.fit[(4*nbasis+1):(5*nbasis)])
    idw.temp.est=cbind(alpha0.est,alpha1.est,beta1.est,beta2.est,beta3.est)
  } else {
    idw.temp.est = c(NA,NA,NA)
  }
},error=function(e){
  idw.temp.est = c(NA,NA,NA)
}) 

temp.est
idw.temp.est






##############################################################################
##
## induced covariance estimation
## sd.idw: covariance estimator
##
##############################################################################


tpseq=c(0.8);tau0=0.5;
n0 <- 200
n1 <- 1000
exp0=1.5
eta0=0.9
#C0=4
C0=3
ne=200


if (1*is.na(idw.temp.est[1,1]) == 1) {  # if non-smooth estimator is NA
  sd.idw=c(NA,NA,NA)
} else {
  
  result.ismb=c()
  
  for (j in 1:ne){
    eta_p1 = rexp(n0,1)
    eta_p = rep(eta_p1, each=12)
    
    Gest.sd=SC.func(da$Y,1-da$delta,eta_p1)
    SC.sd=stepfun(Gest.sd$deathtime,c(1,Gest.sd$survp)) #survival function for censoring time (fitting sd)
    
    wt.reg.sd=c()
    
    for(m in 1:n0){
      ni=ni.vec[m]
      index=sum(ni.vec[1:(m-1)])
      if(m==1){index=0}
      tseq=tvec[(index+1):(index+ni)]
      
      invw.curr.sd=pmax(SC.sd(da$Y[m])/SC.sd(tseq),0.0001) #sd
      #wt.curr.sd=(da$delta[m]*(da$Y[m]>tseq)/invw.curr.sd*come.cur)
      wt.curr.sd=(da$delta[m]/invw.curr.sd)
      
      wt.reg.sd=c(wt.reg.sd,wt.curr.sd)#} #induced smoothing방법의 weight
      #wt.reg.idw.sd=c(wt.reg.sd) #idw weight (sd)
    }
    wt.reg.idw.sd=c(wt.reg.sd) #induced smoothing방법의 weight
    
    #U=obj$Umat
    #I=I.reg.idw
    #W=wt.reg.idw
    #logY=obj$Yvec
    #H = diag(1/n0, nc, nc)
    
    #result = (1/n0)*t(eta*come*U*I) %*%
    result = (1/n0)*t(eta_p*U*come) %*%
      #{eta*wt.reg.idw.sd*(pnorm((U%*%idw.fit-logY)/
      #                            sqrt(matrix(fatdiag(U %*% H %*% t(U), steps=12), 
      #                                                         byrow=T, ncol=12))))-tau0}
      {wt.reg.idw.sd*(pnorm((U%*%idw.fit-logY)/sqrt(diag(U %*% H %*% t(U)))))-tau0}
    result.ismb = cbind(result.ismb,result)
  }
  
  v = cov(t(result.ismb))
  #a.beta = (1/n0)*t(come*U*I*W*
  #                    as.vector(dnorm((U%*%idw.fit-logY)/
  #                      sqrt(matrix(fatdiag(U %*% H %*% t(U), steps=12), 
  #                          byrow=T, ncol=12)))))%*%
  #  (U/sqrt(matrix(fatdiag(U %*% H %*% t(U), steps=12), 
  #                 byrow=T, ncol=12)))
  #a.beta = (1/n0)*t(come*U*I*W*as.vector(dnorm((U%*%idw.fit-logY)/sqrt(diag(U %*% H %*% t(U))))))%*%(U/sqrt(diag(U %*% H %*% t(U))))
  ##a.beta = (1/n0)*t(U*I*W*as.vector(dnorm((U%*%idw.fit-logY)/sqrt(diag(U %*% H %*% t(U))))))%*%(U/sqrt(diag(U %*% H %*% t(U))))
  

  a.beta = (1/n0)*t(U*come*W*
                      #a.beta = (1/n0)*t(U*I_c*W*
                      #a.beta = (1/n0)*t(come*U*I*W*
                      as.vector(dnorm((U%*%idw.fit-logY)/
                                        sqrt(diag(U %*% H %*% t(U))))))%*%
    (U/sqrt(diag(U %*% H %*% t(U))))
  
  #inva.beta = qr.solve(a.beta)
  inva.beta <- solve(a.beta) 
  sigma = t(inva.beta) %*% v %*% inva.beta
  base.sigma = sigma[1:4,1:4]
  X.sigma = sigma[5:8,5:8]
  Z.t.sigma = sigma[9:12,9:12]
  pbase=poly.func(tpseq)
  v.base = pbase%*%base.sigma%*%t(pbase)
  v.X = pbase%*%X.sigma%*%t(pbase)
  v.Z.t = pbase%*%Z.t.sigma%*%t(pbase)
  V.idw = cbind(v.base,v.X,v.Z.t)

  return(list(gamma=gamma.fit, idw=idw.fit, samp=samp, est=temp.est, 
              idw.est=idw.temp.est))
}

################################################################################
##
## S2 simulation
## 
################################################################################
tpseq=c(0.8);tau0=0.5;
#n0 <- 200
n1 <- 1000
exp0=1.5
eta0=0.9
#C0=4
C0=3
ne=200

#results=matrix(0, nrow=n1, ncol=11)
#results=as.data.frame(results)
#index=1

# pararell computing -----------------------------;


#numCores <- detectCores(logical=FALSE) - 1
numCores <- detectCores() - 1
# numCores <- 40
myCluster <- makeCluster(numCores)
registerDoParallel(myCluster)

#sim.func.S1=function(n1, exp0, eta0, C0, ne){
#set.seed(114)

system.time(
  sim_func_S2 <-
    foreach (i=1:n1, .options.RNG = 1234, .packages = c('survival','BB','quantreg','emplik','prodlim')) %dorng% {
      #foreach (i=1:n1, .packages = c('survival','BB','quantreg','emplik','prodlim')) %dopar% {
      #    for (i in 1:n1) {
      obj=da.genrt.weibull.S2(n0, eta0, C0, tvec0, tau0)
      #da=obj$da
      fit=main.func.S2(obj,rep(1,n0),ne,tpseq, tau0)
      
      #readr::write_csv(obj$da, path = paste0(i,'_sim.csv'))
      #results$V1[i]=as.vector(t(fit$est))[1]
      #results$V2[i]=as.vector(t(fit$est))[2]
      #results$V3[i]=as.vector(t(fit$est))[3]
      #results$V4[i]=as.vector(t(fit$idw.est))[1]
      #results$V5[i]=as.vector(t(fit$idw.est))[2]
      #results$V6[i]=as.vector(t(fit$idw.est))[3]
      #results$V7[i]=1-sum(da$delta)/n0
      #results$V8[i]=fit$samp
      #results$V9[i]=as.vector(t(fit$sd.idw))[1]
      #results$V10[i]=as.vector(t(fit$sd.idw))[2]
      #results$V11[i]=as.vector(t(fit$sd.idw))[3]
      #cat("\n", paste0( "Iteration =  ", index, "END" ), "\n" )
      #write.csv(results, file=paste0("sim " ,i,".csv"))
      #index = index + 1
      
      #my_list <- list("dif" = c(pred_logT_p1,pred_logT_g1,pred_logT_r1) - log(pbc_test1$time_yr2), 
      #                "BS" = BS_pbc_pgr)
      
      my_list <- list("li" = fit$est,"idw" = fit$idw.est, 
                      #                "cens_rate"=1-sum(da$delta)/n0, "sample" = fit$samp,
                      "sample" = fit$samp,"idw.sd" = fit$sd.idw)
      
      return(my_list)
    })
stopCluster(myCluster)
closeAllConnections()
#warnings()

# initialization
par_est_li <- par_est_is <- sd_est_is <- matrix(rep(0,n1*3), nrow=n1)
#cens_rate <- samp_size <- c(rep(0,n1))
samp_size <- c(rep(0,n1))

for(i in 1:n1){
  par_est_li[i,1:3] <- sim_func_S2[[i]]$li
  par_est_is[i,1:3] <- sim_func_S2[[i]]$idw
  #cens_rate[i] <- sim_func_S2[[i]]$cens_rate
  samp_size[i] <- sim_func_S2[[i]]$sample
  sd_est_is[i,1:3] <- sim_func_S2[[i]]$idw.sd
}

(is_est_ese <- sqrt(diag(var(par_est_is))))
(li_est_ese <- sqrt(diag(var(par_est_li))))
colMeans(sd_est_is)
colMeans(par_est_li)
colMeans(par_est_is)

par(mfrow=c(1,3))
plot(par_est_li[,1],par_est_is[,1])
plot(par_est_li[,2],par_est_is[,2])
plot(par_est_li[,3],par_est_is[,3])

#res_sim: 시뮬레이션 결과
#res_sim <- cbind(par_est_li,par_est_is,sd_est_is,cens_rate,samp_size)
res_sim <- cbind(par_est_li,par_est_is,sd_est_is,samp_size)

#시뮬레이션 결과를 csv파일로 저장한다.
write.csv(res_sim, file=paste0("S2_sim re H12_0328_", n1,"_tau",tau0, "_t",tpseq,"_ne", ne,".csv"))
  