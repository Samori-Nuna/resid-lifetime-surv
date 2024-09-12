
##########################################################################
 ##
 ## This version contains R code when there is L=6 trunction for trueTimes
 ##
 ## Only one Difference from Lin et al.:
 ## (1) when there is a uniroot error: 50 -> changed to NA
 ## 
 ## Updated: 2024.09.12
 ##
##########################################################################


rm(list=ls()); par(mfrow=c(1,1))

set.seed(1) # seed number setting

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
# other options for FPCA
# library(fda) # more general
# library(fdapace)

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

  L=6

# da.genrt=function(N){
    N=400 #original code: N=400 #2024.09.01
  n=N ##number of sample size
  K=length(c(0,seq(0.5,10,0.5))) #number of repetition
  
  # parameters for the f(t)
  betas <- cbind(-2, 0.2) # original code
    # betas <- cbind(3, -2) # *paper: cbind(3, -2)
  
  sigma.y <- 0.3 # measurement error(epsilon_i) standard deviation
  
  D=matrix(c(0.1,0.04,0.04,0.1),nrow=2) # random effects bi = (bi1,bi2)^T ~ N(0,D]) # original code
    # D=matrix(c(0.4,0.1,0.1,0.2),nrow=2) # *paper: (.4,.1,.1,.2)
  
  bb <- mvrnorm(n, rep(0, nrow(D)), D) # bi
  randomeffect<-bb
  
  # parameters for the survival model
  betasm=0.4
  sigma.t=1.2
  
  # parameters for censoring
  # censoring=0.12 ## small -> censoring rate small
  
  Xcoef=1 # original code
    # Xcoef=0.5 # *paper: Xcoef=0.5
  
  X.id=runif(n,0,1)
  X=rep(X.id,each=K)
  
  # design matrices for the longitudinal measurement model
  times <- c(replicate(n, c(0,seq(0.5,10,0.5)))) # length(c(0,seq(0.5,10,0.5))) = K
  times[which(times!=0)] <- times[which(times!=0)]+ rnorm(n*(K-1),0,0.1)
  # summary(times)
  # hist(times, breaks = 200)
  Ai=betas[,1]
  Bi=betas[,2]
  
  id <- rep(1:n, each = K)
  eta.y <-as.vector(Ai+Bi*times+randomeffect[id,1]+randomeffect[id,2]*times) # Bi(t)
  y <- rnorm(n * K, eta.y, sigma.y) # Zi(t) = Bi(t) + epsilon_i
  
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
  for (i in 1:n){
    Up <- 50
    tries <- 5
    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    while(inherits(Root, "try-error") && tries > 0) # The use of && ensures that if inherits(Root, "try-error") is FALSE, the tries > 0 condition is not even checked, which is the desired behavior. This ensures both efficiency and correctness in the logic. The && operator is appropriate here because it ensures short-circuit evaluation, meaning that if inherits(Root, "try-error") is FALSE, it does not evaluate tries > 0. This is efficient and prevents unnecessary evaluation of the second condition if the first one is already FALSE.
    {
      tries <- tries - 1
      Up <- Up + 50
      Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    }
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA # # original code: else 50
  }
    sum(is.na(trueTimes)) # average 4.5(original code number) vs. 40(paper number)
    # ?uniroot
    # ?try
    # summary(trueTimes)
  na.ind <- !is.na(trueTimes)
    # sum(na.ind)
    # length(na.ind) - sum(na.ind) # 98(/1000) 50's (set.seed(1)) -> 40(/400) 50's..!
  trueTimes <- trueTimes[na.ind]
  n <- length(trueTimes)
  trueTimes <- rep(trueTimes,each=K)
  
  #Ctimes <- -log(runif(n, 0, 1))/(censoring*exp(X.id[na.ind]*Xcoef))
  Ctimes <- runif(n,0,18) # original code: not 50, but 18
  Ctimes<-rep(Ctimes,each=K)
    # par(mfrow=c(2,1))
    # hist(Ctimes)
    # hist(trueTimes)
  
  ### Observed time calculation
  # ### Synario (I) Original code: truncation at L(=6)
  Time=pmin(L,trueTimes)
  # summary(Time)
  #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 0.00766 0.11537 0.27332 0.98903 0.69952 6.00000
  Time=pmin(Time,Ctimes)
  # summary(Time)
  #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 0.00766 0.11082 0.26789 0.91689 0.65854 6.00000
  event=as.numeric(pmin(L,trueTimes) <= Ctimes)
  1 - sum(event)/length(Time) # censoring rate: 0.1016 (10.16%)
  
  ### Synario (II) No truncation
  ## (24/05/23) 위 세 줄 대신에 이렇게 바꾸어보기
  # Time = pmin(trueTimes, Ctimes) # changed!
  # # Time = trueTimes ################################################### when censoring = 0
  # # summary(Time)
  # #.   Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
  # # 0.00766  0.11082  0.26789  1.30096  0.65854 16.69550
  # event = as.numeric(trueTimes <= Ctimes) # event indicator generation(changed)
    # event = as.numeric(rep(1, n*K)) #################################### when censoring = 0
    # sum(event)/length(trueTimes)
    # 1 - sum(event)/length(Time) # censoring rate: 0.5665 (56.65%) 
    # par(mfrow=c(1,2))
    # hist(trueTimes, xlim = c(0,51), ylim = c(0, 6000), breaks = 100)
    # hist(Ctimes, xlim = c(0,51), ylim = c(0, 6000), breaks = 100)
  
  #Time=pmin(trueTimes,Ctimes)
  #Time=pmin(L,Time)
  #event=as.numeric(trueTimes <= Ctimes)
  
  long.na.ind <- rep(na.ind, each = K)
  y <- y[long.na.ind] # y = Zi(t)
  times<-times[long.na.ind] 
  X<-X[long.na.ind]
  
  ind <- times <= Time # sum(ind) = (I) 982; (II) 1181
  id <- id[long.na.ind][ind]
  y <- y[ind] # y = Zi(t)
    # summary(y)
    # par(mfrow=c(1,1))
    # hist(y, breaks = 400)
    # boxplot(y)
  X <- X[ind]
  times=times[ind]
  Time=Time[ind]
  trueTimes=trueTimes[ind]
  Ctimes=Ctimes[ind]
  event=event[ind]
  data=cbind(id,y,X,event,trueTimes,Ctimes,Time,times)
  data.id=data[!duplicated(data[,1]),]
  data.id=as.data.frame(data.id)
  plot(survfit(Surv(Time,event)~1,data=data.id))
  max(data.id[,7]) # (I) 6(=L) -> (II) 16.6955
  1 - sum(data.id[,4])/dim(data.id)[1] # censoring rate = 0.1335 (when line 141, then 0) ############# correction made: 2024.09.01.
  # return(list(data=data,X.id=X.id,Xcoef=Xcoef,Ai=Ai,Bi=Bi,bb=bb,betasm=betasm,sigma.t=sigma.t))
# }
  
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

  # L=6
  t0vec = seq(2, L, 1)
  
  # fpca.score function with original scale (visit times)
  fpca.score2 <- function(data.m, grids.u, muhat, eigenvals, eigenfuncs, sig2hat, K) {
    temp <- table(data.m[, 1])
    n <- length(temp)
    m.l <- as.vector(temp)
    result <- matrix(0, n, K)
    N <- length(grids.u)
    evalmat <- diag(eigenvals[1:K])
    current <- 0
    eigenfuncs.u <- t(eigenfuncs)
    data.u <- matrix(as.numeric(as.vector(data.m[, -1])),
                     nrow = nrow(data.m[, -1]), ncol = ncol(data.m[, -1]))
    
    cat("N:", N, "\n")
    cat("evalmat:\n", evalmat, "\n")
    
    for (i in 1:n) {
      Y <- as.vector(data.u[(current + 1):(current + m.l[i]), 1])
      meastime <- data.u[(current + 1):(current + m.l[i]), 2]
      gridtime <- ceiling(N * (meastime / 10))  # 시간 변수를 다시 10으로 나누어 인덱싱 문제 해결
      
      # gridtime이 올바른 범위 내에 있는지 확인
      gridtime[gridtime < 1] <- 1
      gridtime[gridtime > N] <- N
      
      muy <- muhat[gridtime]
      Phiy <- matrix(eigenfuncs.u[gridtime, 1:K], ncol = K)
      Sigy <- Phiy %*% evalmat %*% t(Phiy) + sig2hat * diag(m.l[i])
      temp.y <- matrix(Y - muy)
      result[i, ] <- evalmat %*% t(Phiy) %*% solve(Sigy, temp.y)
      current <- current + m.l[i]
    }
    return(result)
  }
  
# fpca.fit <- function(data){
  
  fdata=as.data.frame(data) # data=cbind(id,y,X,event,trueTimes,Ctimes,Time,times)
    # dim(fdata[which(fdata$times<=10),]) # (I) (ind <- times <= Time)에서 이미 6 이하로 걸러졌음(변화 X) 982x8; (II) 1170x8
    # dim(fdata[which(fdata$times>0),]) # (I) N(=n)개 줄어듦 582 x 8; (II) 781x8
  fdata=fdata[which(fdata$times<=10 & fdata$times>0),]   # times가 0 초과 10 이하
    # 1 - sum(fdata$event)/dim(fdata)[1] # censoring rate goes up to 79.35%
    # summary(fdata$times)
    # (I)  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    #   0.2203  0.9239  1.9979  2.4160  3.9107  5.9856 
    # (II) Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    #   0.2203  1.0858  3.1079  3.7406  5.9687  9.9866
    # dim(fdata)
    # (I) 582x8; (II) 770x8
  td=fdata
  ID=td[,1]   #id
  measurement=td[,2]  # y: Zi(t) = Bi(t) + epsilon_i
  obstime=td[,8]/10  # (visit) times/10, originally: times/10
  obstime2=td[,8]
  tempdata=cbind(ID,measurement,obstime)
  tempdata2=cbind(ID,measurement,obstime2)
    # summary(measurement)
    #      Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    # -27.1988  -9.9993  -3.9067  -5.3585   0.5269   3.7815
    # summary(obstime)
    # summary(obstime2)
  
  ### candidate models for fitting
  M.set<-c(4,5,6,7)  # number of basis functions
  r.set<-c(3)    # dimension of the process: originally set only to perform when r=3
  
  ### parameters for fpca.mle (default)
  ini.method="EM"    # initial method for Newton or "loc":locally linear method
  basis.method="bs"   # basis functions / bs = cubic B-splines
  sl.v=rep(0.5,10)   # step length for first K steps (length(sl.v)? 10)
  max.step=50        # maximum number of iterations Newton will run 
  grid.l=seq(0,1,0.01) # grid used by the local linear method
  grids=seq(0,1,0.001) # grid used by EM: diff from default (seq(0, 1, 0.002))
  
   # nsim = 1
  
  ### fit candidate models by fpca.mle
  tryCatch({
    result <- fpca.mle(tempdata, M.set, r.set, ini.method, basis.method, sl.v, max.step, grid.l, grids) # 별도 파일
    cat("fpca.mle: successfully done.\n")
  }, error = function(e){
    # error_log[[length(error_log) + 1]] <- paste("nsim =", nsim, "step1(fpca.mle)-scaled:", e$message)
    cat("Error occurred in fpca.mle(scaled), step1, for nsim =", nsim, ": ", e$message, "\n")
    result <- fpca.mle(tempdata, M.set, r.set, ini.method, basis.method, sl.v, max.step, grid.l, grids)
  })

  tryCatch({
    result2 <- fpca.mle(tempdata2, M.set, r.set, ini.method, basis.method, sl.v, max.step, grid.l, grids) # 별도 파일
    cat("fpca.mle: successfully done.\n")
  }, error = function(e){
    # error_log[[length(error_log) + 1]] <- paste("nsim =", nsim, "step1(fpca.mle)-unscaled:", e$message)
    cat("Error occurred in fpca.mle(unscaled), step1, for nsim =", nsim, ": ", e$message, "\n")
    result2 <- fpca.mle(tempdata2, M.set, r.set, ini.method, basis.method, sl.v, max.step, grid.l, grids)
  })
  
  ### re-scaled grid: different
  grids.new <- result$grid
    # summary(grids)
    # summary(grids.new)
  grids.new2 <- result2$grid
    # summary(grids.new2)
    
  ### model selection result: the true model M=4, r=3 is selected with the smallest CV score among all converged models
  M <- result$selected_model[1]
  r <- result$selected_model[2]
    nbasis <- M
    pc <- r
    
  M2 <- result2$selected_model[1]
  r2 <- result2$selected_model[2]
    nbasis2 <- M2
    pc2 <- r2
  
  ## estimated eigenvalues 
  evalest<-result$eigenvalues ## estimated
    # evalest
  evalest2<-result2$eigenvalues ## estimated
    # evalest2
  # round(evalest * 10, 6) == round(evalest2, 6)
  
  ## estimated error variance 
  sig2est<-result$error_var ## estimated
    # sig2est
  sig2est2<-result2$error_var ## estimated
    # sig2est2
    # round(sig2est, 6) == round(sig2est2, 6)
    # View(t(result$eigenfunctions))
  ## estimated eigenfunctions: different..
  eigenfest<-result$eigenfunctions
  eigenfest2<-result2$eigenfunctions
    # sum(round(eigenfest, 4) == round(eigenfest2, 4))
  
  ########################################################## plots: eigenfunction comparison
  # par(mfrow=c(1,2))
  # eigf <- result$eigenfunctions[1,]
  # eigf2 <- result2$eigenfunctions[1,]
  # plot.ts(eigf, col = "blue",
  #         ylim = c(-2.5, 2.5),
  #         main = "Eigenfunc Comparison: first pc")
  # lines(eigf2, col = "red")
  # legend("bottomleft", legend = c("Scaled", "Unscaled"), col = c("blue", "red"), lty = 1)
  # 
  # eigf_2 <- result$eigenfunctions[2,]
  # eigf2_2 <- result2$eigenfunctions[2,]
  # plot.ts(eigf_2, col = "blue",
  #         ylim = c(-2.5, 2.5),
  #         main = "Eigenfunc Comparison: second pc")
  # lines(eigf2_2, col = "red")
  # legend("bottomleft", legend = c("Scaled", "Unscaled"), col = c("blue", "red"), lty = 1)
  ##########################################################################################

  ## estimated mean curve: the same!
  muest<-result$fitted_mean
  muest2<-result2$fitted_mean
    # sum(round(muest, 4) == round(muest2, 4))
  
  ## CV scores and convergence for each model: the same!
  cv_scores <- result$cv_scores ##CV
    # cv_scores
  cv_scores2 <- result2$cv_scores ##CV
    # cv_scores2
    # round(cv_scores, 4) == round(cv_scores2, 4)
  
  convergence <- result$converge ##convergence
    # convergence
  convergence2 <- result2$converge ##convergence
    # convergence2
    # round(convergence, 4) == round(convergence2, 4)
  
  ## construct super dataset with extracted FPCA scores
  splinef=cbind(1,ns(t0vec,nbasis-1))      # ns ; Generate a Basis Matrix for Natural Cubic Splines
    # dim(splinef) # 5x6
    # dim(ns(t0vec,nbasis-1)) # 5x5: 5 from t0vec, df=nbasis-1
    # ?ns
  
  tempdata=0
  tempfdata=tempfdata2=0
  mwfdata=mwfdata2=c() 
  
  # t0vec=c(2,3,4,5)
  for(t in t0vec){
    # t=6
    tempdata=fdata[which(fdata$Time>t & fdata$times<t),] # planned visit time (=t) < observed time (Time)
                                                         # & irregular visit time (=times) < planned visit (=t)
    if( dim(tempdata)[1]==0 || dim(tempdata)[1]==1 ) { next } # if 데이터가 없거나 반복측정 아니면 skip됨
    tempfdata = cbind(tempdata$id,tempdata$y,tempdata$times/10) # unscaled 
    tempfdata2 = cbind(tempdata$id,tempdata$y,tempdata$times) # originally: times/10
    fpcs <- fpca.score(tempfdata, grids.new, muest, evalest, eigenfest, sig2est, r) ### SCORE! (times/10)
    fpcs2 <- fpca.score2(tempfdata2, grids.new2, muest2, evalest2, eigenfest2, sig2est2, r2) ### SCORE! (times)
    # fpcs2 <- fpcs2 * 10 # scaling back??? (multiplying by 10): 24/07/10
    # ?fpca.score
    tempdata.id=tempdata[!duplicated(tempdata$id),] # dim(tempdata.id)=47x8
    mwfdata=eval(parse(text=paste0( # id / Z1~Z3 / X / event / trueTimes / Ctimes / Time / Shat / tij
      "rbind(mwfdata,cbind(id=tempdata.id$id,",
      paste(sprintf("Z%d=fpcs[,%d]",1:pc,1:pc),collapse=','),
      ",tempdata.id[,3:6],Time=tempdata.id[,7],Shat=rep(0,dim(tempdata.id)[1]),tij=rep(t,dim(tempdata.id)[1])))")))
    mwfdata2=eval(parse(text=paste0(
      "rbind(mwfdata2,cbind(id=tempdata.id$id,",
      paste(sprintf("Z%d=fpcs2[,%d]",1:pc2,1:pc2),collapse=','),
      ",tempdata.id[,3:6],Time=tempdata.id[,7],Shat=rep(0,dim(tempdata.id)[1]),tij=rep(t,dim(tempdata.id)[1])))")))
  }
  
  mwfdata=mwfdata[order(mwfdata$id),] # id 순서대로 정렬
  mwfdata.id=mwfdata[!duplicated(mwfdata$id),] # 모든 t(in t0vec)에 대한 id 리스트업
    dim(mwfdata) # 184 x 11
    dim(mwfdata.id) # 47 x 11
  
  mwfdata2=mwfdata2[order(mwfdata2$id),] # id 순서대로 정렬
  mwfdata2.id=mwfdata2[!duplicated(mwfdata2$id),] # 모든 t(in t0vec)에 대한 id 리스트업
    dim(mwfdata2) # 184 x 11
    dim(mwfdata2.id) # 47 x 11
  
  # # PC score comparison
  # sum(round(mwfdata["Z1"],4) == round(mwfdata2["Z1"],4)) # different!
  # sum(round(mwfdata["Z2"],4) == round(mwfdata2["Z2"],4)) # diff
  # summary(mwfdata[c("Z1", "Z2")])
  # summary(mwfdata2[c("Z1", "Z2")])
  
  ## arrange the super data set into the forms used in estimation
  Umat=Umat2=c()
  Yvec=Yvec2=c() # Yvec == Yvec2 (TRUE)
  tvec=tvec2=c() # tvec == tvec2 (TRUE)
  visit.mat=visit.mat2=c() # visit.mat == visit.mat2 (TRUE)
  dataall=eval(parse(text=paste0( # id / Z1 ~ Z3 / X / Time / event / tij(=t0vec; 2, 3, 4, 5, 6)
    "data.frame(id=rep(mwfdata.id[,1],each=length(t0vec)),",
    paste(sprintf("Z%d=999",1:pc),collapse=','),
    ",X=999,Time=999,event=999,tij=rep(t0vec,dim(mwfdata.id)[1]))")))
  
  dataall2=eval(parse(text=paste0( # id / Z1 ~ Z3 / X / Time / event / tij(=t0vec; 2, 3, 4, 5, 6)
    "data.frame(id=rep(mwfdata2.id[,1],each=length(t0vec)),",
    paste(sprintf("Z%d=999",1:pc2),collapse=','),
    ",X=999,Time=999,event=999,tij=rep(t0vec,dim(mwfdata2.id)[1]))")))
  
  for(k in mwfdata.id[,1]){ # for each id (no duplicates)
    for(t in t0vec){
      if(sum(mwfdata$id==k & mwfdata$tij==t)){     ## mwfdata$id==k & mwfdata$tij==t일때만 dataall 데이터를 mwfdata 데이터로 대체함
        dataall[dataall$id==k & dataall$tij==t,2:(pc+1)]= mwfdata[mwfdata$id==k & mwfdata$tij==t, 2:(pc+1)]}
    }
    dataall[dataall$id==k,(pc+2)]= mwfdata[which(mwfdata$id==k)[1] ,]$X # time-indep covariate  ##다른 변수들은 mwfdata$id==k인 경우 mwfdata로 대체함
    dataall[dataall$id==k,(pc+3)]= mwfdata[which(mwfdata$id==k)[1] ,]$Time # Time = observed time
    dataall[dataall$id==k,(pc+4)]= mwfdata[which(mwfdata$id==k)[1] ,]$event # event indicator      
    Umat=eval(parse(text=paste0( #dim(Umat)=235x30 : why? 235=47(unique id)*5(t0vec), 30 = 6(*1) + 6(*X) + 6(*Z1~Z3), and each 6 comes from dim(splinef)=5x6(=col!)
      "rbind(Umat,cbind(splinef,splinef*dataall[dataall$id==k,]$X[1],",
      paste(sprintf("splinef*dataall[dataall$id==k,]$Z%d",1:pc),collapse=','),
      "))")))
    Yvec=c(Yvec,log(pmax(mwfdata.id[mwfdata.id[,1]==k,]$Time-t0vec,0.00001))) # Yvec = log(residual time), length(Yvec)=235=47(unique)x5(t0vec)
    tvec=c(tvec,t0vec) # length(tvec) == 47(unique id #) * 5 =length(t0vec)
    visit.mat=rbind(visit.mat, t0vec %in% mwfdata[mwfdata$id==k,]$tij*1) # visit indicators for each subject(row) x at each t0vec(column)
  }
  
  da=data.frame(mwfdata.id$Time,mwfdata.id$event,mwfdata.id$X)
  names(da)=c("Y","delta","X")
  ni.vec=rep(length(t0vec),dim(mwfdata.id)[1])
  
  for(k in mwfdata2.id[,1]){ # for each id (no duplicates)
    for(t in t0vec){
      if(sum(mwfdata2$id==k & mwfdata2$tij==t)){     ## mwfdata$id==k & mwfdata$tij==t일때만 dataall 데이터를 mwfdata 데이터로 대체함
        dataall2[dataall2$id==k & dataall2$tij==t,2:(pc2+1)]= mwfdata2[mwfdata2$id==k & mwfdata2$tij==t, 2:(pc2+1)]}
    }
    dataall2[dataall2$id==k,(pc2+2)]= mwfdata2[which(mwfdata2$id==k)[1] ,]$X # time-indep covariate  ##다른 변수들은 mwfdata$id==k인 경우 mwfdata로 대체함
    dataall2[dataall2$id==k,(pc2+3)]= mwfdata2[which(mwfdata2$id==k)[1] ,]$Time # Time = observed time
    dataall2[dataall2$id==k,(pc2+4)]= mwfdata2[which(mwfdata2$id==k)[1] ,]$event # event indicator      
    Umat2=eval(parse(text=paste0( #dim(Umat)=235x30 : why? 235=47(unique id)*5(t0vec), 30 = 6(*1) + 6(*X) + 6(*Z1~Z3), and each 6 comes from dim(splinef)=5x6(=col!)
      "rbind(Umat2,cbind(splinef,splinef*dataall2[dataall2$id==k,]$X[1],",
      paste(sprintf("splinef*dataall2[dataall2$id==k,]$Z%d",1:pc2),collapse=','),
      "))")))
    Yvec2=c(Yvec2,log(pmax(mwfdata2.id[mwfdata2.id[,1]==k,]$Time-t0vec,0.00001))) # Yvec = log(residual time), length(Yvec)=235=47(unique)x5(t0vec)
    tvec2=c(tvec2,t0vec) # length(tvec) == 47(unique id #) * 5 =length(t0vec)
    visit.mat2=rbind(visit.mat2, t0vec %in% mwfdata2[mwfdata2$id==k,]$tij*1) # visit indicators for each subject(row) x at each t0vec(column)
  }
  
  da2=data.frame(mwfdata2.id$Time,mwfdata2.id$event,mwfdata2.id$X) ## da == da2 (TRUE)
  names(da2)=c("Y","delta","X") 
  ni.vec2=rep(length(t0vec),dim(mwfdata2.id)[1]) ## ni.vec == ni.vec2 (TRUE)
  
  obj=list(da=da, ni.vec=ni.vec, Umat=Umat, Umat2=Umat2, 
           Yvec=Yvec, tvec=tvec, visit.mat=visit.mat, dataall=dataall, dataall2=dataall2, 
           nbasis=nbasis, nbasis2=nbasis2, pc=pc, pc2=pc2,
           cv=cv_scores, cv2=cv_scores2, converge=convergence, converge2=convergence2)
  # return(obj)
# }

# fpca.analysis <- function(obj)
# {
  ########################
  ## estimation
  ########################
  #wt.p is weight used for re-sampling; for fitting only, set wt.p to 1
  #obj is a list generated by da.genrt() -- ????
    tau0 = 0.5
    tpseq = c(2)
    # tpseq = c(2, 2.5, 3, 3.5, 4) #2024.09.03
    nsim = 1
    M = 1E8
    
  
  
  # scaled
  
  da=obj$da; 
  Umat=obj$Umat;  
  tvec=obj$tvec; 
  ni.vec=obj$ni.vec
  Y.reg=obj$Yvec; 
  U.reg=obj$Umat; 
  visit.mat=obj$visit.mat
  dataall=obj$dataall
  wt.p=rep(1,dim(da)[1])
  nbasis=obj$nbasis # should follow the result from previous function fit.fpca
  pc=obj$pc # should follow the result from previous function fit.fpca
  
  if(dim(da)[1] != sum(da$delta)){
    # censoring != 0,
    Gest <- SC.func(da$Y, 1-da$delta, wt.p) # data frame: n obs. of 4 variables(deathtime, ndeath, nrisk, survp)
    SC=stepfun(Gest$deathtime,c(1,Gest$survp)) #survival function for censoring time
  }

  # Gest=SC.func(da$Y,1-da$delta,wt.p)
  # SC=stepfun(Gest$deathtime,c(1,Gest$survp)) #survival function for censoring time
  wt.reg.1=wt.reg.2=wt.reg.add=c()
  n0=dim(da)[1]
  
  for(i in 1:n0){
    ni=ni.vec[i]
    index=sum(ni.vec[1:(i-1)])
    if(i==1){index=0}      # index -> 0 5 10 15 . . 
    U.mat.curr=matrix(Umat[(index+1):(index+ni),],ni)
    tseq=tvec[(index+1):(index+ni)]
    come.cur=visit.mat[i,] #the visit indicators for the ith subject
    
    if(dim(da)[1] != sum(da$delta)){
      invw.curr=pmax(SC(da$Y[i])/pmax(SC(tseq),0.000001),0.0001) # censoring exists
    } else {
      invw.curr=rep(1,dim(visit.mat)[2]) # no censoring
    }
    
    wt.curr=(da$delta[i]*(da$Y[i]>tseq)/invw.curr*come.cur)*wt.p[i] # weight-dimension check!!!
    wt.curr2=(da$delta[i]/invw.curr)
    pseudo1=-apply(U.mat.curr*wt.curr,2,sum)
    pseudo2=2*apply(U.mat.curr*(da$Y[i]>tseq)*come.cur*wt.p[i]*tau0,2,sum)
    
    Y.reg=c(Y.reg,M,M)
    U.reg=rbind(U.reg,rbind(pseudo1,pseudo2))
    wt.reg.1=c(wt.reg.1,wt.curr) #lin et al방법의 점추정을 위한 weight
    wt.reg.add=c(wt.reg.add,1,1)
    wt.reg.2=c(wt.reg.2,wt.curr2) } #induced smoothing방법의 점추정을 위한 weight
  
  wt.reg = c(wt.reg.1,wt.reg.add)
    wt.reg.sc <- wt.reg
    
    # ?rq.wfit
    # dim(U.reg)
    # View(U.reg)
    # NROW(Y.reg); NCOL(Y.reg)
    # NROW(wt.reg.sc); NCOL(wt.reg.sc)
    # ?as.matrix
    # wt.reg.sc <- as.matrix(wt.reg.sc, nrow = NROW(Y.reg), ncol = NCOL(Y.reg))
    # dim(wt.reg.sc)
    # dim(wt.reg.sc*U.reg)
    
  #evaluate the quantile model at tpseq
  try=try(rq.wfit(U.reg,Y.reg,weights=wt.reg.sc)$coef)
  if (inherits(try, "try-error")) {gamma.fit=rep(NA,nbasis*(pc+2)) }else {
    gamma.fit=rq.wfit(U.reg,Y.reg,weights=wt.reg.sc)$coef
    non.conv=(sum(abs(gamma.fit))>=200)|(is.na(sum(gamma.fit)))
    if(sum(abs(gamma.fit))>=5000){gamma.fit=gamma.fit*NA}}
  
  # #evaluate the quantile model at tpseq
  # try=try(rq.wfit(U.reg,Y.reg,weights=wt.reg.sc)$coef)
  # if (inherits(try, "try-error")) {
  #   gamma.fit=rep(NA,nbasis*(pc+2))
  #     error_msg <- attr(try, "condition")$message
  #     # error_log[[length(error_log) + 1]] <- paste("nsim =", nsim, "step2(rq.wfit) error", error_msg)
  #     cat("Error occurred in rq.wfit(scaled), step2, for nsim =", nsim, ": ", error_msg, "\n")
  #   } else {
  #   gamma.fit=rq.wfit(U.reg,Y.reg,weights=wt.reg.sc)$coef
  #   non.conv=(sum(abs(gamma.fit))>=200)|(is.na(sum(gamma.fit)))
  #   if(sum(abs(gamma.fit))>=5000){gamma.fit=gamma.fit*NA}}
  
  
  # singular design matrix? check
  wx <- U.reg * wt.reg.sc
  print(wx)
  dim(wx) # 826 x 30
  qr(wx)$rank # 20
  
  # ## Penalized quantile regression?
  # ?rq.fit.lasso
  # model_lasso <- rq.fit.lasso(x = wx, y = Y.reg, tau = 0.5, lambda = 0.1)
  
    # tryCatch({
  #   gamma.fit <- rq.wfit(U.reg, Y.reg, weights = wt.reg)$coef # rq.wfit 함수 실행 시도
  # }, error = function(e) {
  #   # 에러 발생 시 에러 메시지를 리스트에 저장
  #   error_log[[length(error_log) + 1]] <- paste("nsim =", nsim, "step2(rq.wfit) error:", e$message)
  #   cat("Error occurred in rq.wfit(scaled), step2, for nsim =", nsim, ": ", e$message, "\n")
  #   gamma.fit <- rep(NA,nbasis*(pc + 2)) # 에러 발생 시 gamma.fit 초기화: 20?=M(4)*5(number of coefficients: alpha0,1,beta1,2,3)
  # })
  # non.conv <- (sum(abs(gamma.fit)) >= 200) | (is.na(sum(gamma.fit))) # non.conv 및 gamma.fit 값 확인
  # if(sum(abs(gamma.fit))>=5000){gamma.fit=gamma.fit*NA}
  
  splinef=cbind(1,ns(t0vec,nbasis-1))[which(abs(t0vec-tpseq)==min(abs(t0vec-tpseq)))[1],]
  alpha0.est=as.vector(splinef%*%gamma.fit[1:nbasis])
  alpha1.est=as.vector(splinef%*%gamma.fit[(nbasis+1):(2*nbasis)])
  beta=c()
  for(i in 1:pc){
    beta=cbind(beta, as.vector(splinef%*%gamma.fit[((i+1)*nbasis+1):((i+2)*nbasis)]))}
  temp.est=data.frame(alpha0.est,alpha1.est,beta)
  names(temp.est)=c("alpha0.est","alpha1.est",sprintf("beta%d.est",1:pc))
  
  temp.est
  
  # unscaled:
  
  da=obj$da;
  Umat=obj$Umat2;  
  tvec=obj$tvec;
  ni.vec=obj$ni.vec
  Y.reg=obj$Yvec;
  U.reg=obj$Umat2; 
  visit.mat=obj$visit.mat
  dataall2=obj$dataall2
  wt.p=rep(1,dim(da)[1])
  nbasis2=obj$nbasis2 # should follow the result from previous function fit.fpca
  pc2=obj$pc2 # should follow the result from previous function fit.fpca
  
  if(dim(da)[1] != sum(da$delta)){
    # censoring != 0,
    Gest <- SC.func(da$Y, 1-da$delta, wt.p) # data frame: n obs. of 4 variables(deathtime, ndeath, nrisk, survp)
    SC=stepfun(Gest$deathtime,c(1,Gest$survp)) #survival function for censoring time
  }
  
  # Gest=SC.func(da$Y,1-da$delta,wt.p)
  # SC=stepfun(Gest$deathtime,c(1,Gest$survp)) #survival function for censoring time
  wt.reg.1=wt.reg.2=wt.reg.add=c()
  n0=dim(da)[1]
  for(i in 1:n0){
    ni=ni.vec[i]
    index=sum(ni.vec[1:(i-1)])
    if(i==1){index=0}      # index -> 0 5 10 15 . . 
    U.mat.curr=matrix(Umat[(index+1):(index+ni),],ni)
    tseq=tvec[(index+1):(index+ni)]
    come.cur=visit.mat[i,] #the visit indicators for the ith subject
    
    if(dim(da)[1] != sum(da$delta)){
      invw.curr=pmax(SC(da$Y[i])/pmax(SC(tseq),0.000001),0.0001)
    } else {
      invw.curr=rep(1,dim(visit.mat)[2])
    }
    
    wt.curr=(da$delta[i]*(da$Y[i]>tseq)/invw.curr*come.cur)*wt.p[i]
    wt.curr2=(da$delta[i]/invw.curr)
    pseudo1=-apply(U.mat.curr*wt.curr,2,sum)
    pseudo2=2*apply(U.mat.curr*(da$Y[i]>tseq)*come.cur*wt.p[i]*tau0,2,sum)
    
    Y.reg=c(Y.reg,M,M)
    U.reg=rbind(U.reg,rbind(pseudo1,pseudo2))
    wt.reg.1=c(wt.reg.1,wt.curr) #lin et al방법의 점추정을 위한 weight
    wt.reg.add=c(wt.reg.add,1,1)
    wt.reg.2=c(wt.reg.2,wt.curr2) } #induced smoothing방법의 점추정을 위한 weight
  
  wt.reg = c(wt.reg.1,wt.reg.add)
    wt.reg.usc <- wt.reg
    
  # #evaluate the quantile model at tpseq
  # try=try(rq.wfit(U.reg,Y.reg,weights=wt.reg)$coef)
  # if (inherits(try, "try-error")) {gamma.fit=rep(NA,nbasis*(pc+2)) } else { # 20?=M(4)*5(number of coefficients: alpha0,1,beta1,2,3)
  #   gamma.fit=rq.wfit(U.reg,Y.reg,weights=wt.reg)$coef
  #   non.conv=(sum(abs(gamma.fit))>=200)|(is.na(sum(gamma.fit)))
  #   if(sum(abs(gamma.fit))>=5000){gamma.fit=gamma.fit*NA}}
  
  try2=try(rq.wfit(U.reg,Y.reg,weights=wt.reg.usc)$coef)
  if (inherits(try2, "try-error")) {
    gamma.fit2=rep(NA,nbasis2*(pc2+2))
      error_msg <- attr(try2, "condition")$message
      # error_log[[length(error_log) + 1]] <- paste("nsim =", nsim, "step2(rq.wfit) error", error_msg)
      cat("Error occurred in rq.wfit(unscaled), step2, for nsim =", nsim, ": ", error_msg, "\n")
    } else { # 20?=M(4)*5(number of coefficients: alpha0,1,beta1,2,3)
    gamma.fit2=rq.wfit(U.reg,Y.reg,weights=wt.reg.usc)$coef
    non.conv=(sum(abs(gamma.fit2))>=200)|(is.na(sum(gamma.fit2)))
    if(sum(abs(gamma.fit2))>=5000){gamma.fit2=gamma.fit2*NA}}
  
  # singular design matrix? check
  wx2 <- U.reg * sqrt(wt.reg.usc)
  print(wx2)
  dim(wx2) # 826 x 30
  qr(wx2)$rank # 20
  
    # tryCatch({
  #   gamma.fit2 <- rq.wfit(U.reg2, Y.reg, weights = wt.reg)$coef # rq.wfit 함수 실행 시도
  #   non.conv <- (sum(abs(gamma.fit2)) >= 200) | (is.na(sum(gamma.fit2))) # non.conv 및 gamma.fit 값 확인
  #   if (sum(abs(gamma.fit2)) >= 5000) {
  #     gamma.fit2 <- gamma.fit2*NA
  #   }
  # }, error = function(e) {
  #   # 에러 발생 시 에러 메시지를 리스트에 저장
  #   error_log[[length(error_log) + 1]] <- paste("nsim =", nsim, "step2(rq.wfit) error:", e$message)
  #   cat("Error occurred in rq.wfit(unscaled), step2, for nsim =", nsim, ": ", e$message, "\n")
  #   gamma.fit2 <- rep(NA, nbasis2 * (pc2 + 2)) # 에러 발생 시 gamma.fit 초기화: 20?=M(4)*5(number of coefficients: alpha0,1,beta1,2,3)
  # })
  
  splinef2=cbind(1,ns(t0vec,nbasis2-1))[which(abs(t0vec-tpseq)==min(abs(t0vec-tpseq)))[1],]
  alpha0.est2=as.vector(splinef2%*%gamma.fit2[1:nbasis2])
  alpha1.est2=as.vector(splinef2%*%gamma.fit2[(nbasis2+1):(2*nbasis2)])
  beta2=c()
  for(i in 1:pc2){
    beta2=cbind(beta2, as.vector(splinef2%*%gamma.fit2[((i+1)*nbasis2+1):((i+2)*nbasis2)]))}
  temp.est2=data.frame(alpha0.est2,alpha1.est2,beta2)
  names(temp.est2)=c("alpha0.est","alpha1.est",sprintf("beta%d.est",1:pc2))
  
    # # checking the values
    # length(gamma.fit)
    # length(gamma.fit2)
    dim(Umat%*%gamma.fit) == dim(Umat2%*%gamma.fit2)
    View(round(Umat%*%gamma.fit,4))
    View(round(Umat2%*%gamma.fit2,4))
    
    # visit.vec <- as.vector(t(visit.mat))
    # 
    # cat(ifelse(sum(round(Umat%*%gamma.fit*visit.vec,4) == round(Umat2%*%gamma.fit2*visit.vec,4)) == length(visit.vec),
    #        "Eq_9 the same between scaled and unscaled \n",
    #        "err in Eq_9 \n")) # Eq_9 <- Umat%*%gamma.fit # 5.16 추가
    # # Eq_9 (compared only "visited" in visit.mat) should turn out to be the same (YES!)
    
    fitted_qt_fpca = list(coef=temp.est, gamma.fit = gamma.fit, data=dataall, wt.reg.scaled = wt.reg.sc, 
                          coef2=temp.est2, gamma.fit2 = gamma.fit2, data2=dataall2, wt.reg.unscaled = wt.reg.usc, 
                          da=da, come.cur = come.cur, wt.reg2 = wt.reg.2)
    
    
    dim(wx) # 826 x 30
    qr(wx)$rank # 20
    dim(wx2) # 826 x 30
    qr(wx2)$rank # 20
    
    range(wx)
    range(wx2)
    
    # return(fitted_qt_fpca)
  # return(list(coef=temp.est, gamma.fit = gamma.fit, data=dataall,
  #             coef2=temp.est2, gamma.fit2 = gamma.fit2, data2=dataall2,
  #             da=da, come.cur = come.cur, wt.reg = wt.reg, wt.reg.2 = wt.reg.2))

# }


###########################################################################
## For induced smoothing
###########################################################################

objectF = function(theta){
  
  theta = as.matrix(theta)
  
  result = (1/N)*t(U*come) %*% (W*(pnorm((U%*%theta - logY)/sqrt(diag(U%*%H%*%t(U)))))-tau0)
  result = as.vector(result)
  
}

###########################################################################
## simulation
###########################################################################

error_log = list()
set.seed(1)
NN = 1
# nbasis = 4
M = 1E8
L = 6
N = 400
tau0 = 0.5
t0vec = seq(2, L, 1)

# vauc.fpca = c()
# vmean.fpca = c()

coef_lin=coef_lin2=c()
coef_IS=coef_IS2=c()

start_time <- Sys.time()

qrl=c()

for (nsim in 1:NN) {
  
  #############################
  ## estimate
  #############################
  ## generate data
  dat = da.genrt(N)$data
  
  ## fpca
  fpcaobj = fpca.fit(dat)
  pc = fpcaobj$pc
  
  for (tpseq in c(2)){ ## c(2, 2.5, 3, 3.5, 4)  
      
      # tpseq = 2
    
      # fpca
      fit_qt_fpca = fpca.analysis(fpcaobj)
      # print(fit_qt_fpca$gamma.fit)
      coef_lin[[nsim]] <- fit_qt_fpca$coef
      print(fit_qt_fpca$coef)
      # print(fit_qt_fpca$gamma.fit2)
      coef_lin2[[nsim]] <- fit_qt_fpca$coef2
      print(fit_qt_fpca$coef2)
      
      #################################
      ## quantile residual lifetime
      #################################
      IDall=unique(fit_qt_fpca$data[which(fit_qt_fpca$data$Time>=tpseq),1]) # ID of subjects: tpseq <= observed time (zero???)
      qt.fpca=qt.fpca2=c()
      qt.true=c()
      
      ##true (true residual lifetime)
      data.id=dat[!duplicated(dat[,1]),] # unique subject ID list 
      qt.true=cbind(IDall,pmin(data.id[data.id[,1]%in%IDall,5]-tpseq,L-tpseq)) # min(trueTimes-tpseq,L-tpseq) : trueTimes=potential event times
      
      for(k in IDall){
        
        ## fpca
        ND.fpca=fit_qt_fpca$data[fit_qt_fpca$data[,1]==k,]
        ND.fpca=ND.fpca[which(abs(ND.fpca$tij-tpseq)==min(abs(ND.fpca$tij-tpseq)))[1],] # (???)
        # qt=exp(fit_qt_fpca$coef[1]+fit_qt_fpca$coef[2]*ND.fpca$X
        #        +fit_qt_fpca$coef[3]*ND.fpca$Z1+fit_qt_fpca$coef[4]*ND.fpca$Z2+fit_qt_fpca$coef[5]*ND.fpca$Z3)
        qt = eval(parse(text=paste0(
          "exp(fit_qt_fpca$coef[1] + fit_qt_fpca$coef[2] * ND.fpca$X",
          paste(sprintf("+ fit_qt_fpca$coef[%d] * ND.fpca$Z%d", 3:(2 + pc), 1:pc), collapse = ''),
          ")"
        )))
        #qt=ifelse(qt<range(qt.true[,2])[2]&qt>range(qt.true[,2])[1],qt,0)
        # qt=exp(fit_qt_fpca$coef[1]+fit_qt_fpca$coef[2]*ND.fpca$X
        #        +fit_qt_fpca$coef[3]*ND.fpca$Z1+fit_qt_fpca$coef[4]*ND.fpca$Z2+fit_qt_fpca$coef[5]*ND.fpca$Z3) #(why same???)
        qt=ifelse(qt < range(qt.true[,2])[2] & qt > range(qt.true[,2])[1], qt, NA) # range(qt.true[,2])[2] = max of true residual time, range(qt.true[,2])[1] = min of true residual time
        qt.fpca=rbind(qt.fpca,c(k,qt))
        # qt.fpcaauc=rbind(qt.fpcaauc,c(k,qt))
        
        
        ND.fpca2=fit_qt_fpca$data2[fit_qt_fpca$data2[,1]==k,]
        ND.fpca2=ND.fpca2[which(abs(ND.fpca2$tij-tpseq)==min(abs(ND.fpca2$tij-tpseq)))[1],] # (???)
        # qt2=exp(fit_qt_fpca$coef2[1]+fit_qt_fpca$coef2[2]*ND.fpca$X
        #        +fit_qt_fpca$coef[3]*ND.fpca$Z1+fit_qt_fpca$coef[4]*ND.fpca$Z2+fit_qt_fpca$coef[5]*ND.fpca$Z3)
        qt2 = eval(parse(text=paste0(
          "exp(fit_qt_fpca$coef2[1] + fit_qt_fpca$coef2[2] * ND.fpca2$X",
          paste(sprintf("+ fit_qt_fpca$coef2[%d] * ND.fpca2$Z%d", 3:(2 + pc), 1:pc), collapse = ''),
          ")"
        )))
        qt2=ifelse(qt2 < range(qt.true[,2])[2] & qt2 > range(qt.true[,2])[1], qt2, NA) # range(qt.true[,2])[2] = max of true residual time, range(qt.true[,2])[1] = min of true residual time
        qt.fpca2=rbind(qt.fpca2,c(k,qt2))
        
      }
      
      qrl[[nsim]] <- cbind(qt.true, qt.fpca, qt.fpca2)

      ### Induced smoothing(scaled):
      # obj <- fpcaobj
      U=fpcaobj$Umat # U: covariate W,Z에 fractional basis function을 곱한 형태
      ni.vec = fpcaobj$ni.vec
      t=fpcaobj$tvec
      come=fit_qt_fpca$come.cur #visit indicator
      W=fit_qt_fpca$wt.reg.2 #censoring weight
      logY=fpcaobj$Yvec #log of Y(observed failure time)
      nc=length(fit_qt_fpca$gamma.fit)
      H = diag(1/(nc*N), nc, nc) # 12 개수 변경
      nbasis=fpcaobj$nbasis # should follow the result from previous function fit_qt_fpca
      pc=fpcaobj$pc # should follow the result from previous function fit_qt_fpca
      
      thetastart <- fit_qt_fpca$gamma.fit
      
      if( sum(is.na(thetastart)) >= 1){
        
        coef_IS[[nsim]] <- NA
        
      } else {
        
        is.fit = dfsane(thetastart, objectF, control = list(tol=5*1e-03))
        idw.fit = is.fit$par  ## 최적화된 변수 값들
        
        tryCatch({
          if (is.fit$convergence == 0){  ## 수렴했는지
            idw.splinef=cbind(1,ns(t0vec,nbasis-1))[which(abs(t0vec-tpseq)==min(abs(t0vec-tpseq)))[1],]
            idw.alpha0.est=as.vector(idw.splinef%*%idw.fit[1:nbasis])
            idw.alpha1.est=as.vector(idw.splinef%*%idw.fit[(nbasis+1):(2*nbasis)])
            beta=c()
            for(i in 1:pc){
              beta=cbind(beta, as.vector(idw.splinef%*%idw.fit[((i+1)*nbasis+1):((i+2)*nbasis)]))}
            idw.temp.est=data.frame(idw.alpha0.est,idw.alpha1.est,beta)
            names(idw.temp.est)=c("idw.alpha0.est","idw.alpha1.est",sprintf("idw.beta%d.est",1:pc))
          } else {   
            idw.temp.est = c(NA,NA,NA,NA,NA) #idw.temp.est: induced smoothing 방법으로 구한 점추정치
          }
        },error=function(e){
          idw.temp.est = c(NA,NA,NA,NA,NA)
          # error_log[[length(error_log) + 1]] <- paste("nsim =", nsim, "step 3(idw.temp.est):", e$message)
          cat("Error occurred in is.fit(scaled), step3, for nsim =", nsim, ": ", e$message, "\n")
        })
        
        coef_IS[[nsim]] <- idw.temp.est
        print(idw.temp.est)
        
      }
  
  
      ### Induced smoothing (unscaled):
      # obj <- fpcaobj
      U=fpcaobj$Umat2 # U: covariate W,Z에 fractional basis function을 곱한 형태
      ni.vec = fpcaobj$ni.vec
      t=fpcaobj$tvec
      come=fit_qt_fpca$come.cur #visit indicator
      W=fit_qt_fpca$wt.reg.2 #censoring weight
      logY=fpcaobj$Yvec #log of Y(observed failure time)
      nc2=length(fit_qt_fpca$gamma.fit2)
      H = diag(1/(nc2*N), nc2, nc2) # 12 개수 변경
      nbasis=fpcaobj$nbasis # should follow the result from previous function fit_qt_fpca
      pc=fpcaobj$pc # should follow the result from previous function fit_qt_fpca
      
      thetastart2 <- fit_qt_fpca$gamma.fit2
      
      if( sum(is.na(thetastart2)) >= 1){
        
        coef_IS2[[nsim]] <- NA
        
      } else {
        
        is.fit = dfsane(thetastart2, objectF, control = list(tol=5*1e-03))
        idw.fit = is.fit$par  ## 최적화된 변수 값들
        
        tryCatch({
          if (is.fit$convergence == 0){  ## 수렴했는지
            idw.splinef=cbind(1,ns(t0vec,nbasis-1))[which(abs(t0vec-tpseq)==min(abs(t0vec-tpseq)))[1],]
            idw.alpha0.est=as.vector(idw.splinef%*%idw.fit[1:nbasis])
            idw.alpha1.est=as.vector(idw.splinef%*%idw.fit[(nbasis+1):(2*nbasis)])
            beta=c()
            for(i in 1:pc){
              beta=cbind(beta, as.vector(idw.splinef%*%idw.fit[((i+1)*nbasis+1):((i+2)*nbasis)]))}
            idw.temp.est2=data.frame(idw.alpha0.est,idw.alpha1.est,beta)
            names(idw.temp.est2)=c("idw.alpha0.est","idw.alpha1.est",sprintf("idw.beta%d.est",1:pc))
          } else {   
            idw.temp.est2 = c(NA,NA,NA,NA,NA) #idw.temp.est: induced smoothing 방법으로 구한 점추정치
          }
        },error=function(e){
          idw.temp.est2 = c(NA,NA,NA,NA,NA)
          # error_log[[length(error_log) + 1]] <- paste("nsim =", nsim, "step 3(idw.temp.est):", e$message)
          cat("Error occurred in is.fit(unscaled), step3, for nsim =", nsim, ": ", e$message, "\n")
        })
        
        coef_IS2[[nsim]] <- idw.temp.est2
        print(idw.temp.est2)
        
      }
      
      }
  
}

end_time <- Sys.time()
end_time - start_time



coef_lin; coef_lin2
coef_IS; coef_IS2
qrl

# 리스트의 각 데이터 프레임에 대해 모든 요소의 개수를 세고, 합산합니다.
total_elements_coef_lin <- sum(sapply(coef_lin, function(df) nrow(df) * ncol(df)))
total_elements_coef_lin2 <- sum(sapply(coef_lin2, function(df) nrow(df) * ncol(df)))
sum_na_coef_lin <- sum(sapply(coef_lin, function(x) sum(is.na(x))))
sum_na_coef_lin2 <- sum(sapply(coef_lin2, function(x) sum(is.na(x))))

total_elements_coef_lin - sum_na_coef_lin;
total_elements_coef_lin2 - sum_na_coef_lin2

length(coef_IS) - sum(is.na(coef_IS));
length(coef_IS2) - sum(is.na(coef_IS2))

# temp.fpca$pc; temp.fpca$nbasis
# temp.fpca$cv

# setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2/Results-IS")
# write.csv(cbind(fit.fpca$coef, idw.temp.est), "coefs.csv")

# temp.fpca2 = fpca.fit2(data=data, pc=4)
#   for (tpseq in c(2, 3, 4)) { ## c(2, 2.5, 3, 3.5, 4)
#     fit.fpca2 = fpca.analysis(temp.fpca2)
#     print(fit.fpca2$gamma.fit) # all the same for different tpseq's
#     print(tpseq)
#     print(fit.fpca2$coef) # different
#   }

# eq9_m <- fit.fpca$Eq_9
# eq9 <- fit.fpca2$Eq_9
# 
# # coefficient x each covariate
# n <- length(eq9_m); vec1 <- as.vector(rep(1, n))
# data1 <- eval(parse(text=paste0("data.frame(vec1, fit.fpca$data$X, ", # model selected by fpca.fit function
#                                 paste(sprintf("fit.fpca$data$Z%d",1:pc), collapse=', '), ")")))
# eq7_m <- as.matrix(data1) %*% t(as.matrix(fit.fpca$coef)) 
# sum(eq9_m == eq7_m) == n/length(t0vec)
# 
# pc = 4
# data2 <- eval(parse(text=paste0("data.frame(vec1, fit.fpca2$data$X, ", # when the number of pc is set different
#                                 paste(sprintf("fit.fpca2$data$Z%d",1:pc), collapse=', '), ")")))
# eq7 <- as.matrix(data2) %*% t(as.matrix(fit.fpca2$coef))
# sum(eq9 == eq7) == n/length(t0vec) 
# 
# # View(eq9_m); View(eq9)
# View(cbind(round(eq9_m,6), round(eq9,6)))
# sum(round(eq9_m,6) == round(eq9,6))
# which(round(eq9_m,6) == round(eq9,6))
# View(temp.fpca$Umat); View(temp.fpca2$Umat)


# ### Induced smoothing:
# obj = temp.fpca2
# U=obj$Umat # U: covariate W,Z에 fractional basis function을 곱한 형태
# ni.vec = obj$ni.vec
# t=obj$tvec
# come=fit.fpca2$come.cur #visit indicator
# W=fit.fpca2$wt.reg.2 #censoring weight
# logY=obj$Yvec #log of Y(observed failure time)
# nc=length(fit.fpca2$gamma.fit)
# H = diag(1/(nc*N), nc, nc) # 12 개수 변경
# nbasis=obj$nbasis # should follow the result from previous function fit.fpca
# pc=obj$pc # should follow the result from previous function fit.fpca
# 
# coef_lin[[nsim]] <- fit.fpca2$coef
# thetastart <- fit.fpca2$gamma.fit
# 
# if( sum(is.na(thetastart)) >= 1){
#   
#   coef_IS[[nsim]] <- NA
#   
# } else {
#   
#   is.fit = dfsane(thetastart, objectF, control = list(tol=5*1e-03))
#   idw.fit = is.fit$par  ## 최적화된 변수 값들
#   
#   tryCatch({
#     if (is.fit$convergence == 0){  ## 수렴했는지
#       idw.splinef=cbind(1,ns(t0vec,nbasis-1))[which(abs(t0vec-tpseq)==min(abs(t0vec-tpseq)))[1],]
#       idw.alpha0.est=as.vector(idw.splinef%*%idw.fit[1:nbasis])
#       idw.alpha1.est=as.vector(idw.splinef%*%idw.fit[(nbasis+1):(2*nbasis)])
#       beta=c()
#       for(i in 1:pc){
#         beta=cbind(beta, as.vector(idw.splinef%*%idw.fit[((i+1)*nbasis+1):((i+2)*nbasis)]))}
#       idw.temp.est=data.frame(idw.alpha0.est,idw.alpha1.est,beta)
#       names(idw.temp.est)=c("idw.alpha0.est","idw.alpha1.est",sprintf("idw.beta%d.est",1:pc))
#     } else {   
#       idw.temp.est = c(NA,NA,NA,NA,NA) #idw.temp.est: induced smoothing 방법으로 구한 점추정치
#     }
#   },error=function(e){
#     idw.temp.est = c(NA,NA,NA,NA,NA)
#     
#   })
#   
#   coef_IS[[nsim]] <- idw.temp.est
#   print(idw.temp.est)
#   
# }
# 
# setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2/Results-IS")
# write.csv(cbind(fit.fpca2$coef, idw.temp.est), "coefs.csv")
# 
# coef_lin
# coef_IS

# write.csv(coef_lin, "coef_lin_0508.csv")
# write.csv(coef_IS, "coef_IS_0508.csv")

# my_list <- coef_lin
# 
# max_length <- max(sapply(my_list, length)) # 가장 긴 벡터의 길이 계산
# 
# coef.lin = matrix(NA, ncol = max_length, nrow = length(my_list))
# coef.lin = data.frame()
# # 리스트의 각 요소에 대해 반복
# for (i in 1:length(my_list)) {
#   current_length <- length(my_list[[i]])
#   print(paste("Current length for element", i, ":", current_length))  # 현재 요소의 길이 출력
#   
#   # 현재 요소 길이가 max_length보다 짧으면 NA로 채워넣기
#   if (current_length < max_length) {
#     padded_vector <- c(my_list[[i]], rep(NA, max_length - current_length))
#   } else {
#     padded_vector <- my_list[[i]]
#   }
#   
#   print(paste("Length of padded vector for element", i, ":", length(padded_vector)))  # 패딩된 벡터의 길이 출력
#   
#   # 행렬의 i번째 행에 padded_vector 할당
#   coef.lin[i, ] <- padded_vector  # 명시적으로 전체 열 범위 지정
# }
# 
# coef.lin <- as.data.frame(coef.lin)
# 
# # 각 요소를 길이에 맞춰 NA로 확장
# expanded_list <- lapply(my_list, function(x) {
#   length(x) <- max_length
#   x
# })
# 
# # 전체 데이터 프레임을 한 번에 생성
# rows <- vector("list", max_length)
# for (i in 1:max_length) {
#   # 각 행에 대해, 리스트의 각 벡터에서 i번째 값을 추출
#   rows[[i]] <- sapply(expanded_list, `[`, i)
# }
# 
# # 모든 행을 하나의 데이터 프레임으로 결합
# combined_df <- data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
# 
# # 데이터 프레임의 컬럼 이름 설정
# names(combined_df) <- names(my_list)
# 
# # # 데이터 프레임 확인
# # print(combined_df)
# 
# # CSV 파일로 저장
# write.csv(df, "coef_lin_0508.csv", row.names = FALSE)
# 
# ##############
# my_list <- coef_IS
# 
# # 가장 긴 벡터의 길이 계산
# max_length <- max(sapply(my_list, length))
# 
# # 각 요소를 길이에 맞춰 NA로 확장
# expanded_list <- lapply(my_list, function(x) {
#   length(x) <- max_length
#   x
# })
# 
# # 전체 데이터 프레임을 한 번에 생성
# rows <- vector("list", max_length)
# for (i in 1:max_length) {
#   # 각 행에 대해, 리스트의 각 벡터에서 i번째 값을 추출
#   rows[[i]] <- sapply(expanded_list, `[`, i)
# }
# 
# # 모든 행을 하나의 데이터 프레임으로 결합
# combined_df <- data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
# 
# # 데이터 프레임의 컬럼 이름 설정
# names(combined_df) <- names(my_list)
# 
# # # 데이터 프레임 확인
# # print(combined_df)
# 
# # CSV 파일로 저장
# write.csv(df, "coef_IS_0508.csv", row.names = FALSE)
