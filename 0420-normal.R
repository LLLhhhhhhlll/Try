rm(list = ls())
library(MASS)
library(data.table)
library(compiler)
library(condMVNorm)
library(pbapply)
library(tmvtnorm)
library(doParallel)
library(foreach)
library(doSNOW)
library(tcltk)

da_gen <- function(sim,setting,J){
  set.seed(20250420+sim)
  if(setting=="NULL"){mu_treatment <- c(0,0,0,0)
  }else if(setting=="H1"){mu_treatment <- c(0.50,0.50,0.50,0.50)
  }else if(setting=="peak"){mu_treatment <- c(0.25,0.25,0.25,0.50)
  }else if(setting=="line"){mu_treatment <- c(0.20,0.30,0.40,0.50)}
  
  if(J == 3) {
    t <- c(1,2,3)/3
    low_bround <- c(0.000000, 1.038459, 2.290364)
    up_bround <- c(Inf, 3.301585, 2.290364)
    mu_control <- 0
    sigma_control <- 1
    sigma_treatment <- c(1,1,1,1)
    unit_n <- 30
  }else if(J == 4) {
    t <- c(1,2,3,4)/4
    low_bround <- c(0.0000000, 0.1267247, 1.4593383, 2.251983)
    up_bround <- c(Inf, 4.013738, 2.826973, 2.251983)
    mu_control <- 0
    sigma_control <- 1
    sigma_treatment <- c(1,1,1,1)
    unit_n <- 25
  }else if(J == 5) {
    t <- c(1,2,3,4,5)/5
    low_bround <- c( 0.0000000, -0.6733603,  0.7312345,  1.6214603, 2.228605)
    up_bround <- c(Inf, 4.609266, 3.263919, 2.639797, 2.228605)
    mu_control <- 0
    sigma_control <- 1
    sigma_treatment <- c(1,1,1,1)
    unit_n <- 22
  }
  mu_control <- 0
  sigma_control <- 1
  sigma_treatment <- c(1,1,1,1)
  K <- length(mu_treatment)
  J <- length(low_bround)
  n <- t*length(t)*unit_n
  I <- function(j,index)  {1/(sigma_treatment[index]**2/n[j]+sigma_control**2/n[j])}
  Z <- function(da_treat,da_control,j,index) {(mean(da_treat)-mean(da_control))*sqrt(I(j,index))}
  all_da <- mvrnorm(n[J],mu=c(mu_control,mu_treatment),
                    Sigma = diag(c(sigma_control,sigma_treatment)**2,nrow = K+1,ncol = K+1))
  select_index <-  which.max(mapply(function(x) Z(all_da[c(1:n[1]),x+1],all_da[c(1:n[1]),1],j=1,index=x),x=1:K))
  
  inter_decide <- function(j,select_index,da){
    mean_temp <-  Z(da[c(1:n[j]),select_index+1],da[c(1:n[j]),1],j,index=select_index)
    if(mean_temp<low_bround[j]){"f"}else if(mean_temp>up_bround[j]){"e"} else{"c"}
  }
  inter_outcome <- as.character()
  for (i in 1:J) {
    inter_outcome[i] <- inter_decide(j=i,select_index,da=all_da)
  }
  end_stage <- first(which(inter_outcome %in% c("e","f")))
  stage2_est <- mean(all_da[c((n[1]+1):n[2]),select_index+1])-mean(all_da[c((n[1]+1):n[2]),1])
  
  naive <- mapply(function(x) mean(all_da[c(1:n[1]),x+1])-mean(all_da[c(1:n[1]),1]), 1:K)
  naive[select_index] <- mean(all_da[c(1:n[end_stage]),select_index+1])-mean(all_da[c(1:n[end_stage]),1])
  
  x_bar <- apply(all_da[1:n[1],], 2, mean)
  x_bar[1] <- mean(all_da[1:n[end_stage],1])
  x_bar[select_index+1] <- mean(all_da[1:n[end_stage],select_index+1])
  
  sigma_hat <- sqrt(apply(all_da[1:n[1],], 2, var))
  sigma_hat[1] <- sqrt(var(all_da[1:n[end_stage],1]))
  sigma_hat[select_index+1] <- sqrt(var(all_da[1:n[end_stage],select_index+1]))
  
  S <- function(i,j){(mean(all_da[1:n[j],i+1])-mean(all_da[1:n[j],1]))*I(j,i)}
  s_all <- c(mapply(S, i=1:K,1),S(select_index,end_stage))
  
  c(end_stage,select_index,naive,x_bar,s_all,sigma_hat,stage2_est,naive[select_index])
}
da_genc <- cmpfun(da_gen)

##### MI & SI #####
mu_cov_f <- function(control_eff,treat_eff,select_index,
                     control_sig,treat_sig,K,J,n){
  V <- function(sig,j) sig**2/n[j]
  I <- function(i,j){1/(V(treat_sig[i],j)+V(control_sig,j))} 
  A <- diag(1,K+J-1)
  diag(A[1:(K-1),1:(K-1)]) <- -1 
  A[1:(K-1),K] <- 1
  
  theta_v <- diag(nrow = K+J-1 )
  diag(theta_v[1:(K-1) ,1:(K-1)]) <- sqrt(I((1:K)[-select_index],1)) 
  diag(theta_v[K:(K+J-1),K:(K+J-1)]) <- sqrt(I(select_index,1:J)) 
  
  cov <-cbind(matrix(V(control_sig,1),nrow = K+J-1,ncol = K-1),
              mapply(function(x) matrix(V(control_sig,x),nrow = K+J-1,ncol = 1),x=1:J))
  diag(cov[1:(K-1),1:(K-1)]) <- diag(cov[1:(K-1),1:(K-1)])+V(treat_sig[-select_index],1)
  for (i in K:(K+J-1)) {
    for (j in K:(K+J-1)) {
      cov[i,j] <- cov[i,j]+V(treat_sig[select_index],j-(K-1))
    }
  }
  for (i in 1:(K+J-1)) {
    for (j in 1:(K+J-1)) {
      if(i<=K-1 & j<=K-1) {
        cov[i,j] <- cov[i,j]*sqrt(I((1:K)[-select_index][i],1)*I((1:K)[-select_index][j],1))
      }else if(i<=K-1 & between(j,K,K+J-1)){
        cov[i,j] <-cov[i,j]*sqrt(I((1:K)[-select_index][i],1)*I(select_index,j-(K-1)))
      }else if(between(i,K,K+J-1) & between(j,K,K+J-1)){
        cov[i,j] <-cov[i,j]*sqrt(I(select_index,i-(K-1))*I(select_index,j-(K-1)))
      }
    }}
  cov[lower.tri(cov)] <- t(cov)[lower.tri(cov)]
  mu= A %*% theta_v %*% as.matrix(c(treat_eff[-select_index]-control_eff,
                                    rep(treat_eff[select_index]-control_eff,J)))
  cov=A %*% cov %*% t(A)
  if(any(diag(cov)<0)){diag(cov)[which(diag(cov)<0)] <- runif(sum(diag(cov)<0))}
  list(mu=mu,cov=cov)
}
bias_theta  <- function(select_index,theta,sigma_control,sigma_treatment,K,J){
  if(J == 3) {
    t <- c(1,2,3)/3
    low_bround <- c(0.000000, 1.038459, 2.290364)
    up_bround <- c(Inf, 3.301585, 2.290364)
    mu_control <- 0
    unit_n <- 30
  }else if(J == 4) {
    t <- c(1,2,3,4)/4
    low_bround <- c(0.0000000, 0.1267247, 1.4593383, 2.3056892)
    up_bround <- c(Inf, 4.013738, 2.826973, 2.251983)
    mu_control <- 0
     unit_n <- 25
  }else if(J == 5) {
    t <- c(1,2,3,4,5)/5
    low_bround <- c( 0.0000000, -0.6733603,  0.7312345,  1.6214603, 2.228605)
    up_bround <- c(Inf, 4.609266, 3.263919, 2.639797, 2.228605)
    mu_control <- 0
    unit_n <- 22
  }
  n <- t*length(t)*unit_n
  I <- function(j,index)  {1/(sigma_treatment[index]**2/n[j]+sigma_control**2/n[j])}
  mu_cov <- mu_cov_f(control_eff=0,treat_eff=theta,
                     control_sig=sigma_control,treat_sig=sigma_treatment,
                     select_index,K,J,n)
  mu <- mu_cov$mu; cov <- mu_cov$cov
  prob_Q <- pmvnorm(lower = c(rep(0,K-1),low_bround[1]),upper = c(rep(Inf,K-1),up_bround[1]),
                    mean = mu[1:K,],sigma = cov[1:K,1:K])[1]
  
  bias <- rep(NA,K)
  if(!is.na(prob_Q)){
    mtm <- mtmvnorm(mean = mu[1:K,],sigma = cov[1:K,1:K],
                    lower = c(rep(0,K-1),low_bround[1]),upper = c(rep(Inf,K-1),up_bround[1]),
                    doComputeVariance=FALSE)$tmean
    bias[-select_index] <- (mtm[K]-mtm[-K])/sqrt(I(1,c(1:K)[-select_index]))-theta[-select_index]
    E_select <- function(stage){
      mu_t <- mu[1:(K+stage-1),]
      cov_t <- cov[1:(K+stage-1),1:(K+stage-1)]
      lower1 <- c(rep(0,K-1),low_bround[1:(stage-1)],-Inf); upper1 <- c(rep(Inf,K-1),up_bround[1:(stage-1)],low_bround[stage])
      lower2 <- c(rep(0,K-1),low_bround[1:(stage-1)],up_bround[stage]);upper2 <- c(rep(Inf,K-1),up_bround[1:(stage-1)],Inf)
      (pmvnorm(lower = lower1,upper = upper1,mean = mu_t,sigma = cov_t)[1]*
          mtmvnorm(mean = mu_t,sigma = cov_t,lower = lower1,upper = upper1,doComputeVariance=FALSE)$tmean[K+stage-1]+
          pmvnorm(lower = lower2,upper = upper2,mean = mu_t,sigma = cov_t)[1]*
          mtmvnorm(mean = mu_t,sigma = cov_t,lower = lower2,upper = upper2,doComputeVariance=FALSE)$tmean[K+stage-1])/
        sqrt(I(stage,select_index))
    }
    bias[select_index] <- sum(mapply(E_select, stage=2:J))/prob_Q-theta[select_index]
  }
  bias
}
MI <- function(initial,select_index,max_iterations=20,tol=0.001,sigma_control,sigma_treatment,K,J){
  conver <- 0
  theta_hat <- initial
  for (i in 1:max_iterations) {
    bias_temp <- bias_theta(select_index,theta=theta_hat,sigma_control,sigma_treatment,K,J)
    if(is.na(sum(bias_temp)) | is.infinite(sum(bias_temp))){
      theta_hat <- initial; conver <- -1; break()
    }else{
      theta_hat_new <- initial-bias_temp
    }
    euc.dis <- sqrt(sum((theta_hat - theta_hat_new)^2))
    if(euc.dis<=tol){
      theta_hat <- theta_hat_new; conver <- 1; break()
    }else{
      theta_hat <- theta_hat_new
    }
  }
  c(conver,theta_hat[select_index])
}
SI <- function(initial,select_index,sigma_control,sigma_treatment,K,J){
  (initial-bias_theta(select_index,initial,sigma_control,sigma_treatment,K,J))[select_index]
}  
MIc <- cmpfun(MI)

##### MUE #####
est_MUE <- function(MLE,end_stage,select_index,sub,sigma_control,sigma_treatment,K,J){
  if(J == 3) {
    t <- c(1,2,3)/3
    low_bround <- c(0.000000, 1.038459, 2.290364)
    up_bround <- c(Inf, 3.301585, 2.290364)
    mu_control <- 0
    unit_n <- 30
  }else if(J == 4) {
    t <- c(1,2,3,4)/4
    low_bround <- c(0.0000000, 0.1267247, 1.4593383, 2.3056892)
    up_bround <- c(Inf, 4.013738, 2.826973, 2.251983)
    mu_control <- 0
    unit_n <- 25
  }else if(J == 5) {
    t <- c(1,2,3,4,5)/5
    low_bround <- c( 0.0000000, -0.6733603,  0.7312345,  1.6214603, 2.228605)
    up_bround <- c(Inf, 4.609266, 3.263919, 2.639797, 2.228605)
    mu_control <- 0
    unit_n <- 22
  }
  
  n <- t*length(t)*unit_n
  I <- function(j,index)  {1/(sigma_treatment[index]**2/n[j]+sigma_control**2/n[j])}
  Z_end <- MLE[select_index]*sqrt(I(end_stage,select_index))
  cond_pvalue <- function(theta,select_index,end_stage,Z_end){
    mu_cov <- mu_cov_f(control_eff=0,treat_eff=theta,select_index,
                       control_sig=sigma_control,treat_sig=sigma_treatment,K,J,n)
    mu <- mu_cov$mu; cov <- mu_cov$cov
    
    prob_Q <- mvtnorm::pmvnorm(lower = c(rep(0,K-1),low_bround[1]),
                               upper = c(rep(Inf,K-1),up_bround[1]),
                               mean = mu[1:K,],sigma = cov[1:K,1:K])[1]
    p_value <- function(end_stage,select_index,Z_end){
      P_stage <- function(stage,stop,Z_end){
        new_mu_s <- mu[1:(K+stage-1),]; new_cov_s=cov[1:(K+stage-1),1:(K+stage-1)]
        if(stop==T){
          P <-  mvtnorm::pmvnorm(lower = c(rep(0,K-1),low_bround[1:(stage-1)],Z_end),
                                 upper = c(rep(Inf,K-1),up_bround[1:(stage-1)],Inf),
                                 mean = new_mu_s,sigma= new_cov_s)[1]      
        }else if(stop==F){
          P <-   mvtnorm::pmvnorm(lower = c(rep(0,K-1),low_bround[1:(stage-1)],up_bround[stage]),
                                  upper = c(rep(Inf,K-1),up_bround[1:(stage-1)],Inf),
                                  mean = new_mu_s,sigma= new_cov_s)[1]
        }
        if(is.na(P)){0}else{P}
      }  
      if(end_stage==2){
        numerator <- P_stage(stage=end_stage,stop=T,Z_end)
      }else if(end_stage > 2){
        numerator <- sum(mapply(P_stage, stage=2:(end_stage-1),stop=F,Z_end),
                         P_stage(stage=end_stage,stop=T,Z_end)) 
      }
      value <- 1
      if((!is.na(numerator)) & (!is.na(prob_Q))){
        if(prob_Q > 0){
          if(numerator/prob_Q >=1){ value <- 1
          }else{value <- numerator/prob_Q}
        }
      }
      value
    }
    p_value(end_stage,select_index,Z_end)
  }  
  cond_pvaluec <- cmpfun(cond_pvalue)
  P_value_sub <- function(select_index,sub="MLE",x){
    theta <- rep(x,K)
    if(sub=="MLE"){
      theta[-select_index] <- MLE[-select_index]
    }else{
      theta[-select_index] <- rep(0,K-1)
    }
    cond_pvaluec(theta,select_index,end_stage,Z_end)
  }
  a <- optimise(function(x) (P_value_sub(select_index,sub,x)-1/2)^2,c(-1,1),tol=0.001)
  c(a$minimum,a$objective)
}

##### UMVCUE #####
mu_cov_UMVCUE <- function(control_eff,treat_eff,control_sig,treat_sig,
                          select_index,end_stage,K,J,n){
  V <- function(sig,j) sig**2/n[j]
  I <- function(i,j){1/(V(treat_sig[i],j)+V(control_sig,j))} 
  mu <- rep(0,end_stage+2)
  cov <-cbind(matrix(nrow = end_stage+2,ncol = end_stage-1),
              matrix(sqrt(I(select_index,1)),nrow = end_stage+2,ncol = 1),
              matrix(sqrt(I(select_index,end_stage)),nrow = end_stage+2,ncol = 2))
  for(i in 1:(end_stage-1)){
    for(j in i:(end_stage-1)){
      cov[i,j] <- sqrt(I(select_index,i)/I(select_index,j))
    }
  }
  for(i in 1:(end_stage+2)){
    for(j in (end_stage:(end_stage+2))){
      if(i <= end_stage-1){cov[i,j] <- cov[i,j] * sqrt(I(select_index,i))}
      if(i == end_stage){cov[i,j] <- cov[i,j] * sqrt(I(select_index,1))}
      if(i >= end_stage+1){cov[i,j] <- cov[i,j] * sqrt(I(select_index,end_stage))}
    }
  }
  for(i in 1:(end_stage-1)){
    for(j in (end_stage:(end_stage+2))){
      if(j == end_stage){cov[i,j] <- -cov[i,j] * V(control_sig,i)}
      if(j == end_stage+1){cov[i,j] <- -cov[i,j] * V(control_sig,end_stage)}
      if(j == end_stage+2){cov[i,j] <- cov[i,j] * V(treat_sig[select_index],end_stage)}
    }
  }
  diag(cov[(end_stage):(end_stage+2),(end_stage):(end_stage+2)]) <- diag(
    cov[(end_stage):(end_stage+2),(end_stage):(end_stage+2)])*
    c(V(control_sig,1),V(control_sig,end_stage),V(treat_sig[select_index],end_stage))
  cov[end_stage,(end_stage+1):(end_stage+2)] <- cov[end_stage,(end_stage+1):(end_stage+2)]*
    c(V(control_sig,end_stage),0)
  cov[end_stage+1,end_stage+2] <- 0
  cov[lower.tri(cov)] <- t(cov)[lower.tri(cov)]
  list(mu=mu,cov=cov)
}
UMVCUE_mc <- function(X,select_index,end_stage,mc_n,alpha,sigma_control,sigma_treatment,K,J){
  if(J == 3) {
    t <- c(1,2,3)/3
    low_bround <- c(0.000000, 1.038459, 2.290364)
    up_bround <- c(Inf, 3.301585, 2.290364)
    mu_control <- 0
    unit_n <- 30
  }else if(J == 4) {
    t <- c(1,2,3,4)/4
    low_bround <- c(0.0000000, 0.1267247, 1.4593383, 2.3056892)
    up_bround <- c(Inf, 4.013738, 2.826973, 2.251983)
    mu_control <- 0
    unit_n <- 25
  }else if(J == 5) {
    t <- c(1,2,3,4,5)/5
    low_bround <- c( 0.0000000, -0.6733603,  0.7312345,  1.6214603, 2.228605)
    up_bround <- c(Inf, 4.609266, 3.263919, 2.639797, 2.228605)
    mu_control <- 0
    unit_n <- 22
  }
  n <- t*length(t)*unit_n
  I <- function(j,index)  {1/(sigma_treatment[index]**2/n[j]+sigma_control**2/n[j])}
  mu_cov <- mu_cov_UMVCUE(control_eff=0,treat_eff=rep(0,K),
                          control_sig=sigma_control,treat_sig=sigma_treatment,
                          select_index,end_stage,K,J,n)
  mc <- rcmvnorm(mc_n,mean=mu_cov$mu,sigma = mu_cov$cov,
                 dependent=1:end_stage,given=c((end_stage+1):(end_stage+2)),
                 X.given=c(X[1],X[select_index+1])*sqrt(I(end_stage,select_index)))
  if(end_stage > 2){
    mc_stage2 <- (sqrt(I(2,select_index))*mc[,2]-sqrt(I(1,select_index))*mc[,1])/
      (I(2,select_index)-I(1,select_index))
  }else{
    mc_stage2 <- (I(2,select_index)*(X[select_index+1]-X[1])-sqrt(I(1,select_index))*mc[,1])/
      (I(2,select_index)-I(1,select_index))
  }
  low_bround_mc <- matrix(low_bround[1:(end_stage-1)],nrow = mc_n, ncol = end_stage-1,byrow = T)
  up_bround_mc <- matrix(up_bround[1:(end_stage-1)],nrow = mc_n, ncol = end_stage-1,byrow = T)
  low_bround_mc[,1] <- apply(rbind(matrix(low_bround[1],ncol = mc_n),
                                   (matrix(X[-c(1,select_index+1)],nrow = K-1, ncol = mc_n)-
                                      matrix(mc[,end_stage]/sqrt(I(1,select_index)),nrow = K-1, ncol = mc_n, byrow = T))*
                                     sqrt(I(1,c(1:K)[-select_index]))),2,max)
  mc_re <- apply(cbind(mc[,1:(end_stage-1)] >= low_bround_mc,
                       mc[,1:(end_stage-1)] <= up_bround_mc), 1, all)
  re <- rep(NA,4)
  if(sum(mc_re) > 1){
    mc_est <- mean(mc_stage2[mc_re])
    SE <- sqrt((sigma_treatment[select_index]**2+sigma_control**2)/(n[2]-n[1])-
                 var(mc_stage2[mc_re]))
    re <- c(mean(mc_re),mc_est,mc_est+qnorm(alpha/2)*SE,mc_est-qnorm(alpha/2)*SE)
  }
  re
}

##### RB #####
RB_mc <- function(s_all,select_index,end_stage,mc_num,alpha,sigma_control,sigma_treatment,K,J){
  if(J == 3) {
    t <- c(1,2,3)/3
    low_bround <- c(0.000000, 1.038459, 2.290364)
    up_bround <- c(Inf, 3.301585, 2.290364)
    mu_control <- 0
    unit_n <- 30
  }else if(J == 4) {
    t <- c(1,2,3,4)/4
    low_bround <- c(0.0000000, 0.1267247, 1.4593383, 2.3056892)
    up_bround <- c(Inf, 4.013738, 2.826973, 2.251983)
    mu_control <- 0
    unit_n <- 25
  }else if(J == 5) {
    t <- c(1,2,3,4,5)/5
    low_bround <- c( 0.0000000, -0.6733603,  0.7312345,  1.6214603, 2.228605)
    up_bround <- c(Inf, 4.609266, 3.263919, 2.639797, 2.228605)
    mu_control <- 0
    unit_n <- 22
  }
  n <- t*length(t)*unit_n
  theta <- rep(0,4)
  I <- function(i,j) n[j]/(sigma_control^2+sigma_treatment[i]^2)
  I_C <- function(i,ii,j) n[j]*sigma_control^2/(sigma_control^2+sigma_treatment[i]^2)/
    (sigma_control^2+sigma_treatment[ii]^2)
  
  mu <- c(theta*I(c(1:K),1),rep(theta[select_index]*I(select_index,1),end_stage-1))
  cov <- matrix(0,K+end_stage-1,K+end_stage-1,byrow = T)
  for(i in 1:K){
    for(ii in 1:K){
      cov[i,ii] <- I_C((1:K)[i],c(1:K)[ii],1)
    }
  }
  diag(cov[1:K,1:K]) <- I(1:K,1)
  if(end_stage>2){
    diag(cov[(K+1):(K+end_stage-1),(K+1):(K+end_stage-1)]) <- I(select_index,1)
  }else{
    cov[K+end_stage-1,K+end_stage-1] <- I(select_index,1)
  }
  cov[lower.tri(cov)] <- t(cov)[lower.tri(cov)]
  
  trans_T <- diag(1,K+end_stage-1,K+end_stage-1)
  trans_T[select_index,] <- 1
  trans_T[1:K,1:K] <- solve(cov[1:K,1:K])*I(c(1:K),1)
  mu_T <- trans_T %*% mu
  cov_T <- trans_T %*% cov %*% t(trans_T)
  
  T_star <- (I(c(1:K),1)*solve(cov[1:K,1:K])) %*% as.matrix(s_all[1:K])
  T_star[select_index] <- T_star[select_index]+(s_all[K+1]-s_all[select_index])
  mc_s <- rcmvnorm(mc_num,mu_T,cov_T,given.ind = 1:K,dependent.ind = (K+1):(K+end_stage-1),X.given = T_star)
  s1 <- matrix(NA, nrow = K, ncol = mc_num)
  for(i in 1:mc_num){
    T_star_s <- T_star
    T_star_s[select_index] <- T_star[select_index]-sum(mc_s[i,])
    s1[,i] <- solve((I(c(1:K),1)*solve(cov[1:K,1:K]))) %*% T_star_s
  }
  re <- rep(NA,4)
  z_stage1 <- apply(s1, 2, function(x) x/sqrt(I(1:K,1)))
  mc_s_all <- cbind(s1[select_index,],mc_s)
  mc_z <- t(apply(mc_s_all, 1, function(x) cumsum(x)/sqrt(I(select_index,1:end_stage))))
  bround_L <- matrix(low_bround[1:(end_stage-1)],mc_num,(end_stage-1),byrow = T)
  bround_U <- matrix(up_bround[1:(end_stage-1)],mc_num,(end_stage-1),byrow = T)
  
  index <- apply(cbind(apply(mc_z[,1:(end_stage-1)] >= bround_L,1,function(x)identical(x,y=rep(TRUE,end_stage-1))),
                       apply(mc_z[,1:(end_stage-1)] <= bround_U,1,function(x)identical(x,y=rep(TRUE,end_stage-1))),
                       apply(z_stage1,2,which.max)==select_index), 1, all)
  if(sum(index)>1){
    mc_est <- mean(mc_s_all[index,2]/(I(select_index,2)-I(select_index,1)))
    SE <- sqrt(1/(I(select_index,2)-I(select_index,1))-
                 var(mc_s_all[index,2]/(I(select_index,2)-I(select_index,1))))
    re <- c(sum(index)/mc_num,mc_est,mc_est+qnorm(alpha/2)*SE,mc_est-qnorm(alpha/2)*SE)
  }
  re
}

setting <- "peak"
sim <- 10
J=4
setwd("C:/Users/13379/Desktop/essay/submit-SIM/res/4")
time1 <- Sys.time()
cls <- makeSOCKcluster(2)
registerDoSNOW(cls)

pb <- txtProgressBar(max=sim, style=3, char = "*",)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

res  <- foreach(i=c(1:sim),
                  .options.snow=opts,
                  .combine = cbind,
                  .packages = c("MASS","data.table","compiler","condMVNorm","pbapply","tmvtnorm","BSDA","pracma")) %dopar% {
                  da_genc(i,setting,J)
                  }
end_stage <- res[1,]
select_index <- res[2,]
naive <- res[3:6,]
x_bar <- res[7:11,]
s_all <- res[12:16,]
# sigma_hat <- res[17:21,]
sigma_hat <- matrix(1,5,sim)
bias <- res[22:23,]

re_SI  <- foreach(i=which(end_stage>=2),
                  .options.snow=opts,
                  .combine = cbind,
                  .packages = c("MASS","data.table","compiler","condMVNorm","pbapply","tmvtnorm","BSDA","pracma")) %dopar% {
                    SI(naive[,i],select_index[i],sigma_hat[1,i],sigma_hat[-1,i],K=4,J)
                  }
write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],as.numeric(re_SI)),paste(setting,"SI.csv",sep = "-"))
re_MI  <- foreach(i=which(end_stage>=2),
                  .options.snow=opts,
                  .combine = cbind,
                  .packages = c("MASS","data.table","compiler","condMVNorm","pbapply","tmvtnorm","BSDA","pracma")) %dopar% {
                    MIc(naive[,i],select_index[i],max_iterations=20,tol=0.001,
                        sigma_hat[1,i],sigma_hat[-1,i],K=4,J)
                  }
write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],t(re_MI)),paste(setting,"MI.csv",sep = "-"))

re_MUE_MLE  <- foreach(i=which(end_stage>=2),
                       .options.snow=opts,
                       .combine = cbind,
                       .packages = c("MASS","data.table","compiler","condMVNorm","pbapply","tmvtnorm","BSDA","pracma")) %dopar% {
                         est_MUE(naive[,i],end_stage[i],select_index[i],sub="MLE",
                                 sigma_hat[1,i],sigma_hat[-1,i],K=4,J)
                       }
write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],t(re_MUE_MLE)),paste(setting,"MUE_MLE.csv",sep = "-"))
re_MUE_ZERO  <- foreach(i=which(end_stage>=2),
                        .options.snow=opts,
                        .combine = cbind,
                        .packages = c("MASS","data.table","compiler","condMVNorm","pbapply","tmvtnorm","BSDA","pracma")) %dopar% {
                          est_MUE(naive[,i],end_stage[i],select_index[i],sub="ZERO",
                                  sigma_hat[1,i],sigma_hat[-1,i],K=4,J)
                        }
write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],t(re_MUE_ZERO)),paste(setting,"MUE_ZERO.csv",sep = "-"))

re_UMVCUE  <- foreach(i=which(end_stage>=2),
                        .options.snow=opts,
                        .combine = cbind,
                        .packages = c("MASS","data.table","compiler","condMVNorm","pbapply","tmvtnorm","BSDA","pracma")) %dopar% {
                          UMVCUE_mc(x_bar[,i],select_index[i],end_stage[i],
                                    mc_n=10000,alpha=0.05,
                                    sigma_hat[1,i],sigma_hat[-1,i],K=4,J)
                        }
write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],t(re_UMVCUE)),paste(setting,"UMVCUE.csv",sep = "-"))

re_RB  <- foreach(i=which(end_stage>=2),
                      .options.snow=opts,
                      .combine = cbind,
                      .packages = c("MASS","data.table","compiler","condMVNorm","pbapply","tmvtnorm","BSDA","pracma")) %dopar% {
                        RB_mc(s_all[,i],select_index[i],end_stage[i],
                              mc_num=10000,alpha=0.05,
                              sigma_hat[1,i],sigma_hat[-1,i],K=4,J)
                      }
write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],t(re_RB)),paste(setting,"RB.csv",sep = "-"))
close(pb)
stopCluster(cls)
time2 <- Sys.time()
print(time2-time1)
