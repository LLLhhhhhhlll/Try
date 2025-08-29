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
mu_cov_f <- function(p_control,p_treatment,select_index,theta,K,J,n,I_type="p_bar"){
  if(is.null(theta)){theta <- log(p_treatment/(1-p_treatment))-log(p_control/(1-p_control))}
  if(I_type=="p_bar"){
    I <- function(i,j) {p_bar <- (p_treatment+p_control)/2; n[j]/2*p_bar[i]*(1-p_bar[i])}
  }else{
    I <- function(i,j) {n[j]/4*p_treatment[i]*(1-p_treatment[i])+n[j]/4*p_control*(1-p_control)}
  }
  
  cov <- cbind(matrix(n[1]/4*p_control*(1-p_control),K+J-1,K,byrow = T),
               matrix(0,K+J-1,J-1,byrow = T))
  diag(cov[1:K,1:K]) <- I(c(c(1:K)[-select_index],select_index),1)
  diag(cov[(K+1):(K+J-1),(K+1):(K+J-1)]) <- diff(I(select_index,1:J))
  cov[lower.tri(cov)] <- t(cov)[lower.tri(cov)]
  
  mu <- c(theta[-select_index]*I(c(1:K)[-select_index],1),
          theta[select_index]*c(I(select_index,1),diff(I(select_index,1:J))))
  
  trans_s <- diag(1,K+J-1,K+J-1)
  trans_s[K:(K+J-1),K:(K+J-1)][lower.tri(trans_s[K:(K+J-1),K:(K+J-1)])] <- 1
  
  theta_v <- matrix(0,K+J-1, K+J-1)
  diag(theta_v[1:(K-1),1:(K-1)]) <- 1/sqrt(I(c(1:K)[-select_index],1))
  diag(theta_v[K:(K+J-1),K:(K+J-1)]) <- 1/sqrt(I(select_index,1:J))
  
  mu_Z= theta_v %*% trans_s %*% mu 
  cov_Z = theta_v %*% (trans_s %*% cov %*% t(trans_s)) %*% t(theta_v)
  
  A <- diag(1,K+J-1)
  diag(A[1:(K-1),1:(K-1)]) <- -1
  A[1:(K-1),K] <- 1
  
  list(mu=A %*% mu_Z,cov=A %*% cov_Z %*% t(A))
}
n_find <- function(target_power,p_control,p_treatment,
                   t,aftype="sfLDOF",alpha=0.025,
                   select_bround=NULL,low_unit_n,up_unit_n){
  rr <- function(unit_n,t){
    library(gsDesign)
    library(mvtnorm)
    K <- length(p_treatment)
    J <- length(t)
    n <- t*length(t)*ceiling(unit_n)
    prob_H1 <- function(select_index,lower,upper,stage,n){
      H1 <- mu_cov_f(p_control,p_treatment,select_index,theta=NULL,K,J,n,I_type="non-p_bar")
      mu_H1 <- H1$mu; cov_H1 <- H1$cov
      pmvnorm(lower = lower,
              upper = upper,
              mean = mu_H1[1:(K+stage-1),],
              sigma=cov_H1[1:(K+stage-1),1:(K+stage-1)])[1]
    }
    H0 <- mu_cov_f(p_control,rep(p_control,K),
                   select_index=1,theta=NULL,K,J,n,I_type="non-p_bar")
    mu_H0 <- H0$mu; cov_H0 <- H0$cov
    if (all(p_treatment > p_control)){
      prob_select_wrong <- 0
    }else{
      select_wrong_index <- which(p_treatment <= p_control)  
      prob_select_wrong <- sum(sapply(X=select_wrong_index,prob_H1,
                                      lower = rep(0,K-1),upper = rep(Inf,K-1),stage=0,n))
    }
    if(prob_select_wrong >= 1-target_power){
      print("enter a bigger unit_n"); break
    }else{
      # the remaining beta for spending at following stage
      beta <- 1-target_power-prob_select_wrong
    }
    
    low_bround <- up_bround <- alpha_spend <- beta_spend <- numeric(length=J)
    select_index <- which(p_treatment > p_control)  
    # select_bround refers to the pre-defined boundary for stage 1
    if(is.null(select_bround)){
      alpha_spend <- eval(parse(text=paste(aftype,"(alpha = alpha, t)$spend")))
      beta_spend <- eval(parse(text=paste(aftype,"(alpha = beta, t)$spend")))
      up_bround[1] <- uniroot(function(x) K*pmvnorm(
        lower = c(rep(0,K-1),x),upper = rep(Inf,K),
        mean=mu_H0[1:K,],sigma=cov_H0[1:K,1:K])[1]-
          alpha_spend[1], c(-10, 10), tol = 0.0001)$root
      low_bround[1] <- uniroot(function(x)  sum(sapply(X=select_index,prob_H1,
                                                       lower = c(rep(0,K-1),-Inf),
                                                       upper = c(rep(Inf,K-1),x),stage=1,unit_n))-
                                 beta_spend[1], c(-10, 10), tol = 0.0001)$root
    }else{
      low_bround[1] <- select_bround[1]
      up_bround[1] <- select_bround[2]
      alpha_spend[1] <- K*pmvnorm(
        lower = c(rep(0,K-1),up_bround[1]),upper = rep(Inf,K),
        mean=mu_H0[1:K,],sigma=cov_H0[1:K,1:K])
      beta_spend[1] <- sum(sapply(X=select_index,prob_H1,
                                  lower = c(rep(0,K-1),-Inf),
                                  upper = c(rep(Inf,K-1),low_bround[1]),stage=1,n))
      if(beta_spend[1] >= beta){
        print("enter a bigger unit_n"); break
      }else{
        alpha_spend[-1] <- eval(parse(text=paste(aftype,"(alpha = alpha-alpha_spend[1], 
                                         (t[-1]-t[1])/(1-t[1]))$spend")))+alpha_spend[1]
        beta_spend[-1] <- eval(parse(text=paste(aftype,"(alpha = beta-beta_spend[1],
                                        (t[-1]-t[1])/(1-t[1]))$spend")))+beta_spend[1]
        
      }
    }
    for (i in 2:J) {
      if(any(low_bround > up_bround)){
        print("enter a lower unit_n"); break
      }else{     
        up_bround[i] <- uniroot(function(x) K*pmvnorm(
          lower = c(rep(0,K-1),low_bround[1:(i-1)],x),
          upper = c(rep(Inf,K-1),up_bround[1:(i-1)],Inf),
          mean = mu_H0[1:(K+i-1),],sigma=cov_H0[1:(K+i-1),1:(K+i-1)])[1]-
            (alpha_spend[i]-alpha_spend[i-1]), c(-10, 10), tol = 0.0001)$root
        low_bround[i] <- uniroot(function(x)  sum(sapply(X=select_index,prob_H1,
                                                         lower = c(rep(0,K-1),low_bround[1:(i-1)],-Inf),
                                                         upper =c(rep(Inf,K-1),up_bround[1:(i-1)],x),stage=i,n))-
                                   (beta_spend[i]-beta_spend[i-1]), c(-10, 10), tol = 0.0001)$root
      }
    }
    # calculate the power under a specified stage 
    power_cal <- function(stage,select_index){
      if(stage==1){
        sum(sapply(X=select_index,prob_H1,
                   lower = c(rep(0,K-1),up_bround[1]),
                   upper = c(rep(Inf,K-1),Inf),stage=1,n))
      }else{
        sum(sapply(X=select_index,prob_H1,
                   lower = c(rep(0,K-1),low_bround[1:(stage-1)],up_bround[stage]),
                   upper =c(rep(Inf,K-1),up_bround[1:(stage-1)],Inf),stage,n))
      }
    }
    # calculate the power for the entire trial 
    power <- cumsum(sapply(X=1:J,power_cal,select_index))
    list(samplesize_onegroup=t*length(t)*(unit_n),
         low_bround=low_bround[-J],
         up_bround=up_bround,power=power,
         c=up_bround[J]-low_bround[J])
  }
  # find the root to make U_J=L_J
  rr(unit_n=ceiling(uniroot(function(x,t) 
    rr(x,t)$c, c(low_unit_n,up_unit_n), tol = 0.0001,t=t)$root),t)
}
# set.seed(20250425)
# br = n_find(target_power = 0.8,p_control=0.4,p_treatment=c(0.6,0.4,0.4,0.4),
#             t=c(1,2,3,4)/4,aftype="sfLDOF",alpha=0.025,select_bround=c(0,Inf),
#             low_unit_n=30,up_unit_n=40)
# br
# 
# br = n_find(target_power = 0.8,p_control=0.4,p_treatment=c(0.6,0.4,0.4,0.4),
#             t=c(1,2,3)/3,aftype="sfLDOF",alpha=0.025,select_bround=c(0,Inf),
#             low_unit_n=40,up_unit_n=50)
# br
# br = n_find(target_power = 0.8,p_control=0.4,p_treatment=c(0.6,0.4,0.4,0.4),
#             t=c(1,2,3,4,5)/5,aftype="sfLDOF",alpha=0.025,select_bround=c(0,Inf),
#             low_unit_n=30,up_unit_n=40)
# br
##### sim #####
da_gen <- function(sim,setting,J){
  set.seed(20250425+sim)
  p_control <- 0.4
  if(setting=="NULL"){p_treatment <- rep(0.4,4)
  }else if(setting=="H1"){p_treatment <- c(0.6,0.6,0.6,0.6)
  }else if(setting=="peak"){p_treatment <- c(0.5,0.5,0.5,0.6)
  }else if(setting=="line"){p_treatment <- c(0.45,0.50,0.55,0.60)}
  theta <- log(p_treatment/(1-p_treatment))-log(p_control/(1-p_control))
  
  if(J==4){t <- c(1,2,3,4)/4
  low_bround <- c(0.00000000, 0.08743326, 1.41906914, 2.256008)
  up_bround <- c(Inf, 4.013693, 2.826977, 2.256008)
  unit_n <- 39    
  n <- t*length(t)*unit_n
  }else if (J==3){t <- c(1,2,3)/3
  low_bround <- c(0.000000, 1.015946, 2.291169)
  up_bround <- c(Inf, 3.301556, 2.291169)
  unit_n <- 47    
  n <- t*length(t)*unit_n
  }else if (J==5){t <- c(1,2,3,4,5)/5
  low_bround <- c(0.0000000, -0.7601533,  0.6465035,  1.5388462, 2.240586)
  up_bround <- c(Inf, 4.609263, 3.263943, 2.639832, 2.240586)
  unit_n <- 34    
  n <- t*length(t)*unit_n
  }
  
  K <- length(p_treatment)
  all_da <- apply(replicate(J,mapply(function(i) rbinom(1,unit_n, c(p_control,p_treatment)[i]),i=1:(K+1))),1,cumsum)
  S <- function(E,C,N) (N*E-N*C)/(2*N)
  I <- function(E,C,N) (E+C)*(2*N-E-C)*N^2/(2*N)^3
  Z <- function(E,C,N) S(E,C,N)/sqrt(I(E,C,N))
  select_index <- which.is.max(Z(all_da[1,-1],all_da[1,1],n[1]))
  inter_decide <- function(j,select_index,da){
    Z_temp <-  Z(all_da[j,select_index+1],all_da[j,1],n[j])
    if(Z_temp < low_bround[j]){"f"}else if(Z_temp > up_bround[j]){"e"} else{"c"}
  }
  inter_outcome <- as.character()
  for (i in 1:J) {
    inter_outcome[i] <- inter_decide(j=i,select_index,da=all_da)
  }
  end_stage <- first(which(inter_outcome %in% c("e","f")))
  
  stage2_est <- S(all_da[2,select_index+1]-all_da[1,select_index+1],all_da[2,1]-all_da[1,1],unit_n)/
    I(all_da[2,select_index+1]-all_da[1,select_index+1],all_da[2,1]-all_da[1,1],unit_n)
  
  naive_P <- all_da[1,]/unit_n
  naive_P[1] <- all_da[end_stage,1]/n[end_stage]
  naive_P[select_index+1] <- all_da[end_stage,select_index+1]/n[end_stage]
  
  naive_est <- S(all_da[1,-1],all_da[1,1],n[1])/I(all_da[1,-1],all_da[1,1],n[1])
  naive_est[select_index] <- S(all_da[end_stage,select_index+1],all_da[end_stage,1],n[end_stage])/
    I(all_da[end_stage,select_index+1],all_da[end_stage,1],n[end_stage])
  
  x_bar <- all_da[1,]
  x_bar[1] <- all_da[end_stage,1]
  x_bar[select_index+1] <- all_da[end_stage,select_index+1]
  
  s_all <- c(S(all_da[1,-1],all_da[1,1],n[1]),
             S(all_da[end_stage,select_index+1],all_da[end_stage,1],n[end_stage]))
  
  c(end_stage,select_index,naive_P,naive_est,x_bar,s_all,theta[select_index],
    stage2_est,naive_est[select_index])
}
da_genc <- cmpfun(da_gen)

##### MI & SI #####
bias_theta  <- function(select_index,theta,K,J,p_control,p_treatment){
  if(J==4){t <- c(1,2,3,4)/4
  low_bround <- c(0.00000000, 0.08743326, 1.41906914, 2.256008)
  up_bround <- c(Inf, 4.013693, 2.826977, 2.256008)
  unit_n <- 39    
  n <- t*length(t)*unit_n
  }else if (J==3){t <- c(1,2,3)/3
  low_bround <- c(0.000000, 1.015946, 2.291169)
  up_bround <- c(Inf, 3.301556, 2.291169)
  unit_n <- 47    
  n <- t*length(t)*unit_n
  }else if (J==5){t <- c(1,2,3,4,5)/5
  low_bround <- c(0.0000000, -0.7601533,  0.6465035,  1.5388462, 2.240586)
  up_bround <- c(Inf, 4.609263, 3.263943, 2.639832, 2.240586)
  unit_n <- 34    
  n <- t*length(t)*unit_n
  }
  
  I <- function(p_E,p_C,N) {p_bar <- (p_E+p_C)/2; N/2*p_bar*(1-p_bar)}
  mu_cov <- mu_cov_f(p_control,p_treatment,select_index,theta,K,J,n)
  mu <- mu_cov$mu; cov <- mu_cov$cov
  prob_Q <- pmvnorm(lower = c(rep(0,K-1),low_bround[1]),upper = c(rep(Inf,K-1),up_bround[1]),
                    mean = mu[1:K,],sigma = cov[1:K,1:K])[1]
  
  bias <- rep(NA,K)
  if(!is.na(prob_Q)){
    mtm <- mtmvnorm(mean = mu[1:K,],sigma = cov[1:K,1:K],
                    lower = c(rep(0,K-1),low_bround[1]),upper = c(rep(Inf,K-1),up_bround[1]),
                    doComputeVariance=FALSE)$tmean
    bias[-select_index] <- (mtm[K]-mtm[-K])/sqrt(I(p_treatment[-select_index],p_control,n[1]))-theta[-select_index]
    E_select <- function(stage){
      mu_t <- mu[1:(K+stage-1),]
      cov_t <- cov[1:(K+stage-1),1:(K+stage-1)]
      lower1 <- c(rep(0,K-1),low_bround[1:(stage-1)],-Inf); upper1 <- c(rep(Inf,K-1),up_bround[1:(stage-1)],low_bround[stage])
      lower2 <- c(rep(0,K-1),low_bround[1:(stage-1)],up_bround[stage]);upper2 <- c(rep(Inf,K-1),up_bround[1:(stage-1)],Inf)
      (pmvnorm(lower = lower1,upper = upper1,mean = mu_t,sigma = cov_t)[1]*
          mtmvnorm(mean = mu_t,sigma = cov_t,lower = lower1,upper = upper1,doComputeVariance=FALSE)$tmean[K+stage-1]+
          pmvnorm(lower = lower2,upper = upper2,mean = mu_t,sigma = cov_t)[1]*
          mtmvnorm(mean = mu_t,sigma = cov_t,lower = lower2,upper = upper2,doComputeVariance=FALSE)$tmean[K+stage-1])/
        sqrt(I(p_treatment[select_index],p_control,n[stage]))
    }
    bias[select_index] <- sum(mapply(E_select, stage=2:J))/prob_Q-theta[select_index]
  }
  bias
}
MI <- function(initial,select_index,max_iterations=20,tol=0.001,K,J,p_control,p_treatment){
  conver <- 0
  theta_hat <- initial
  for (i in 1:max_iterations) {
    bias_temp <- bias_theta(select_index,theta=theta_hat,K,J,p_control,p_treatment)
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
SI <- function(initial,select_index,K,J,p_control,p_treatment){
  (initial-bias_theta(select_index,initial,K,J,p_control,p_treatment))[select_index]
}  
MIc <- cmpfun(MI)
##### MUE #####
est_MUE <- function(naive,end_stage,select_index,sub,K,J,p_control,p_treatment){
  if(J==4){t <- c(1,2,3,4)/4
  low_bround <- c(0.00000000, 0.08743326, 1.41906914, 2.256008)
  up_bround <- c(Inf, 4.013693, 2.826977, 2.256008)
  unit_n <- 39    
  n <- t*length(t)*unit_n
  }else if (J==3){t <- c(1,2,3)/3
  low_bround <- c(0.000000, 1.015946, 2.291169)
  up_bround <- c(Inf, 3.301556, 2.291169)
  unit_n <- 47    
  n <- t*length(t)*unit_n
  }else if (J==5){t <- c(1,2,3,4,5)/5
  low_bround <- c(0.0000000, -0.7601533,  0.6465035,  1.5388462, 2.240586)
  up_bround <- c(Inf, 4.609263, 3.263943, 2.639832, 2.240586)
  unit_n <- 34    
  n <- t*length(t)*unit_n
  }
  
  I <- function(p_E,p_C,N) {p_bar <- (p_E+p_C)/2; N/2*p_bar*(1-p_bar)}
  Z_end <- sqrt(n[end_stage]/2)*(p_treatment[select_index]-p_control)/
    sqrt((p_treatment[select_index]+p_control)/2*(1-(p_treatment[select_index]+p_control)/2))
  
  cond_pvalue <- function(theta,select_index,end_stage,Z_end,p_control,p_treatment){
    mu_cov <- mu_cov_f(p_control,p_treatment,select_index,theta,K,J,n)
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
      theta[-select_index] <- naive[-select_index]
    }else{
      theta[-select_index] <- rep(0,K-1)
    }
    cond_pvaluec(theta,select_index,end_stage,Z_end,p_control,p_treatment)
  }
  a <- optimise(function(x) (P_value_sub(select_index,sub,x)-1/2)^2,c(-1,1),tol=0.001)
  c(a$minimum,a$objective)
}
##### UMVCUE #####
UMVCUE_mc <- function(x_bar,select_index,end_stage,mc_num,alpha,K,J){
  if(J==4){t <- c(1,2,3,4)/4
  low_bround <- c(0.00000000, 0.08743326, 1.41906914, 2.256008)
  up_bround <- c(Inf, 4.013693, 2.826977, 2.256008)
  unit_n <- 39    
  n <- t*length(t)*unit_n
  }else if (J==3){t <- c(1,2,3)/3
  low_bround <- c(0.000000, 1.015946, 2.291169)
  up_bround <- c(Inf, 3.301556, 2.291169)
  unit_n <- 47    
  n <- t*length(t)*unit_n
  }else if (J==5){t <- c(1,2,3,4,5)/5
  low_bround <- c(0.0000000, -0.7601533,  0.6465035,  1.5388462, 2.240586)
  up_bround <- c(Inf, 4.609263, 3.263943, 2.639832, 2.240586)
  unit_n <- 34    
  n <- t*length(t)*unit_n
  }
  
  S <- function(E,C,N) (N*E-N*C)/(2*N)
  I <- function(E,C,N) (E+C)*(2*N-E-C)*N^2/(2*N)^3
  Z <- function(E,C,N) S(E,C,N)/sqrt(I(E,C,N))
  mc_C <- matrix(x_bar[1],nrow=mc_num,ncol=end_stage-1)
  mc_E <- matrix(x_bar[select_index+1],nrow=mc_num,ncol=end_stage-1)
  mc_C[,1] <- rhyper(mc_num, x_bar[1], n[end_stage]-x_bar[1], n[1])
  mc_E[,1] <- rhyper(mc_num, x_bar[select_index+1], n[end_stage]-x_bar[select_index+1], n[1])
  if(end_stage > 2){
    for(j in 2:(end_stage-1)){
      mc_C[,j] <- mapply(function(x) rhyper(1, x_bar[1]-x, n[end_stage-j+1]-
                                              (x_bar[1]-x), unit_n),x=mc_C[,j-1])+mc_C[,j-1]
      mc_E[,j] <-  mapply(function(x) rhyper(1, x_bar[select_index+1]-x, n[end_stage-j+1]-
                                               (x_bar[select_index+1]-x), unit_n),x=mc_E[,j-1])+mc_E[,j-1]
    }
  }
  Z_stage1 <- matrix(nrow = K,ncol = mc_num)
  Z_stage1[select_index,] <- Z(mc_E[,1],mc_C[,1],n[1])
  Z_stage1[-select_index,] <- mapply(function(i) Z(x_bar[-c(1,select_index+1)],mc_C[i,1],n[1]), i=1:mc_num)
  re <- rep(NA,2)
  
  Z_matirx <- t(mapply(function(i) Z(mc_E[i,],mc_C[i,],n[1:(end_stage-1)]), i=1:mc_num))
  bround_L <- matrix(low_bround[1:(end_stage-1)],mc_num,(end_stage-1),byrow = T)
  bround_U <- matrix(up_bround[1:(end_stage-1)],mc_num,(end_stage-1),byrow = T)
  if(end_stage > 2){
    index <- apply(cbind(apply(Z_matirx >= bround_L,1,function(x)identical(x,y=rep(TRUE,end_stage-1))),
                         apply(Z_matirx <= bround_U,1,function(x)identical(x,y=rep(TRUE,end_stage-1))),
                         apply(Z_stage1,2,which.is.max)==select_index), 1, all)
  }else{
    index <- apply(cbind(t(Z_matirx >= as.numeric(bround_L)),
                         t(Z_matirx <= as.numeric(bround_U)),
                         apply(Z_stage1,2,which.is.max)==select_index), 1, all)
  } 
  if(sum(index) > 1){
    if(end_stage > 2){
      mc_est <- mean((S(mc_E[,2]-mc_E[,1],mc_C[,2]-mc_C[,1],unit_n)/I(mc_E[,2]-mc_E[,1],mc_C[,2]-mc_C[,1],unit_n))[index])
      SE <- mean(sqrt(1/I(mc_E[,2]-mc_E[,1],mc_C[,2]-mc_C[,1],unit_n)-
                        var((S(mc_E[,2]-mc_E[,1],mc_C[,2]-mc_C[,1],unit_n)/I(mc_E[,2]-mc_E[,1],mc_C[,2]-mc_C[,1],unit_n))[index])))
      re <- c(mean(index),mc_est,mc_est+qnorm(alpha/2)*SE,mc_est-qnorm(alpha/2)*SE)
    }else{
      mc_est <- mean((S(x_bar[select_index+1]-mc_E[,1],x_bar[1]-mc_C[,1],unit_n)/I(x_bar[select_index+1]-mc_E[,1],x_bar[1]-mc_C[,1],unit_n))[index])
      SE <- mean(sqrt(1/I(x_bar[select_index+1]-mc_E[,1],x_bar[1]-mc_C[,1],unit_n)-
                        var((S(x_bar[select_index+1]-mc_E[,1],x_bar[1]-mc_C[,1],unit_n)/
                               I(x_bar[select_index+1]-mc_E[,1],x_bar[1]-mc_C[,1],unit_n))[index])),na.rm=T)
      re <- c(mean(index),mc_est,mc_est+qnorm(alpha/2)*SE,mc_est-qnorm(alpha/2)*SE)
    }
  }
  re
}

##### RB #####
RB_mc <- function(s_all,p_control,p_treatment,select_index,end_stage,mc_num,alpha,K,J){
  if(J==4){t <- c(1,2,3,4)/4
  low_bround <- c(0.00000000, 0.08743326, 1.41906914, 2.256008)
  up_bround <- c(Inf, 4.013693, 2.826977, 2.256008)
  unit_n <- 39    
  n <- t*length(t)*unit_n
  }else if (J==3){t <- c(1,2,3)/3
  low_bround <- c(0.000000, 1.015946, 2.291169)
  up_bround <- c(Inf, 3.301556, 2.291169)
  unit_n <- 47    
  n <- t*length(t)*unit_n
  }else if (J==5){t <- c(1,2,3,4,5)/5
  low_bround <- c(0.0000000, -0.7601533,  0.6465035,  1.5388462, 2.240586)
  up_bround <- c(Inf, 4.609263, 3.263943, 2.639832, 2.240586)
  unit_n <- 34    
  n <- t*length(t)*unit_n
  }
  
  theta <- rep(0,K)
  # I <- function(i,j) {n[j]/2*((all_da[j,i+1]+all_da[j,1])/(2*n[j]))*(1-(all_da[j,i+1]+all_da[j,1])/(2*n[j]))}
  I <- function(i,j) {p_bar <- (p_treatment+p_control)/2; n[j]/2*p_bar[i]*(1-p_bar[i])}
  cov <- cbind(matrix(n[1]/4*p_control*(1-p_control),K+end_stage-1,K,byrow = T),
               matrix(0,K+end_stage-1,end_stage-1,byrow = T))
  diag(cov[1:K,1:K]) <- I(1:K,1)
  if(end_stage>2){
    diag(cov[(K+1):(K+end_stage-1),(K+1):(K+end_stage-1)]) <- diff(I(select_index,1:end_stage))
  }else{
    cov[K+end_stage-1,K+end_stage-1] <- diff(I(select_index,1:end_stage))
  }
  cov[lower.tri(cov)] <- t(cov)[lower.tri(cov)]
  mu <- c(theta*I(1:K,1),theta[select_index]*diff(I(select_index,1:end_stage)))
  
  trans_T <- diag(1,K+end_stage-1,K+end_stage-1)
  trans_T[select_index,] <- 1
  trans_T[1:K,1:K] <- solve(cov[1:K,1:K])*I(1:K,1)
  mu_T <- trans_T %*% mu
  cov_T <- trans_T %*% cov %*% t(trans_T)
  
  T_star <- (I(1:K,1)*solve(cov[1:K,1:K])) %*% as.matrix(s_all[1:K])
  T_star[select_index] <- T_star[select_index]+(s_all[K+1]-s_all[select_index])
  mc_s <- rcmvnorm(mc_num,mu_T,cov_T,given.ind = 1:K,dependent.ind = (K+1):(K+end_stage-1),X.given = T_star)
  s1 <- matrix(NA, nrow = K, ncol = mc_num)
  for(i in 1:mc_num){
    T_star_s <- T_star
    T_star_s[select_index] <- T_star[select_index]-sum(mc_s[i,])
    s1[,i] <- solve((I(1:K,1)*solve(cov[1:K,1:K]))) %*% T_star_s
  }
  re <- rep(NA,4)
  z_stage1 <- apply(s1, 2, function(x) x/sqrt(I(1:K,1)))
  mc_s_all <- cbind(s1[select_index,],mc_s)
  mc_z <- t(apply(mc_s_all, 1, function(x) cumsum(x)/sqrt(I(select_index,1:end_stage))))
  bround_L <- matrix(low_bround[1:(end_stage-1)],mc_num,(end_stage-1),byrow = T)
  bround_U <- matrix(up_bround[1:(end_stage-1)],mc_num,(end_stage-1),byrow = T)
  
  index <- apply(cbind(apply(mc_z[,1:(end_stage-1)] >= bround_L,1,function(x)identical(x,y=rep(TRUE,end_stage-1))),
                       apply(mc_z[,1:(end_stage-1)] <= bround_U,1,function(x)identical(x,y=rep(TRUE,end_stage-1))),
                       apply(z_stage1,2,which.is.max)==select_index), 1, all)
  if(sum(index)>1){
    mc_est <- mean(mc_s_all[index,2]/(I(select_index,2)-I(select_index,1)))
    SE <- sqrt(1/(I(select_index,2)-I(select_index,1))-
                 var(mc_s_all[index,2]/(I(select_index,2)-I(select_index,1))))
    re <- c(sum(index)/mc_num,mc_est,mc_est+qnorm(alpha/2)*SE,mc_est-qnorm(alpha/2)*SE)
  }
  re
}

setwd("C:/Users/13379/Desktop/essay/submit-SIM/binary-logHR/5")
# setting <- "peak"
for(setting in c("NULL","H1","peak","line")){
  J <- 5
  sim <- 100000
  
  time1 <- Sys.time()
  cls <- makeSOCKcluster(18)
  registerDoSNOW(cls)
  
  pb <- txtProgressBar(max=sim, style=3, char = "*",)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  res  <- foreach(i=c(1:sim),
                  .options.snow=opts,
                  .combine = cbind,
                  .packages = c("MASS","data.table","compiler","condMVNorm",
                                "pbapply","tmvtnorm","BSDA","pracma","nnet")) %dopar% {
                                  da_genc(i,setting,J)
                                }
  end_stage <- res[1,]
  select_index <- res[2,]
  naive_P <- res[3:7,]
  naive <- res[8:11,]
  x_bar <- res[12:16,]
  s_all <- res[17:21,]
  bias <- res[22:24,]
  write.csv(cbind(end_stage,select_index,t(bias))[which(end_stage>=2),],paste(setting,"stage2-naive.csv",sep = "-"))
  
  re_SI <- foreach(i=which(end_stage>=2),
                   .options.snow=opts,
                   .combine = cbind,
                   .packages = c("MASS","data.table","compiler","condMVNorm",
                                 "pbapply","tmvtnorm","BSDA","pracma","nnet")) %dopar% {
                                   SI(naive[,i],select_index[i],K=4,J,naive_P[1,i],naive_P[-1,i])
                                 }
  write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],as.numeric(re_SI)),paste(setting,"SI.csv",sep = "-"))
  
  re_MI <- foreach(i=which(end_stage>=2),
                   .options.snow=opts,
                   .combine = cbind,
                   .packages = c("MASS","data.table","compiler","condMVNorm",
                                 "pbapply","tmvtnorm","BSDA","pracma","nnet")) %dopar% {
                                   MIc(naive[,i],select_index[i],max_iterations=20,tol=0.001,K=4,J,
                                       naive_P[1,i],naive_P[-1,i])
                                 }
  write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],t(re_MI)),paste(setting,"MI.csv",sep = "-"))
  
  re_MUE_MLE <- foreach(i=which(end_stage>=2),
                        .options.snow=opts,
                        .combine = cbind,
                        .packages = c("MASS","data.table","compiler","condMVNorm",
                                      "pbapply","tmvtnorm","BSDA","pracma","nnet")) %dopar% {
                                        est_MUE(naive[,i],end_stage[i],select_index[i],sub="MLE",
                                                K=4,J,naive_P[1,i],naive_P[-1,i])
                                      }
  write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],t(re_MUE_MLE)),paste(setting,"MUE_MLE.csv",sep = "-"))
  
  re_MUE_ZERO <- foreach(i=which(end_stage>=2),
                         .options.snow=opts,
                         .combine = cbind,
                         .packages = c("MASS","data.table","compiler","condMVNorm",
                                       "pbapply","tmvtnorm","BSDA","pracma","nnet")) %dopar% {
                                         est_MUE(naive[,i],end_stage[i],select_index[i],sub="ZERO",
                                                 K=4,J,naive_P[1,i],naive_P[-1,i])
                                       }
  write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],t(re_MUE_ZERO)),paste(setting,"MUE_ZERO.csv",sep = "-"))
  
  re_UMVCUE <- foreach(i=which(end_stage>=2),
                       .options.snow=opts,
                       .combine = cbind,
                       .packages = c("MASS","data.table","compiler","condMVNorm",
                                     "pbapply","tmvtnorm","BSDA","pracma","nnet")) %dopar% {
                                       UMVCUE_mc(x_bar[,i],select_index[i],end_stage[i],
                                                 mc_num=10000,alpha=0.05,K=4,J)
                                     }
  write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],t(re_UMVCUE)),paste(setting,"UMVCUE.csv",sep = "-"))
  
  re_RB <- foreach(i=which(end_stage>=2),
                   .options.snow=opts,
                   .combine = cbind,
                   .packages = c("MASS","data.table","compiler","condMVNorm",
                                 "pbapply","tmvtnorm","BSDA","pracma","nnet")) %dopar% {
                                   RB_mc(s_all[,i],naive_P[1,i],naive_P[-1,i],select_index[i],end_stage[i],
                                         mc_num=10000,alpha=0.05,K=4,J)
                                 }
  write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],t(re_RB)),paste(setting,"RB.csv",sep = "-"))
  
  close(pb)
  stopCluster(cls)
  time2 <- Sys.time()
  print(time2-time1)
}

