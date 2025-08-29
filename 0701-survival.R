rm(list = ls())
mu_cov_f <- function(lambda_C,lambda_E,shape,select_index,theta,D_n,K,J){
  if(is.null(theta)){theta <- log(lambda_E/lambda_C)}
  I <- function(i,j){1/4*(D_n[j,1]+D_n[j,i+1])}
  
  cov <- cbind(matrix(1/4*D_n[1,1],K+J-1,K,byrow = T),
               matrix(0,K+J-1,J-1,byrow = T))
  diag(cov[1:K,1:K]) <- mapply(I,i=c(c(1:K)[-select_index],select_index),j=1)
  diag(cov[(K+1):(K+J-1),(K+1):(K+J-1)]) <- diff(mapply(I,i=select_index,j=1:J))
  cov[lower.tri(cov)] <- t(cov)[lower.tri(cov)]
  
  mu <- c(c(theta[-select_index],theta[select_index])*mapply(I,i=c(c(1:K)[-select_index],select_index),j=1),
          theta[select_index]*diff(mapply(I,i=select_index,j=1:J)))
  
  trans_s <- diag(1,K+J-1,K+J-1)
  trans_s[K:(K+J-1),K:(K+J-1)][lower.tri(trans_s[K:(K+J-1),K:(K+J-1)])] <- 1
  
  theta_v <- matrix(0,K+J-1, K+J-1)
  diag(theta_v[1:(K-1),1:(K-1)]) <- 1/sqrt(mapply(I,i=c(c(1:K)[-select_index]),j=1))
  diag(theta_v[K:(K+J-1),K:(K+J-1)]) <- 1/sqrt(mapply(I,i=select_index,j=1:J))
  
  mu_Z= theta_v %*% trans_s %*% mu 
  cov_Z = theta_v %*% (trans_s %*% cov %*% t(trans_s)) %*% t(theta_v)
  
  A <- diag(1,K+J-1)
  diag(A[1:(K-1),1:(K-1)]) <- -1
  A[1:(K-1),K] <- 1
  
  list(mu=A %*% mu_Z,cov=A %*% cov_Z %*% t(A))
}
bround_f <- function(lambda_C,shape,unit_D,aftype,alpha,select_bround,K,J,t){
  H0 <- mu_cov_f(lambda_C,rep(lambda_C,K),shape,select_index=1,theta=NULL,
                 matrix(t*length(t)*unit_D,nrow = J, ncol = K+1, byrow = F),K,J)
  mu_H0 <- H0$mu; cov_H0 <- H0$cov
  futility_bround <- efficacy_bround <- alpha_spend <- numeric(length=J)
  if(is.null(select_bround)){
    alpha_spend <- eval(parse(text=paste(aftype,"(alpha = alpha, t)$spend")))
    efficacy_bround[1] <- uniroot(function(x) K*pmvnorm(
      lower = rep(-Inf,K),upper = c(rep(0,K-1),x),
      mean = mu_H0[1:K,],sigma = cov_H0[1:K,1:K])[1]-
        alpha_spend[1], c(-10, 10), tol = 0.0001)$root
  }else{
    efficacy_bround[1] <- select_bround[1]
    futility_bround[1] <- select_bround[2]
    alpha_spend[1] <- K*pmvnorm(
      lower = rep(-Inf,K),upper = c(rep(0,K-1),efficacy_bround[1]),
      mean=mu_H0[1:K,],sigma=cov_H0[1:K,1:K])
    alpha_spend[-1] <- eval(parse(text=paste(aftype,"(alpha = alpha-alpha_spend[1], 
                                         (t[-1]-t[1])/(1-t[1]))$spend")))+alpha_spend[1]
  }
  for (i in 2:J) {
    efficacy_bround[i] <- uniroot(function(x) K*pmvnorm(
      lower = c(rep(-Inf,K-1),efficacy_bround[1:(i-1)],-Inf),
      upper = c(rep(0,K-1),rep(Inf,i-1),x),
      mean = mu_H0[1:(K+i-1),],sigma=cov_H0[1:(K+i-1),1:(K+i-1)])[1]-
        (alpha_spend[i]-alpha_spend[i-1]), c(-10, 10), tol = 0.0001)$root
  }
  futility_bround[-1] <- -efficacy_bround[-1] 
  futility_bround[J] <- efficacy_bround[J]
  list(futility_bround=futility_bround,efficacy_bround=efficacy_bround)
}
n_find <-function(target_power,lambda_C,lambda_E,shape,aftype="sfLDOF",alpha=0.025,
                  select_bround=NULL,low_unit_D,up_unit_D,C,a_t,t,K,J){
  power_f <- function(unit_D){
    library(gsDesign)
    bround <- bround_f(lambda_C,shape,unit_D,aftype,alpha,select_bround,K,J,t)
    efficacy_bround <- bround$efficacy_bround
    futility_bround <- bround$futility_bround 
    
    inter_t <- length(J)
    D_n <- matrix(nrow = J, ncol = K+1)
    S_f <- function(t,lambda,shape) {exp(-lambda*t^shape)}
    D_f <- function(f_t,a_t,lambda,shape,C) {integrate(function(u) 1-S_f(f_t-u,lambda,shape), 0, min(a_t,f_t))$value*C}
    inter_t[1] <- uniroot(function(f_t) sum(mapply(D_f,f_t,a_t,lambda=c(lambda_C,lambda_E),shape,C))-unit_D*(K+1),
                          lower = 0.001, upper = 100)$root
    D_n[1,] <- mapply(D_f,inter_t[1],a_t,lambda=c(lambda_C,lambda_E),shape,C)
    select_index <- which.min(lambda_E)
    for(j in 2:J){
      inter_t[j] <- uniroot(function(f_t) sum(mapply(D_f,f_t,a_t,lambda=c(lambda_C,lambda_E[select_index]),shape,C))-
                              (sum(D_n[1,c(1,select_index+1)])+(t[j]-t[1])*length(t)*unit_D*2),
                            lower = inter_t[j-1], upper = 100)$root
      D_n[j,c(1,select_index+1)] <- mapply(D_f,inter_t[j],a_t,lambda=c(lambda_C,lambda_E[select_index]),shape,C)
    }
    
    prob_H1 <- function(select_index,lower,upper,stage,N){
      H1 <- mu_cov_f(lambda_C,lambda_E,shape,select_index,theta=NULL,D_n,K,J)
      mu_H1 <- H1$mu; cov_H1 <- H1$cov
      pmvnorm(lower = lower,
              upper = upper,
              mean = mu_H1[1:(K+stage-1),],
              sigma = cov_H1[1:(K+stage-1),1:(K+stage-1)])[1]
    }
    
    power_cal <- function(stage,select_index){
      if(stage==1){
        sum(sapply(X=select_index,prob_H1,
                   lower = c(rep(-Inf,K-1),-Inf),
                   upper = c(rep(0,K-1),efficacy_bround[1]),stage=1,N))
      }else{
        sum(sapply(X=select_index,prob_H1,
                   lower = c(rep(-Inf,K-1),efficacy_bround[1:(stage-1)],-Inf),
                   upper =c(rep(0,K-1),futility_bround[1:(stage-1)],efficacy_bround[stage]),stage,N))
      }
    }
    power <- sum(sapply(X=1:J,power_cal,select_index))
    list(power=power, futility_bround=futility_bround,efficacy_bround=efficacy_bround, 
         D_n = unit_D*(K+1)+(t-t[1])*length(t)*unit_D*2,inter_t=inter_t)
  }
  power_f(ceiling(uniroot(function(unit_D) power_f(unit_D)$power-target_power, c(low_unit_D,up_unit_D), tol = 0.0001)$root))
}
# br <- n_find(target_power = 0.8,lambda_C=log(2),lambda_E= exp(c(log(0.8),0,0,0))*log(2),
#        shape=0.5,aftype="sfLDOF",alpha=0.025,select_bround=c(-Inf,0),
#        low_unit_D=100,up_unit_D=200,C=1000,a_t=1,t=c(1,2,3,4)/4,K=4,J=4)
# br <- n_find(target_power = 0.8,lambda_C=log(2),lambda_E= exp(c(log(0.8),0,0,0))*log(2),
#        shape=0.5,aftype="sfLDOF",alpha=0.025,select_bround=c(-Inf,0),
#        low_unit_D=100,up_unit_D=200,C=1000,a_t=1,t=c(1,2,3)/3,K=4,J=3)
# br <- n_find(target_power = 0.8,lambda_C=log(2),lambda_E= exp(c(log(0.8),0,0,0))*log(2),
#              shape=0.5,aftype="sfLDOF",alpha=0.025,select_bround=c(-Inf,0),
#              low_unit_D=100,up_unit_D=200,C=1000,a_t=1,t=c(1,2,3,4,5)/5,K=4,J=5)

# unit_D <- 100
# t <- c(1:4)/4
# J <- 4
# K <- 4
# select_index <- 1
# theta <- rep(0,K)
# select_bround <- c(-Inf,0)
# target_power = 0.8;
# C <- 1000; a_t <- 1
# lambda_C=log(2);
# lambda_E= exp(c(log(0.8),0,0,0))*log(2);
# shape=0.5;a_t=1;f_t=c(1,2,3,4);aftype="sfLDOF";alpha=0.025;
# low_unit_D <- 100
# up_unit_D <- 200

##### sim #####
rm(list = ls())
# sim <- 1
# setting <- "peak"
# J <- 4
# shape <- 0.5
# a_t <- 1
# N <- 1000
da_gen <- function(sim,setting,J){
  set.seed(20250425+sim)
  lambda_C <- log(2)
  if(setting=="NULL"){lambda_E <- rep(lambda_C,4)
  }else if(setting=="H1"){lambda_E <- exp(rep(log(0.8),4))*lambda_C
  }else if(setting=="peak"){lambda_E <- exp(c(log(0.8)/2,log(0.8)/2,log(0.8)/2,log(0.8)))*lambda_C
  }else if(setting=="line"){lambda_E <- exp(c(log(0.8)/4,log(0.8)/2,log(0.8)*3/4,log(0.8)))*lambda_C}
  if(J==4){
    low_bround <- c(0.000000,  4.013691,  2.828157, -2.285393)
    up_bround <- c(-Inf, -4.013691, -2.828157, -2.285393)
    D_n <- c(625,  875, 1125, 1375) 
  }else if(J==3){
    low_bround <- c(0.000000,  3.301594, -2.302597)
    up_bround <- c(-Inf, -3.301594, -2.302597)
    D_n <- c(765, 1071, 1377) 
  }else if(J==5){
    low_bround <- c(0.000000,  4.609252,  3.264385,  2.644414, -2.278346)
    up_bround <- c(-Inf, -4.609252, -3.264385, -2.644414, -2.278346)
    D_n <- c(555,  777,  999, 1221, 1443) 
  }
  K <- length(lambda_E)
  theta <- log(lambda_E/lambda_C)
  stage_da_f <- function(f_t){
    status <- da$D_time <= f_t
    time <- f_t - da$a_t
    time[status] <- da$y[status]
    data.frame(treatment = da$treatment,status,time)
  }
  Z_sel_f <- function(i,j){
    stage_da <- stage_da_f(f_t=inter_t[j])
    summary(coxph(Surv(time, status) ~ factor(treatment), data= stage_da[stage_da$treatment %in% c(0,i),]))$coef
  }
  I <- function(i,j){ 1/(Z_sel_f(i,j)[3])^2 }
  S <- function(i,j){ Z_sel_f(i,j)[1]*I(i,j) }
  da <- data.frame(cbind(reshape2::melt(mapply(function(scale) 
    rweibullPH(N,shape,scale), scale=c(lambda_C,lambda_E)),value.name="y"),
    a_t=runif((K+1)*N,0,a_t)))
  da$treatment <- rep(c(0:K),each=N)
  da$D_time <- da$y+da$a
  inter_t <- length(J)
  inter_t[1] <- sort(da$D_time)[D_n[1]]
  stage1_cox <- coxph(Surv(time, status) ~ factor(treatment), data= stage_da_f(inter_t[1]))
  select_index <- which.min(stage1_cox$coefficients)
  further_time <- da[(da$treatment %in% c(0,select_index)) & da$D_time > inter_t[1],]$D_time
  inter_t[2:J] <- sort(further_time)[D_n[2:J]-D_n[1]]
  
  inter_decide <- function(j,select_index){
    Z_temp <-  Z_sel_f(select_index,j)[4]
    if(Z_temp > low_bround[j]){"f"}else if(Z_temp < up_bround[j]){"e"} else{"c"}
  }
  inter_outcome <- as.character()
  for (j in 1:J) {
    inter_outcome[j] <- inter_decide(j,select_index)
  }
  end_stage <- first(which(inter_outcome %in% c("e","f")))
  
  naive <- stage1_cox$coefficients
  naive[select_index] <- S(select_index,end_stage)/I(select_index,end_stage)
  stage2 <- (S(select_index,2)-S(select_index,1))/(I(select_index,2)-I(select_index,1))
  
  cov_S <- matrix(0, nrow = K+J-1, ncol = K+J-1)
  I_stage1 <- diag(mapply(I, 1:K,1))
  cov_stage1 <- stage1_cox$var
  diag(cov_stage1) <- 1/mapply(I, 1:K,1)
  cov_S[1:K,1:K] <- I_stage1 %*% cov_stage1 %*% t(I_stage1)
  if(end_stage > 1){
    if(end_stage>2){
      diag(cov_S[(K+1):(K+end_stage-1),(K+1):(K+end_stage-1)]) <- diff(mapply(I,select_index,1:end_stage))
    }else{
      cov_S[(K+1):(K+end_stage-1),(K+1):(K+end_stage-1)] <- diff(mapply(I,select_index,1:end_stage))
    }}
  if(J > end_stage){
    if(J-end_stage > 1){
      diag(cov_S[(K+end_stage):(K+J-1),(K+end_stage):(K+J-1)]) <- 1/4*diff(D_n)[end_stage:(J-1)]
    }else{
      cov_S[(K+end_stage):(K+J-1),(K+end_stage):(K+J-1)] <- 1/4*diff(D_n)[end_stage:(J-1)]
    }}
  S_all=c(mapply(S,1:K,1),S(select_index,end_stage))
  list(end_stage=end_stage,select_index=select_index,test=inter_outcome[end_stage],naive=naive,
       bias=c(stage2,naive[select_index]),cov_S=cov_S,S_all=S_all)
  # c(end_stage=end_stage,select_index=select_index,naive,
  #      bias=c(stage2,naive[select_index])-theta[select_index])
}
da_genc <- cmpfun(da_gen)

# for(setting in c("H1","NULL","peak","line")){
#   aa <- pbmapply(da_gen, sim=1:sim,setting,J)
#   print(table(aa[3,]))
# }
# aa <- pbmapply(da_gen, sim=1:10000,"LFC",J)
# table(aa[2,])

S_cum_f <- function(cov_S,select_index,K,J){
  trans_t <- diag(1,K+J-1,K+J-1)
  trans_t[1:K,1:K][-K,] <- diag(nrow = K, ncol = K)[-select_index,]
  trans_t[1:K,1:K][K,]  <- diag(nrow = K, ncol = K)[select_index,]
  
  trans_c <- diag(1,K+J-1,K+J-1)
  trans_c[K:(K+J-1),K:(K+J-1)][lower.tri(trans_c[K:(K+J-1),K:(K+J-1)])] <- 1
  list(S_sum = trans_c %*% (trans_t %*% cov_S %*% t(trans_t)) %*% t(trans_c), 
       I_stage = diag(trans_c %*% (trans_t %*% cov_S %*% t(trans_t)) %*% t(trans_c)))
}
S_Z_diff_f <- function(cov_S,theta,select_index,K,J){
  I_stage <- S_cum_f(cov_S,select_index,K,J)$I_stage
  theta_v <- diag(1/sqrt(I_stage),K+J-1,K+J-1)
  cov_Z = theta_v %*% S_cum_f(cov_S,select_index,K,J)$S_sum %*%t(theta_v) 
  
  A <- diag(1,K+J-1)
  diag(A[1:(K-1),1:(K-1)]) <- -1
  A[1:(K-1),K] <- 1
  
  mu_Z <- c(theta[-select_index]*sqrt(I_stage[(1:(K-1))]),
            theta[select_index]*sqrt(I_stage[K:(K+J-1)]))
  list(mu=A %*% mu_Z,cov=A %*% cov_Z %*% t(A))
}
##### MI & SI #####
bias_theta  <- function(select_index,theta,K,J,cov_S){
  if(J==4){
    low_bround <- c(0.000000,  4.013691,  2.828157, -2.285393)
    up_bround <- c(-Inf, -4.013691, -2.828157, -2.285393)
  }else if(J==3){
    low_bround <- c(0.000000,  3.301594, -2.302597)
    up_bround <- c(-Inf, -3.301594, -2.302597)
  }else if(J==5){
    low_bround <- c(0.000000,  4.609252,  3.264385,  2.644414, -2.278346)
    up_bround <- c(-Inf, -4.609252, -3.264385, -2.644414, -2.278346)
  }
  I_stage <- S_cum_f(cov_S,select_index,K,J)$I_stage
  mu_cov <- S_Z_diff_f(cov_S,theta,select_index,K,J)
  mu <- mu_cov$mu; cov <- mu_cov$cov
  
  prob_Q <- pmvnorm(lower = c(rep(-Inf,K-1),up_bround[1]),upper = c(rep(0,K-1),low_bround[1]),
                    mean = mu[1:K,],sigma = cov[1:K,1:K])[1]
  
  bias <- rep(NA,K)
  if(!is.na(prob_Q)){
    mtm <- mtmvnorm(mean = mu[1:K,],sigma = cov[1:K,1:K],
                    lower = c(rep(-Inf,K-1),up_bround[1]),upper = c(rep(0,K-1),low_bround[1]),
                    doComputeVariance=FALSE)$tmean
    bias[-select_index] <- (mtm[K]-mtm[-K])/sqrt(I_stage[1:(K-1)])-theta[-select_index]
    E_select <- function(stage){
      mu_t <- mu[1:(K+stage-1),]
      cov_t <- cov[1:(K+stage-1),1:(K+stage-1)]
      lower1 <- c(rep(-Inf,K-1),up_bround[1:(stage-1)],-Inf); upper1 <- c(rep(0,K-1),low_bround[1:(stage-1)],up_bround[stage])
      lower2 <- c(rep(-Inf,K-1),up_bround[1:(stage-1)],low_bround[stage]);upper2 <- c(rep(0,K-1),low_bround[1:(stage-1)],Inf)
      sum(pmvnorm(lower = lower1,upper = upper1,mean = mu_t,sigma = cov_t)[1]*
            mtmvnorm(mean = mu_t,sigma = cov_t,lower = lower1,upper = upper1,doComputeVariance=FALSE)$tmean[K+stage-1],
          pmvnorm(lower = lower2,upper = upper2,mean = mu_t,sigma = cov_t)[1]*
            mtmvnorm(mean = mu_t,sigma = cov_t,lower = lower2,upper = upper2,doComputeVariance=FALSE)$tmean[K+stage-1],na.rm = T)/
        sqrt(I_stage[K:(K+J-1)][stage])
    }
    bias[select_index] <- sum(mapply(E_select, stage=2:J))/prob_Q-theta[select_index]
  }
  bias
}
MI <- function(initial,select_index,max_iterations=20,tol=0.001,K,J,cov_S){
  conver <- 0
  theta_hat <- initial
  for (i in 1:max_iterations) {
    bias_temp <- bias_theta(select_index,theta_hat,K,J,cov_S)
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
SI <- function(initial,select_index,K,J,cov_S){
  (initial-bias_theta(select_index,initial,K,J,cov_S))[select_index]
}  
MIc <- cmpfun(MI)

##### MUE #####
est_MUE <- function(naive,end_stage,select_index,sub,K,J,cov_S,S_all){
  if(J==4){
    low_bround <- c(0.000000,  4.013691,  2.828157, -2.285393)
    up_bround <- c(-Inf, -4.013691, -2.828157, -2.285393)
  }else if(J==3){
    low_bround <- c(0.000000,  3.301594, -2.302597)
    up_bround <- c(-Inf, -3.301594, -2.302597)
  }else if(J==5){
    low_bround <- c(0.000000,  4.609252,  3.264385,  2.644414, -2.278346)
    up_bround <- c(-Inf, -4.609252, -3.264385, -2.644414, -2.278346)
  }
  Z_end <- S_all[K+1]/sqrt(S_cum_f(cov_S,select_index,K,J)$I_stage[K+end_stage-1])
  
  cond_pvalue <- function(theta,select_index,end_stage,Z_end,cov_S,K,J){
    mu_cov <- S_Z_diff_f(cov_S,theta,select_index,K,J)
    mu <- mu_cov$mu; cov <- mu_cov$cov
    
    prob_Q <- mvtnorm::pmvnorm(lower = c(rep(-Inf,K-1),up_bround[1]),
                               upper = c(rep(0,K-1),low_bround[1]),
                               mean = mu[1:K,],sigma = cov[1:K,1:K])[1]
    p_value <- function(end_stage,select_index,Z_end){
      P_stage <- function(stage,stop,Z_end){
        new_mu_s <- mu[1:(K+stage-1),]; new_cov_s=cov[1:(K+stage-1),1:(K+stage-1)]
        if(stop==T){
          P <-  mvtnorm::pmvnorm(lower = c(rep(-Inf,K-1),up_bround[1:(stage-1)],-Inf),
                                 upper = c(rep(0,K-1),low_bround[1:(stage-1)],Z_end),
                                 mean = new_mu_s,sigma= new_cov_s)[1]      
        }else if(stop==F){
          P <-   mvtnorm::pmvnorm(lower = c(rep(-Inf,K-1),up_bround[1:(stage-1)],-Inf),
                                  upper = c(rep(0,K-1),low_bround[1:(stage-1)],up_bround[stage]),
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
    cond_pvaluec(theta,select_index,end_stage,Z_end,cov_S,K,J)
  }
  a <- optimise(function(x) (P_value_sub(select_index,sub,x)-1/2)^2,c(-2,2)*naive[select_index],tol=0.001)
  c(a$minimum,a$objective)
}

##### RB #####
RB_mc <- function(S_all,cov_S,select_index,end_stage,mc_num,alpha,K,J){
  if(J==4){
    low_bround <- c(0.000000,  4.013691,  2.828157, -2.285393)
    up_bround <- c(-Inf, -4.013691, -2.828157, -2.285393)
  }else if(J==3){
    low_bround <- c(0.000000,  3.301594, -2.302597)
    up_bround <- c(-Inf, -3.301594, -2.302597)
  }else if(J==5){
    low_bround <- c(0.000000,  4.609252,  3.264385,  2.644414, -2.278346)
    up_bround <- c(-Inf, -4.609252, -3.264385, -2.644414, -2.278346)
  }
  I_stage <- S_cum_f(cov_S,select_index,K,J)$I_stage
  mu <- rep(0,K+end_stage-1)
  I_stage1 <- diag(cov_S)[1:K]
  
  trans_T <- diag(1,K+end_stage-1,K+end_stage-1)
  trans_T[select_index,] <- 1
  trans_T[1:K,1:K] <- solve(cov_S[1:K,1:K])*I_stage1
  mu_T <- trans_T %*% mu
  cov_T <- trans_T %*% (cov_S[1:(K+end_stage-1),1:(K+end_stage-1)]) %*% t(trans_T)
  
  T_star <- (I_stage1*solve(cov_S[1:K,1:K])) %*% as.matrix(S_all[1:K])
  T_star[select_index] <- T_star[select_index]+(S_all[K+1]-S_all[select_index])
  
  mc_s <- rcmvnorm(mc_num,mu_T,cov_T,given.ind = 1:K,dependent.ind = (K+1):(K+end_stage-1),X.given = T_star)
  s1 <- matrix(NA, nrow = K, ncol = mc_num)
  for(i in 1:mc_num){
    T_star_s <- T_star
    T_star_s[select_index] <- T_star[select_index]-sum(mc_s[i,])
    s1[,i] <- solve((I_stage1*solve(cov_S[1:K,1:K]))) %*% T_star_s
  }
  re <- rep(NA,4)
  z_stage1 <- apply(s1, 2, function(x) x/sqrt(I_stage1))
  mc_s_all <- cbind(s1[select_index,],mc_s)
  mc_z <- t(apply(mc_s_all, 1, function(x) cumsum(x)/sqrt(I_stage[K:(K+end_stage-1)])))
  bround_L <- matrix(low_bround[1:(end_stage-1)],mc_num,(end_stage-1),byrow = T)
  bround_U <- matrix(up_bround[1:(end_stage-1)],mc_num,(end_stage-1),byrow = T)
  
  index <- apply(cbind(apply(mc_z[,1:(end_stage-1)] >= bround_U,1,function(x)identical(x,y=rep(TRUE,end_stage-1))),
                       apply(mc_z[,1:(end_stage-1)] <= bround_L,1,function(x)identical(x,y=rep(TRUE,end_stage-1))),
                       apply(z_stage1,2,which.min)==select_index), 1, all)
  if(sum(index)>1){
    mc_est <- mean(mc_s_all[index,2]/(I_stage[K+1]-I_stage[K]))
    SE <- sqrt(1/(I_stage[K+1]-I_stage[K])-var(mc_s_all[index,2]/(I_stage[K+1]-I_stage[K])))
    re <- c(sum(index)/mc_num,mc_est,mc_est+qnorm(alpha/2)*SE,mc_est-qnorm(alpha/2)*SE)
  }
  re
}


setwd("C:/Users/13379/Desktop/essay/submit-SIM/survival-cov/4")
shape <- 0.5
setting <- "line"
a_t <- 1
N <- 1000
J <- 4
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
                              "pbapply","tmvtnorm","BSDA","pracma","nnet","flexsurv")) %dopar% {
                                da_genc(i,setting,J)
                              }
end_stage <- apply(res,2,function(x) x$end_stage)
select_index <- apply(res,2,function(x) x$select_index)
naive <- apply(res,2,function(x) x$naive)
bias <- apply(res,2,function(x) x$bias)
write.csv(cbind(end_stage,select_index,t(bias))[which(end_stage>=2),],paste(setting,"stage2-naive.csv",sep = "-"))

re_SI <- foreach(i=which(end_stage>=2),
                 .options.snow=opts,
                 .combine = cbind,
                 .packages = c("MASS","data.table","compiler","condMVNorm",
                               "pbapply","tmvtnorm","BSDA","pracma","nnet","flexsurv")) %dopar% {
                                 SI(res[,i]$naive,res[,i]$select_index,K=4,J,cov_S=res[,i]$cov_S)
                               }
write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],as.numeric(re_SI)),paste(setting,"SI.csv",sep = "-"))

re_MUE_MLE <- foreach(i=which(end_stage>=2),
                      .options.snow=opts,
                      .combine = cbind,
                      .packages = c("MASS","data.table","compiler","condMVNorm",
                                    "pbapply","tmvtnorm","BSDA","pracma","nnet","flexsurv")) %dopar% {
                                      est_MUE(res[,i]$naive,end_stage[i],res[,i]$select_index,sub="MLE",
                                              K=4,J,res[,i]$cov_S,res[,i]$S_all)
                                    }
write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],t(re_MUE_MLE)),paste(setting,"MUE_MLE.csv",sep = "-"))

re_MUE_ZERO <- foreach(i=which(end_stage>=2),
                       .options.snow=opts,
                       .combine = cbind,
                       .packages = c("MASS","data.table","compiler","condMVNorm",
                                     "pbapply","tmvtnorm","BSDA","pracma","nnet","flexsurv")) %dopar% {
                                       est_MUE(res[,i]$naive,end_stage[i],res[,i]$select_index,sub="ZERO",
                                               K=4,J,res[,i]$cov_S,res[,i]$S_all)
                                     }
write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],t(re_MUE_ZERO)),paste(setting,"MUE_ZERO.csv",sep = "-"))

re_MI <- foreach(i=which(end_stage>=2),
                 .options.snow=opts,
                 .combine = cbind,
                 .packages = c("MASS","data.table","compiler","condMVNorm",
                               "pbapply","tmvtnorm","BSDA","pracma","nnet","flexsurv")) %dopar% {
                                 MIc(res[,i]$naive,res[,i]$select_index,max_iterations=20,tol=0.001,
                                     K=4,J,cov_S=res[,i]$cov_S)
                               }
write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],t(re_MI)),paste(setting,"MI.csv",sep = "-"))

re_RB <- foreach(i=which(end_stage>=2),
                 .options.snow=opts,
                 .combine = cbind,
                 .packages = c("MASS","data.table","compiler","condMVNorm",
                               "pbapply","tmvtnorm","BSDA","pracma","nnet")) %dopar% {
                                 RB_mc(res[,i]$S_all,res[,i]$cov_S,select_index[i],end_stage[i],
                                       mc_num=10000,alpha=0.05,K=4,J)
                               }
write.csv(cbind(end_stage[which(end_stage>=2)],select_index[which(end_stage>=2)],t(re_RB)),paste(setting,"RB.csv",sep = "-"))


close(pb)
stopCluster(cls)
time2 <- Sys.time()
print(time2-time1)
