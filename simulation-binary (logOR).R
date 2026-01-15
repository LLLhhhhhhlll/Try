###############################################################################
## Seamless phase II/III simulation (multi-arm, multi-stage early stopping)
##
## This script:
##   1) Simulates raw trial data under different effect-size configurations
##   2) Implements several estimators (Naive, Stage2, CBAE-SI/MI, CMUE, UMVCUE, RB)
##   3) Runs Monte Carlo simulation in parallel and saves results to CSV files
##
## Main inputs (global):
##   setting: "NULL", "H1", "peak", "line"  (effect-size configurations)
##   J:       total number of stages (3/4/5)
##   sim:     number of simulation replicates
###############################################################################

rm(list = ls())


###############################################################################
## 0) Load required packages
###############################################################################
library(MASS)         # mvrnorm (not heavily used in this binary part, but kept)
library(data.table)   # first()
library(compiler)     # cmpfun() for speed
library(condMVNorm)   # rcmvnorm() conditional MVN sampling (RB method)
library(pbapply)      # progress helpers (optional)
library(tmvtnorm)     # mtmvnorm() truncated MVN mean (CBAE bias)
library(doParallel)   # parallel backend
library(foreach)      # foreach
library(doSNOW)       # foreach + progress bar
library(tcltk)        # txtProgressBar
library(mvtnorm)      # pmvnorm() MVN probabilities (used extensively)

###############################################################################
## 1) Data generator for one simulation replicate
##    - Generate patient-level binary outcomes for control + K experimental arms
##    - Select best arm based on Z-statistic at end of phase II (stage 1)
##    - Apply group-sequential boundaries in phase III (stages 2..J)
##    - Output sufficient components used by subsequent estimators
###############################################################################
da_gen <- function(sim, setting, J){
  
  set.seed(20250425 + sim)
  
  ## Fixed control response rate
  p_control <- 0.4
  
  ## Scenario: treatment response rates for K=4 arms
  if(setting=="NULL"){
    p_treatment <- rep(0.4, 4)
  }else if(setting=="H1"){
    p_treatment <- c(0.6, 0.6, 0.6, 0.6)
  }else if(setting=="peak"){
    p_treatment <- c(0.5, 0.5, 0.5, 0.6)
  }else if(setting=="line"){
    p_treatment <- c(0.45, 0.50, 0.55, 0.60)
  }else{
    stop("Unknown setting.")
  }
  
  ## True log-OR vs control
  theta <- log(p_treatment/(1-p_treatment)) - log(p_control/(1-p_control))
  
  if(J==4){
    t <- c(1,2,3,4)/4
    low_bround <- c(0.00000000, 0.08743326, 1.41906914, 2.256008)
    up_bround  <- c(Inf, 4.013693, 2.826977, 2.256008)
    unit_n <- 39
  }else if(J==3){
    t <- c(1,2,3)/3
    low_bround <- c(0.000000, 1.015946, 2.291169)
    up_bround  <- c(Inf, 3.301556, 2.291169)
    unit_n <- 47
  }else if(J==5){
    t <- c(1,2,3,4,5)/5
    low_bround <- c(0.0000000, -0.7601533, 0.6465035, 1.5388462, 2.240586)
    up_bround  <- c(Inf, 4.609263, 3.263943, 2.639832, 2.240586)
    unit_n <- 34
  }
  
  ## Cumulative sample size at each stage
  n <- t * length(t) * unit_n
  
  K <- length(p_treatment)
  
  ## all_da is cumulative event counts at each stage:
  ## rows = stages 1..J
  ## cols = (control, arm1..armK)
  
  all_da <- apply(
    replicate(J, mapply(function(i) rbinom(1, unit_n, c(p_control, p_treatment)[i]),
                        i = 1:(K+1))),
    1, cumsum
  )
  
  ## Define score-like statistic S and information I for difference in proportions
  S <- function(E, C, N) (N*E - N*C) / (2*N)
  I <- function(E, C, N) (E + C) * (2*N - E - C) * N^2 / (2*N)^3
  Z <- function(E, C, N) S(E, C, N) / sqrt(I(E, C, N))
  
  ## Phase II selection at stage 1: choose arm with largest Z
  select_index <- which.is.max(Z(all_da[1, -1], all_da[1, 1], n[1]))
  
  ## Interim decision at stage j based on selected arm Z
  inter_decide <- function(j, select_index, da){
    Z_temp <- Z(all_da[j, select_index+1], all_da[j, 1], n[j])
    if(Z_temp < low_bround[j]) "f" else if(Z_temp > up_bround[j]) "e" else "c"
  }
  
  ## Apply decision rule sequentially
  inter_outcome <- character(J)
  for(i in 1:J){
    inter_outcome[i] <- inter_decide(j=i, select_index=select_index, da=all_da)
  }
  
  ## Trial stops at first efficacy/futility crossing
  end_stage <- data.table::first(which(inter_outcome %in% c("e","f")))
  
  ## Stage-2 increment estimator (post-selection), expressed on S/I scale:
  ## uses increments between stage 1 and 2 counts (unit_n patients)
  stage2_est <- S(all_da[2, select_index+1]-all_da[1, select_index+1],
                  all_da[2, 1]-all_da[1, 1], unit_n) /
    I(all_da[2, select_index+1]-all_da[1, select_index+1],
      all_da[2, 1]-all_da[1, 1], unit_n)
  
  ## naive_P: plug-in probabilities used later as "p_control" and "p_treatment" inputs
  ## - start from stage-1 sample proportions
  ## - replace control and selected arm with cumulative proportions at end_stage
  naive_P <- all_da[1, ] / unit_n
  naive_P[1] <- all_da[end_stage, 1] / n[end_stage]
  naive_P[select_index+1] <- all_da[end_stage, select_index+1] / n[end_stage]
  
  ## naive_est: stage-1 estimator for each arm, except selected uses cumulative to end_stage
  naive_est <- S(all_da[1, -1], all_da[1, 1], n[1]) / I(all_da[1, -1], all_da[1, 1], n[1])
  naive_est[select_index] <- S(all_da[end_stage, select_index+1], all_da[end_stage, 1], n[end_stage]) /
    I(all_da[end_stage, select_index+1], all_da[end_stage, 1], n[end_stage])
  
  ## x_bar: “observed sufficient summaries” used in UMVCUE (cumulative counts)
  x_bar <- all_da[1, ]
  x_bar[1] <- all_da[end_stage, 1]
  x_bar[select_index+1] <- all_da[end_stage, select_index+1]
  
  ## s_all: score-like summaries for RB
  ## here: (K arm stage-1 scores) + (selected arm cumulative score at end_stage)
  s_all <- c(S(all_da[1, -1], all_da[1, 1], n[1]),
             S(all_da[end_stage, select_index+1], all_da[end_stage, 1], n[end_stage]))
  
  ## Return vector; downstream slices it into blocks
  c(end_stage, select_index,
    naive_P,          # length K+1 (=5)
    naive_est,        # length K (=4)
    x_bar,            # length K+1 (=5)
    s_all,            # length K+1 (=5)
    theta[select_index],  # true theta of selected arm (for bias checking)
    stage2_est,
    naive_est[select_index])
}
da_genc <- cmpfun(da_gen)

###############################################################################
## 2) Estimator implementations (CBAE-SI/MI, CMUE, UMVCUE, RB)
###############################################################################

## --- 2.1 Mean vector and covariance matrix for CBAE/CMUE under MVN framework ---
mu_cov_f <- function(p_control, p_treatment, select_index, theta,
                     K, J, n, I_type = "p_bar"){
  
  ## If theta is not given, use log-odds ratio between treatment and control
  if(is.null(theta)){
    theta <- log(p_treatment/(1-p_treatment)) - log(p_control/(1-p_control))
  }
  
  ## Fisher information choice:
  ##   - "p_bar": use pooled p_bar = (pE+pC)/2 for information
  ##   - otherwise: use separate pE and pC components
  if(I_type == "p_bar"){
    I <- function(i, j){
      p_bar <- (p_treatment + p_control)/2
      n[j]/2 * p_bar[i] * (1 - p_bar[i])
    }
  }else{
    I <- function(i, j){
      n[j]/4 * p_treatment[i] * (1 - p_treatment[i]) +
        n[j]/4 * p_control * (1 - p_control)
    }
  }
  
  ## Build covariance matrix for score-like statistics:
  ## Dimension: (K + J - 1) x (K + J - 1)
  ##   - first K components: stage-1 arm-wise statistics
  ##   - remaining (J-1) components: increments for selected arm across stages
  cov <- cbind(
    matrix(n[1]/4 * p_control * (1 - p_control), K+J-1, K, byrow=TRUE),
    matrix(0, K+J-1, J-1, byrow=TRUE)
  )
  
  ## Stage-1 variance terms:
  ## Here the diagonal is set to I(...) at stage 1
  diag(cov[1:K, 1:K]) <- I(c((1:K)[-select_index], select_index), 1)
  
  ## Increment variances for selected arm: diff(I) across cumulative stages
  diag(cov[(K+1):(K+J-1), (K+1):(K+J-1)]) <- diff(I(select_index, 1:J))
  
  ## Symmetrize
  cov[lower.tri(cov)] <- t(cov)[lower.tri(cov)]
  
  ## Mean vector for the score-like statistics
  mu <- c(
    theta[-select_index] * I((1:K)[-select_index], 1),
    theta[select_index] * c(I(select_index, 1), diff(I(select_index, 1:J)))
  )
  
  ## Transform increments into cumulative sums (so you can express Z boundaries)
  trans_s <- diag(1, K+J-1)
  trans_s[K:(K+J-1), K:(K+J-1)][lower.tri(trans_s[K:(K+J-1), K:(K+J-1)])] <- 1
  
  ## Standardization to Z-scale:
  ## - first (K-1) are compared at stage 1
  ## - remaining are cumulative Z’s for selected arm across stages
  theta_v <- matrix(0, K+J-1, K+J-1)
  diag(theta_v[1:(K-1), 1:(K-1)]) <- 1/sqrt(I((1:K)[-select_index], 1))
  diag(theta_v[K:(K+J-1), K:(K+J-1)]) <- 1/sqrt(I(select_index, 1:J))
  
  mu_Z  <- theta_v %*% trans_s %*% mu
  cov_Z <- theta_v %*% (trans_s %*% cov %*% t(trans_s)) %*% t(theta_v)
  
  ## A: contrast transform for selection event:
  ## for non-selected arms, use (Z_selected - Z_other) >= 0
  A <- diag(1, K+J-1)
  diag(A[1:(K-1), 1:(K-1)]) <- -1
  A[1:(K-1), K] <- 1
  
  list(mu = A %*% mu_Z,
       cov = A %*% cov_Z %*% t(A))
}

###############################################################################
## MI / SI / CMUE / UMVCUE / RB estimators under seamless phase II/III design
## (multi-arm selection in phase II + group-sequential stopping in phase III)
###############################################################################

##### -------------------------------------------------------------------------
##### (A) CBAE: conditional bias-adjusted estimator
#####     - SI: single-step bias correction
#####     - MI: multi-step (iterative) bias correction
##### -------------------------------------------------------------------------

## (A1) bias_theta(): compute conditional bias term b_theta( theta_hat_MLE | Q )
##      where Q denotes the selection event (best arm selected + continuation).
##      This is the key ingredient for CBAE-SI and CBAE-MI.
bias_theta <- function(select_index, theta, K, J, p_control, p_treatment){
  
  ## Design parameters (same as simulation)
  if(J==4){
    t <- c(1,2,3,4)/4
    low_bround <- c(0.00000000, 0.08743326, 1.41906914, 2.256008)
    up_bround  <- c(Inf, 4.013693, 2.826977, 2.256008)
    unit_n <- 39
  }else if(J==3){
    t <- c(1,2,3)/3
    low_bround <- c(0.000000, 1.015946, 2.291169)
    up_bround  <- c(Inf, 3.301556, 2.291169)
    unit_n <- 47
  }else if(J==5){
    t <- c(1,2,3,4,5)/5
    low_bround <- c(0.0000000, -0.7601533, 0.6465035, 1.5388462, 2.240586)
    up_bround  <- c(Inf, 4.609263, 3.263943, 2.639832, 2.240586)
    unit_n <- 34
  }else{
    stop("Only J=3/4/5 supported.")
  }
  n <- t * length(t) * unit_n
  
  ## Information approximation (pooled p_bar) used for Z-scaling inside bias
  I <- function(p_E, p_C, N){
    p_bar <- (p_E + p_C)/2
    N/2 * p_bar * (1 - p_bar)
  }
  
  ## MVN mean/cov of transformed Z vector under theta
  mu_cov <- mu_cov_f(p_control, p_treatment, select_index, theta, K, J, n)
  mu <- mu_cov$mu
  cov <- mu_cov$cov
  
  ## Selection event probability prob_Q
  prob_Q <- pmvnorm(
    lower = c(rep(0, K-1), low_bround[1]),
    upper = c(rep(Inf, K-1), up_bround[1]),
    mean  = mu[1:K,],
    sigma = cov[1:K, 1:K]
  )[1]
  
  bias <- rep(NA, K)
  
  if(!is.na(prob_Q)){
    
    ## Truncated MVN mean for the first K components under selection region
    mtm <- mtmvnorm(
      mean = mu[1:K,], sigma = cov[1:K, 1:K],
      lower = c(rep(0, K-1), low_bround[1]),
      upper = c(rep(Inf, K-1), up_bround[1]),
      doComputeVariance = FALSE
    )$tmean
    
    ## Bias for non-selected arms
    bias[-select_index] <-
      (mtm[K] - mtm[-K]) / sqrt(I(p_treatment[-select_index], p_control, n[1])) -
      theta[-select_index]
    
    ## Bias for selected arm averaged over stopping stage 2..J
    E_select <- function(stage){
      
      mu_t  <- mu[1:(K+stage-1),]
      cov_t <- cov[1:(K+stage-1), 1:(K+stage-1)]
      
      lower1 <- c(rep(0, K-1), low_bround[1:(stage-1)], -Inf)
      upper1 <- c(rep(Inf, K-1), up_bround[1:(stage-1)], low_bround[stage])
      
      lower2 <- c(rep(0, K-1), low_bround[1:(stage-1)], up_bround[stage])
      upper2 <- c(rep(Inf, K-1), up_bround[1:(stage-1)], Inf)
      
      (
        pmvnorm(lower1, upper1, mean=mu_t, sigma=cov_t)[1] *
          mtmvnorm(mean=mu_t, sigma=cov_t, lower=lower1, upper=upper1, doComputeVariance=FALSE)$tmean[K+stage-1] +
          pmvnorm(lower2, upper2, mean=mu_t, sigma=cov_t)[1] *
          mtmvnorm(mean=mu_t, sigma=cov_t, lower=lower2, upper=upper2, doComputeVariance=FALSE)$tmean[K+stage-1]
      ) / sqrt(I(p_treatment[select_index], p_control, n[stage]))
    }
    
    bias[select_index] <- sum(mapply(E_select, stage=2:J)) / prob_Q - theta[select_index]
  }
  
  bias
}

## (A2) MI: iterative bias correction
##      theta_{new} = initial - bias_theta(theta_old) until convergence.
MI <- function(initial, select_index, max_iterations=20, tol=0.001,
               K, J, p_control, p_treatment){
  
  conver <- 0
  theta_hat <- initial
  
  for(i in 1:max_iterations){
    bias_temp <- bias_theta(select_index, theta_hat, K, J, p_control, p_treatment)
    
    if(is.na(sum(bias_temp)) || is.infinite(sum(bias_temp))){
      theta_hat <- initial
      conver <- -1
      break()
    }
    
    theta_hat_new <- initial - bias_temp
    euc.dis <- sqrt(sum((theta_hat - theta_hat_new)^2))
    
    if(euc.dis <= tol){
      theta_hat <- theta_hat_new
      conver <- 1
      break()
    }else{
      theta_hat <- theta_hat_new
    }
  }
  
  c(conver, theta_hat[select_index])
}

## (A3) SI: single-step bias correction
SI <- function(initial, select_index, K, J, p_control, p_treatment){
  (initial - bias_theta(select_index, initial, K, J, p_control, p_treatment))[select_index]
}
MIc <- cmpfun(MI)  # compile MI for speed

##### -------------------------------------------------------------------------
##### (B) CMUE: conditional median unbiased estimator
#####     - compute conditional p-value under candidate theta
#####     - invert p-value to get CMUE point estimate (median) and CI endpoints
##### -------------------------------------------------------------------------
est_MUE <- function(naive, end_stage, select_index, sub,
                    K, J, p_control, p_treatment){
  
  ## Design parameters (same as simulation)
  if(J==4){
    t <- c(1,2,3,4)/4
    low_bround <- c(0.00000000, 0.08743326, 1.41906914, 2.256008)
    up_bround  <- c(Inf, 4.013693, 2.826977, 2.256008)
    unit_n <- 39
  }else if(J==3){
    t <- c(1,2,3)/3
    low_bround <- c(0.000000, 1.015946, 2.291169)
    up_bround  <- c(Inf, 3.301556, 2.291169)
    unit_n <- 47
  }else if(J==5){
    t <- c(1,2,3,4,5)/5
    low_bround <- c(0.0000000, -0.7601533, 0.6465035, 1.5388462, 2.240586)
    up_bround  <- c(Inf, 4.609263, 3.263943, 2.639832, 2.240586)
    unit_n <- 34
  }else{
    stop("Only J=3/4/5 supported.")
  }
  n <- t * length(t) * unit_n
  
  Z_end <- sqrt(n[end_stage]/2) * (p_treatment[select_index] - p_control) /
    sqrt(((p_treatment[select_index] + p_control)/2) * (1 - (p_treatment[select_index] + p_control)/2))
  
  cond_pvalue <- function(theta, select_index, end_stage, Z_end, p_control, p_treatment){
    
    mu_cov <- mu_cov_f(p_control, p_treatment, select_index, theta, K, J, n)
    mu <- mu_cov$mu
    cov <- mu_cov$cov
    
    prob_Q <- pmvnorm(
      lower = c(rep(0, K-1), low_bround[1]),
      upper = c(rep(Inf, K-1), up_bround[1]),
      mean  = mu[1:K,],
      sigma = cov[1:K, 1:K]
    )[1]
    
    p_value <- function(end_stage, select_index, Z_end){
      
      P_stage <- function(stage, stop, Z_end){
        new_mu  <- mu[1:(K+stage-1),]
        new_cov <- cov[1:(K+stage-1), 1:(K+stage-1)]
        
        if(stop){
          P <- pmvnorm(
            lower = c(rep(0,K-1), low_bround[1:(stage-1)], Z_end),
            upper = c(rep(Inf,K-1), up_bround[1:(stage-1)], Inf),
            mean  = new_mu, sigma = new_cov
          )[1]
        }else{
          P <- pmvnorm(
            lower = c(rep(0,K-1), low_bround[1:(stage-1)], up_bround[stage]),
            upper = c(rep(Inf,K-1), up_bround[1:(stage-1)], Inf),
            mean  = new_mu, sigma = new_cov
          )[1]
        }
        if(is.na(P)) 0 else P
      }
      
      if(end_stage == 2){
        numerator <- P_stage(end_stage, TRUE, Z_end)
      }else{
        numerator <- sum(mapply(P_stage, stage=2:(end_stage-1), stop=FALSE, Z_end=Z_end),
                         P_stage(end_stage, TRUE, Z_end))
      }
      
      if(is.na(numerator) || is.na(prob_Q) || prob_Q <= 0) return(1)
      min(1, numerator/prob_Q)
    }
    
    p_value(end_stage, select_index, Z_end)
  }
  
  cond_pvaluec <- cmpfun(cond_pvalue)
  
  ## Plug-in for non-selected arms:
  ##   sub="MLE": fix them at their naive estimates
  ##   sub="ZERO": fix them at 0
  P_value_sub <- function(select_index, sub="MLE", x){
    theta <- rep(x, K)
    if(sub == "MLE"){
      theta[-select_index] <- naive[-select_index]
    }else{
      theta[-select_index] <- rep(0, K-1)
    }
    cond_pvaluec(theta, select_index, end_stage, Z_end, p_control, p_treatment)
  }
  
  ## Solve p(theta)=0.5 (median) using optimization on [-1,1]
  a <- optimise(function(x) (P_value_sub(select_index, sub, x) - 1/2)^2,
                c(-1,1), tol=0.001)
  c(a$minimum, a$objective)
}


##### -------------------------------------------------------------------------
##### (C) RB: Rao–Blackwell estimator via conditional simulation
#####     - condition an unbiased stage-wise estimator on sufficient statistic
#####     - Monte Carlo implementation of conditional expectation
##### -------------------------------------------------------------------------
RB_mc <- function(s_all, p_control, p_treatment, select_index, end_stage,
                  mc_num, alpha, K, J){
  
  ## Design parameters
  if(J==4){
    t <- c(1,2,3,4)/4
    low_bround <- c(0.00000000, 0.08743326, 1.41906914, 2.256008)
    up_bround  <- c(Inf, 4.013693, 2.826977, 2.256008)
    unit_n <- 39
  }else if(J==3){
    t <- c(1,2,3)/3
    low_bround <- c(0.000000, 1.015946, 2.291169)
    up_bround  <- c(Inf, 3.301556, 2.291169)
    unit_n <- 47
  }else if(J==5){
    t <- c(1,2,3,4,5)/5
    low_bround <- c(0.0000000, -0.7601533, 0.6465035, 1.5388462, 2.240586)
    up_bround  <- c(Inf, 4.609263, 3.263943, 2.639832, 2.240586)
    unit_n <- 34
  }else{
    stop("Only J=3/4/5 supported.")
  }
  n <- t * length(t) * unit_n
  
  ## Working model: theta=0 (used for conditional simulation)
  theta <- rep(0, K)
  
  ## Information approximation for stage j under pooled p_bar
  I <- function(i, j){
    p_bar <- (p_treatment + p_control)/2
    n[j]/2 * p_bar[i] * (1 - p_bar[i])
  }
  
  ## Build covariance for stage-1 arms + increments for selected arm
  cov <- cbind(
    matrix(n[1]/4 * p_control * (1 - p_control), K+end_stage-1, K, byrow=TRUE),
    matrix(0, K+end_stage-1, end_stage-1, byrow=TRUE)
  )
  diag(cov[1:K, 1:K]) <- I(1:K, 1)
  
  if(end_stage > 2){
    diag(cov[(K+1):(K+end_stage-1), (K+1):(K+end_stage-1)]) <- diff(I(select_index, 1:end_stage))
  }else{
    cov[K+end_stage-1, K+end_stage-1] <- diff(I(select_index, 1:end_stage))
  }
  cov[lower.tri(cov)] <- t(cov)[lower.tri(cov)]
  
  ## Mean vector (score scale): here all zeros under working model
  mu <- c(theta * I(1:K, 1),
          theta[select_index] * diff(I(select_index, 1:end_stage)))
  
  ## Transform to sufficient statistic T and condition on observed T_star
  trans_T <- diag(1, K+end_stage-1)
  trans_T[select_index,] <- 1
  trans_T[1:K, 1:K] <- solve(cov[1:K, 1:K]) * I(1:K, 1)
  
  mu_T  <- trans_T %*% mu
  cov_T <- trans_T %*% cov %*% t(trans_T)
  
  ## Observed T_star from observed s_all
  T_star <- (I(1:K,1) * solve(cov[1:K,1:K])) %*% as.matrix(s_all[1:K])
  T_star[select_index] <- T_star[select_index] + (s_all[K+1] - s_all[select_index])
  
  ## Conditional simulation for the dependent components given T_star
  mc_s <- rcmvnorm(mc_num, mu_T, cov_T,
                   given.ind=1:K,
                   dependent.ind=(K+1):(K+end_stage-1),
                   X.given=T_star)
  
  ## Recover stage-1 score vector s1 for all arms
  s1 <- matrix(NA, nrow=K, ncol=mc_num)
  for(i in 1:mc_num){
    T_star_s <- T_star
    T_star_s[select_index] <- T_star[select_index] - sum(mc_s[i,])
    s1[,i] <- solve(I(1:K,1) * solve(cov[1:K,1:K])) %*% T_star_s
  }
  
  ## Screen acceptable MC draws:
  ## - continuation region up to end_stage-1
  ## - selection at stage 1
  z_stage1 <- apply(s1, 2, function(x) x / sqrt(I(1:K,1)))
  mc_s_all <- cbind(s1[select_index,], mc_s)
  
  mc_z <- t(apply(mc_s_all, 1, function(x) cumsum(x) / sqrt(I(select_index, 1:end_stage))))
  
  bround_L <- matrix(low_bround[1:(end_stage-1)], mc_num, end_stage-1, byrow=TRUE)
  bround_U <- matrix(up_bround[1:(end_stage-1)],  mc_num, end_stage-1, byrow=TRUE)
  
  index <- apply(cbind(
    apply(mc_z[,1:(end_stage-1)] >= bround_L, 1, function(x) identical(x, rep(TRUE, end_stage-1))),
    apply(mc_z[,1:(end_stage-1)] <= bround_U, 1, function(x) identical(x, rep(TRUE, end_stage-1))),
    apply(z_stage1, 2, which.is.max) == select_index
  ), 1, all)
  
  ## RB estimator and CI from accepted draws
  re <- rep(NA, 4)
  if(sum(index) > 1){
    mc_est <- mean(mc_s_all[index, 2] / (I(select_index,2) - I(select_index,1)))
    SE <- sqrt(1/(I(select_index,2) - I(select_index,1)) -
                 var(mc_s_all[index,2] / (I(select_index,2) - I(select_index,1))))
    re <- c(sum(index)/mc_num,
            mc_est,
            mc_est + qnorm(alpha/2)*SE,
            mc_est - qnorm(alpha/2)*SE)
  }
  re
}

###############################################################################
## 3) Parallel simulation driver
##    - Generate 'sim' replicates
##    - For those with end_stage >= 2 (i.e., phase III occurs), compute estimators
##    - Save outputs to CSV
###############################################################################

# setwd("C:/Users/13379/Desktop/essay/submit-SIM/binary-logHR/5")

for(setting in c("NULL","H1","peak","line")){
  
  J <- 5
  sim <- 10
  
  time1 <- Sys.time()
  
  ## Parallel cluster
  cls <- makeSOCKcluster(18)
  registerDoSNOW(cls)
  
  ## Progress bar for foreach
  pb <- txtProgressBar(max=sim, style=3, char="*")
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  ## --- 3.1 Generate trial-level summaries for all replicates ---
  res <- foreach(i=1:sim,
                 .options.snow=opts,
                 .combine=cbind,
                 .packages=c("data.table","compiler","nnet")) %dopar% {
                   da_genc(i, setting, J)
                 }
  
  ## Parse result blocks (based on da_gen return order)
  end_stage    <- res[1,]
  select_index <- res[2,]
  naive_P      <- res[3:7,]    # 5 rows: control + 4 arms
  naive        <- res[8:11,]   # 4 rows: arms only (theta-like estimates)
  x_bar        <- res[12:16,]  # 5 rows: cumulative counts summaries
  s_all        <- res[17:21,]  # 5 rows: score summaries for RB
  bias         <- res[22:24,]  # true theta_selected, stage2_est, naive_selected (your chosen)
  idx <- which(end_stage >= 2)
  
  ## Save some baseline outputs for trials that reached stage >=2
  write.csv(cbind(end_stage, select_index, t(bias))[idx,],
            paste(setting, "stage2-naive.csv", sep="-"),
            row.names=FALSE)
  
  ###############################################################################
  ## 3.2 Run each estimator for eligible replicates and save to CSV
  ###############################################################################
  
  re_SI <- foreach(i=idx,
                   .options.snow=opts,
                   .combine=cbind,
                   .packages=c("tmvtnorm","mvtnorm","compiler")) %dopar% {
                     SI(naive[,i], select_index[i],
                        K=4, J=J,
                        p_control = naive_P[1,i],
                        p_treatment = naive_P[-1,i])
                   }
  write.csv(cbind(end_stage[idx], select_index[idx], as.numeric(re_SI)),
            paste(setting, "SI.csv", sep="-"),
            row.names=FALSE)
  
  re_MI <- foreach(i=idx,
                   .options.snow=opts,
                   .combine=cbind,
                   .packages=c("tmvtnorm","mvtnorm","compiler")) %dopar% {
                     MIc(naive[,i], select_index[i],
                         max_iterations=20, tol=0.001,
                         K=4, J=J,
                         p_control = naive_P[1,i],
                         p_treatment = naive_P[-1,i])
                   }
  write.csv(cbind(end_stage[idx], select_index[idx], t(re_MI)),
            paste(setting, "MI.csv", sep="-"),
            row.names=FALSE)

  re_MUE_MLE <- foreach(i=idx,
                        .options.snow=opts,
                        .combine=cbind,
                        .packages=c("mvtnorm","compiler")) %dopar% {
                          est_MUE(naive[,i], end_stage[i], select_index[i],
                                  sub="MLE",
                                  K=4, J=J,
                                  p_control = naive_P[1,i],
                                  p_treatment = naive_P[-1,i])
                        }
  write.csv(cbind(end_stage[idx], select_index[idx], t(re_MUE_MLE)),
            paste(setting, "MUE_MLE.csv", sep="-"),
            row.names=FALSE)
  
  re_MUE_ZERO <- foreach(i=idx,
                         .options.snow=opts,
                         .combine=cbind,
                         .packages=c("mvtnorm","compiler")) %dopar% {
                           est_MUE(naive[,i], end_stage[i], select_index[i],
                                   sub="ZERO",
                                   K=4, J=J,
                                   p_control = naive_P[1,i],
                                   p_treatment = naive_P[-1,i])
                         }
  write.csv(cbind(end_stage[idx], select_index[idx], t(re_MUE_ZERO)),
            paste(setting, "MUE_ZERO.csv", sep="-"),
            row.names=FALSE)
  
  re_RB <- foreach(i=idx,
                   .options.snow=opts,
                   .combine=cbind,
                   .packages=c("condMVNorm","nnet","compiler")) %dopar% {
                     RB_mc(s_all[,i],
                           p_control   = naive_P[1,i],
                           p_treatment = naive_P[-1,i],
                           select_index = select_index[i],
                           end_stage    = end_stage[i],
                           mc_num=10000, alpha=0.05,
                           K=4, J=J)
                   }
  write.csv(cbind(end_stage[idx], select_index[idx], t(re_RB)),
            paste(setting, "RB.csv", sep="-"),
            row.names=FALSE)
  
  ## Cleanup
  close(pb)
  stopCluster(cls)
  
  time2 <- Sys.time()
  print(time2 - time1)
}

###############################################################################
## End of script
###############################################################################
