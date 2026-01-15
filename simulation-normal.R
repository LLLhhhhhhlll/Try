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
library(dplyr)

###############################################################################
## 1) Data generator for one simulation replicate
##    - Generate patient-level normal outcomes for control + K experimental arms
##    - Select best arm based on Z-statistic at end of phase II (stage 1)
##    - Apply group-sequential boundaries in phase III (stages 2..J)
##    - Output sufficient components used by subsequent estimators
###############################################################################
da_gen <- function(sim, setting, J){

  ## Fix seed for reproducibility across replicates
  set.seed(20250420 + sim)

  ## --- 1.1 Define treatment effect sizes under different scenarios ---
  if(setting == "NULL"){
    mu_treatment <- c(0, 0, 0, 0)
  }else if(setting == "H1"){
    mu_treatment <- c(0.50, 0.50, 0.50, 0.50)
  }else if(setting == "peak"){
    mu_treatment <- c(0.25, 0.25, 0.25, 0.50)
  }else if(setting == "line"){
    mu_treatment <- c(0.20, 0.30, 0.40, 0.50)
  }

  ## --- 1.2 Specify design parameters depending on number of stages J ---
  ## t: information fraction; low_bround/up_bround: group-seq boundaries;
  ## unit_n: per-stage unit sample size scaling
  if(J == 3){
    t <- c(1,2,3)/3
    low_bround <- c(0.000000, 1.038459, 2.290364)
    up_bround  <- c(Inf, 3.301585, 2.290364)
    unit_n <- 30
  }else if(J == 4){
    t <- c(1,2,3,4)/4
    low_bround <- c(0.0000000, 0.1267247, 1.4593383, 2.251983)
    up_bround  <- c(Inf, 4.013738, 2.826973, 2.251983)
    unit_n <- 25
  }else if(J == 5){
    t <- c(1,2,3,4,5)/5
    low_bround <- c(0.0000000, -0.6733603, 0.7312345, 1.6214603, 2.228605)
    up_bround  <- c(Inf, 4.609266, 3.263919, 2.639797, 2.228605)
    unit_n <- 22
  }

  ## Basic outcome distribution (normal endpoints)
  mu_control <- 0
  sigma_control <- 1
  sigma_treatment <- c(1,1,1,1)

  ## K experimental arms; redefine J from boundary length for safety
  K <- length(mu_treatment)
  J <- length(low_bround)

  ## Cumulative sample size per stage for the control arm
  ## (and assumed same for each treatment arm)
  n <- t * length(t) * unit_n

  ## Fisher information for treatment-control difference at stage j
  I <- function(j,index){
    1 / (sigma_treatment[index]^2/n[j] + sigma_control^2/n[j])
  }

  ## Wald Z-statistic based on mean difference at stage j
  Z <- function(da_treat, da_control, j, index){
    (mean(da_treat) - mean(da_control)) * sqrt(I(j, index))
  }

  ## --- 1.3 Generate patient-level data ---
  ## all_da: n[J] x (K+1) matrix; col1=control, col2..K+1=treatments
  all_da <- mvrnorm(
    n[J],
    mu = c(mu_control, mu_treatment),
    Sigma = diag(c(sigma_control, sigma_treatment)^2, nrow = K+1, ncol = K+1)
  )

  ## --- 1.4 Select best arm at end of phase II (stage 1) ---
  select_index <- which.max(
    mapply(function(x){
      Z(all_da[1:n[1], x+1], all_da[1:n[1], 1], j=1, index=x)
    }, x=1:K)
  )

  ## --- 1.5 Determine stopping stage in phase III using boundaries ---
  ## inter_decide returns:
  ##   "f" = futility stop, "e" = efficacy stop, "c" = continue
  inter_decide <- function(j, select_index, da){
    z_temp <- Z(da[1:n[j], select_index+1], da[1:n[j], 1], j, index=select_index)
    if(z_temp < low_bround[j]) "f" else if(z_temp > up_bround[j]) "e" else "c"
  }

  inter_outcome <- character(J)
  for(i in 1:J){
    inter_outcome[i] <- inter_decide(j=i, select_index=select_index, da=all_da)
  }

  ## First stage where decision is made ("e" or "f")
  end_stage <- data.table::first(which(inter_outcome %in% c("e","f")))

  ## --- 1.6 Construct summary quantities used by estimators ---
  ## Stage-2 incremental estimator (using data between n1+1 .. n2)
  stage2_est <- mean(all_da[(n[1]+1):n[2], select_index+1]) -
                mean(all_da[(n[1]+1):n[2], 1])

  ## Naive MLE for each arm based on stage 1, except selected arm uses cumulative up to end_stage
  naive <- mapply(function(x){
    mean(all_da[1:n[1], x+1]) - mean(all_da[1:n[1], 1])
  }, 1:K)
  naive[select_index] <- mean(all_da[1:n[end_stage], select_index+1]) -
                         mean(all_da[1:n[end_stage], 1])

  ## x_bar: sample means used in UMVCUE (control and selected arm updated to end_stage)
  x_bar <- colMeans(all_da[1:n[1], ])
  x_bar[1] <- mean(all_da[1:n[end_stage], 1])
  x_bar[select_index+1] <- mean(all_da[1:n[end_stage], select_index+1])

  ## Score-like quantities S(i,j) for RB estimator input
  S <- function(i,j){
    (mean(all_da[1:n[j], i+1]) - mean(all_da[1:n[j], 1])) * I(j, i)
  }
  s_all <- c(mapply(S, i=1:K, 1), S(select_index, end_stage))

  ## Return a vector that will be combined across simulation replicates:
  c(end_stage, select_index, naive, x_bar, s_all, sigma_hat, stage2_est, naive[select_index])
}
da_genc <- cmpfun(da_gen)   # compile for speed

###############################################################################
## 2) Estimator implementations (CBAE-SI/MI, CMUE, UMVCUE, RB)
###############################################################################

## --- 2.1 Mean vector and covariance matrix for CBAE/CMUE under MVN framework ---
mu_cov_f <- function(control_eff, treat_eff, select_index,
                     control_sig, treat_sig, K, J, n){

  V <- function(sig,j) sig^2 / n[j]
  I <- function(i,j)  1 / (V(treat_sig[i], j) + V(control_sig, j))

  ## A is the linear transform used to define the selection event and sequential statistics
  A <- diag(1, K+J-1)
  diag(A[1:(K-1), 1:(K-1)]) <- -1
  A[1:(K-1), K] <- 1

  ## theta_v rescales components by sqrt(information)
  theta_v <- diag(nrow = K+J-1)
  diag(theta_v[1:(K-1), 1:(K-1)]) <- sqrt(I((1:K)[-select_index], 1))
  diag(theta_v[K:(K+J-1), K:(K+J-1)]) <- sqrt(I(select_index, 1:J))

  ## Build covariance before applying A
  cov <- cbind(
    matrix(V(control_sig, 1), nrow = K+J-1, ncol = K-1),
    mapply(function(x) matrix(V(control_sig, x), nrow = K+J-1, ncol = 1), x = 1:J)
  )

  ## Add treatment variances
  diag(cov[1:(K-1), 1:(K-1)]) <- diag(cov[1:(K-1), 1:(K-1)]) + V(treat_sig[-select_index], 1)

  for(i in K:(K+J-1)){
    for(j in K:(K+J-1)){
      cov[i,j] <- cov[i,j] + V(treat_sig[select_index], j-(K-1))
    }
  }

  ## Rescale by sqrt(I) to match Z-type statistics
  for(i in 1:(K+J-1)){
    for(j in 1:(K+J-1)){
      if(i <= K-1 && j <= K-1){
        cov[i,j] <- cov[i,j] * sqrt(I((1:K)[-select_index][i],1) * I((1:K)[-select_index][j],1))
      }else if(i <= K-1 && (j >= K && j <= K+J-1)){
        cov[i,j] <- cov[i,j] * sqrt(I((1:K)[-select_index][i],1) * I(select_index, j-(K-1)))
      }else if((i >= K && i <= K+J-1) && (j >= K && j <= K+J-1)){
        cov[i,j] <- cov[i,j] * sqrt(I(select_index, i-(K-1)) * I(select_index, j-(K-1)))
      }
    }
  }

  cov[lower.tri(cov)] <- t(cov)[lower.tri(cov)]

  ## Mean vector under true treatment effects
  mu <- A %*% theta_v %*% as.matrix(
    c(treat_eff[-select_index] - control_eff,
      rep(treat_eff[select_index] - control_eff, J))
  )

  ## Apply linear transform A
  cov <- A %*% cov %*% t(A)

  ## Numerical safety
  if(any(diag(cov) < 0)){
    diag(cov)[diag(cov) < 0] <- runif(sum(diag(cov) < 0))
  }

  list(mu = mu, cov = cov)
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
bias_theta <- function(select_index, theta, sigma_control, sigma_treatment, K, J){

  ## --- Define design parameters depending on J ---
  if(J == 3){
    t <- c(1,2,3)/3
    low_bround <- c(0.000000, 1.038459, 2.290364)
    up_bround  <- c(Inf, 3.301585, 2.290364)
    unit_n <- 30
  }else if(J == 4){
    t <- c(1,2,3,4)/4
    low_bround <- c(0.0000000, 0.1267247, 1.4593383, 2.3056892)
    up_bround  <- c(Inf, 4.013738, 2.826973, 2.251983)
    unit_n <- 25
  }else if(J == 5){
    t <- c(1,2,3,4,5)/5
    low_bround <- c(0.0000000, -0.6733603, 0.7312345, 1.6214603, 2.228605)
    up_bround  <- c(Inf, 4.609266, 3.263919, 2.639797, 2.228605)
    unit_n <- 22
  }

  ## Cumulative sample sizes per stage
  n <- t * length(t) * unit_n

  ## Fisher information for treatment-control difference at stage j
  I <- function(j, index){
    1 / (sigma_treatment[index]^2/n[j] + sigma_control^2/n[j])
  }

  ## --- Build mean vector and covariance matrix for the score/Z-statistics ---
  mu_cov <- mu_cov_f(control_eff = 0, treat_eff = theta,
                     control_sig = sigma_control, treat_sig = sigma_treatment,
                     select_index = select_index, K = K, J = J, n = n)
  mu  <- mu_cov$mu
  cov <- mu_cov$cov

  ## --- Compute prob(Q): probability of selection/continuation event ---
  ## Q = { selected arm is best at stage 1 } ∩ { continuation after stage 1 }
  ## This uses MVN probability over first K components.
  prob_Q <- pmvnorm(
    lower = c(rep(0, K-1), low_bround[1]),
    upper = c(rep(Inf, K-1), up_bround[1]),
    mean  = mu[1:K,],
    sigma = cov[1:K, 1:K]
  )[1]

  ## Initialize bias vector (length K arms)
  bias <- rep(NA, K)

  if(!is.na(prob_Q)){

    ## --- Compute truncated mean E[ S | Q ] over first K components ---
    ## mtm is the conditional expectation of the vector under truncation region Q
    mtm <- mtmvnorm(
      mean  = mu[1:K,],
      sigma = cov[1:K, 1:K],
      lower = c(rep(0, K-1), low_bround[1]),
      upper = c(rep(Inf, K-1), up_bround[1]),
      doComputeVariance = FALSE
    )$tmean

    ## --- Bias for NON-selected arms ---
    ## This corresponds to the conditional bias adjustment for arms not selected
    ## (matches your main manuscript formula for non-selected arm bias).
    bias[-select_index] <-
      (mtm[K] - mtm[-K]) / sqrt(I(1, (1:K)[-select_index])) - theta[-select_index]

    ## --- Bias for SELECTED arm ---
    ## Need integrate over possible stopping stage (2..J) in phase III.
    ## E_select(stage) computes contribution at a given stage.
    E_select <- function(stage){

      ## Use first (K + stage - 1) dimensions for event up to that stage
      mu_t  <- mu[1:(K + stage - 1),]
      cov_t <- cov[1:(K + stage - 1), 1:(K + stage - 1)]

      ## Split event at stage into:
      ##  (1) cross lower boundary at 'stage'  OR
      ##  (2) cross upper boundary at 'stage'
      ## while staying within boundaries at earlier stages.
      lower1 <- c(rep(0, K-1), low_bround[1:(stage-1)], -Inf)
      upper1 <- c(rep(Inf, K-1), up_bround[1:(stage-1)], low_bround[stage])

      lower2 <- c(rep(0, K-1), low_bround[1:(stage-1)], up_bround[stage])
      upper2 <- c(rep(Inf, K-1), up_bround[1:(stage-1)], Inf)

      ## Weighted truncated mean at that stage (dimension K+stage-1)
      ## Divide by sqrt(I) to map back to theta scale.
      (
        pmvnorm(lower = lower1, upper = upper1, mean = mu_t, sigma = cov_t)[1] *
          mtmvnorm(mean = mu_t, sigma = cov_t, lower = lower1, upper = upper1,
                   doComputeVariance = FALSE)$tmean[K + stage - 1] +
        pmvnorm(lower = lower2, upper = upper2, mean = mu_t, sigma = cov_t)[1] *
          mtmvnorm(mean = mu_t, sigma = cov_t, lower = lower2, upper = upper2,
                   doComputeVariance = FALSE)$tmean[K + stage - 1]
      ) / sqrt(I(stage, select_index))
    }

    ## Aggregate over stages (2..J) and normalize by prob(Q)
    bias[select_index] <- sum(mapply(E_select, stage = 2:J)) / prob_Q - theta[select_index]
  }

  bias
}

## (A2) MI: iterative bias correction
##      theta_{new} = initial - bias_theta(theta_old) until convergence.
MI <- function(initial, select_index,
               max_iterations = 20, tol = 0.001,
               sigma_control, sigma_treatment, K, J){

  conver <- 0
  theta_hat <- initial

  for(i in 1:max_iterations){

    ## Compute conditional bias at current iterate
    bias_temp <- bias_theta(select_index, theta = theta_hat,
                            sigma_control, sigma_treatment, K, J)

    ## If numerical failure, stop and mark as not converged
    if(is.na(sum(bias_temp)) || is.infinite(sum(bias_temp))){
      theta_hat <- initial
      conver <- -1
      break()
    }

    ## One-step update
    theta_hat_new <- initial - bias_temp

    ## Euclidean convergence criterion
    euc.dis <- sqrt(sum((theta_hat - theta_hat_new)^2))

    if(euc.dis <= tol){
      theta_hat <- theta_hat_new
      conver <- 1
      break()
    }else{
      theta_hat <- theta_hat_new
    }
  }

  ## Return: (convergence flag, corrected estimate for selected arm)
  c(conver, theta_hat[select_index])
}

## (A3) SI: single-step bias correction
SI <- function(initial, select_index, sigma_control, sigma_treatment, K, J){
  (initial - bias_theta(select_index, initial, sigma_control, sigma_treatment, K, J))[select_index]
}

MIc <- cmpfun(MI)   # compile MI for speed


##### -------------------------------------------------------------------------
##### (B) CMUE: conditional median unbiased estimator
#####     - compute conditional p-value under candidate theta
#####     - invert p-value to get CMUE point estimate (median) and CI endpoints
##### -------------------------------------------------------------------------

## est_MUE(): returns CMUE point estimate (median) by solving p(theta)=0.5
est_MUE <- function(MLE, end_stage, select_index, sub,
                    sigma_control, sigma_treatment, K, J){

  ## --- design parameters depend on J ---
  if(J == 3){
    t <- c(1,2,3)/3
    low_bround <- c(0.000000, 1.038459, 2.290364)
    up_bround  <- c(Inf, 3.301585, 2.290364)
    unit_n <- 30
  }else if(J == 4){
    t <- c(1,2,3,4)/4
    low_bround <- c(0.0000000, 0.1267247, 1.4593383, 2.3056892)
    up_bround  <- c(Inf, 4.013738, 2.826973, 2.251983)
    unit_n <- 25
  }else if(J == 5){
    t <- c(1,2,3,4,5)/5
    low_bround <- c(0.0000000, -0.6733603, 0.7312345, 1.6214603, 2.228605)
    up_bround  <- c(Inf, 4.609266, 3.263919, 2.639797, 2.228605)
    unit_n <- 22
  }

  n <- t * length(t) * unit_n

  ## Fisher information and observed Z at stopping stage
  I <- function(j, index){
    1 / (sigma_treatment[index]^2/n[j] + sigma_control^2/n[j])
  }
  Z_end <- MLE[select_index] * sqrt(I(end_stage, select_index))

  ## cond_pvalue(): conditional p-value under candidate theta, given Q
  ## where Q is the selection event (best arm at stage 1 and continuation).
  cond_pvalue <- function(theta, select_index, end_stage, Z_end){

    ## Build MVN mean/cov under theta
    mu_cov <- mu_cov_f(control_eff = 0, treat_eff = theta,
                       select_index = select_index,
                       control_sig  = sigma_control,
                       treat_sig    = sigma_treatment,
                       K = K, J = J, n = n)
    mu  <- mu_cov$mu
    cov <- mu_cov$cov

    ## P(Q): selection/continuation probability at stage 1
    prob_Q <- pmvnorm(
      lower = c(rep(0, K-1), low_bround[1]),
      upper = c(rep(Inf, K-1), up_bround[1]),
      mean  = mu[1:K,],
      sigma = cov[1:K, 1:K]
    )[1]

    ## Numerator: probability of crossing efficacy boundary at end_stage with Z >= Z_end
    ## while respecting boundaries at earlier stages (continue region)
    P_stage <- function(stage, stop, Z_end){
      new_mu  <- mu[1:(K + stage - 1),]
      new_cov <- cov[1:(K + stage - 1), 1:(K + stage - 1)]

      if(stop){
        ## crossing upper boundary at 'stage' with Z >= Z_end
        P <- pmvnorm(
          lower = c(rep(0, K-1), low_bround[1:(stage-1)], Z_end),
          upper = c(rep(Inf, K-1), up_bround[1:(stage-1)], Inf),
          mean  = new_mu,
          sigma = new_cov
        )[1]
      }else{
        ## crossing upper boundary at earlier stages (used to build up numerator)
        P <- pmvnorm(
          lower = c(rep(0, K-1), low_bround[1:(stage-1)], up_bround[stage]),
          upper = c(rep(Inf, K-1), up_bround[1:(stage-1)], Inf),
          mean  = new_mu,
          sigma = new_cov
        )[1]
      }

      if(is.na(P)) 0 else P
    }

    ## Build numerator depending on end_stage
    if(end_stage == 2){
      numerator <- P_stage(stage = end_stage, stop = TRUE, Z_end = Z_end)
    }else{
      numerator <- sum(
        mapply(P_stage, stage = 2:(end_stage-1), stop = FALSE, Z_end = Z_end),
        P_stage(stage = end_stage, stop = TRUE, Z_end = Z_end)
      )
    }

    ## Conditional p-value = numerator / prob_Q
    if(is.na(numerator) || is.na(prob_Q) || prob_Q <= 0){
      return(1)
    }
    min(1, numerator / prob_Q)
  }

  cond_pvaluec <- cmpfun(cond_pvalue)

  ## Plug-in strategy for non-selected arms:
  ## sub="MLE": fix non-selected arms at their MLEs
  ## sub="ZERO": fix non-selected arms at 0
  P_value_sub <- function(select_index, sub, x){
    theta <- rep(x, K)
    if(sub == "MLE"){
      theta[-select_index] <- MLE[-select_index]
    }else{
      theta[-select_index] <- rep(0, K-1)
    }
    cond_pvaluec(theta, select_index, end_stage, Z_end)
  }

  ## Solve for CMUE median: p(theta) = 0.5
  ## NOTE: bracket range (-1,1) may need to widen for larger effects.
  a <- optimise(function(x) (P_value_sub(select_index, sub, x) - 1/2)^2,
                interval = c(-1, 1), tol = 0.001)

  ## Return: (median estimate, optimization objective)
  c(a$minimum, a$objective)
}


##### -------------------------------------------------------------------------
##### (C) RB: Rao–Blackwell estimator via conditional simulation
#####     - condition an unbiased stage-wise estimator on sufficient statistic
#####     - Monte Carlo implementation of conditional expectation
##### -------------------------------------------------------------------------

RB_mc <- function(s_all, select_index, end_stage, mc_num, alpha,
                  sigma_control, sigma_treatment, K, J){

  ## --- design parameters depend on J ---
  if(J == 3){
    t <- c(1,2,3)/3
    low_bround <- c(0.000000, 1.038459, 2.290364)
    up_bround  <- c(Inf, 3.301585, 2.290364)
    unit_n <- 30
  }else if(J == 4){
    t <- c(1,2,3,4)/4
    low_bround <- c(0.0000000, 0.1267247, 1.4593383, 2.3056892)
    up_bround  <- c(Inf, 4.013738, 2.826973, 2.251983)
    unit_n <- 25
  }else if(J == 5){
    t <- c(1,2,3,4,5)/5
    low_bround <- c(0.0000000, -0.6733603, 0.7312345, 1.6214603, 2.228605)
    up_bround  <- c(Inf, 4.609266, 3.263919, 2.639797, 2.228605)
    unit_n <- 22
  }

  n <- t * length(t) * unit_n

  ## Working model used to derive conditional distribution (theta=0)
  theta <- rep(0, K)

  ## Information terms for score statistics
  I <- function(i,j) n[j] / (sigma_control^2 + sigma_treatment[i]^2)

  ## Cov between treatment-control scores due to shared control
  I_C <- function(i, ii, j){
    n[j] * sigma_control^2 / (sigma_control^2 + sigma_treatment[i]^2) /
      (sigma_control^2 + sigma_treatment[ii]^2)
  }

  ## Mean vector of score statistics up to end_stage
  mu <- c(theta * I(1:K, 1),
          rep(theta[select_index] * I(select_index, 1), end_stage-1))

  ## Covariance matrix of score statistics
  cov <- matrix(0, K + end_stage - 1, K + end_stage - 1)

  for(i in 1:K){
    for(ii in 1:K){
      cov[i,ii] <- I_C(i, ii, 1)
    }
  }
  diag(cov[1:K, 1:K]) <- I(1:K, 1)

  if(end_stage > 2){
    diag(cov[(K+1):(K+end_stage-1), (K+1):(K+end_stage-1)]) <- I(select_index, 1)
  }else{
    cov[K+end_stage-1, K+end_stage-1] <- I(select_index, 1)
  }
  cov[lower.tri(cov)] <- t(cov)[lower.tri(cov)]

  ## Construct transformation to sufficient statistic T (as in manuscript)
  trans_T <- diag(1, K + end_stage - 1)
  trans_T[select_index,] <- 1
  trans_T[1:K, 1:K] <- solve(cov[1:K, 1:K]) * I(1:K, 1)

  mu_T  <- trans_T %*% mu
  cov_T <- trans_T %*% cov %*% t(trans_T)

  ## Observed sufficient statistic T_star computed from observed s_all
  T_star <- (I(1:K,1) * solve(cov[1:K,1:K])) %*% as.matrix(s_all[1:K])
  T_star[select_index] <- T_star[select_index] + (s_all[K+1] - s_all[select_index])

  ## Step 1: conditional simulation of remaining components given T_star
  mc_s <- rcmvnorm(mc_num, mu_T, cov_T,
                   given.ind = 1:K,
                   dependent.ind = (K+1):(K+end_stage-1),
                   X.given = T_star)

  ## Recover stage-wise scores s1 (stage 1 scores for all arms) from simulated T’s
  s1 <- matrix(NA, nrow = K, ncol = mc_num)
  for(i in 1:mc_num){
    T_star_s <- T_star
    T_star_s[select_index] <- T_star[select_index] - sum(mc_s[i,])
    s1[,i] <- solve(I(1:K,1) * solve(cov[1:K,1:K])) %*% T_star_s
  }

  ## Step 2: accept only samples consistent with selection event and boundaries
  z_stage1 <- apply(s1, 2, function(x) x / sqrt(I(1:K, 1)))

  ## Assemble full stage-wise scores for selected arm (stage 1 + increments)
  mc_s_all <- cbind(s1[select_index,], mc_s)

  ## Convert to cumulative Z’s for the selected arm
  mc_z <- t(apply(mc_s_all, 1, function(x) cumsum(x) / sqrt(I(select_index, 1:end_stage))))

  bround_L <- matrix(low_bround[1:(end_stage-1)], mc_num, end_stage-1, byrow = TRUE)
  bround_U <- matrix(up_bround[1:(end_stage-1)],  mc_num, end_stage-1, byrow = TRUE)

  ## Accept if:
  ##   (i) selected-arm Z stays within boundaries up to end_stage-1
  ##   (ii) stage-1 selection chooses select_index as the max Z
  index <- apply(
    cbind(
      apply(mc_z[,1:(end_stage-1)] >= bround_L, 1, function(x) identical(x, rep(TRUE, end_stage-1))),
      apply(mc_z[,1:(end_stage-1)] <= bround_U, 1, function(x) identical(x, rep(TRUE, end_stage-1))),
      apply(z_stage1, 2, which.max) == select_index
    ),
    1, all
  )

  ## Step 3: RB estimate = conditional expectation of stage2 unbiased increment
  re <- rep(NA, 4)
  if(sum(index) > 1){

    ## stage2 unbiased increment estimator based on simulated increments
    mc_est <- mean(mc_s_all[index, 2] / (I(select_index,2) - I(select_index,1)))

    ## SE uses law of total variance: Var(U) - Var(E[U|T])
    SE <- sqrt(1/(I(select_index,2) - I(select_index,1)) -
                 var(mc_s_all[index,2] / (I(select_index,2) - I(select_index,1))))

    re <- c(sum(index)/mc_num,
            mc_est,
            mc_est + qnorm(alpha/2) * SE,
            mc_est - qnorm(alpha/2) * SE)
  }

  re
}

###############################################################################
## 3) Parallel simulation driver
##    - Generate 'sim' replicates
##    - For those with end_stage >= 2 (i.e., phase III occurs), compute estimators
##    - Save outputs to CSV
###############################################################################

setting <- "peak"
sim <- 10
J <- 4

setwd("C:/Users/13379/Desktop/essay/submit-SIM/res/4")

time1 <- Sys.time()

## Create PSOCK cluster and register for foreach
cls <- makeSOCKcluster(2)
registerDoSNOW(cls)

## Progress bar in parallel foreach
pb <- txtProgressBar(max = sim, style = 3, char = "*")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## --- 3.1 Generate trial-level summaries for all replicates ---
res <- foreach(i = 1:sim,
               .options.snow = opts,
               .combine = cbind,
               .packages = c("MASS","data.table","compiler","condMVNorm","pbapply",
                             "tmvtnorm","BSDA","pracma")) %dopar% {
  da_genc(i, setting, J)
}

## Unpack res matrix into components
end_stage    <- res[1,]
select_index <- res[2,]
naive        <- res[3:6,]
x_bar        <- res[7:11,]
s_all        <- res[12:16,]

stage2_bias_related <- res[22:23,] 

## Only analyze replicates where trial reaches phase III interim(s)
idx <- which(end_stage >= 2)

###############################################################################
## 3.2 Run each estimator for eligible replicates and save to CSV
###############################################################################

## --- CBAE-SI ---
re_SI <- foreach(i = idx,
                 .options.snow = opts,
                 .combine = cbind,
                 .packages = c("MASS","data.table","compiler","condMVNorm","pbapply",
                               "tmvtnorm","BSDA","pracma")) %dopar% {
  SI(naive[,i], select_index[i], sigma_hat[1,i], sigma_hat[-1,i], K=4, J=J)
}
write.csv(cbind(end_stage[idx], select_index[idx], as.numeric(re_SI)),
          paste(setting, "SI.csv", sep = "-"))

## --- CBAE-MI ---
re_MI <- foreach(i = idx,
                 .options.snow = opts,
                 .combine = cbind,
                 .packages = c("MASS","data.table","compiler","condMVNorm","pbapply",
                               "tmvtnorm","BSDA","pracma")) %dopar% {
  MIc(naive[,i], select_index[i],
      max_iterations = 20, tol = 0.001,
      sigma_hat[1,i], sigma_hat[-1,i], K=4, J=J)
}
write.csv(cbind(end_stage[idx], select_index[idx], t(re_MI)),
          paste(setting, "MI.csv", sep = "-"))

## --- CMUE-MLE ---
re_MUE_MLE <- foreach(i = idx,
                      .options.snow = opts,
                      .combine = cbind,
                      .packages = c("MASS","data.table","compiler","condMVNorm","pbapply",
                                    "tmvtnorm","BSDA","pracma")) %dopar% {
  est_MUE(naive[,i], end_stage[i], select_index[i], sub="MLE",
          sigma_hat[1,i], sigma_hat[-1,i], K=4, J=J)
}
write.csv(cbind(end_stage[idx], select_index[idx], t(re_MUE_MLE)),
          paste(setting, "MUE_MLE.csv", sep = "-"))

## --- CMUE-ZERO ---
re_MUE_ZERO <- foreach(i = idx,
                       .options.snow = opts,
                       .combine = cbind,
                       .packages = c("MASS","data.table","compiler","condMVNorm","pbapply",
                                     "tmvtnorm","BSDA","pracma")) %dopar% {
  est_MUE(naive[,i], end_stage[i], select_index[i], sub="ZERO",
          sigma_hat[1,i], sigma_hat[-1,i], K=4, J=J)
}
write.csv(cbind(end_stage[idx], select_index[idx], t(re_MUE_ZERO)),
          paste(setting, "MUE_ZERO.csv", sep = "-"))

## --- UMVCUE ---
re_UMVCUE <- foreach(i = idx,
                     .options.snow = opts,
                     .combine = cbind,
                     .packages = c("MASS","data.table","compiler","condMVNorm","pbapply",
                                   "tmvtnorm","BSDA","pracma")) %dopar% {
  UMVCUE_mc(x_bar[,i], select_index[i], end_stage[i],
            mc_n = 10000, alpha = 0.05,
            sigma_hat[1,i], sigma_hat[-1,i], K=4, J=J)
}
write.csv(cbind(end_stage[idx], select_index[idx], t(re_UMVCUE)),
          paste(setting, "UMVCUE.csv", sep = "-"))

## --- RB ---
re_RB <- foreach(i = idx,
                 .options.snow = opts,
                 .combine = cbind,
                 .packages = c("MASS","data.table","compiler","condMVNorm","pbapply",
                               "tmvtnorm","BSDA","pracma")) %dopar% {
  RB_mc(s_all[,i], select_index[i], end_stage[i],
        mc_num = 10000, alpha = 0.05,
        sigma_hat[1,i], sigma_hat[-1,i], K=4, J=J)
}
write.csv(cbind(end_stage[idx], select_index[idx], t(re_RB)),
          paste(setting, "RB.csv", sep = "-"))

###############################################################################
## 4) Clean up parallel backend and report elapsed time
###############################################################################
close(pb)
stopCluster(cls)

time2 <- Sys.time()
print(time2 - time1)
