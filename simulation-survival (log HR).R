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
## (1) Mean/Covariance under MVN approximation for Z-vector (survival/cox score)
###############################################################################
mu_cov_f <- function(lambda_C, lambda_E, shape, select_index, theta,
                     D_n, K, J){
  
  ## theta default: log hazard ratio (treatment vs control) for each arm
  if(is.null(theta)){
    theta <- log(lambda_E / lambda_C)
  }
  
  ## I(i,j): information approximation at stage j for arm i
  ## Here you use event counts confirmed by D_n:
  ##   D_n[j,1] = cumulative events in control by stage j
  ##   D_n[j,i+1] = cumulative events in arm i by stage j
  ## Factor 1/4 corresponds to approx Var(logHR_hat) ~ 4 / total events (rough scale)
  I <- function(i, j){
    1/4 * (D_n[j, 1] + D_n[j, i+1])
  }
  
  ## Build covariance for "S" vector:
  ## dimension (K + J - 1):
  ##   - first K are stage-1 score-like stats for each arm
  ##   - remaining J-1 are increments for selected arm across later stages
  cov <- cbind(
    matrix(1/4 * D_n[1,1], K+J-1, K,   byrow=TRUE),
    matrix(0,             K+J-1, J-1, byrow=TRUE)
  )
  
  ## Stage-1 variances for arms
  diag(cov[1:K, 1:K]) <- mapply(I, i=c((1:K)[-select_index], select_index), j=1)
  
  ## Increment variances for selected arm across looks: diff(I)
  diag(cov[(K+1):(K+J-1), (K+1):(K+J-1)]) <- diff(mapply(I, i=select_index, j=1:J))
  
  ## Symmetrize
  cov[lower.tri(cov)] <- t(cov)[lower.tri(cov)]
  
  ## Mean vector on "S" scale:
  ## - non-selected arms: theta_i * I(i,1)
  ## - selected arm: theta_sel * I(sel,1) plus increments
  mu <- c(
    c(theta[-select_index], theta[select_index]) *
      mapply(I, i=c((1:K)[-select_index], select_index), j=1),
    theta[select_index] * diff(mapply(I, i=select_index, j=1:J))
  )
  
  ## Transform increments to cumulative sums for selected arm across looks
  trans_s <- diag(1, K+J-1)
  trans_s[K:(K+J-1), K:(K+J-1)][lower.tri(trans_s[K:(K+J-1), K:(K+J-1)])] <- 1
  
  ## Standardize to Z scale using sqrt(I)
  theta_v <- matrix(0, K+J-1, K+J-1)
  diag(theta_v[1:(K-1), 1:(K-1)]) <- 1/sqrt(mapply(I, i=(1:K)[-select_index], j=1))
  diag(theta_v[K:(K+J-1), K:(K+J-1)]) <- 1/sqrt(mapply(I, i=select_index, j=1:J))
  
  mu_Z  <- theta_v %*% trans_s %*% mu
  cov_Z <- theta_v %*% (trans_s %*% cov %*% t(trans_s)) %*% t(theta_v)
  
  ## A: contrast matrix encoding selection event
  ## selection means (Z_sel - Z_other) <= 0 depending on sign conventions;
  ## here you use A with -1 on others and +1 on selected
  A <- diag(1, K+J-1)
  diag(A[1:(K-1), 1:(K-1)]) <- -1
  A[1:(K-1), K] <- 1
  
  list(mu = A %*% mu_Z,
       cov = A %*% cov_Z %*% t(A))
}


###############################################################################
## 1) Data generator for one simulation replicate
##    - Generate patient-level normal outcomes for control + K experimental arms
##    - Select best arm based on Z-statistic at end of phase II (stage 1)
##    - Apply group-sequential boundaries in phase III (stages 2..J)
##    - Output sufficient components used by subsequent estimators
###############################################################################
da_gen <- function(sim, setting, J){
  
  set.seed(20250425 + sim)
  
  ## Baseline hazard parameter (control)
  lambda_C <- log(2)
  
  ## Scenario: hazards for 4 arms (lambda_E < lambda_C implies benefit)
  if(setting=="NULL"){
    lambda_E <- rep(lambda_C, 4)
  }else if(setting=="H1"){
    lambda_E <- exp(rep(log(0.8), 4)) * lambda_C
  }else if(setting=="peak"){
    lambda_E <- exp(c(log(0.8)/2, log(0.8)/2, log(0.8)/2, log(0.8))) * lambda_C
  }else if(setting=="line"){
    lambda_E <- exp(c(log(0.8)/4, log(0.8)/2, log(0.8)*3/4, log(0.8))) * lambda_C
  }
  
  if(J==4){
    low_bround <- c(0.000000,  4.013691,  2.828157, -2.285393)
    up_bround  <- c(-Inf,     -4.013691, -2.828157, -2.285393)
    D_n <- c(625, 875, 1125, 1375)
  }else if(J==3){
    low_bround <- c(0.000000,  3.301594, -2.302597)
    up_bround  <- c(-Inf,     -3.301594, -2.302597)
    D_n <- c(765, 1071, 1377)
  }else if(J==5){
    low_bround <- c(0.000000,  4.609252,  3.264385,  2.644414, -2.278346)
    up_bround  <- c(-Inf,     -4.609252, -3.264385, -2.644414, -2.278346)
    D_n <- c(555, 777, 999, 1221, 1443)
  }
  
  K <- length(lambda_E)
  theta <- log(lambda_E / lambda_C)
  
  ## -----------------------------
  ## Individual-level data sim
  ## -----------------------------
  ## The code uses:
  ##   rweibullPH(N, shape, scale) for each arm’s event time
  ##   runif for entry times (accrual) in [0, a_t]
  ## Generate event times per arm (control + K experimental)
  da <- data.frame(
    cbind(
      reshape2::melt(
        mapply(function(scale) flexsurv::rweibullPH(N, shape, scale),
               scale = c(lambda_C, lambda_E)),
        value.name = "y"
      ),
      a_t = runif((K+1)*N, 0, a_t)
    )
  )
  
  ## Treatment indicator (0=control, 1..K=experimental arms)
  da$treatment <- rep(0:K, each=N)
  
  ## Calendar event time = entry time + event time
  da$D_time <- da$y + da$a_t
  
  ## stage_da_f: construct data at calendar analysis time f_t:
  ## - status indicates whether event occurred by f_t
  ## - time is follow-up time: min(y, f_t - entry)
  stage_da_f <- function(f_t){
    status <- da$D_time <= f_t
    time <- f_t - da$a_t
    time[status] <- da$y[status]
    data.frame(treatment = da$treatment, status = status, time = time)
  }
  
  ## Cox model summary function for arm i vs control at interim look j
  ## Returns coefficient table from summary(coxph)$coef (matrix)
  Z_sel_f <- function(i, j){
    stage_da <- stage_da_f(f_t = inter_t[j])
    summary(survival::coxph(survival::Surv(time, status) ~ factor(treatment),
                            data = stage_da[stage_da$treatment %in% c(0, i), ]))$coef
  }
  
  ## I(i,j): approximate information from Cox model SE:
  ## Here using 1 / (SE^2), where SE is in coef table col 3 
  I <- function(i, j){
    1 / (Z_sel_f(i, j)[3])^2
  }
  
  ## S(i,j): score-like quantity = beta_hat * I
  S <- function(i, j){
    Z_sel_f(i, j)[1] * I(i, j)
  }
  
  ## -----------------------------
  ## Determine interim calendar times inter_t by planned event numbers D_n
  ## -----------------------------
  ## NOTE: your raw code had inter_t <- length(J) (wrong). Use numeric(J).
  inter_t <- numeric(J)
  
  ## First look: calendar time at which total pooled D_time reaches D_n[1]
  inter_t[1] <- sort(da$D_time)[D_n[1]]
  
  ## Stage-1 Cox across all arms vs control
  stage1_cox <- survival::coxph(survival::Surv(time, status) ~ factor(treatment),
                                data = stage_da_f(inter_t[1]))
  
  ## Selection rule: choose arm with minimum coefficient (most negative logHR)
  select_index <- which.min(stage1_cox$coefficients)
  
  ## For later looks, only consider patients in control + selected arm not yet in data cutoff
  further_time <- da[(da$treatment %in% c(0, select_index)) & (da$D_time > inter_t[1]), ]$D_time
  
  ## Define later looks so that additional events among (control + selected arm)
  ## reach the planned event totals D_n[2:J]
  inter_t[2:J] <- sort(further_time)[D_n[2:J] - D_n[1]]
  
  ## -----------------------------
  ## Stopping rule at each stage
  ## -----------------------------
  inter_decide <- function(j, select_index){
    Z_temp <- Z_sel_f(select_index, j)[4]
    
    ## Your sign convention:
    ##   - If Z_temp > low_bround[j] => futility ("f")
    ##   - Else if Z_temp < up_bround[j] => efficacy ("e")
    ##   - Else continue ("c")
    if(Z_temp > low_bround[j]) "f" else if(Z_temp < up_bround[j]) "e" else "c"
  }
  
  inter_outcome <- character(J)
  for(j in 1:J){
    inter_outcome[j] <- inter_decide(j, select_index)
  }
  
  ## Stop at first e/f
  end_stage <- data.table::first(which(inter_outcome %in% c("e","f")))
  
  ## -----------------------------
  ## Estimators
  ## -----------------------------
  ## naive: stage-1 cox coefficients for K arms
  naive <- stage1_cox$coefficients
  
  ## Replace selected arm's naive by cumulative estimator at end_stage
  naive[select_index] <- S(select_index, end_stage) / I(select_index, end_stage)
  
  ## Stage-2 increment estimator for selected arm (between look1 and look2)
  stage2 <- (S(select_index, 2) - S(select_index, 1)) / (I(select_index, 2) - I(select_index, 1))
  
  ## -----------------------------
  ## Build covariance matrix cov_S for RB / conditional methods
  ## -----------------------------
  cov_S <- matrix(0, nrow=K+J-1, ncol=K+J-1)
  
  ## Covariance of stage-1 coefficients from coxph
  ## stage1_cox$var is Var(beta_hat). Convert to Var(S) = I * Var(beta) * I
  I_stage1 <- diag(mapply(I, 1:K, 1))
  
  cov_stage1 <- stage1_cox$var
  diag(cov_stage1) <- 1 / mapply(I, 1:K, 1)   ## enforce diagonal consistent with I
  
  cov_S[1:K, 1:K] <- I_stage1 %*% cov_stage1 %*% t(I_stage1)
  
  ## Selected-arm increment variances up to end_stage
  if(end_stage > 1){
    if(end_stage > 2){
      diag(cov_S[(K+1):(K+end_stage-1), (K+1):(K+end_stage-1)]) <- diff(mapply(I, select_index, 1:end_stage))
    }else{
      cov_S[(K+1):(K+end_stage-1), (K+1):(K+end_stage-1)] <- diff(mapply(I, select_index, 1:end_stage))
    }
  }
  
  ## For remaining looks after stopping, you fill with a 1/4 * diff(D_n) approximation.
  ## This is a modeling choice for the future increments.
  if(J > end_stage){
    if(J - end_stage > 1){
      diag(cov_S[(K+end_stage):(K+J-1), (K+end_stage):(K+J-1)]) <- 1/4 * diff(D_n)[end_stage:(J-1)]
    }else{
      cov_S[(K+end_stage):(K+J-1), (K+end_stage):(K+J-1)] <- 1/4 * diff(D_n)[end_stage:(J-1)]
    }
  }
  
  ## S_all: stage-1 S for all arms + selected cumulative S at end_stage
  S_all <- c(mapply(S, 1:K, 1), S(select_index, end_stage))
  
  list(end_stage=end_stage,
       select_index=select_index,
       test=inter_outcome[end_stage],
       naive=naive,
       bias=c(stage2, naive[select_index]),  ## your stored "bias" pieces
       cov_S=cov_S,
       S_all=S_all)
}
da_genc <- cmpfun(da_gen)

## (5) Helper transforms: cumulative S and Z differences for selection event
S_cum_f <- function(cov_S, select_index, K, J){
  
  ## trans_t: reorders the first K components so that selected arm is at row K
  trans_t <- diag(1, K+J-1)
  trans_t[1:K, 1:K][-K,] <- diag(K)[-select_index,]
  trans_t[1:K, 1:K][K, ] <- diag(K)[ select_index,]
  
  ## trans_c: converts increments into cumulative sums (for selected arm across looks)
  trans_c <- diag(1, K+J-1)
  trans_c[K:(K+J-1), K:(K+J-1)][lower.tri(trans_c[K:(K+J-1), K:(K+J-1)])] <- 1
  
  S_sum <- trans_c %*% (trans_t %*% cov_S %*% t(trans_t)) %*% t(trans_c)
  
  list(S_sum  = S_sum,
       I_stage = diag(S_sum))  ## "information" is the diagonal variances on S_sum scale
}

S_Z_diff_f <- function(cov_S, theta, select_index, K, J){
  
  I_stage <- S_cum_f(cov_S, select_index, K, J)$I_stage
  
  ## Standardize to Z-scale
  theta_v <- diag(1/sqrt(I_stage), K+J-1)
  cov_Z <- theta_v %*% S_cum_f(cov_S, select_index, K, J)$S_sum %*% t(theta_v)
  
  ## Contrast matrix A encodes selection event comparisons (selected vs others)
  A <- diag(1, K+J-1)
  diag(A[1:(K-1), 1:(K-1)]) <- -1
  A[1:(K-1), K] <- 1
  
  ## Mean on Z-scale:
  ## - non-selected arms use first (K-1) entries of I_stage after reordering
  ## - selected arm across looks use entries K..K+J-1
  mu_Z <- c(theta[-select_index] * sqrt(I_stage[1:(K-1)]),
            theta[select_index] * sqrt(I_stage[K:(K+J-1)]))
  
  list(mu  = A %*% mu_Z,
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
bias_theta <- function(select_index, theta, K, J, cov_S){
  
  ## Boundaries (same sign convention as da_gen)
  if(J==4){
    low_bround <- c(0.000000,  4.013691,  2.828157, -2.285393)
    up_bround  <- c(-Inf,     -4.013691, -2.828157, -2.285393)
  }else if(J==3){
    low_bround <- c(0.000000,  3.301594, -2.302597)
    up_bround  <- c(-Inf,     -3.301594, -2.302597)
  }else if(J==5){
    low_bround <- c(0.000000,  4.609252,  3.264385,  2.644414, -2.278346)
    up_bround  <- c(-Inf,     -4.609252, -3.264385, -2.644414, -2.278346)
    }
  
  I_stage <- S_cum_f(cov_S, select_index, K, J)$I_stage
  mu_cov <- S_Z_diff_f(cov_S, theta, select_index, K, J)
  mu  <- mu_cov$mu
  cov <- mu_cov$cov
  
  ## --- Compute prob(Q): probability of selection/continuation event ---
  ## Q = { selected arm is best at stage 1 } ∩ { continuation after stage 1 }
  ## This uses MVN probability over first K components.
  prob_Q <- mvtnorm::pmvnorm(
    lower = c(rep(-Inf, K-1), up_bround[1]),
    upper = c(rep(0,   K-1), low_bround[1]),
    mean  = mu[1:K,],
    sigma = cov[1:K, 1:K]
  )[1]
  
  bias <- rep(NA, K)
  
  if(!is.na(prob_Q)){
    
    ## Truncated MVN mean under selection region
    mtm <- tmvtnorm::mtmvnorm(
      mean = mu[1:K,], sigma = cov[1:K, 1:K],
      lower = c(rep(-Inf, K-1), up_bround[1]),
      upper = c(rep(0,   K-1), low_bround[1]),
      doComputeVariance = FALSE
    )$tmean
    
    ## Bias for non-selected arms
    bias[-select_index] <- (mtm[K] - mtm[-K]) / sqrt(I_stage[1:(K-1)]) - theta[-select_index]
    
    ## Bias contribution for selected arm averaged over possible stopping stages 2..J
    E_select <- function(stage){
      
      mu_t  <- mu[1:(K+stage-1),]
      cov_t <- cov[1:(K+stage-1), 1:(K+stage-1)]
      
      lower1 <- c(rep(-Inf, K-1), up_bround[1:(stage-1)], -Inf)
      upper1 <- c(rep(0,   K-1), low_bround[1:(stage-1)], up_bround[stage])
      
      lower2 <- c(rep(-Inf, K-1), up_bround[1:(stage-1)], low_bround[stage])
      upper2 <- c(rep(0,   K-1), low_bround[1:(stage-1)], Inf)
      
      ## Weighted truncated means in the two terminal regions
      ( sum(
        mvtnorm::pmvnorm(lower1, upper1, mean=mu_t, sigma=cov_t)[1] *
          tmvtnorm::mtmvnorm(mean=mu_t, sigma=cov_t, lower=lower1, upper=upper1,
                             doComputeVariance=FALSE)$tmean[K+stage-1],
        mvtnorm::pmvnorm(lower2, upper2, mean=mu_t, sigma=cov_t)[1] *
          tmvtnorm::mtmvnorm(mean=mu_t, sigma=cov_t, lower=lower2, upper=upper2,
                             doComputeVariance=FALSE)$tmean[K+stage-1],
        na.rm=TRUE
      )
      ) / sqrt(I_stage[K:(K+J-1)][stage])
    }
    
    bias[select_index] <- sum(mapply(E_select, stage=2:J)) / prob_Q - theta[select_index]
  }
  
  bias
}

## (A2) MI: iterative bias correction
##      theta_{new} = initial - bias_theta(theta_old) until convergence.
MI <- function(initial, select_index, max_iterations=20, tol=0.001, K, J, cov_S){
  
  conver <- 0
  theta_hat <- initial
  
  for(i in 1:max_iterations){
    
    bias_temp <- bias_theta(select_index, theta_hat, K, J, cov_S)
    
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
SI <- function(initial, select_index, K, J, cov_S){
  (initial - bias_theta(select_index, initial, K, J, cov_S))[select_index]
}
MIc <- cmpfun(MI)

##### -------------------------------------------------------------------------
##### (B) CMUE: conditional median unbiased estimator
#####     - compute conditional p-value under candidate theta
#####     - invert p-value to get CMUE point estimate (median) and CI endpoints
##### -------------------------------------------------------------------------
est_MUE <- function(naive, end_stage, select_index, sub,
                    K, J, cov_S, S_all){
  
  if(J==4){
    low_bround <- c(0.000000,  4.013691,  2.828157, -2.285393)
    up_bround  <- c(-Inf,     -4.013691, -2.828157, -2.285393)
  }else if(J==3){
    low_bround <- c(0.000000,  3.301594, -2.302597)
    up_bround  <- c(-Inf,     -3.301594, -2.302597)
  }else if(J==5){
    low_bround <- c(0.000000,  4.609252,  3.264385,  2.644414, -2.278346)
    up_bround  <- c(-Inf,     -4.609252, -3.264385, -2.644414, -2.278346)
  }
  
  ## Observed Z at end_stage for selected arm (from cumulative S)
  Z_end <- S_all[K+1] / sqrt(S_cum_f(cov_S, select_index, K, J)$I_stage[K+end_stage-1])
  
  ## Conditional p-value function under candidate theta
  cond_pvalue <- function(theta, select_index, end_stage, Z_end, cov_S, K, J){
    
    mu_cov <- S_Z_diff_f(cov_S, theta, select_index, K, J)
    mu  <- mu_cov$mu
    cov <- mu_cov$cov
    
    prob_Q <- mvtnorm::pmvnorm(
      lower = c(rep(-Inf, K-1), up_bround[1]),
      upper = c(rep(0,   K-1), low_bround[1]),
      mean  = mu[1:K,],
      sigma = cov[1:K, 1:K]
    )[1]
    
    p_value <- function(end_stage, Z_end){
      
      P_stage <- function(stage, stop, Z_end){
        
        new_mu  <- mu[1:(K+stage-1),]
        new_cov <- cov[1:(K+stage-1), 1:(K+stage-1)]
        
        if(stop){
          P <- mvtnorm::pmvnorm(
            lower = c(rep(-Inf, K-1), up_bround[1:(stage-1)], -Inf),
            upper = c(rep(0,   K-1), low_bround[1:(stage-1)], Z_end),
            mean  = new_mu, sigma=new_cov
          )[1]
        }else{
          P <- mvtnorm::pmvnorm(
            lower = c(rep(-Inf, K-1), up_bround[1:(stage-1)], -Inf),
            upper = c(rep(0,   K-1), low_bround[1:(stage-1)], up_bround[stage]),
            mean  = new_mu, sigma=new_cov
          )[1]
        }
        
        if(is.na(P)) 0 else P
      }
      
      if(end_stage == 2){
        numerator <- P_stage(2, TRUE, Z_end)
      }else{
        numerator <- sum(mapply(P_stage, stage=2:(end_stage-1), stop=FALSE, Z_end=Z_end),
                         P_stage(end_stage, TRUE, Z_end))
      }
      
      if(is.na(numerator) || is.na(prob_Q) || prob_Q <= 0) return(1)
      min(1, numerator/prob_Q)
    }
    
    p_value(end_stage, Z_end)
  }
  
  cond_pvaluec <- cmpfun(cond_pvalue)
  
  ## Plug-in strategy for non-selected arms:
  ## sub="MLE": fix non-selected arms at their MLEs
  ## sub="ZERO": fix non-selected arms at 0
  P_value_sub <- function(x){
    theta <- rep(x, K)
    if(sub == "MLE"){
      theta[-select_index] <- naive[-select_index]
    }else{
      theta[-select_index] <- rep(0, K-1)
    }
    cond_pvaluec(theta, select_index, end_stage, Z_end, cov_S, K, J)
  }
  
  ## Solve for CMUE median: p(theta) = 0.5
  ## NOTE: bracket range (-1,1) may need to widen for larger effects.
  a <- optimise(function(x) (P_value_sub(x) - 1/2)^2,
                interval = c(-2, 2) * naive[select_index],
                tol = 0.001)
  
  c(a$minimum, a$objective)
}

##### -------------------------------------------------------------------------
##### (C) RB: Rao–Blackwell estimator via conditional simulation
#####     - condition an unbiased stage-wise estimator on sufficient statistic
#####     - Monte Carlo implementation of conditional expectation
##### -------------------------------------------------------------------------
RB_mc <- function(S_all, cov_S, select_index, end_stage, mc_num, alpha, K, J){
  
  ## Boundaries
  if(J==4){
    low_bround <- c(0.000000,  4.013691,  2.828157, -2.285393)
    up_bround  <- c(-Inf,     -4.013691, -2.828157, -2.285393)
  }else if(J==3){
    low_bround <- c(0.000000,  3.301594, -2.302597)
    up_bround  <- c(-Inf,     -3.301594, -2.302597)
  }else if(J==5){
    low_bround <- c(0.000000,  4.609252,  3.264385,  2.644414, -2.278346)
    up_bround  <- c(-Inf,     -4.609252, -3.264385, -2.644414, -2.278346)
  }else{
    stop("Only J=3/4/5 supported.")
  }
  
  ## “Information” from diagonal of S_sum after transforms
  I_stage <- S_cum_f(cov_S, select_index, K, J)$I_stage
  
  ## Working MVN model with mean zero
  mu <- rep(0, K+end_stage-1)
  
  ## Stage-1 diagonal variances (used in transform)
  I_stage1 <- diag(cov_S)[1:K]
  
  ## Transform to T statistics and condition on observed T_star
  trans_T <- diag(1, K+end_stage-1)
  trans_T[select_index,] <- 1
  trans_T[1:K, 1:K] <- solve(cov_S[1:K, 1:K]) * I_stage1
  
  mu_T  <- trans_T %*% mu
  cov_T <- trans_T %*% cov_S[1:(K+end_stage-1), 1:(K+end_stage-1)] %*% t(trans_T)
  
  ## Observed T_star from observed S_all
  T_star <- (I_stage1 * solve(cov_S[1:K, 1:K])) %*% as.matrix(S_all[1:K])
  T_star[select_index] <- T_star[select_index] + (S_all[K+1] - S_all[select_index])
  
  ## Step 1: conditional simulation of remaining components given T_star
  mc_s <- condMVNorm::rcmvnorm(mc_num, mu_T, cov_T,
                               given.ind=1:K,
                               dependent.ind=(K+1):(K+end_stage-1),
                               X.given=T_star)
  
  s1 <- matrix(NA, nrow=K, ncol=mc_num)
  for(i in 1:mc_num){
    T_star_s <- T_star
    T_star_s[select_index] <- T_star[select_index] - sum(mc_s[i,])
    s1[,i] <- solve(I_stage1 * solve(cov_S[1:K, 1:K])) %*% T_star_s
  }
  
  z_stage1 <- apply(s1, 2, function(x) x / sqrt(I_stage1))
  mc_s_all <- cbind(s1[select_index,], mc_s)
  
  ## Step 2: accept only samples consistent with selection event and boundaries
  mc_z <- t(apply(mc_s_all, 1, function(x) cumsum(x) / sqrt(I_stage[K:(K+end_stage-1)])))
  
  bround_L <- matrix(low_bround[1:(end_stage-1)], mc_num, end_stage-1, byrow=TRUE)
  bround_U <- matrix(up_bround[1:(end_stage-1)],  mc_num, end_stage-1, byrow=TRUE)
  
   index <- apply(cbind(
    apply(mc_z[,1:(end_stage-1)] >= bround_U, 1, function(x) identical(x, rep(TRUE, end_stage-1))),
    apply(mc_z[,1:(end_stage-1)] <= bround_L, 1, function(x) identical(x, rep(TRUE, end_stage-1))),
    apply(z_stage1, 2, which.min) == select_index
  ), 1, all)
  
  re <- rep(NA, 4)
  if(sum(index) > 1){
    ## Step 3: RB estimate = conditional expectation of stage2 unbiased increment
    mc_est <- mean(mc_s_all[index, 2] / (I_stage[K+1] - I_stage[K]))
    
    SE <- sqrt(1/(I_stage[K+1] - I_stage[K]) -
                 var(mc_s_all[index, 2] / (I_stage[K+1] - I_stage[K])))
    
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
# setwd("C:/Users/13379/Desktop/essay/submit-SIM/survival-cov/4")

## Global simulation settings
shape   <- 0.5
setting <- "line"
a_t     <- 1
N       <- 1000
J       <- 4
sim     <- 10

time1 <- Sys.time()

## Parallel cluster
cls <- makeSOCKcluster(18)
registerDoSNOW(cls)

## Progress bar
pb <- txtProgressBar(max=sim, style=3, char="*")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

res <- foreach(i=1:sim,
               .options.snow=opts,
               .combine=cbind,
               .packages=c("data.table","compiler","condMVNorm","tmvtnorm",
                           "nnet","flexsurv","reshape2","survival","mvtnorm")) %dopar% {
                             da_genc(i, setting, J)
                           }

## Extract components from each replicate list
end_stage    <- apply(res, 2, function(x) x$end_stage)
select_index <- apply(res, 2, function(x) x$select_index)
naive        <- apply(res, 2, function(x) x$naive)
bias         <- apply(res, 2, function(x) x$bias)

idx <- which(end_stage >= 2)

write.csv(cbind(end_stage, select_index, t(bias))[idx,],
          paste(setting, "stage2-naive.csv", sep="-"),
          row.names=FALSE)

###############################################################################
## 3.2 Run each estimator for eligible replicates and save to CSV
###############################################################################
## --- CBAE-SI ---
re_SI <- foreach(i=idx,
                 .options.snow=opts,
                 .combine=cbind,
                 .packages=c("tmvtnorm","mvtnorm","compiler")) %dopar% {
                   SI(res[,i]$naive, res[,i]$select_index, K=4, J=J, cov_S=res[,i]$cov_S)
                 }

write.csv(cbind(end_stage[idx], select_index[idx], as.numeric(re_SI)),
          paste(setting, "SI.csv", sep="-"),
          row.names=FALSE)
## --- CBAE-MI ---
re_MI <- foreach(i=idx,
                 .options.snow=opts,
                 .combine=cbind,
                 .packages=c("tmvtnorm","mvtnorm","compiler")) %dopar% {
                   MIc(res[,i]$naive, res[,i]$select_index,
                       max_iterations=20, tol=0.001,
                       K=4, J=J, cov_S=res[,i]$cov_S)
                 }

write.csv(cbind(end_stage[idx], select_index[idx], t(re_MI)),
          paste(setting, "MI.csv", sep="-"),
          row.names=FALSE)


## --- CMUE-MLE ---
re_MUE_MLE <- foreach(i=idx,
                      .options.snow=opts,
                      .combine=cbind,
                      .packages=c("mvtnorm","tmvtnorm","compiler")) %dopar% {
                        est_MUE(res[,i]$naive, end_stage[i], res[,i]$select_index,
                                sub="MLE", K=4, J=J,
                                cov_S=res[,i]$cov_S, S_all=res[,i]$S_all)
                      }

write.csv(cbind(end_stage[idx], select_index[idx], t(re_MUE_MLE)),
          paste(setting, "MUE_MLE.csv", sep="-"),
          row.names=FALSE)

## --- CMUE-ZERO ---
re_MUE_ZERO <- foreach(i=idx,
                       .options.snow=opts,
                       .combine=cbind,
                       .packages=c("mvtnorm","tmvtnorm","compiler")) %dopar% {
                         est_MUE(res[,i]$naive, end_stage[i], res[,i]$select_index,
                                 sub="ZERO", K=4, J=J,
                                 cov_S=res[,i]$cov_S, S_all=res[,i]$S_all)
                       }

write.csv(cbind(end_stage[idx], select_index[idx], t(re_MUE_ZERO)),
          paste(setting, "MUE_ZERO.csv", sep="-"),
          row.names=FALSE)

## --- RB ---
re_RB <- foreach(i=idx,
                 .options.snow=opts,
                 .combine=cbind,
                 .packages=c("condMVNorm","compiler")) %dopar% {
                   RB_mc(S_all=res[,i]$S_all,
                         cov_S=res[,i]$cov_S,
                         select_index=select_index[i],
                         end_stage=end_stage[i],
                         mc_num=10000, alpha=0.05,
                         K=4, J=J)
                 }

write.csv(cbind(end_stage[idx], select_index[idx], t(re_RB)),
          paste(setting, "RB.csv", sep="-"),
          row.names=FALSE)


###############################################################################
## 4) Clean up parallel backend and report elapsed time
###############################################################################
close(pb)
stopCluster(cls)

time2 <- Sys.time()
print(time2 - time1)
