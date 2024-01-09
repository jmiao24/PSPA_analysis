library(stats)
library(rootSolve)
library(parallel)
library(purrr)

### OLS ###
ols <- function(features, outcome) {
  features <- as.matrix(features)
  outcome <- as.matrix(outcome)
  ols_coeffs <- solve(t(features) %*% features) %*% t(features) %*% outcome
  return(ols_coeffs)
}

# Sampling variance of ols estimator
## Y1 and Y2 must have the same length
ols_var_est <- function(X, Y1, Y2, sandwich=TRUE) {
  n <- nrow(X)
  beta_Y1 <- ols(X, Y1)
  beta_Y2 <- ols(X, Y2)
  Sigmainv <- solve(1/n * t(X) %*% X)
  if (sandwich) {
    M <- 1/n * t(X) %*% (X * ((Y1 - X %*% beta_Y1) * (Y2 - X %*% beta_Y2)))
  } else {
    M <- 1/n * mean((Y1 - X %*% beta_Y1) * (Y2 - X %*% beta_Y2)) * t(X) %*% X
  }
  V <- Sigmainv %*% M %*% Sigmainv
  var_est <- diag(V)
  return(var_est)
}


classic_ols_asymptotic <- function(X, Y, alpha, sandwich=TRUE) {
  X <- cbind(1, as.matrix(X))
  Y <- as.matrix(Y)
  n <- nrow(X)
  beta_classic <- ols(X, Y)
  stderr_classic <- sqrt(ols_var_est(X, Y, Y, sandwich = sandwich)/n)
  halfwidth <- qnorm(1 - alpha/2) * stderr_classic
  lower = beta_classic - halfwidth
  upper = beta_classic + halfwidth
  P_classic <- 2*pnorm(abs(beta_classic/stderr_classic), 0, 1, lower.tail = F)
  return(list(theta = beta_classic, se = stderr_classic, lower_ci = lower, upper_ci = upper, P = P_classic))
}

pp_ols_asymptotic <- function(X_labeled, X_unlabeled, Y_labeled, Yhat_labeled, Yhat_unlabeled, alpha, sandwich=TRUE) {
  X_labeled <- cbind(1, as.matrix(X_labeled))
  X_unlabeled <- cbind(1, as.matrix(X_unlabeled))
  Y_labeled <- as.matrix(Y_labeled)
  Yhat_labeled <- as.matrix(Yhat_labeled)
  Yhat_unlabeled <- as.matrix(Yhat_unlabeled)
  n <- nrow(X_labeled)
  N <- nrow(X_unlabeled)
  thetatildef <- ols(X_unlabeled, Yhat_unlabeled)
  rectifierhat <- ols(X_labeled, Y_labeled - Yhat_labeled)
  pp_thetahat <- thetatildef + rectifierhat
  stderr_tildef <- sqrt(ols_var_est(X_unlabeled, Yhat_unlabeled, Yhat_unlabeled, sandwich=sandwich)/N)
  stderr_rec <- sqrt(ols_var_est(X_labeled, Y_labeled - Yhat_labeled, Y_labeled - Yhat_labeled, sandwich=sandwich)/n)
  # Calculations
  beta_PP <- thetatildef + rectifierhat
  stderr_PP <- sqrt(stderr_tildef^2 + stderr_rec^2)
  halfwidth <- qnorm(1 - alpha/2) * stderr_PP
  lower <- beta_PP - halfwidth
  upper <- beta_PP + halfwidth
  P_PP <- 2*pnorm(abs(beta_PP/stderr_PP), 0, 1, lower.tail = F) 
  return(list(theta = beta_PP, se = stderr_PP, lower_ci = lower, upper_ci = upper, P = P_PP))
}

eff_ols_asymptotic <- function(X_labeled, X_unlabeled, 
                               Y_labeled, Yhat_labeled, Yhat_unlabeled, 
                               alpha, sandwich = TRUE) {
  X_labeled <- cbind(1, as.matrix(X_labeled))
  X_unlabeled <- cbind(1, as.matrix(X_unlabeled))
  X <- rbind(X_labeled, X_unlabeled)
  
  Y_labeled <- as.matrix(Y_labeled)
  Yhat_labeled <- as.matrix(Yhat_labeled)
  Yhat_unlabeled <-as.matrix(Yhat_unlabeled)
  Yhat <- rbind(Yhat_labeled, Yhat_unlabeled)
  
  n <- nrow(X_labeled)
  N <- nrow(X_unlabeled)
  
  beta_resid <- ols(X_labeled, Y_labeled - Yhat_labeled)
  beta_imp <- ols(X, Yhat)
  beta_eff <- beta_resid + beta_imp
  
  stderr_imp <- sqrt(ols_var_est(X, Yhat, Yhat, sandwich=sandwich)/(n + N))
  stderr_resid <- sqrt(ols_var_est(X_labeled, Y_labeled - Yhat_labeled, 
                                   Y_labeled - Yhat_labeled, sandwich=sandwich)/n)
  # Calculations
  stderr_eff <- sqrt(stderr_imp ^ 2 + stderr_resid ^ 2)
  halfwidth <- qnorm(1 - alpha/2) * stderr_eff
  lower <- beta_eff - halfwidth
  upper <- beta_eff + halfwidth
  P_eff <- 2 * pnorm(abs(beta_eff/stderr_eff), 0, 1, lower.tail = F) 
  return(list(theta = beta_eff, se = stderr_eff, 
              lower_ci = lower, upper_ci = upper, P = P_eff))
}

oracle_ols_asymptotic <- function(X_labeled, X_unlabeled, 
                                  Y_labeled, Yhat_labeled, Yhat_unlabeled, 
                                  alpha, sandwich = TRUE) {
  X_labeled <- cbind(1, as.matrix(X_labeled))
  X_unlabeled <- cbind(1, as.matrix(X_unlabeled))
  X <- rbind(X_labeled, X_unlabeled)
  
  Y_labeled <- as.matrix(Y_labeled)
  Yhat_labeled <- as.matrix(Yhat_labeled)
  Yhat_unlabeled <- as.matrix(Yhat_unlabeled)
  Yhat <- rbind(Yhat_labeled, Yhat_unlabeled)
  
  n <- nrow(X_labeled)
  N <- nrow(X_unlabeled)
  
  beta_resid <- ols(X_labeled, Y_labeled - Yhat_labeled)
  beta_imp <- ols(X, Yhat)
  beta_eff <- beta_resid + beta_imp
  
  stderr_imp <- sqrt(ols_var_est(X, Yhat, Yhat, sandwich=sandwich)/(n + N))
  stderr_resid <- sqrt(ols_var_est(X_labeled, Y_labeled - Yhat_labeled, 
                                   Y_labeled - Yhat_labeled, sandwich=sandwich)/n)
  # Calculations
  stderr_eff <- sqrt(stderr_imp ^ 2 + stderr_resid ^ 2)
  halfwidth <- qnorm(1 - alpha/2) * stderr_eff
  lower <- beta_eff - halfwidth
  upper <- beta_eff + halfwidth
  P_eff <- 2 * pnorm(abs(beta_eff/stderr_eff), 0, 1, lower.tail = F) 
  return(list(theta = beta_eff, se = stderr_eff, 
              lower_ci = lower, upper_ci = upper, P = P_eff))
}

PopInf_ols_asymptotic <- function(X_labeled, X_unlabeled, Y_labeled, Yhat_labeled, Yhat_unlabeled, alpha, w = NA, sandwich=TRUE) {
  
  X_labeled <- cbind(1, as.matrix(X_labeled))
  X_unlabeled <- cbind(1, as.matrix(X_unlabeled))
  Y_labeled <- as.matrix(Y_labeled)
  Yhat_labeled <- as.matrix(Yhat_labeled)
  Yhat_unlabeled <- as.matrix(Yhat_unlabeled)
  n <- nrow(X_labeled)
  N <- nrow(X_unlabeled)
  
  # Point estimates
  beta_Yhat_unlabeled <- ols(X_unlabeled, Yhat_unlabeled)
  beta_Y_labeled <- ols(X_labeled, Y_labeled)
  beta_Yhat_labeled <- ols(X_labeled, Yhat_labeled)
  
  # Variance of point estimates
  var_Yhat_unlabeled <- ols_var_est(X_unlabeled, Yhat_unlabeled, Yhat_unlabeled, sandwich=sandwich)
  var_Y_labeled <- ols_var_est(X_labeled, Y_labeled, Y_labeled, sandwich=sandwich)
  var_Yhat_labeled <- ols_var_est(X_labeled, Yhat_labeled, Yhat_labeled, sandwich=sandwich)
  cov_Y_Yhat_labeled <- ols_var_est(X_labeled, Yhat_labeled, Y_labeled, sandwich=sandwich)
  
  if (is.na(w)){
    w <- (cov_Y_Yhat_labeled/n) / (var_Yhat_unlabeled/N + var_Yhat_labeled/n)
  }
  
  # Results
  beta_PopInf <- beta_Y_labeled - w * beta_Yhat_labeled + w * beta_Yhat_unlabeled
  stderr_PopInf <- sqrt(var_Y_labeled/n + w^2 * var_Yhat_labeled/n - 2 * w * cov_Y_Yhat_labeled/n + w^2 * var_Yhat_unlabeled/N)
  P_PopInf <- 2*pnorm(abs(beta_PopInf/stderr_PopInf), 0, 1, lower.tail = F)
  halfwidth <- qnorm(1 - alpha/2) * stderr_PopInf
  upper <- beta_PopInf + halfwidth
  lower <- beta_PopInf - halfwidth
  return(list(theta = beta_PopInf, se = stderr_PopInf, lower_ci = lower, upper_ci = upper, P = P_PopInf, w = w))
}

