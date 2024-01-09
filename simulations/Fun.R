library(stats)
library(rootSolve)
library(parallel)
library(purrr)

### Mean estimations ###
classic_mean_asymptotic <- function(Y_labeled, Yhat_labeled, Yhat_unlabeled, alpha) {
  n <- length(Y_labeled)
  beta_classic <- mean(Y_labeled)
  var_Y_labeled <- var(Y_labeled)
  # Calculations
  beta_classic <-mean(Y_labeled)
  se_classic <- sqrt((var_Y_labeled / n))
  P_classic <- 2*pnorm(abs(beta_classic/se_classic), 0, 1, lower.tail = F)
  upper <- beta_classic + qnorm(1 - alpha/2) * se_classic
  lower <- beta_classic - qnorm(1 - alpha/2) * se_classic
  return(list(theta = beta_classic, se = se_classic, lower_ci = lower, upper_ci = upper, P = P_classic))
}

pp_mean_asymptotic <- function(Y_labeled, Yhat_labeled, Yhat_unlabeled, alpha) {
  n <- length(Y_labeled)
  N <- length(Yhat_unlabeled)
  tildethetaf <- mean(Yhat_unlabeled)
  rechat <- mean(Yhat_labeled - Y_labeled)
  thetahatPP <- tildethetaf - rechat
  sigmaftilde <- sd(Yhat_unlabeled)
  sigmarec <- sd(Yhat_labeled - Y_labeled)
  # Calculations
  beta_PP <- tildethetaf - rechat
  stderr_PP <- sqrt((sigmaftilde^2 / N) + (sigmarec^2 / n))
  P_PP <- 2*pnorm(abs(beta_PP/stderr_PP), 0, 1, lower.tail = F)
  upper <- beta_PP + qnorm(1 - alpha/2) * stderr_PP
  lower <- beta_PP - qnorm(1 - alpha/2) * stderr_PP
  return(list(theta = beta_PP, se = stderr_PP, lower_ci = lower, upper_ci = upper, P = P_PP))
}

PopInf_mean_asymptotic <- function(Y_labeled, Yhat_labeled, Yhat_unlabeled, alpha, w = NA) {
  n <- length(Y_labeled)
  N <- length(Yhat_unlabeled)
  beta_Yhat_unlabeled <- mean(Yhat_unlabeled)
  beta_Y_labeled <- mean(Y_labeled)
  beta_Yhat_labeled <- mean(Yhat_labeled)
  var_Yhat_unlabeled <- var(Yhat_unlabeled)
  var_Y_labeled <- var(Y_labeled)
  var_Yhat_labeled <- var(Yhat_labeled)
  cov_Y_Yhat_labeled <- cov(Y_labeled, Yhat_labeled)
  if (is.na(w)){
      w <- (cov_Y_Yhat_labeled/n) / (var_Yhat_unlabeled/N + var_Yhat_labeled/n)
  }
  # Calculations
  beta_PopInf <- beta_Y_labeled - w * beta_Yhat_labeled + w * beta_Yhat_unlabeled
  stderr_PopInf <- sqrt(var_Y_labeled/n + w^2 * var_Yhat_labeled/n - 2 * w * cov_Y_Yhat_labeled/n + w^2 * var_Yhat_unlabeled/N)
  P_PopInf <- 2*pnorm(abs(beta_PopInf/stderr_PopInf), 0, 1, lower.tail = F)
  upper <- beta_PopInf + qnorm(1 - alpha/2) * stderr_PopInf
  lower <- beta_PopInf - qnorm(1 - alpha/2) * stderr_PopInf
  return(list(theta = beta_PopInf, se = stderr_PopInf, lower_ci = lower, upper_ci = upper, P = P_PopInf, w = w))
}

# Calculate efficient influence function based estimator
eff_mean_asymptotic <- function(Y_labeled, Yhat_labeled, Yhat_unlabeled, alpha) {
  n <- length(Y_labeled)
  N <- length(Yhat_unlabeled)
  
  # (weighted) average of residuals on labeled data
  mean_resid <- mean(Y_labeled - Yhat_labeled)
  # average of imputations on all data
  mean_imp <- mean(c(Yhat_labeled, Yhat_unlabeled))
  
  beta_eff <- mean_resid + mean_imp
  
  r. <- c(rep(1, n), rep(0, N))
  pii <- n / (n + N)
  y <- c(Y_labeled, rep(0, N))
  yhat <- c(Yhat_labeled, Yhat_unlabeled)
  # influence function
  phi <- r. * (y - yhat) / pii + yhat - beta_eff
  stderr_eff <- sqrt(mean(phi ^ 2)) / sqrt(n + N)
  
  P_eff <- 2 * pnorm(abs(beta_eff / stderr_eff), 0, 1, lower.tail = FALSE)
  upper <- beta_eff + qnorm(1 - alpha/2) * stderr_eff
  lower <- beta_eff - qnorm(1 - alpha/2) * stderr_eff
  
  return(list(theta = beta_eff, se = stderr_eff, 
              lower_ci = lower, upper_ci = upper, P = P_eff))
}


# Oracle estimator guided by EIF with true E(Y|x)
oracle_mean_asymptotic <- function(Y_labeled = lab_data$Y,
                       EYx = EYx,
                       alpha = 0.05) {
  n <- length(Y_labeled)
  N <- length(EYx) - n
  mean_resid <- mean(Y_labeled - EYx[1:n])
  # average of imputations on all data
  mean_imp <- mean(EYx)
  
  beta_oracle <- mean_resid + mean_imp
  
  r. <- c(rep(1, n), rep(0, N))
  pii <- n / (n + N)
  y <- c(Y_labeled, rep(0, N))
  # influence function
  phi <- r. * (y - EYx) / pii + EYx - beta_oracle
  stderr_oracle <- sqrt(mean(phi ^ 2)) / sqrt(n + N)
  
  P_oracle <- 2 * pnorm(abs(beta_oracle / stderr_oracle), 0, 1, lower.tail = FALSE)
  upper <- beta_oracle + qnorm(1 - alpha/2) * stderr_oracle
  lower <- beta_oracle - qnorm(1 - alpha/2) * stderr_oracle
  
  return(list(theta = beta_oracle, se = stderr_oracle, 
              lower_ci = lower, upper_ci = upper, P = P_oracle))
}

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
  X <- as.matrix(X)
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
  X_labeled <- as.matrix(X_labeled)
  X_unlabeled <- as.matrix(X_unlabeled)
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
  X_labeled <- as.matrix(X_labeled)
  X_unlabeled <- as.matrix(X_unlabeled)
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

oracle_ols_asymptotic <- function(X_labeled, X_unlabeled, 
                               Y_labeled, Yhat_labeled, Yhat_unlabeled, 
                               alpha, sandwich = TRUE) {
  X_labeled <- as.matrix(X_labeled)
  X_unlabeled <- as.matrix(X_unlabeled)
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

  X_labeled <- as.matrix(X_labeled)
  X_unlabeled <- as.matrix(X_unlabeled)
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