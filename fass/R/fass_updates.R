################################################################################
### FASS Support and Log-Odds Update Functions (Internal)
################################################################################

#' Sample from Inverse Gaussian Distribution
#'
#' Helper function to draw random variates from an Inverse-Gaussian distribution
#' using Michael, Schucany, and Haas (1976) transformation.
#'
#' @param n Integer, number of observations to draw
#' @param mu Mean parameter (positive numeric or vector)
#' @param lambda Shape parameter (positive numeric or vector)
#' @return Numeric vector of length n of Inverse-Gaussian draws
#' @noRd
rinversegauss <- function(n, mu, lambda) {
  if (length(mu) != n) mu <- rep(mu, length.out = n)
  if (length(lambda) != n) lambda <- rep(lambda, length.out = n)
  mu <- pmin(pmax(mu, 1e-12), 1e7)
  z <- stats::rnorm(n, 0, 1)
  u <- z^2
  x <- mu + mu^2 * u / (2 * lambda) - mu / (2 * lambda) * sqrt(4 * mu * lambda * u + mu^2 * u^2)
  x[is.na(x) | x <= 0] <- lambda[is.na(x) | x <= 0] / pmax(1e-12, u[is.na(x) | x <= 0])
  
  v <- stats::runif(n, 0, 1)
  idx <- which(v > mu / (mu + x))
  if (length(idx) != 0) x[idx] <- mu[idx]^2 / x[idx]
  return(x)
}

#' Log of the First-Order Correction Factor with Zeroth-Order Fallback
#'
#' The first-order surrogate multiplies the prior ratio by a correction factor
#' \code{1 + E(sigma^2 | y) / (2n) * A}, which is not guaranteed to be positive
#' for the Gaussian and Student-t families. Following the paper, a non-positive
#' (or undefined) factor is replaced by the zeroth-order approximation, i.e. the
#' factor is set to 1 and its log contribution to 0. The same rule is used when
#' tracking the running model odds in \code{FASS.sampler}, so the proposal
#' probabilities and the recorded odds always agree.
#'
#' @param x Numeric vector of correction factors
#' @return Numeric vector of log correction factors
#' @noRd
log_correction <- function(x) {
  out <- numeric(length(x))
  ok <- !is.na(x) & x > 0
  out[ok] <- log(x[ok])
  out
}

#' Log-Odds Update Step for Gaussian Spike-and-Slab
#'
#' Computes the vector of candidate model log-odds for adding or deleting
#' each predictor under Gaussian priors.
#'
#' @param n Integer, sample size
#' @param m Integer, number of candidate covariates
#' @param a_pi Prior Beta shape parameter for inclusion
#' @param b_pi Prior Beta shape parameter for inclusion
#' @param tau0 Spike prior precision (Gaussian/Student-t) or rate (Laplace)
#' @param tau1 Slab prior precision (Gaussian/Student-t) or rate (Laplace)
#' @param gamma Binary indicator vector of length m
#' @param gamma.sum Current number of active features (sum(gamma))
#' @param beta.hat OLS coefficient estimate vector (1/n * Q %*% X^T y)
#' @param sigma.hat Residual variance estimate
#' @param Q Gram inverse matrix n * (X^T X)^{-1}
#' @param Q.tau1 Precomputed vector tau.star %*% Q
#' @param Q.tau Surrogate term A = z' Q z - U for the current model
#' @param tau.star Current z vector (e.g. tau * beta.hat for the Gaussian family)
#' @return Numeric vector of candidate log-odds (length m + 1, terminating with 0 for no change)
#' @noRd
update.gauss <- function(n, m, a_pi, b_pi, tau0, tau1, gamma, gamma.sum, beta.hat, sigma.hat, Q, Q.tau1, Q.tau, tau.star,
                         Q_diag = NULL, beta2 = NULL, const_add = NULL, const_rem = NULL, scale_factor = NULL) {
  Q.tau <- as.numeric(Q.tau)
  if (is.null(scale_factor)) scale_factor <- sigma.hat / (2 * n)
  log_temp2 <- log_correction(1 + Q.tau * scale_factor)
  if (is.null(Q_diag)) Q_diag <- diag(Q)
  Q_tau1_vec <- as.vector(Q.tau1)
  if (is.null(beta2)) beta2 <- beta.hat^2
  
  temp1 <- log(a_pi + (gamma.sum - gamma)) - log(b_pi + m - (gamma.sum - gamma) - 1) + 
    log(tau1) / 2 - log(tau0) / 2 - (tau1 - tau0) * beta2 / 2
  
  if (is.null(const_add)) const_add <- (tau1^2 - tau0^2) * Q_diag * beta2 - (tau1 - tau0) * Q_diag
  if (is.null(const_rem)) const_rem <- (tau0^2 - tau1^2) * Q_diag * beta2 - (tau0 - tau1) * Q_diag
  
  # For gamma == 0: adding feature
  t3_add <- 1 + (Q.tau + 2 * (Q_tau1_vec - tau.star * Q_diag) * (tau1 - tau0) * beta.hat + const_add) * scale_factor
  odds_add <- temp1 + log_correction(t3_add) - log_temp2
  
  # For gamma == 1: removing feature
  t3_rem <- 1 + (Q.tau + 2 * (Q_tau1_vec - tau.star * Q_diag) * (tau0 - tau1) * beta.hat + const_rem) * scale_factor
  odds_rem <- -temp1 - log_temp2 + log_correction(t3_rem)
  
  odds <- ifelse(gamma == 0, odds_add, odds_rem)
  odds[is.na(odds)] <- temp1[is.na(odds)] * ifelse(gamma[is.na(odds)] == 0, 1, -1)
  return(c(odds, 0))
}

#' Log-Odds Update Step for Laplace Spike-and-Slab
#'
#' Computes the vector of candidate model log-odds for adding or deleting
#' each predictor under Laplace (L1) priors.
#'
#' @inheritParams update.gauss
#' @return Numeric vector of candidate log-odds (length m + 1, terminating with 0 for no change)
#' @noRd
update.laplace <- function(n, m, a_pi, b_pi, tau0, tau1, gamma, gamma.sum, beta.hat, sigma.hat, Q, Q.tau1, Q.tau, tau.star,
                           Q_diag = NULL, sign_beta = NULL, abs_beta = NULL, const_add = NULL, const_rem = NULL, scale_factor = NULL) {
  Q.tau <- as.numeric(Q.tau)
  if (is.null(scale_factor)) scale_factor <- sigma.hat / (2 * n)
  log_temp2 <- log_correction(1 + Q.tau * scale_factor)
  if (is.null(Q_diag)) Q_diag <- diag(Q)
  Q_tau1_vec <- as.vector(Q.tau1)
  if (is.null(sign_beta)) sign_beta <- sign(beta.hat)
  if (is.null(abs_beta)) abs_beta <- abs(beta.hat)
  
  temp1 <- log(a_pi + (gamma.sum - gamma)) - log(b_pi + m - (gamma.sum - gamma) - 1) + 
    log(tau1) - log(tau0) - (tau1 - tau0) * abs_beta
  
  if (is.null(const_add)) const_add <- (tau1^2 - tau0^2) * Q_diag
  if (is.null(const_rem)) const_rem <- (tau0^2 - tau1^2) * Q_diag
  
  # For gamma == 0: adding feature
  t3_add <- 1 + (Q.tau + 2 * (Q_tau1_vec - tau.star * Q_diag) * (tau1 - tau0) * sign_beta + const_add) * scale_factor
  odds_add <- temp1 + log_correction(t3_add) - log_temp2
  
  # For gamma == 1: removing feature
  t3_rem <- 1 + (Q.tau + 2 * (Q_tau1_vec - tau.star * Q_diag) * (tau0 - tau1) * sign_beta + const_rem) * scale_factor
  odds_rem <- -temp1 - log_temp2 + log_correction(t3_rem)
  
  odds <- ifelse(gamma == 0, odds_add, odds_rem)
  return(c(odds, 0))
}

#' Log-Odds Update Step for Student-t Spike-and-Slab
#'
#' Computes the vector of candidate model log-odds for adding or deleting
#' each predictor under Student-t priors with nu degrees of freedom.
#'
#' @inheritParams update.gauss
#' @param nu Degrees of freedom for Student-t prior
#' @return Numeric vector of candidate log-odds (length m + 1, terminating with 0 for no change)
#' @noRd
update.student <- function(n, m, a_pi, b_pi, tau0, tau1, gamma, gamma.sum, beta.hat, sigma.hat, Q, Q.tau1, Q.tau, tau.star, nu,
                           Q_diag = NULL, beta2 = NULL, tau_star_add = NULL, tau_star_rem = NULL,
                           u_diff_add = NULL, u_diff_rem = NULL, temp1_prior_diff = NULL, scale_factor = NULL) {
  Q.tau <- as.numeric(Q.tau)
  if (is.null(scale_factor)) scale_factor <- sigma.hat / (2 * n)
  temp <- -log_correction(1 + Q.tau * scale_factor)
  
  if (is.null(Q_diag)) Q_diag <- diag(Q)
  Q_tau1_vec <- as.vector(Q.tau1)
  if (is.null(beta2)) beta2 <- beta.hat^2
  
  if (is.null(temp1_prior_diff)) {
    temp1_prior_diff <- log(tau1) / 2 - log(tau0) / 2 - (nu + 1) / 2 * (log(1 + tau1 / nu * beta2) - log(1 + tau0 / nu * beta2))
  }
  temp1 <- log(a_pi + (gamma.sum - gamma)) - log(b_pi + m - (gamma.sum - gamma) - 1) + temp1_prior_diff
  
  if (is.null(tau_star_add)) tau_star_add <- (nu + 1) * beta.hat / (beta2 + nu / tau1)
  if (is.null(tau_star_rem)) tau_star_rem <- (nu + 1) * beta.hat / (beta2 + nu / tau0)
  
  if (is.null(u_diff_add)) {
    u_diff_add <- (nu + 1) * (nu / tau1 - beta2) / ((beta2 + nu / tau1)^2) - 
      (nu + 1) * (nu / tau0 - beta2) / ((beta2 + nu / tau0)^2)
  }
  if (is.null(u_diff_rem)) {
    u_diff_rem <- (nu + 1) * (nu / tau0 - beta2) / ((beta2 + nu / tau0)^2) - 
      (nu + 1) * (nu / tau1 - beta2) / ((beta2 + nu / tau1)^2)
  }
  
  t2_add <- 1 + (Q.tau + 2 * (Q_tau1_vec - tau.star * Q_diag) * (tau_star_add - tau.star) + 
                   (tau_star_add^2 - tau.star^2) * Q_diag - u_diff_add * Q_diag) * scale_factor
  t2_add <- log_correction(t2_add)
  odds_add <- temp + temp1 + t2_add
  
  t2_rem <- 1 + (Q.tau + 2 * (Q_tau1_vec - tau.star * Q_diag) * (tau_star_rem - tau.star) + 
                   (tau_star_rem^2 - tau.star^2) * Q_diag - u_diff_rem * Q_diag) * scale_factor
  t2_rem <- log_correction(t2_rem)
  odds_rem <- temp - temp1 + t2_rem
  
  odds <- ifelse(gamma == 0, odds_add, odds_rem)
  odds[is.na(odds)] <- temp1[is.na(odds)] * ifelse(gamma[is.na(odds)] == 0, 1, -1)
  return(c(odds, 0))
}

#' Multinomial Proposal Step
#'
#' Samples the next candidate move from normalized exponential model odds.
#'
#' @param odds Vector of relative candidate log-odds
#' @return Integer index of the selected candidate move
#' @noRd
sampler <- function(odds) {
  prob <- exp(odds - max(odds))
  prob <- prob / sum(prob)
  sample.int(length(prob), 1L, prob = prob)
}

#' Post-Hoc Gibbs Sampler for Null Model
#'
#' Samples intercept and residual variance when no predictors are selected.
#'
#' @param y Numeric response vector
#' @param n Sample size
#' @param niter2 Post-burn-in iterations
#' @param nburn Burn-in iterations
#' @param intercept Logical, whether to include intercept
#' @return List with outsigma, outalpha, times
#' @noRd
Gibbs_posthoc_null <- function(y, n, niter2, nburn, intercept = TRUE) {
  outsigma <- c()
  outalpha <- c()
  sigma2 <- 1
  times2 <- c()
  
  for (itr in 1:(nburn + niter2)) {
    start <- Sys.time()
    
    if (intercept == TRUE) alpha <- stats::rnorm(1, mean(y), sqrt(sigma2 / n))
    
    if (intercept == TRUE) {
      sigma2 <- 1 / stats::rgamma(1, n / 2, sum((y - alpha)^2) / 2)  
    } else {
      sigma2 <- 1 / stats::rgamma(1, n / 2, sum(y^2) / 2)
    }
    
    times2 <- c(times2, difftime(Sys.time(), start, units = "secs"))
    
    if (itr > nburn) {
      outsigma <- c(outsigma, sigma2)
      if (intercept == TRUE) outalpha <- c(outalpha, alpha)
    }
  }
  
  posterior <- list(outsigma = outsigma, outalpha = outalpha, times = times2)
  return(posterior)
}

#' Post-Hoc Gibbs Sampler for Selected Submodel
#'
#' Runs MCMC on the submodel chosen by shotgun search to draw posterior coefficients.
#'
#' @param y Numeric response vector
#' @param X Submodel design matrix (selected columns)
#' @param X.y Crossproduct X^T y
#' @param gram Crossproduct matrix X^T X
#' @param n Sample size
#' @param p Number of active predictors
#' @param niter2 Post-burn-in iterations
#' @param nburn Burn-in iterations
#' @param intercept Logical, whether model includes intercept
#' @param tau1 Slab precision
#' @param family Slab prior family
#' @param nu Degrees of freedom when family is "Student"
#' @return List with outbeta, outW, outsigma, outalpha, times
#' @noRd
Gibbs_posthoc <- function(y, X, X.y, gram, n, p, niter2, nburn, intercept = TRUE, tau1 = 1, family = "Gaussian", nu = NA) {
  outbeta <- matrix(nrow = niter2, ncol = p)
  out.W <- matrix(nrow = niter2, ncol = p)
  outsigma <- numeric(niter2)
  outalpha <- if (intercept) numeric(niter2) else NULL
  beta <- rep(0, p)
  sigma2 <- 1
  W <- rep(1, p)
  times2 <- c()
  
  for (itr in 1:(nburn + niter2)) {
    start <- Sys.time()
    
    if (intercept == TRUE) alpha <- stats::rnorm(1, mean(y), sqrt(sigma2 / n))
    
    u <- stats::rnorm(p, 0, 1)
    prec_diag <- if (family == "Laplace") (tau1^2 / W) else if (family == "Gaussian") rep(tau1, p) else (tau1 / W)
    L.t <- chol(gram / sigma2 + diag(prec_diag, nrow = p))
    v <- forwardsolve(t(L.t), X.y / sigma2)
    beta <- as.vector(rowSums(backsolve(L.t, cbind(v, u))))
    
    if (intercept == TRUE) {
      sigma2 <- 1 / stats::rgamma(1, n / 2, sum((y - X %*% beta - alpha)^2) / 2)  
    } else {
      sigma2 <- 1 / stats::rgamma(1, n / 2, sum((y - X %*% beta)^2) / 2)
    }
    
    if (family == 'Laplace') W <- 1 / rinversegauss(p, 1 / (tau1 * pmax(1e-12, abs(beta))), 1)
    if (family == 'Student') W <- 1 / stats::rgamma(p, (nu + 1) / 2, (nu + tau1 * beta^2) / 2)
    
    times2 <- c(times2, difftime(Sys.time(), start, units = "secs"))
    
    if (itr > nburn) {
      out.W[(itr - nburn), ] <- W
      outbeta[(itr - nburn), ] <- beta
      outsigma[itr - nburn] <- sigma2
      if (intercept == TRUE) outalpha[itr - nburn] <- alpha
    }
  }
  
  posterior <- list(outbeta = outbeta, outW = out.W, outsigma = outsigma, outalpha = outalpha, times = times2)
  return(posterior)
}
