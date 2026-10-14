################################################################################
### Fast Approximate Spike and Slab (FASS) Sampler
################################################################################

#' Fast Approximate Stochastic Search (FASS) for Variable Selection
#'
#' Runs the FASS stochastic model search of Fontana and Sanso for linear
#' regression under continuous spike-and-slab priors. Candidate moves (adding or
#' removing one covariate) are scored with a first-order surrogate for the
#' posterior odds, built from the OLS fit of the full model. After a one-time
#' pre-computation of \eqn{Q = n (X^T X)^{-1}}, \eqn{\hat\beta} and the residual
#' variance, every search iteration costs \eqn{O(p)}. The highest-odds model
#' visited can then be passed to a post-hoc posterior approximation.
#'
#' @section Prior parameterisation:
#' \code{tau0} and \code{tau1} are inverse-scale parameters of the spike and
#' slab respectively:
#' \itemize{
#'   \item \code{"Gaussian"}: precision, \eqn{\beta_j \sim N(0, 1/\tau)}.
#'   \item \code{"Student"}: \eqn{t_\nu} with scale \eqn{1/\sqrt{\tau}}.
#'   \item \code{"Laplace"}: rate, density \eqn{(\tau/2) \exp(-\tau|\beta_j|)}.
#' }
#' In the notation of the paper (where the spike and slab have inverse scales
#' \eqn{\tau_0, \tau_1}), \code{tau} corresponds to \eqn{\tau^2} for the Gaussian
#' and Student-t families and to \eqn{\tau} for the Laplace family. The settings
#' used in the paper, \eqn{\tau_0 = n^{1/2}, \tau_1 = 1}, are therefore
#' \code{tau0 = n, tau1 = 1} (Gaussian, Student) and \code{tau0 = sqrt(n), tau1 = 1}
#' (Laplace). The coefficient priors are not scaled by \eqn{\sigma^2}.
#'
#' @param y Numeric response vector of length n
#' @param X Numeric n x p covariate matrix (n > p, full column rank). Do not
#'   include a column of ones; use \code{intercept = TRUE} instead.
#' @param niter Integer, number of stochastic search iterations
#' @param a_pi,b_pi Shape parameters of the Beta(a_pi, b_pi) prior on the
#'   inclusion probability (default 1, 1)
#' @param family Spike-and-slab family: \code{"Laplace"}, \code{"Gaussian"}, or
#'   \code{"Student"}. Student-t priors are robust to model superinduction.
#' @param tau1 Slab inverse-scale parameter (see Prior parameterisation)
#' @param tau0 Spike inverse-scale parameter (see Prior parameterisation)
#' @param nu Degrees of freedom when \code{family = "Student"}; values of 5 or
#'   more are recommended
#' @param full_posterior Post-hoc inference for the selected model: \code{NA}
#'   (selection only, the default), \code{"Gibbs"}, \code{"EP"},
#'   \code{"Block-VB"}, \code{"MMFVB"}, \code{"MFVB-CAVI"}, or \code{"LA"}
#' @param nburn Burn-in iterations for \code{full_posterior = "Gibbs"} (default 1000)
#' @param niter2 Retained iterations for \code{full_posterior = "Gibbs"} (default 1000)
#' @param verbose Logical, whether to print progress every 500 iterations
#' @param intercept Logical, whether to include an (unpenalised) intercept (default TRUE)
#' @param rescale_y Logical, whether to divide y by sd(y) internally and
#'   rescale the posterior output back (default FALSE)
#' @return A list with
#'   \item{gamma}{\code{niter} x p 0/1 matrix; row i is the model after iteration i}
#'   \item{odds}{Vector of length \code{niter + 1} of approximate log posterior
#'     odds against the null model; \code{odds[1] = 0} is the null model and
#'     \code{odds[i + 1]} corresponds to \code{gamma[i, ]}}
#'   \item{posterior}{Post-hoc fit of the highest-odds model (\code{NA} if
#'     \code{full_posterior = NA}). For the approximation methods this is a list
#'     with \code{mu}, \code{Sigma}, \code{sigma2} and \code{alpha}; for
#'     \code{"Gibbs"} it holds draws \code{outbeta}, \code{outsigma} and
#'     \code{outalpha}. Coefficients refer to the selected columns only and the
#'     intercept is on the scale of the original (uncentred) \code{X}.}
#' @seealso \code{\link{fit_approx_posterior}}, \code{\link{predict_approx}}
#' @examples
#' set.seed(1)
#' n <- 2000; p <- 20
#' X <- matrix(rnorm(n * p), n, p)
#' y <- drop(1 + X[, 1:3] %*% c(2, -1.5, 1) + rnorm(n))
#' fit <- FASS.sampler(y, X, niter = 500, family = "Student", nu = 5,
#'                     tau1 = 1, tau0 = n, full_posterior = "Block-VB")
#' best <- which.max(fit$odds)
#' if (best > 1) which(fit$gamma[best - 1, ] == 1)
#' @export
FASS.sampler <- function(y, X, niter, a_pi = 1, b_pi = 1, family = c("Laplace", "Gaussian", "Student"),
                         tau1, tau0, nu = NA, full_posterior = NA, nburn = NA, niter2 = NA,
                         verbose = FALSE, intercept = TRUE, rescale_y = FALSE) {
  
  # 1. Input Validation
  data_val <- validate_regression_data(y, X, check_invertibility = TRUE)
  y <- data_val$y
  X <- data_val$X
  
  params_val <- validate_fass_params(niter, family, tau1, tau0, nu, a_pi, b_pi, 
                                     full_posterior, nburn, niter2)
  niter <- params_val$niter
  family <- params_val$family
  tau1 <- params_val$tau1
  tau0 <- params_val$tau0
  nu <- params_val$nu
  a_pi <- params_val$a_pi
  b_pi <- params_val$b_pi
  full_posterior <- params_val$full_posterior
  nburn <- params_val$nburn
  niter2 <- params_val$niter2
  
  if (rescale_y) {
    scale_y <- stats::sd(y)
    y <- y / scale_y
  } else {
    scale_y <- 1
  }
  
  start0 <- Sys.time()
  n <- length(y)
  m <- ncol(X)
  
  if (intercept) {
    x_mean <- colMeans(X)
    y_mean <- mean(y)
    X_work <- sweep(X, 2, x_mean, "-")
    y_work <- y - y_mean
  } else {
    x_mean <- rep(0, m)
    y_mean <- 0
    X_work <- X
    y_work <- y
  }
  
  ### pre-computes Q (n * (X^T X)^{-1}), beta.hat, and RSS
  X2 <- crossprod(X_work)
  R <- tryCatch(chol(X2), error = function(e) {
    stop("Cross-product matrix X^T X is not positive definite. Candidate features in 'X' appear collinear or rank-deficient.")
  })
  diag_R <- diag(R)
  rel_pivots <- diag_R / sqrt(diag(X2))
  if (any(rel_pivots < 1e-7 | is.na(rel_pivots))) {
    stop("Cross-product matrix X^T X is computationally singular (near-zero pivot). Candidate features in 'X' appear collinear or rank-deficient.")
  }
  Q <- n * chol2inv(R)
  X.y <- crossprod(X_work, y_work)
  beta.hat <- as.vector((1 / n) * Q %*% X.y)
  
  resid0 <- y_work - X_work %*% beta.hat
  df_sigma <- max(1, n - m - (if (intercept) 3 else 2))
  sigma.hat <- sum(resid0^2) / df_sigma
  alpha.hat <- if (intercept) (y_mean - sum(x_mean * beta.hat)) else 0
  rm(R)
  if (is.na(full_posterior)) {
    # Selection only: release the n x p copies before the search loop
    rm(X, X2, X.y, y, X_work, y_work)
  }
  
  gamma.store <- matrix(nrow = niter, ncol = m)
  gamma <- rep(0, m)
  
  if (verbose) {
    message(paste("Initialized in:", format(Sys.time() - start0)))
  }
  
  model.odds <- numeric(niter + 1)
  model.odds[1] <- 0
  
  tau <- rep(tau0, m)
  gamma.sum <- sum(gamma == 1)
  tau[which(gamma == 1)] <- tau1
  
  Q_diag <- diag(Q)
  scale_factor <- sigma.hat / (2 * n)
  
  # --- 1. Laplace Family ---
  if (family == "Laplace") {
    tau.star <- tau * sign(beta.hat)
    Q.tau1 <- tau.star %*% Q
    Q.tau <- as.numeric(Q.tau1 %*% tau.star)
    
    sign_beta <- sign(beta.hat)
    abs_beta <- abs(beta.hat)
    const_add <- (tau1^2 - tau0^2) * Q_diag
    const_rem <- (tau0^2 - tau1^2) * Q_diag
    d_tau_add <- tau1 - tau0
    d_tau_rem <- tau0 - tau1
    d_tau2_add <- tau1^2 - tau0^2
    d_tau2_rem <- tau0^2 - tau1^2
    
    for (i in 1:niter) {
      odds <- update.laplace(n, m, a_pi, b_pi, tau0, tau1, gamma, gamma.sum, beta.hat, sigma.hat, Q, Q.tau1, Q.tau, tau.star,
                             Q_diag = Q_diag, sign_beta = sign_beta, abs_beta = abs_beta,
                             const_add = const_add, const_rem = const_rem, scale_factor = scale_factor)
      update <- sampler(odds)
      
      if (update != length(odds)) {
        s_beta <- sign_beta[update]
        Q_uu <- Q_diag[update]
        
        if (gamma[update] == 0) {
          gamma[update] <- 1
          gamma.sum <- gamma.sum + 1
          diff_tau <- d_tau_add * s_beta
          Q.tau <- (Q.tau + 2 * Q.tau1[update] * diff_tau + 
                      d_tau2_add * Q_uu - 2 * tau.star[update] * Q_uu * diff_tau)
          tau.star[update] <- tau.star[update] + diff_tau 
          Q.tau1 <- Q.tau1 + diff_tau * Q[update, ]
        } else {
          gamma[update] <- 0
          gamma.sum <- gamma.sum - 1
          diff_tau <- d_tau_rem * s_beta
          Q.tau <- (Q.tau + 2 * Q.tau1[update] * diff_tau + 
                      d_tau2_rem * Q_uu - 2 * tau.star[update] * Q_uu * diff_tau)
          tau.star[update] <- tau.star[update] + diff_tau 
          Q.tau1 <- Q.tau1 + diff_tau * Q[update, ]
        }
      }
      
      model.odds[i + 1] <- model.odds[i] + odds[update]
      gamma.store[i, ] <- gamma
      
      if (verbose && i %% 500 == 0) {
        message(paste0(100 * round(i / niter, 4), "% | Active: ", gamma.sum))
      }
    }
  }
  
  # --- 2. Gaussian Family ---
  if (family == "Gaussian") {
    tau.star <- tau * beta.hat
    Q.tau1 <- tau.star %*% Q
    Q.tau <- as.numeric(Q.tau1 %*% tau.star - sum(tau * Q_diag))
    
    log_pi_ratio <- 0
    R <- 1 + Q.tau * scale_factor
    log.R <- log_correction(R)
    
    beta2 <- beta.hat^2
    const_add <- (tau1^2 - tau0^2) * Q_diag * beta2 - (tau1 - tau0) * Q_diag
    const_rem <- (tau0^2 - tau1^2) * Q_diag * beta2 - (tau0 - tau1) * Q_diag
    d_tau_add <- tau1 - tau0
    d_tau_rem <- tau0 - tau1
    d_tau2_add <- tau1^2 - tau0^2
    d_tau2_rem <- tau0^2 - tau1^2
    log_tau_ratio <- (log(tau1) - log(tau0)) / 2
    
    for (i in 1:niter) {
      odds <- update.gauss(n, m, a_pi, b_pi, tau0, tau1, gamma, gamma.sum, beta.hat, sigma.hat, Q, Q.tau1, Q.tau, tau.star,
                           Q_diag = Q_diag, beta2 = beta2, const_add = const_add, const_rem = const_rem, scale_factor = scale_factor)
      update <- sampler(odds)
      
      if (update != length(odds)) {
        b_up <- beta.hat[update]
        b2_up <- beta2[update]
        Q_uu <- Q_diag[update]
        log_pi_update <- log(a_pi + (gamma.sum - gamma[update])) - log(b_pi + m - (gamma.sum - gamma[update]) - 1) + 
          log_tau_ratio - d_tau_add * b2_up / 2
        
        if (gamma[update] == 0) {
          log_pi_ratio <- log_pi_ratio + log_pi_update
          gamma[update] <- 1
          gamma.sum <- gamma.sum + 1
          diff_tau_b <- d_tau_add * b_up
          Q.tau <- (Q.tau + 2 * Q.tau1[update] * diff_tau_b + 
                      d_tau2_add * Q_uu * b2_up - 
                      2 * tau.star[update] * Q_uu * diff_tau_b - d_tau_add * Q_uu)
          tau.star[update] <- tau.star[update] + diff_tau_b 
          Q.tau1 <- Q.tau1 + diff_tau_b * Q[update, ]
        } else {
          log_pi_ratio <- log_pi_ratio - log_pi_update
          gamma[update] <- 0
          gamma.sum <- gamma.sum - 1
          diff_tau_b <- d_tau_rem * b_up
          Q.tau <- (Q.tau + 2 * Q.tau1[update] * diff_tau_b + 
                      d_tau2_rem * Q_uu * b2_up - 
                      2 * tau.star[update] * Q_uu * diff_tau_b - d_tau_rem * Q_uu)
          tau.star[update] <- tau.star[update] + diff_tau_b 
          Q.tau1 <- Q.tau1 + diff_tau_b * Q[update, ]
        }
      }
      
      R.new <- 1 + Q.tau * scale_factor
      log.R.new <- log_correction(R.new)
      
      model.odds[i + 1] <- log_pi_ratio + log.R.new - log.R
      gamma.store[i, ] <- gamma
      
      if (verbose && i %% 500 == 0) {
        message(paste0(100 * round(i / niter, 4), "% | Active: ", gamma.sum))
      }
    }
  }
  
  # --- 3. Student-t Family ---
  if (family == "Student") {
    beta2 <- beta.hat^2
    tau.star <- as.vector((nu + 1) * beta.hat / (beta2 + nu / tau))
    Q.tau1 <- tau.star %*% Q
    Q.tau <- Q.tau1 %*% tau.star - sum(Q_diag * (nu + 1) * (nu / tau - beta2) / ((beta2 + nu / tau)^2))
    Q.tau <- as.numeric(Q.tau)
    
    log_pi_ratio <- 0
    R <- 1 + Q.tau * scale_factor
    log.R <- log_correction(R)
    
    temp1_prior_diff <- log(tau1) / 2 - log(tau0) / 2 - (nu + 1) / 2 * (log(1 + tau1 / nu * beta2) - log(1 + tau0 / nu * beta2))
    tau_star_add <- (nu + 1) * beta.hat / (beta2 + nu / tau1)
    tau_star_rem <- (nu + 1) * beta.hat / (beta2 + nu / tau0)
    
    u_diff_add <- (nu + 1) * (nu / tau1 - beta2) / ((beta2 + nu / tau1)^2) - 
      (nu + 1) * (nu / tau0 - beta2) / ((beta2 + nu / tau0)^2)
    u_diff_rem <- -u_diff_add
    
    for (i in 1:niter) {
      odds <- update.student(n, m, a_pi, b_pi, tau0, tau1, gamma, gamma.sum, beta.hat, sigma.hat, Q, Q.tau1, Q.tau, tau.star, nu,
                             Q_diag = Q_diag, beta2 = beta2, tau_star_add = tau_star_add, tau_star_rem = tau_star_rem,
                             u_diff_add = u_diff_add, u_diff_rem = u_diff_rem,
                             temp1_prior_diff = temp1_prior_diff, scale_factor = scale_factor)
      update <- sampler(odds)
      
      if (update != length(odds)) {
        log_pi_update <- log(a_pi + (gamma.sum - gamma[update])) - log(b_pi + m - (gamma.sum - gamma[update]) - 1) + 
          temp1_prior_diff[update]
        
        Q_uu <- Q_diag[update]
        if (gamma[update] == 0) {
          log_pi_ratio <- log_pi_ratio + log_pi_update
          gamma[update] <- 1
          gamma.sum <- gamma.sum + 1
          tau.star.new <- tau_star_add[update]
          u.diff <- u_diff_add[update]
        } else {
          log_pi_ratio <- log_pi_ratio - log_pi_update
          gamma[update] <- 0
          gamma.sum <- gamma.sum - 1
          tau.star.new <- tau_star_rem[update]
          u.diff <- u_diff_rem[update]
        }
        
        d_tau_s <- tau.star.new - tau.star[update]
        Q.tau <- (Q.tau + 2 * (Q.tau1[update] - tau.star[update] * Q_uu) * d_tau_s + 
                    (tau.star.new^2 - tau.star[update]^2) * Q_uu - u.diff * Q_uu)
        Q.tau1 <- Q.tau1 + d_tau_s * Q[update, ]
        tau.star[update] <- tau.star.new
      }
      
      R.new <- 1 + Q.tau * scale_factor
      log.R.new <- log_correction(R.new)
      
      model.odds[i + 1] <- log_pi_ratio + log.R.new - log.R
      gamma.store[i, ] <- gamma
      
      if (verbose && i %% 500 == 0) {
        message(paste0(100 * round(i / niter, 4), "% | Active: ", gamma.sum))
      }
    }
  }
  
  if (verbose) {
    message(paste0("Total Selection Time: ", format(Sys.time() - start0)))
  }
  
  # --- Post-Hoc Posterior Estimation ---
  posterior <- NA
  
  if (!is.na(full_posterior)) {
    best_idx <- which.max(model.odds)
    is_null_model <- (best_idx == 1)
    
    if (full_posterior == 'Gibbs') {
      if (is_null_model) {
        posterior <- Gibbs_posthoc_null(y, n, niter2, nburn, intercept)
      } else {
        gamma_selected <- as.logical(as.vector(gamma.store[best_idx - 1, ]))
        if (sum(gamma_selected) == 0) {
          posterior <- Gibbs_posthoc_null(y, n, niter2, nburn, intercept)
        } else {
          X_sub_work <- X_work[, gamma_selected, drop = FALSE]
          p_sub <- ncol(X_sub_work)
          gram_sub <- X2[gamma_selected, gamma_selected, drop = FALSE]
          X.y_sub <- X.y[gamma_selected]
          posterior <- Gibbs_posthoc(y, X_sub_work, X.y_sub, gram_sub, n, p_sub, niter2, nburn, 
                                     intercept = intercept, tau1 = tau1, family = family, nu = nu)
          # Gibbs_posthoc samples the intercept of the centred design; convert
          # it to the intercept for the original (uncentred) covariates
          if (intercept) {
            posterior$outalpha <- posterior$outalpha -
              as.vector(posterior$outbeta %*% x_mean[gamma_selected])
          }
        }
      }
    } else if (full_posterior %in% c('EP', 'Block-VB', 'MMFVB', 'MFVB-CAVI', 'LA')) {
      if (is_null_model) {
        posterior <- list(mu = numeric(0), Sigma = matrix(0, 0, 0), sigma2 = sum((y - y_mean)^2) / max(1, n - 1),
                          alpha = if (intercept) y_mean else 0, alpha.out = if (intercept) y_mean else 0,
                          sigma2.out = sum((y - y_mean)^2) / max(1, n - 1))
      } else {
        gamma_selected <- as.logical(as.vector(gamma.store[best_idx - 1, ]))
        if (sum(gamma_selected) == 0) {
          posterior <- list(mu = numeric(0), Sigma = matrix(0, 0, 0), sigma2 = sum((y - y_mean)^2) / max(1, n - 1),
                            alpha = if (intercept) y_mean else 0, alpha.out = if (intercept) y_mean else 0,
                            sigma2.out = sum((y - y_mean)^2) / max(1, n - 1))
        } else {
          X_sub_work <- X_work[, gamma_selected, drop = FALSE]
          p_sub <- ncol(X_sub_work)
          gram_sub <- X2[gamma_selected, gamma_selected, drop = FALSE]
          X.y_sub <- X.y[gamma_selected]
          beta_sub <- as.vector(chol2inv(chol(gram_sub)) %*% X.y_sub)
          
          resid <- y_work - X_sub_work %*% beta_sub
          sigma2.hat <- sum(resid^2) / max(1, n - p_sub - (if (intercept) 1 else 0))
          
          prior_map <- c("Gaussian" = "gaussian", "Laplace" = "laplace", "Student" = "student_t")
          posterior <- fit_approx_posterior(y_work, X_sub_work, method = full_posterior,
                                            prior_type = prior_map[[family]],
                                            hyper_params = list(tau = tau1, nu = nu),
                                            intercept = FALSE,
                                            sigma2.hat = sigma2.hat,
                                            XtX = gram_sub,
                                            X.y = X.y_sub)
          if (intercept) {
            posterior$alpha <- y_mean - sum(x_mean[gamma_selected] * posterior$mu)
          } else {
            posterior$alpha <- 0
          }
          posterior$alpha.out <- posterior$alpha
        }
      }
    }
  }
  
  # Rescale back to original y scale if rescaled
  if (rescale_y && is.list(posterior)) {
    if (identical(full_posterior, 'Gibbs')) {
      if (!is.null(posterior$outbeta)) posterior$outbeta <- posterior$outbeta * scale_y
      if (!is.null(posterior$outalpha)) posterior$outalpha <- posterior$outalpha * scale_y
      if (!is.null(posterior$outsigma)) posterior$outsigma <- posterior$outsigma * (scale_y^2)
    } else {
      if (!is.null(posterior$mu)) posterior$mu <- posterior$mu * scale_y
      if (!is.null(posterior$Sigma)) posterior$Sigma <- posterior$Sigma * (scale_y^2)
      if (!is.null(posterior$sigma2)) posterior$sigma2 <- posterior$sigma2 * (scale_y^2)
      if (!is.null(posterior$alpha)) posterior$alpha <- posterior$alpha * scale_y
      if (!is.null(posterior$alpha.out)) posterior$alpha.out <- posterior$alpha.out * scale_y
      if (!is.null(posterior$sigma2.out)) posterior$sigma2.out <- posterior$sigma2.out * (scale_y^2)
    }
  }
  
  return(list(gamma = gamma.store, odds = model.odds, posterior = posterior))
}
