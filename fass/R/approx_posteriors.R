################################################################################
### Unified Approximate Posteriors (EP, Variational Bayes, Laplace Approximation)
################################################################################

# ==============================================================================
# Shared Helpers (Internal)
# ==============================================================================

#' Gauss-Hermite Quadrature Nodes and Weights
#'
#' Computes nodes and weights via the Golub-Welsch tridiagonal eigenvalue algorithm.
#'
#' @param n Integer, number of quadrature nodes
#' @return List with nodes (numeric vector) and weights (numeric vector)
#' @noRd
gauss_hermite_nodes_weights <- function(n) {
  i <- 1:(n - 1)
  b <- sqrt(i / 2)
  J <- matrix(0, nrow = n, ncol = n)
  diag(J[-1, -n]) <- b
  diag(J[-n, -1]) <- b
  
  eig <- eigen(J, symmetric = TRUE)
  idx <- order(eig$values)
  nodes <- eig$values[idx]
  weights <- (eig$vectors[1, idx]^2) * sqrt(pi)
  
  return(list(nodes = nodes, weights = weights))
}

#' Shared EM Residual Variance Update
#'
#' Calculates the updated residual variance estimate from residuals and trace term.
#'
#' @param y Numeric response vector
#' @param X Numeric design matrix
#' @param mu_beta Posterior mean coefficient vector
#' @param trace_term Numeric scalar trace term Tr(X^T X Sigma)
#' @param sigma2 Current residual variance
#' @param n Sample size
#' @param intercept Logical, whether model includes intercept
#' @param include_prior_sigma2 Logical, whether to include sigma2 in numerator (default TRUE)
#' @return Updated sigma2 estimate
#' @noRd
calc_sigma2_em_update <- function(y, X, mu_beta, trace_term, sigma2, n, intercept, include_prior_sigma2 = TRUE) {
  resid <- if (intercept) (y - X %*% mu_beta - mean(y)) else (y - X %*% mu_beta)
  rss <- sum(resid^2)
  if (intercept && include_prior_sigma2) {
    (rss + trace_term + sigma2) / (n + 2)
  } else {
    (rss + trace_term) / (n + 2)
  }
}

# ==============================================================================
# EP Moment Matching Helpers (Internal)
# ==============================================================================

#' @noRd
match_moments_gaussian <- function(mu_cav, v_cav, hyper) {
  tau <- hyper$tau
  v_prior <- 1 / tau
  prec_new <- 1 / v_cav + 1 / v_prior
  v_new <- 1 / prec_new
  mu_new <- v_new * (mu_cav / v_cav)
  return(list(mu = mu_new, var = v_new))
}

#' @noRd
match_moments_laplace <- function(mu_cav, v_cav, hyper) {
  tau <- hyper$tau
  lam <- tau
  sigma <- sqrt(v_cav)
  
  mu_plus <- mu_cav - v_cav * lam
  mu_minus <- mu_cav + v_cav * lam
  
  alpha_plus <- (0 - mu_plus) / sigma
  alpha_minus <- (0 - mu_minus) / sigma 
  
  const_plus <- -lam * mu_cav + 0.5 * v_cav * lam^2
  const_minus <- lam * mu_cav + 0.5 * v_cav * lam^2
  
  log_Z_plus <- const_plus + stats::pnorm(alpha_plus, lower.tail = FALSE, log.p = TRUE)
  log_Z_minus <- const_minus + stats::pnorm(alpha_minus, lower.tail = TRUE, log.p = TRUE)
  
  max_log_Z <- max(log_Z_plus, log_Z_minus)
  log_Z_total <- max_log_Z + log(exp(log_Z_plus - max_log_Z) + exp(log_Z_minus - max_log_Z))
  
  w_plus <- exp(log_Z_plus - log_Z_total)
  w_minus <- exp(log_Z_minus - log_Z_total)
  
  ratio_plus <- exp(stats::dnorm(alpha_plus, log = TRUE) - stats::pnorm(alpha_plus, lower.tail = FALSE, log.p = TRUE))
  E_plus <- mu_plus + sigma * ratio_plus
  var_trunc_plus <- v_cav * (1 + alpha_plus * ratio_plus - ratio_plus^2)
  E2_plus <- var_trunc_plus + E_plus^2
  
  ratio_minus <- -exp(stats::dnorm(alpha_minus, log = TRUE) - stats::pnorm(alpha_minus, lower.tail = TRUE, log.p = TRUE))
  E_minus <- mu_minus + sigma * ratio_minus
  var_trunc_minus <- v_cav * (1 + alpha_minus * ratio_minus - ratio_minus^2)
  E2_minus <- var_trunc_minus + E_minus^2
  
  mu_new <- w_plus * E_plus + w_minus * E_minus
  E2_new <- w_plus * E2_plus + w_minus * E2_minus
  v_new <- E2_new - mu_new^2
  
  if (v_new < 1e-10) v_new <- 1e-10
  return(list(mu = mu_new, var = v_new))
}

#' @noRd
match_moments_student_AGHQ <- function(mu_cav, v_cav, hyper) {
  gh_nodes <- hyper$gh_nodes
  gh_weights <- hyper$gh_weights
  
  nu_tau2 <- hyper$nu / hyper$tau
  nu_plus_1 <- hyper$nu + 1
  
  log_f <- function(x) {
    -0.5 * (x - mu_cav)^2 / v_cav - (nu_plus_1 / 2) * log(1 + x^2 / nu_tau2)
  }
  
  x_mode <- mu_cav 
  log_val <- log_f(x_mode)
  converged <- FALSE
  
  for (i in 1:15) {
    grad <- -(x_mode - mu_cav) / v_cav - (nu_plus_1 * x_mode) / (nu_tau2 + x_mode^2)
    hess <- -1 / v_cav - nu_plus_1 * (nu_tau2 - x_mode^2) / (nu_tau2 + x_mode^2)^2
    
    if (hess >= 0) hess <- -1 / v_cav 
    
    step <- grad / hess
    x_new <- x_mode - step
    
    for (j in 1:5) {
      log_val_new <- log_f(x_new)
      if (log_val_new >= log_val) break
      step <- step / 2
      x_new <- x_mode - step
    }
    
    if (abs(x_new - x_mode) < 1e-6) {
      x_mode <- x_new
      converged <- TRUE
      break
    }
    
    x_mode <- x_new
    log_val <- log_val_new
  }
  
  hess_final <- -1 / v_cav - nu_plus_1 * (nu_tau2 - x_mode^2) / (nu_tau2 + x_mode^2)^2
  
  if (!converged || hess_final >= 0) {
    sd_cav <- sqrt(v_cav)
    lower_b <- min(0, mu_cav) - 3 * sd_cav
    upper_b <- max(0, mu_cav) + 3 * sd_cav
    opt <- stats::optimize(log_f, interval = c(lower_b, upper_b), maximum = TRUE)
    x_mode <- opt$maximum
    hess_final <- -1 / v_cav - nu_plus_1 * (nu_tau2 - x_mode^2) / (nu_tau2 + x_mode^2)^2
    if (hess_final >= 0) hess_final <- -1 / v_cav
  }
  
  var_mode <- -1 / hess_final
  if (var_mode > 100 * v_cav) var_mode <- 100 * v_cav
  
  x_eval <- x_mode + sqrt(2 * var_mode) * gh_nodes
  log_tilted <- -0.5 * (x_eval - mu_cav)^2 / v_cav - (nu_plus_1 / 2) * log(1 + (x_eval^2) / nu_tau2)
  log_integrand <- log_tilted + gh_nodes^2
  max_log <- max(log_integrand)
  w_tilde <- gh_weights * exp(log_integrand - max_log)
  
  Z <- sum(w_tilde)
  mu_new <- sum(w_tilde * x_eval) / Z
  E_x2 <- sum(w_tilde * (x_eval^2)) / Z
  v_new <- E_x2 - mu_new^2
  
  if (is.na(v_new) || v_new <= 1e-10) v_new <- 1e-10
  if (is.na(mu_new)) mu_new <- mu_cav
  
  return(list(mu = mu_new, var = v_new))
}

# ==============================================================================
# 1. Expectation Propagation
# ==============================================================================

#' Expectation Propagation (EP) for Linear Regression
#'
#' Fits posterior distribution using Expectation Propagation with Gaussian,
#' Laplace, or Student-t priors and EM for residual variance estimation.
#' Most users should call \code{\link{fit_approx_posterior}}, which handles the
#' intercept; these lower-level functions assume \code{X} and \code{y} are
#' column-centred when \code{intercept = TRUE}.
#'
#' @param y Numeric response vector
#' @param X Numeric design matrix (column-centred if \code{intercept = TRUE})
#' @param sigma2.hat Initial residual variance estimate
#' @param XtX Crossproduct matrix X^T X
#' @param X.y Crossproduct vector X^T y
#' @param prior_type One of "gaussian", "laplace", "student_t"
#' @param hyper_params List with slab parameter \code{tau} (precision for
#'   Gaussian/Student-t, rate for Laplace) and Student-t degrees of freedom \code{nu}
#' @param opts Control parameters: max_iter, tol, damping, em_iter, tol_em, n_nodes
#' @param intercept Logical, whether model includes an intercept
#' @return List with mu, Sigma, sigma2, alpha, sites
#' @export
EP <- function(y, X, sigma2.hat, XtX, X.y, prior_type = c("gaussian", "laplace", "student_t"), 
               hyper_params = list(tau = 1, nu = 3), 
               opts = list(max_iter = 50, tol = 1e-4, damping = 0.8, em_iter = 10, tol_em = 1e-10, n_nodes = 21), 
               intercept = TRUE) {
  
  prior_type <- match.arg(prior_type)
  default_opts <- list(max_iter = 50, tol = 1e-4, damping = 0.8, em_iter = 10, tol_em = 1e-10, n_nodes = 21)
  opts <- if (missing(opts) || is.null(opts)) default_opts else utils::modifyList(default_opts, as.list(opts))
  
  matcher_func <- switch(prior_type,
                         "gaussian"  = match_moments_gaussian,
                         "laplace"   = match_moments_laplace,
                         "student_t" = match_moments_student_AGHQ)
  
  n <- nrow(X)
  p <- ncol(XtX)
  sigma2 <- sigma2.hat
  
  if (prior_type == "student_t") {
    gh_cache <- gauss_hermite_nodes_weights(opts$n_nodes)
    hyper_params$gh_nodes <- gh_cache$nodes
    hyper_params$gh_weights <- gh_cache$weights
  }
  
  site_prec <- rep(1e-6, p)
  site_mean_prec <- rep(0, p)
  mu_beta <- rep(0, p)
  Sigma_beta <- diag(p)
  
  for (outer in 1:opts$em_iter) { 
    for (iter in 1:opts$max_iter) {
      params_old <- site_prec
      Prec_Global <- XtX / sigma2 + diag(x = site_prec, nrow = p, ncol = p)
      
      Sigma_beta <- tryCatch(chol2inv(chol(Prec_Global)), 
                             error = function(e) chol2inv(chol(Prec_Global + diag(1e-6, p))))
      mu_beta <- Sigma_beta %*% (X.y / sigma2 + site_mean_prec)
      
      for (j in 1:p) {
        v_marg <- Sigma_beta[j, j]
        mu_marg <- mu_beta[j]
        denom <- 1 - v_marg * site_prec[j]
        if (denom < 1e-10) next
        
        v_cav <- v_marg / denom
        mu_cav <- (mu_marg - v_marg * site_mean_prec[j]) / denom
        
        moments <- matcher_func(mu_cav, v_cav, hyper_params)
        mu_new <- moments$mu
        v_new <- moments$var
        
        d_prec <- 1 / v_new - 1 / v_cav
        d_mean_prec <- mu_new / v_new - mu_cav / v_cav
        
        delta <- opts$damping
        site_prec[j] <- delta * d_prec + (1 - delta) * site_prec[j]
        site_mean_prec[j] <- delta * d_mean_prec + (1 - delta) * site_mean_prec[j]
        
        if (site_prec[j] < 1e-8) site_prec[j] <- 1e-8
      }
      
      if (max(abs(site_prec - params_old)) < opts$tol) break
    }
    
    trace_term <- sum(XtX * Sigma_beta)
    sigma2_new <- calc_sigma2_em_update(y, X, mu_beta, trace_term, sigma2, n, intercept, include_prior_sigma2 = TRUE)
    
    if (abs(sigma2_new - sigma2) < opts$tol_em) {
      sigma2 <- sigma2_new
      break
    }
    sigma2 <- sigma2_new
  }
  
  alpha_est <- if (intercept) mean(y) else 0
  return(list(mu = mu_beta, Sigma = Sigma_beta, sigma2 = sigma2, alpha = alpha_est,
              sites = list(prec = site_prec, mean_prec = site_mean_prec)))
}

# ==============================================================================
# 2. Variational Bayes Methods
# ==============================================================================

#' Block-Coordinate Variational Bayes (EM)
#'
#' Coordinate-ascent VB with a full-covariance Gaussian factor for \eqn{\beta}
#' and independent factors for the latent scales; \eqn{\sigma^2} by EM.
#'
#' @inheritParams EP
#' @return List with mu, Sigma, sigma2, alpha, omega
#' @export
VB <- function(y, X, sigma2.hat, XtX, X.y, prior_type = c("laplace", "student_t", "gaussian"), 
               hyper_params = list(tau = 1, nu = 3), 
               opts = list(max_iter = 500, tol = 1e-4, em_iter = 100, tol_em = 1e-10), 
               intercept = TRUE) {
  
  prior_type <- match.arg(prior_type)
  default_opts <- list(max_iter = 500, tol = 1e-4, em_iter = 100, tol_em = 1e-10)
  opts <- if (missing(opts) || is.null(opts)) default_opts else utils::modifyList(default_opts, as.list(opts))
  
  n <- nrow(X)
  p <- ncol(XtX)
  
  sigma2 <- sigma2.hat
  omega <- rep(1, p) 
  mu_beta <- rep(0, p)
  Sigma_beta <- diag(p)
  
  for (outer in 1:opts$em_iter) { 
    for (iter in 1:opts$max_iter) {
      omega_old <- omega
      W <- diag(x = omega, nrow = p, ncol = p)
      Prec_q <- XtX / sigma2 + W
      
      Sigma_beta <- tryCatch(chol2inv(chol(Prec_q)), 
                             error = function(e) chol2inv(chol(Prec_q + diag(1e-6, nrow = p, ncol = p))))
      mu_beta <- Sigma_beta %*% (X.y / sigma2)
      E_beta2 <- as.vector(mu_beta)^2 + diag(Sigma_beta)
      
      if (prior_type == "laplace") {
        tau1 <- hyper_params$tau
        omega <- tau1 / sqrt(E_beta2)
      } else if (prior_type == "student_t") {
        tau1 <- hyper_params$tau
        nu <- hyper_params$nu
        omega <- (nu + 1) / (nu / tau1 + E_beta2)
      } else if (prior_type == "gaussian") {
        tau1 <- hyper_params$tau
        omega <- rep(tau1, p)
      }
      
      omega[omega < 1e-10] <- 1e-10
      if (max(abs(omega - omega_old)) < opts$tol) break
    }
    
    trace_term <- sum(XtX * Sigma_beta)
    sigma2_new <- calc_sigma2_em_update(y, X, mu_beta, trace_term, sigma2, n, intercept, include_prior_sigma2 = TRUE)
    
    if (abs(sigma2_new - sigma2) < opts$tol_em) {
      sigma2 <- sigma2_new
      break
    }
    sigma2 <- sigma2_new
  }
  
  alpha_est <- if (intercept) mean(y) else 0
  return(list(mu = mu_beta, Sigma = Sigma_beta, sigma2 = sigma2, alpha = alpha_est, omega = omega))
}

#' Mean-Field Variational Bayes (CAVI)
#'
#' Mean-field coordinate-ascent VB with independent Gaussian factors for each
#' coefficient; \eqn{\sigma^2} by EM.
#'
#' @inheritParams EP
#' @return List with mu, Sigma (diagonal), sigma2, alpha, omega
#' @export
VB_MF <- function(y, X, sigma2.hat, XtX, X.y, prior_type = c("laplace", "student_t", "gaussian"), 
                  hyper_params = list(tau = 1, nu = 3), 
                  opts = list(max_iter = 500, tol = 1e-4, em_iter = 100, tol_em = 1e-10), 
                  intercept = TRUE) {
  
  prior_type <- match.arg(prior_type)
  default_opts <- list(max_iter = 500, tol = 1e-4, em_iter = 100, tol_em = 1e-10)
  opts <- if (missing(opts) || is.null(opts)) default_opts else utils::modifyList(default_opts, as.list(opts))
  
  n <- nrow(X)
  p <- ncol(XtX)
  
  sigma2 <- sigma2.hat
  omega <- rep(1, p) 
  mu_beta <- rep(0, p)
  s2_beta <- rep(1, p)
  XtX_diag <- diag(XtX)
  
  for (outer in 1:opts$em_iter) { 
    for (iter in 1:opts$max_iter) {
      omega_old <- omega
      
      for (j in 1:p) {
        s2_beta[j] <- 1.0 / (XtX_diag[j] / sigma2 + omega[j])
        partial_resid_j <- X.y[j] - sum(XtX[j, ] * mu_beta) + XtX_diag[j] * mu_beta[j]
        mu_beta[j] <- s2_beta[j] * (partial_resid_j / sigma2)
      }
      
      E_beta2 <- mu_beta^2 + s2_beta
      
      if (prior_type == "laplace") {
        tau1 <- hyper_params$tau
        omega <- tau1 / sqrt(E_beta2)
      } else if (prior_type == "student_t") {
        tau1 <- hyper_params$tau
        nu <- hyper_params$nu
        omega <- (nu + 1) / (nu / tau1 + E_beta2)
      } else if (prior_type == "gaussian") {
        tau1 <- hyper_params$tau
        omega <- rep(tau1, p)
      }
      
      omega[omega < 1e-10] <- 1e-10
      if (max(abs(omega - omega_old)) < opts$tol) break
    }
    
    trace_term <- sum(XtX_diag * s2_beta)
    sigma2_new <- calc_sigma2_em_update(y, X, mu_beta, trace_term, sigma2, n, intercept, include_prior_sigma2 = TRUE)
    
    if (abs(sigma2_new - sigma2) < opts$tol_em) {
      sigma2 <- sigma2_new
      break
    }
    sigma2 <- sigma2_new
  }
  
  alpha_est <- if (intercept) mean(y) else 0
  return(list(mu = mu_beta, Sigma = diag(x = s2_beta, nrow = p, ncol = p), sigma2 = sigma2, alpha = alpha_est, omega = omega))
}

#' Direct ELBO Variational Bayes
#'
#' Mean-field Gaussian VB that maximises the ELBO of the marginal (scale-mixture)
#' prior directly; \eqn{\sigma^2} by EM.
#'
#' @inheritParams EP
#' @return List with mu, Sigma (diagonal), sigma2, alpha
#' @export
VB_Direct <- function(y, X, sigma2.hat, XtX, X.y, prior_type = c("laplace", "student_t", "gaussian"), 
                      hyper_params = list(tau = 1, nu = 3), 
                      opts = list(max_iter = 500, tol = 1e-4, em_iter = 100, tol_em = 1e-10), 
                      intercept = TRUE) {
  
  prior_type <- match.arg(prior_type)
  default_opts <- list(max_iter = 500, tol = 1e-4, em_iter = 100, tol_em = 1e-10)
  opts <- if (missing(opts) || is.null(opts)) default_opts else utils::modifyList(default_opts, as.list(opts))
  
  n <- nrow(X)
  p <- ncol(XtX)
  
  sigma2 <- sigma2.hat
  mu_beta <- rep(0, p)
  sigma_beta <- rep(1, p)
  XtX_diag <- diag(XtX)
  
  if (prior_type == "laplace") {
    tau1 <- hyper_params$tau
    for (outer in 1:opts$em_iter) { 
      for (iter in 1:opts$max_iter) {
        mu_old <- mu_beta
        sigma_old <- sigma_beta
        
        for (j in 1:p) {
          a <- XtX_diag[j] / sigma2
          b <- (2 * tau1) * stats::dnorm(mu_beta[j] / sigma_beta[j])
          sigma_beta[j] <- (-b + sqrt(b^2 + 4 * a)) / (2 * a)
          
          sum_others <- sum(XtX[j, ] * mu_beta) - XtX_diag[j] * mu_beta[j]
          mu_beta[j] <- (X.y[j] - sum_others - (sigma2 * tau1) * (2 * stats::pnorm(mu_beta[j] / sigma_beta[j]) - 1)) / XtX_diag[j]
        }
        
        diff <- max(abs(mu_beta - mu_old), abs(sigma_beta - sigma_old))
        if (diff < opts$tol) break
      }
      
      trace_term <- sum(XtX_diag * sigma_beta^2) 
      sigma2_new <- calc_sigma2_em_update(y, X, mu_beta, trace_term, sigma2, n, intercept, include_prior_sigma2 = TRUE)
      
      if (abs(sigma2_new - sigma2) < opts$tol_em) {
        sigma2 <- sigma2_new
        break
      }
      sigma2 <- sigma2_new
    }
    
  } else if (prior_type == "gaussian") {
    tau1 <- hyper_params$tau
    for (outer in 1:opts$em_iter) { 
      for (iter in 1:opts$max_iter) {
        mu_old <- mu_beta
        for (j in 1:p) {
          sigma_beta[j] <- 1 / sqrt(XtX_diag[j] / sigma2 + tau1)
          sum_others <- sum(XtX[j, ] * mu_beta) - XtX_diag[j] * mu_beta[j]
          mu_beta[j] <- (sigma_beta[j]^2) * ((X.y[j] - sum_others) / sigma2)
        }
        diff <- max(abs(mu_beta - mu_old))
        if (diff < opts$tol) break
      }
      
      trace_term <- sum(XtX_diag * sigma_beta^2) 
      sigma2_new <- calc_sigma2_em_update(y, X, mu_beta, trace_term, sigma2, n, intercept, include_prior_sigma2 = TRUE)
      
      if (abs(sigma2_new - sigma2) < opts$tol_em) {
        sigma2 <- sigma2_new
        break
      }
      sigma2 <- sigma2_new
    }
    
  } else if (prior_type == "student_t") {
    tau1 <- hyper_params$tau
    nu <- hyper_params$nu
    
    gh <- gauss_hermite_nodes_weights(15)
    nodes <- gh$nodes
    weights <- gh$weights
    
    obj_fun <- function(par) {
      mu <- par[1:p]
      log_sig <- par[(p + 1):(2 * p)]
      sig <- exp(log_sig)
      
      resid <- y - X %*% mu
      if (intercept) resid <- resid - mean(y)
      ll <- -0.5 * sum(resid^2) / sigma2 - 0.5 * sum(XtX_diag * sig^2) / sigma2
      ent <- sum(log_sig)
      
      prior_val <- 0
      for (j in 1:p) {
        beta_jk <- mu[j] + sqrt(2) * sig[j] * nodes
        h_vals <- -(nu + 1) / 2 * log1p(beta_jk^2 * tau1 / nu)
        prior_val <- prior_val + (1 / sqrt(pi)) * sum(weights * h_vals)
      }
      return(-(ll + ent + prior_val))
    }
    
    grad_fun <- function(par) {
      mu <- par[1:p]
      log_sig <- par[(p + 1):(2 * p)]
      sig <- exp(log_sig)
      
      grad_mu <- numeric(p)
      grad_log_sig <- numeric(p)
      grad_ll_mu <- (X.y - XtX %*% mu) / sigma2
      
      for (j in 1:p) {
        beta_jk <- mu[j] + sqrt(2) * sig[j] * nodes
        h_prime <- -(nu + 1) * beta_jk / (nu / tau1 + beta_jk^2)
        grad_prior_mu <- (1 / sqrt(pi)) * sum(weights * h_prime)
        grad_prior_sig <- (1 / sqrt(pi)) * sum(weights * h_prime * sqrt(2) * nodes)
        
        grad_mu[j] <- grad_ll_mu[j] + grad_prior_mu
        grad_sig <- - (XtX_diag[j] * sig[j]) / sigma2 + 1 / sig[j] + grad_prior_sig
        grad_log_sig[j] <- grad_sig * sig[j]
      }
      return(-c(grad_mu, grad_log_sig))
    }
    
    for (outer in 1:opts$em_iter) { 
      init_par <- c(mu_beta, log(sigma_beta))
      opt_res <- stats::optim(init_par, obj_fun, grad_fun, method = "L-BFGS-B", control = list(maxit = opts$max_iter, factr = 1e7))
      
      mu_beta <- opt_res$par[1:p]
      sigma_beta <- exp(opt_res$par[(p + 1):(2 * p)])
      
      trace_term <- sum(XtX_diag * sigma_beta^2) 
      sigma2_new <- calc_sigma2_em_update(y, X, mu_beta, trace_term, sigma2, n, intercept, include_prior_sigma2 = TRUE)
      
      if (abs(sigma2_new - sigma2) < opts$tol_em) {
        sigma2 <- sigma2_new
        break
      }
      sigma2 <- sigma2_new
    }
  }
  
  alpha_est <- if (intercept) mean(y) else 0
  return(list(mu = mu_beta, Sigma = diag(x = sigma_beta^2, nrow = p, ncol = p), sigma2 = sigma2, alpha = alpha_est))
}

# ==============================================================================
# 3. Laplace Approximation
# ==============================================================================

#' Laplace Approximation for Spike-and-Slab Linear Models
#'
#' Gaussian approximation centred at the posterior mode with covariance from the
#' Hessian; \eqn{\sigma^2} by EM. For the Laplace prior the mode is found by
#' coordinate descent (soft-thresholding) and the covariance is computed on the
#' non-zero coordinates.
#'
#' @inheritParams EP
#' @return List with mu, Sigma, sigma2, alpha
#' @export
Laplace_Approx <- function(y, X, sigma2.hat, XtX, X.y, prior_type = c("laplace", "student_t", "gaussian"), 
                           hyper_params = list(tau = 1, nu = 3), 
                           opts = list(max_iter = 500, tol = 1e-4, em_iter = 100, tol_em = 1e-6), 
                           intercept = TRUE) {
  
  prior_type <- match.arg(prior_type)
  default_opts <- list(max_iter = 500, tol = 1e-4, em_iter = 100, tol_em = 1e-6)
  opts <- if (missing(opts) || is.null(opts)) default_opts else utils::modifyList(default_opts, as.list(opts))
  
  n <- nrow(X)
  p <- ncol(X)
  
  sigma2 <- sigma2.hat
  mu_beta <- rep(0, p)
  sigma_beta <- rep(1, p)
  XtX_diag <- diag(XtX)
  
  if (prior_type == "laplace") {
    lambda <- hyper_params$tau
    
    for (outer in 1:opts$em_iter) {
      for (iter in 1:opts$max_iter) {
        mu_old <- mu_beta
        
        for (j in 1:p) {
          rho_j <- X.y[j] - sum(XtX[j, -j] * mu_beta[-j])
          c_j <- XtX_diag[j] / sigma2
          tilde_beta <- (rho_j / sigma2) / c_j
          penalty <- lambda / c_j
          
          if (tilde_beta > penalty) {
            mu_beta[j] <- tilde_beta - penalty
          } else if (tilde_beta < -penalty) {
            mu_beta[j] <- tilde_beta + penalty
          } else {
            mu_beta[j] <- 0
          }
        }
        
        if (max(abs(mu_beta - mu_old)) < opts$tol) break
      }
      
      active_idx <- which(mu_beta != 0)
      Sigma <- matrix(0, p, p)
      
      if (length(active_idx) > 0) {
        H_active <- XtX[active_idx, active_idx, drop = FALSE] / sigma2
        Sigma_active <- tryCatch(solve(H_active), error = function(e) MASS::ginv(H_active))
        Sigma[active_idx, active_idx] <- Sigma_active
      }
      
      trace_term <- sum(XtX * Sigma)
      sigma2_new <- calc_sigma2_em_update(y, X, mu_beta, trace_term, sigma2, n, intercept, include_prior_sigma2 = FALSE)
      
      if (abs(sigma2_new - sigma2) < opts$tol_em) {
        sigma2 <- sigma2_new
        break
      }
      sigma2 <- sigma2_new
    }
    
  } else if (prior_type == "student_t") {
    nu <- hyper_params$nu
    tau1 <- hyper_params$tau
    nu_scale2 <- nu / tau1
    
    neg_log_post <- function(beta, sig2) {
      resid <- if (intercept) (y - X %*% beta - mean(y)) else (y - X %*% beta)
      rss <- sum(resid^2) / (2 * sig2)
      prior <- sum(((nu + 1) / 2) * log1p(beta^2 / nu_scale2))
      return(rss + prior)
    }
    
    grad_neg_log_post <- function(beta, sig2) {
      g_lik <- -(X.y - XtX %*% beta) / sig2
      g_prior <- (nu + 1) * beta / (nu_scale2 + beta^2)
      return(as.vector(g_lik + g_prior))
    }
    
    for (outer in 1:opts$em_iter) {
      opt_res <- stats::optim(par = mu_beta, fn = neg_log_post, gr = grad_neg_log_post,
                              sig2 = sigma2, method = "L-BFGS-B",
                              control = list(maxit = opts$max_iter, factr = 1e3))
      mu_beta <- opt_res$par
      
      d2_prior <- (nu + 1) * (nu_scale2 - mu_beta^2) / (nu_scale2 + mu_beta^2)^2
      d2_prior <- pmax(d2_prior, 0)
      
      H <- XtX / sigma2 + diag(as.vector(d2_prior), p, p)
      Sigma <- tryCatch(solve(H), error = function(e) MASS::ginv(H))
      
      trace_term <- sum(XtX * Sigma)
      sigma2_new <- calc_sigma2_em_update(y, X, mu_beta, trace_term, sigma2, n, intercept, include_prior_sigma2 = FALSE)
      
      if (abs(sigma2_new - sigma2) < opts$tol_em) {
        sigma2 <- sigma2_new
        break
      }
      sigma2 <- sigma2_new
    }
    
  } else if (prior_type == "gaussian") {
    tau1 <- hyper_params$tau
    
    for (outer in 1:opts$em_iter) {
      H <- XtX / sigma2 + diag(tau1, p, p)
      Sigma <- tryCatch(solve(H), error = function(e) MASS::ginv(H))
      mu_beta <- as.vector(Sigma %*% (X.y / sigma2))
      
      trace_term <- sum(XtX * Sigma)
      sigma2_new <- calc_sigma2_em_update(y, X, mu_beta, trace_term, sigma2, n, intercept, include_prior_sigma2 = FALSE)
      
      if (abs(sigma2_new - sigma2) < opts$tol_em) {
        sigma2 <- sigma2_new
        break
      }
      sigma2 <- sigma2_new
    }
  }
  
  alpha_est <- if (intercept) mean(y) else 0
  return(list(mu = mu_beta, Sigma = Sigma, sigma2 = sigma2, alpha = alpha_est))
}

# ==============================================================================
# 4. Unified Dispatcher Function
# ==============================================================================

#' Fit Approximate Posterior Distribution
#'
#' Unified interface for approximating \eqn{p(\beta, \sigma^2 \mid \gamma, y)}
#' for a selected model, as described in Section 4.1.1 of the paper. The slab
#' prior is applied to every column of \code{X}; \eqn{\sigma^2} is estimated by
#' EM (a MAP estimate) and \eqn{\beta} is approximated by one of:
#' \describe{
#'   \item{\code{"Block-VB"}}{Block coordinate-ascent VB with a full-covariance
#'     Gaussian \eqn{q(\beta)} (recommended default).}
#'   \item{\code{"MFVB-CAVI"}}{Mean-field CAVI with independent Gaussian factors.}
#'   \item{\code{"MMFVB"}}{Mean-field VB maximising the marginal ELBO directly.}
#'   \item{\code{"LA"}}{Laplace approximation at the posterior mode.}
#'   \item{\code{"EP"}}{Expectation propagation with Gaussian sites.}
#' }
#'
#' When \code{intercept = TRUE}, \code{X} and \code{y} are centred internally
#' (supplied \code{XtX} and \code{X.y} are assumed to be the uncentred
#' cross-products and are adjusted accordingly), and the returned \code{alpha}
#' is the intercept for the original covariates. The lower-level functions
#' (\code{VB}, \code{EP}, ...) instead assume \code{X} is already centred.
#'
#' @param y Numeric response vector
#' @param X Numeric covariate matrix for the selected model (no column of ones)
#' @param method Approximation algorithm: "Block-VB", "MFVB-CAVI", "MMFVB", "LA", or "EP"
#' @param prior_type One of "laplace", "gaussian", "student_t"
#' @param hyper_params List with the slab parameter \code{tau} (precision for
#'   Gaussian/Student-t, rate for Laplace, as in \code{\link{FASS.sampler}}) and
#'   \code{nu} (Student-t degrees of freedom)
#' @param opts Control options list (iterations, tolerances, damping); missing
#'   entries take method-specific defaults
#' @param intercept Logical, whether to include an intercept
#' @param sigma2.hat Initial residual variance (computed from OLS if NULL)
#' @param XtX Crossproduct matrix X^T X (computed if NULL)
#' @param X.y Crossproduct vector X^T y (computed if NULL)
#' @return List containing posterior mean \code{mu}, covariance \code{Sigma},
#'   residual variance \code{sigma2} and intercept \code{alpha}
#' @export
fit_approx_posterior <- function(y, X, method = c("Block-VB", "MFVB-CAVI", "MMFVB", "LA", "EP"),
                                 prior_type = c("laplace", "gaussian", "student_t"),
                                 hyper_params = list(tau = 1, nu = 3),
                                 opts = NULL,
                                 intercept = TRUE,
                                 sigma2.hat = NULL,
                                 XtX = NULL,
                                 X.y = NULL) {
  
  method <- match.arg(method)
  prior_type <- match.arg(prior_type)
  
  data_val <- validate_regression_data(y, X, check_invertibility = FALSE)
  y <- data_val$y
  X <- data_val$X
  
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(XtX)) XtX <- crossprod(X)
  if (is.null(X.y)) X.y <- crossprod(X, y)
  XtX <- matrix(XtX, nrow = p)
  X.y <- as.vector(X.y)
  
  # Profile out the intercept by centring (flat prior on the intercept)
  if (intercept) {
    x_bar <- colMeans(X)
    y_bar <- mean(y)
    X <- sweep(X, 2, x_bar, "-")
    y <- y - y_bar
    XtX <- XtX - n * tcrossprod(x_bar)
    X.y <- X.y - n * x_bar * y_bar
  }
  
  if (is.null(sigma2.hat)) {
    beta_ols <- tryCatch(as.vector(chol2inv(chol(XtX)) %*% X.y),
                         error = function(e) as.vector(MASS::ginv(XtX) %*% X.y))
    resid <- y - X %*% beta_ols
    sigma2.hat <- sum(resid^2) / max(1, n - p - (if (intercept) 1 else 0))
  }
  
  if (is.null(opts)) {
    if (method == "EP") {
      opts <- if (prior_type == "student_t") {
        list(max_iter = 500, tol = 1e-8, damping = 0.2, em_iter = 500, tol_em = 1e-8, n_nodes = 21)
      } else {
        list(max_iter = 500, tol = 1e-8, damping = 0.8, em_iter = 500, tol_em = 1e-8, n_nodes = 21)
      }
    } else {
      opts <- list(max_iter = 500, tol = 1e-4, em_iter = 100, tol_em = if (method == "LA") 1e-6 else 1e-10)
    }
  }
  
  fit <- switch(method,
         "Block-VB"  = VB(y, X, sigma2.hat, XtX, X.y, prior_type = prior_type, hyper_params = hyper_params, opts = opts, intercept = intercept),
         "MFVB-CAVI" = VB_MF(y, X, sigma2.hat, XtX, X.y, prior_type = prior_type, hyper_params = hyper_params, opts = opts, intercept = intercept),
         "MMFVB"     = VB_Direct(y, X, sigma2.hat, XtX, X.y, prior_type = prior_type, hyper_params = hyper_params, opts = opts, intercept = intercept),
         "LA"        = Laplace_Approx(y, X, sigma2.hat, XtX, X.y, prior_type = prior_type, hyper_params = hyper_params, opts = opts, intercept = intercept),
         "EP"        = EP(y, X, sigma2.hat, XtX, X.y, prior_type = prior_type, hyper_params = hyper_params, opts = opts, intercept = intercept))
  
  if (intercept) fit$alpha <- y_bar - sum(x_bar * as.vector(fit$mu))
  fit
}
