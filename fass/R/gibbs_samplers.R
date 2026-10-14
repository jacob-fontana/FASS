################################################################################
### Full MCMC Comparison Samplers (Gibbs, PCG, Skinny Gibbs)
################################################################################

#' Sufficient Statistics with the Intercept Profiled Out
#'
#' With a flat prior on the intercept, the regression on centred covariates
#' decouples the intercept from beta. The centred Gram matrix and cross-product
#' are obtained from the uncentred ones in O(np) without copying X.
#'
#' @param X Numeric n x p design matrix
#' @param y Numeric response vector
#' @param intercept Logical, whether the model includes an intercept
#' @return List with gram, X.y, x_bar, y_bar
#' @noRd
centered_sufficient_stats <- function(X, y, intercept) {
  n <- nrow(X)
  p <- ncol(X)
  gram <- matrix(crossprod(X), nrow = p)
  X.y <- as.vector(crossprod(X, y))
  if (intercept) {
    x_bar <- colMeans(X)
    y_bar <- mean(y)
    gram <- gram - n * tcrossprod(x_bar)
    X.y <- X.y - n * x_bar * y_bar
  } else {
    x_bar <- rep(0, p)
    y_bar <- 0
  }
  list(gram = gram, X.y = X.y, x_bar = x_bar, y_bar = y_bar)
}

#' Standard Gibbs Sampler for Spike-and-Slab
#'
#' Uncollapsed Gibbs sampler (stochastic search variable selection) for the
#' continuous spike-and-slab model. Inclusion indicators are updated from the
#' Gaussian conditional odds given the latent scales, which is cheap but can mix
#' poorly; see \code{\link{Gibbs_PC}} for the partially collapsed variant.
#' Each iteration costs \eqn{O(p^3)} after an \eqn{O(np^2)} pre-computation,
#' plus \eqn{O(np)} for the residual variance update.
#'
#' @section Prior parameterisation:
#' Unlike \code{\link{FASS.sampler}}, these samplers take squared \emph{scales}:
#' \code{tau0_2} and \code{tau1_2} are the variance of the Gaussian spike/slab,
#' the squared scale of the Student-t spike/slab, and the squared inverse rate of
#' the Laplace spike/slab. The FASS setting \code{tau0 = n, tau1 = 1} (Gaussian,
#' Student) or \code{tau0 = sqrt(n), tau1 = 1} (Laplace) corresponds to
#' \code{tau0_2 = 1/n, tau1_2 = 1}.
#'
#' Selection is usually based on the median probability model,
#' \code{colMeans(fit$outgamma) > 0.5}.
#'
#' @param X Numeric n x p design matrix
#' @param y Numeric response vector of length n
#' @param tau0_2 Squared scale of the spike (see Prior parameterisation)
#' @param tau1_2 Squared scale of the slab (see Prior parameterisation)
#' @param a_pi Prior Beta shape parameter for inclusion probability
#' @param b_pi Prior Beta shape parameter for inclusion probability
#' @param nburn Number of burn-in iterations
#' @param niter Number of post-burn-in iterations
#' @param family Prior slab family: "Gaussian", "Laplace", or "Student"
#' @param nu Degrees of freedom when family is "Student" (default NA)
#' @param intercept Logical, whether to include an intercept (default TRUE)
#' @param rescale_y Logical, whether to standardize y by sd(y) internally and rescale outputs (default FALSE)
#' @return List of retained draws: outbeta, outgamma, outW (latent scales),
#'   outheta (inclusion probability), outsigma (residual variance), outalpha
#'   (intercept on the original covariate scale, NULL if intercept = FALSE),
#'   plus active_counts and per-iteration times (seconds) for all iterations
#' @export
Gibbs <- function(X, y, tau0_2, tau1_2, a_pi, b_pi, nburn = 2000, niter = 2000, 
                  family = c("Gaussian", "Laplace", "Student"), nu = NA, 
                  intercept = TRUE, rescale_y = FALSE) {
  
  data_val <- validate_regression_data(y, X, check_invertibility = FALSE)
  y <- data_val$y
  X <- data_val$X
  
  params_val <- validate_gibbs_params(tau0_2, tau1_2, a_pi, b_pi, nburn, niter, family, nu)
  tau0_2 <- params_val$tau0_2
  tau1_2 <- params_val$tau1_2
  a_pi <- params_val$a_pi
  b_pi <- params_val$b_pi
  nburn <- params_val$nburn
  niter <- params_val$niter
  family <- params_val$family
  nu <- params_val$nu
  times <- c()
  
  if (rescale_y) {
    scale_y <- stats::sd(y)
    y <- y / scale_y
  } else {
    scale_y <- 1
  }
  
  n <- nrow(X)
  p <- ncol(X)
  ss <- centered_sufficient_stats(X, y, intercept)
  gram <- ss$gram
  X.y <- ss$X.y
  x_bar <- ss$x_bar
  y_bar <- ss$y_bar
  
  outbeta <- matrix(nrow = niter, ncol = p)
  outgamma <- matrix(nrow = niter, ncol = p)
  out.W <- matrix(nrow = niter, ncol = p)
  outtheta <- numeric(niter)
  outalpha <- if (intercept) numeric(niter) else NULL
  outsigma <- numeric(niter)
  
  beta <- rep(0, p)
  gamma <- rep(0, p)
  theta <- 0.01
  sigma2 <- 1
  W <- rep(1, p)
  active.counts <- c()
  
  idx_active <- which(gamma == 1)
  idx_inactive <- which(gamma == 0)
  tau.vec <- rep(tau0_2, p)
  tau.vec[idx_active] <- tau1_2
  n_active <- length(idx_active)
  
  for (itr in 1:(nburn + niter)) {
    start <- Sys.time()
    
    u <- stats::rnorm(p, 0, 1)
    L.t <- chol(gram / sigma2 + diag(1 / (W * tau.vec), nrow = p))
    v <- forwardsolve(t(L.t), X.y / sigma2)
    beta <- as.vector(rowSums(backsolve(L.t, cbind(v, u))))
    
    # Intercept of the centred design is independent of beta given sigma2;
    # map it back to the intercept for the original covariates
    if (intercept == TRUE) alpha <- stats::rnorm(1, y_bar, sqrt(sigma2 / n)) - sum(x_bar * beta)
    
    if (intercept == TRUE) {
      sigma2 <- 1 / stats::rgamma(1, n / 2, sum((y - X %*% beta - alpha)^2) / 2)
    } else {
      sigma2 <- 1 / stats::rgamma(1, n / 2, sum((y - X %*% beta)^2) / 2)
    }
    
    theta <- stats::rbeta(1, a_pi + n_active, b_pi + p - n_active)
    
    temp <- log(theta) - log(1 - theta) + 
      stats::dnorm(beta, 0, sd = sqrt(tau1_2 * W), log = TRUE) - 
      stats::dnorm(beta, 0, sd = sqrt(tau0_2 * W), log = TRUE)
    
    prob <- exp(temp) / (1 + exp(temp))
    prob[is.na(prob)] <- 1
    gamma <- stats::rbinom(length(prob), 1, prob)
    
    idx_active <- which(gamma == 1)
    idx_inactive <- which(gamma == 0)
    tau.vec <- rep(tau0_2, p)
    tau.vec[idx_active] <- tau1_2
    n_active <- length(idx_active)
    active.counts <- c(active.counts, n_active)
    
    if (family == 'Laplace') W <- 1 / rinversegauss(p, sqrt(tau.vec / beta^2), 1)
    if (family == 'Student') W <- 1 / stats::rgamma(p, (nu + 1) / 2, (nu + beta^2 / tau.vec) / 2)
    
    times <- c(times, difftime(Sys.time(), start, units = "secs"))
    
    if (itr > nburn) {
      out.W[(itr - nburn), ] <- W
      outbeta[(itr - nburn), ] <- beta
      outgamma[(itr - nburn), ] <- gamma
      outtheta[itr - nburn] <- theta
      outsigma[itr - nburn] <- sigma2
      if (intercept == TRUE) outalpha[itr - nburn] <- alpha
    }
  }
  
  if (rescale_y) {
    outbeta <- outbeta * scale_y
    outsigma <- outsigma * (scale_y^2)
    if (intercept == TRUE) outalpha <- outalpha * scale_y
  }
  
  return(list(outbeta = outbeta, outgamma = outgamma, outW = out.W, outheta = outtheta, 
              outsigma = outsigma, outalpha = outalpha, active_counts = active.counts, times = times))
}

#' Partially Collapsed Gibbs (PCG) Sampler
#'
#' Gibbs sampler in which the inclusion indicators are updated from the
#' marginal prior odds of the spike-and-slab family (integrating out the latent
#' scales), so selection reflects the Laplace or Student-t prior rather than a
#' conditionally Gaussian one.
#'
#' @inheritParams Gibbs
#' @inheritSection Gibbs Prior parameterisation
#' @return Same structure as \code{\link{Gibbs}}
#' @export
Gibbs_PC <- function(X, y, tau0_2, tau1_2, a_pi, b_pi, nburn = 2000, niter = 2000, 
                     family = c("Gaussian", "Laplace", "Student"), nu = NA, 
                     intercept = TRUE, rescale_y = FALSE) {
  
  data_val <- validate_regression_data(y, X, check_invertibility = FALSE)
  y <- data_val$y
  X <- data_val$X
  
  params_val <- validate_gibbs_params(tau0_2, tau1_2, a_pi, b_pi, nburn, niter, family, nu)
  tau0_2 <- params_val$tau0_2
  tau1_2 <- params_val$tau1_2
  a_pi <- params_val$a_pi
  b_pi <- params_val$b_pi
  nburn <- params_val$nburn
  niter <- params_val$niter
  family <- params_val$family
  nu <- params_val$nu
  times <- c()
  
  if (rescale_y) {
    scale_y <- stats::sd(y)
    y <- y / scale_y
  } else {
    scale_y <- 1
  }
  
  n <- nrow(X)
  p <- ncol(X)
  ss <- centered_sufficient_stats(X, y, intercept)
  gram <- ss$gram
  X.y <- ss$X.y
  x_bar <- ss$x_bar
  y_bar <- ss$y_bar
  
  outbeta <- matrix(nrow = niter, ncol = p)
  outgamma <- matrix(nrow = niter, ncol = p)
  out.W <- matrix(nrow = niter, ncol = p)
  outtheta <- numeric(niter)
  outalpha <- if (intercept) numeric(niter) else NULL
  outsigma <- numeric(niter)
  
  beta <- stats::rnorm(p, 0, 1)
  theta <- 0.5
  sigma2 <- 1
  W <- rep(1, p)
  active.counts <- c()
  
  for (itr in 1:(nburn + niter)) {
    start <- Sys.time()
    
    if (family == 'Gaussian') {
      temp <- log(theta) - log(1 - theta) + 
        stats::dnorm(beta, 0, sd = sqrt(tau1_2 * W), log = TRUE) - 
        stats::dnorm(beta, 0, sd = sqrt(tau0_2 * W), log = TRUE)
    } else if (family == 'Laplace') {
      temp <- log(theta) - log(1 - theta) + 0.5 * log(tau0_2 / tau1_2) + 
        (abs(beta) * (1 / sqrt(tau0_2) - 1 / sqrt(tau1_2)))
    } else if (family == 'Student') {
      temp <- log(theta) - log(1 - theta) + 
        (stats::dt(beta / sqrt(tau1_2), df = nu, log = TRUE) - 0.5 * log(tau1_2)) - 
        (stats::dt(beta / sqrt(tau0_2), df = nu, log = TRUE) - 0.5 * log(tau0_2))
    }
    
    prob <- exp(temp) / (1 + exp(temp))
    prob[is.na(prob)] <- 1
    gamma <- stats::rbinom(length(prob), 1, prob)
    
    idx_active <- which(gamma == 1)
    idx_inactive <- which(gamma == 0)
    tau.vec <- rep(tau0_2, p)
    tau.vec[idx_active] <- tau1_2
    n_active <- length(idx_active)
    active.counts <- c(active.counts, n_active)
    
    if (family == 'Laplace') W <- 1 / rinversegauss(p, sqrt(tau.vec / beta^2), 1)
    if (family == 'Student') W <- 1 / stats::rgamma(p, (nu + 1) / 2, (nu + beta^2 / tau.vec) / 2)
    
    u <- stats::rnorm(p, 0, 1)
    L.t <- chol(gram / sigma2 + diag(1 / (W * tau.vec), nrow = p))
    v <- forwardsolve(t(L.t), X.y / sigma2)
    beta <- as.vector(rowSums(backsolve(L.t, cbind(v, u))))
    
    # Intercept of the centred design is independent of beta given sigma2;
    # map it back to the intercept for the original covariates
    if (intercept == TRUE) alpha <- stats::rnorm(1, y_bar, sqrt(sigma2 / n)) - sum(x_bar * beta)
    
    if (intercept == TRUE) {
      sigma2 <- 1 / stats::rgamma(1, n / 2, sum((y - X %*% beta - alpha)^2) / 2)
    } else {
      sigma2 <- 1 / stats::rgamma(1, n / 2, sum((y - X %*% beta)^2) / 2)
    }
    
    theta <- stats::rbeta(1, a_pi + n_active, b_pi + p - n_active)
    
    times <- c(times, difftime(Sys.time(), start, units = "secs"))
    
    if (itr > nburn) {
      out.W[(itr - nburn), ] <- W
      outbeta[(itr - nburn), ] <- beta
      outgamma[(itr - nburn), ] <- gamma
      outtheta[itr - nburn] <- theta
      outsigma[itr - nburn] <- sigma2
      if (intercept == TRUE) outalpha[itr - nburn] <- alpha
    }
  }
  
  if (rescale_y) {
    outbeta <- outbeta * scale_y
    outsigma <- outsigma * (scale_y^2)
    if (intercept == TRUE) outalpha <- outalpha * scale_y
  }
  
  return(list(outbeta = outbeta, outgamma = outgamma, outW = out.W, outheta = outtheta, 
              outsigma = outsigma, outalpha = outalpha, active_counts = active.counts, times = times))
}

#' Skinny Gibbs Sampler
#'
#' Skinny Gibbs sampler of Narisetty, Shen and He (2019): only the coefficients
#' currently assigned to the slab are drawn jointly, while the remaining ones are
#' drawn independently, reducing the per-iteration cost to roughly
#' \eqn{O(np + |\gamma|^3)}.
#'
#' @inheritParams Gibbs
#' @inheritSection Gibbs Prior parameterisation
#' @return Same structure as \code{\link{Gibbs}}
#' @export
skinnyGibbs <- function(X, y, tau0_2, tau1_2, a_pi, b_pi, nburn = 2000, niter = 2000, 
                        family = c("Gaussian", "Laplace", "Student"), nu = NA, 
                        intercept = TRUE, rescale_y = FALSE) {
  
  data_val <- validate_regression_data(y, X, check_invertibility = FALSE)
  y <- data_val$y
  X <- data_val$X
  
  params_val <- validate_gibbs_params(tau0_2, tau1_2, a_pi, b_pi, nburn, niter, family, nu)
  tau0_2 <- params_val$tau0_2
  tau1_2 <- params_val$tau1_2
  a_pi <- params_val$a_pi
  b_pi <- params_val$b_pi
  nburn <- params_val$nburn
  niter <- params_val$niter
  family <- params_val$family
  nu <- params_val$nu
  times <- c()
  
  if (rescale_y) {
    scale_y <- stats::sd(y)
    y <- y / scale_y
  } else {
    scale_y <- 1
  }
  
  n <- nrow(X)
  p <- ncol(X)
  ss <- centered_sufficient_stats(X, y, intercept)
  gram <- ss$gram
  diags <- diag(gram)
  X.y <- ss$X.y
  x_bar <- ss$x_bar
  y_bar <- ss$y_bar
  
  outbeta <- matrix(nrow = niter, ncol = p)
  outgamma <- matrix(nrow = niter, ncol = p)
  out.W <- matrix(nrow = niter, ncol = p)
  outtheta <- numeric(niter)
  outalpha <- if (intercept) numeric(niter) else NULL
  outsigma <- numeric(niter) 
  
  beta <- rep(0, p)
  gamma <- rep(0, p)
  theta <- 0.01
  sigma2 <- 1
  W <- rep(1, p)
  active.counts <- c()
  
  idx_active <- which(gamma == 1)
  idx_inactive <- which(gamma == 0)
  tau.vec <- rep(tau0_2, p)
  tau.vec[idx_active] <- tau1_2
  n_active <- length(idx_active)
  
  for (itr in 1:(nburn + niter)) {
    start <- Sys.time()
    
    if (n_active == 0) {
      beta <- stats::rnorm(n = p, mean = 0, sd = 1 / sqrt(diags + 1 / (W * tau0_2))) 
    } else {
      u <- stats::rnorm(n_active, 0, 1)
      L.t <- chol(gram[idx_active, idx_active] / sigma2 + diag(1 / (W[idx_active] * tau1_2), nrow = n_active))
      v <- forwardsolve(t(L.t), X.y[idx_active] / sigma2)
      beta[idx_active] <- as.vector(rowSums(backsolve(L.t, cbind(v, u))))
      beta[idx_inactive] <- stats::rnorm(n = p - n_active, mean = 0, sd = 1 / sqrt(diags[idx_inactive] + 1 / (W[idx_inactive] * tau0_2)))
    }
    
    # Skinny Gibbs keeps only the active block in the likelihood, so the
    # intercept is mapped back using the active coefficients only
    if (intercept == TRUE) {
      alpha <- stats::rnorm(1, y_bar, sqrt(sigma2 / n))
      if (n_active > 0) alpha <- alpha - sum(x_bar[idx_active] * beta[idx_active])
    }
    
    if (n_active == 0) {
      if (intercept == TRUE) {
        sigma2 <- 1 / stats::rgamma(1, n / 2, sum((y - alpha)^2) / 2)
      } else {
        sigma2 <- 1 / stats::rgamma(1, n / 2, sum(y^2) / 2)
      }
    } else {
      if (intercept == TRUE) {
        rss <- calc_rss_cpp(X, y - alpha, beta[idx_active], idx_active)
      } else {
        rss <- calc_rss_cpp(X, y, beta[idx_active], idx_active)
      }
      sigma2 <- 1 / stats::rgamma(1, n / 2, rss / 2)
    }
    
    theta <- stats::rbeta(1, a_pi + n_active, b_pi + p - n_active)
    
    temp <- log(theta) - log(1 - theta) + 
      stats::dnorm(beta, 0, sd = sqrt(tau1_2 * W), log = TRUE) - 
      stats::dnorm(beta, 0, sd = sqrt(tau0_2 * W), log = TRUE) + 
      diags * (1 - 1 / sigma2) * beta^2 / 2
    
    if (n_active == 0) {
      log_d_j_star <- temp + beta * X.y / sigma2
      prob <- exp(log_d_j_star) / (1 + exp(log_d_j_star))
      prob[is.na(prob)] <- 1
      gamma <- stats::rbinom(length(prob), 1, prob)
    } else {
      temp2 <- beta * (X.y - as.vector(gram[, idx_active] %*% as.matrix(beta[idx_active]))) / sigma2 
      log_d_j_star <- temp + temp2
      for (j in idx_active) {
        log_d_j_star[j] <- log_d_j_star[j] + beta[j] * gram[j, j] * beta[j] / sigma2
      }
      prob <- exp(log_d_j_star) / (1 + exp(log_d_j_star))
      prob[is.na(prob)] <- 1
      gamma <- stats::rbinom(length(prob), 1, prob)
    }
    
    idx_active <- which(gamma == 1)
    idx_inactive <- which(gamma == 0)
    tau.vec <- rep(tau0_2, p)
    tau.vec[idx_active] <- tau1_2
    n_active <- length(idx_active)
    active.counts <- c(active.counts, n_active)
    
    if (family == 'Laplace') W <- 1 / rinversegauss(p, sqrt(tau.vec / beta^2), 1)
    if (family == 'Student') W <- 1 / stats::rgamma(p, (nu + 1) / 2, (nu + beta^2 / tau.vec) / 2)
    
    times <- c(times, difftime(Sys.time(), start, units = "secs"))
    
    if (itr > nburn) {
      out.W[(itr - nburn), ] <- W
      outbeta[(itr - nburn), ] <- beta
      outgamma[(itr - nburn), ] <- gamma
      outtheta[itr - nburn] <- theta
      outsigma[itr - nburn] <- sigma2
      if (intercept == TRUE) outalpha[itr - nburn] <- alpha
    }
  }
  
  if (rescale_y) {
    outbeta <- outbeta * scale_y
    outsigma <- outsigma * (scale_y^2)
    if (intercept == TRUE) outalpha <- outalpha * scale_y
  }
  
  return(list(outbeta = outbeta, outgamma = outgamma, outW = out.W, outheta = outtheta, 
              outalpha = outalpha, outsigma = outsigma, active_counts = active.counts, times = times))
}
