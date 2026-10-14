################################################################################
### Data Generating Processes (DGPs) for Synthetic Benchmarks
################################################################################

#' Generate Variable Selection Data (Superinduction & M-Closed)
#'
#' Generates synthetic regression datasets for variable selection benchmarks.
#' Supports both the M-Closed setting (epsilon = 0) and the Superinduction
#' setting (epsilon > 0) where noise features leak an unobserved latent effect.
#'
#' @param n Integer, sample size
#' @param p Integer, total number of candidate predictors
#' @param p_true Integer, number of active (true signal) predictors
#' @param epsilon Numeric, contamination weight for superinduction (default 0.01; set 0 for M-Closed)
#' @param seed Integer, random seed (default 4747)
#' @param scale_X Logical, whether to scale columns of X by their standard deviations (default TRUE)
#' @return A list containing:
#'   \item{y}{Response vector of length n}
#'   \item{X}{Matrix of dimension n x p containing active and inactive predictors}
#'   \item{p_true}{Number of true predictors}
#'   \item{true_idx}{Integer indices of active predictors (1:p_true)}
#'   \item{beta}{True coefficient vector for active predictors}
#'   \item{alpha}{True scalar intercept}
#'   \item{scales}{Column scales applied to X}
#' @export
generate_selection_data <- function(n, p, p_true, epsilon = 0.01, seed = 4747, scale_X = TRUE) {
  set.seed(seed)
  
  # 1. True signals
  means <- rep(0, p_true)
  sds <- sqrt(stats::rchisq(p_true, df = 20))
  K <- matrix(0, nrow = n, ncol = p_true)
  for (i in seq_len(p_true)) {
    K[, i] <- stats::rnorm(n, mean = means[i], sd = sds[i])
  }
  
  alpha <- stats::rnorm(1, mean = 0, sd = sqrt(10))
  beta <- stats::rnorm(p_true, mean = 0, sd = 10)
  mu_delta <- stats::rnorm(n, mean = 0, sd = sqrt(10))
  
  # Response vector with latent confounding effect mu_delta
  y <- stats::rnorm(n, mean = alpha + as.vector(K %*% beta) + mu_delta, sd = 1)
  
  # 2. Construct design matrix J
  J <- matrix(0, nrow = n, ncol = p)
  J[, seq_len(p_true)] <- K
  
  if (p > p_true) {
    for (i in (p_true + 1):p) {
      if (epsilon > 0) {
        # Superinduction: inactive features correlate with mu_delta
        J[, i] <- (1 - epsilon) * stats::rnorm(n, 0, sqrt(10)) + epsilon * mu_delta
      } else {
        # M-Closed: pure noise features
        J[, i] <- stats::rnorm(n, 0, sqrt(10))
      }
    }
  }
  
  scales <- NULL
  if (scale_X) {
    J <- scale(J, center = FALSE, scale = TRUE)
    scales <- attr(J, "scaled:scale")
  }
  
  list(
    y = y,
    X = J,
    p_true = p_true,
    true_idx = seq_len(p_true),
    beta = beta,
    alpha = alpha,
    scales = scales
  )
}

#' Generate Post-Hoc Inference Data (Correlated & Uncorrelated)
#'
#' Generates regression datasets with AR(1) or identity correlation structures
#' for testing post-hoc estimation methods (EP, VB, Laplace Approximation).
#'
#' @param n Integer, sample size
#' @param p Integer, number of active predictors in model
#' @param rho Numeric, AR(1) autocorrelation parameter (rho = 0 for uncorrelated)
#' @param signal_sd Numeric, standard deviation of non-zero coefficients (default 10)
#' @param min_signal Numeric, minimum absolute coefficient magnitude (default 0.1)
#' @param seed Integer, random seed (default 4747 + p)
#' @return A list containing y, X, beta_true, XtX, Xy, and sigma2_true
#' @export
generate_posthoc_data <- function(n, p, rho = 0.8, signal_sd = 10, min_signal = 0.1, seed = 4747 + p) {
  set.seed(seed)
  
  # Draw coefficients above threshold
  beta_true <- numeric(0)
  while (length(beta_true) < p) {
    cand <- stats::rnorm(p * 2, 0, signal_sd)
    beta_true <- c(beta_true, cand[abs(cand) >= min_signal])
  }
  beta_true <- beta_true[seq_len(p)]
  
  # Generate predictors with AR(1) structure
  if (abs(rho) > 1e-6) {
    # Cholesky of Toeplitz AR(1)
    Sigma_X <- matrix(0, nrow = p, ncol = p)
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        Sigma_X[i, j] <- rho^abs(i - j)
      }
    }
    L <- chol(Sigma_X)
    Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    X <- Z %*% L
  } else {
    X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  }
  
  alpha_true <- stats::rnorm(1, 0, 1)
  sigma2_true <- 1.0
  y <- stats::rnorm(n, mean = alpha_true + as.vector(X %*% beta_true), sd = sqrt(sigma2_true))
  
  XtX <- crossprod(X)
  Xy <- as.vector(crossprod(X, y))
  
  list(
    y = y,
    X = X,
    beta_true = beta_true,
    alpha_true = alpha_true,
    XtX = XtX,
    Xy = Xy,
    sigma2_true = sigma2_true
  )
}
