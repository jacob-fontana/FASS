################################################################################
### Simulation Evaluation and Metric Helpers
################################################################################

#' Evaluate Variable Selection Performance
#'
#' Evaluates True Positive Rate (TPR), True Negative Rate (TNR), and model size
#' from either a FASS fit or full MCMC Gibbs chain.
#'
#' @param fit Model output list from FASS.sampler, Gibbs, Gibbs_PC, or skinnyGibbs
#' @param true_idx Integer vector of active feature indices (e.g. 1:10)
#' @param p Integer, total number of candidate predictors
#' @param is_gibbs Logical, whether fit originates from a Gibbs sampler (default FALSE)
#' @return A data.frame row with TPR, TNR, active_count, and false_discoveries
#' @export
eval_variable_selection <- function(fit, true_idx, p, is_gibbs = FALSE) {
  p_true <- length(true_idx)
  inactive_idx <- setdiff(seq_len(p), true_idx)
  
  if (is_gibbs) {
    # Median probability model from Gibbs posterior samples
    pip <- colMeans(fit$outgamma)
    selected <- pip > 0.5
  } else {
    # Highest marginal odds model from FASS search trajectory
    best_idx <- which.max(fit$odds)
    if (best_idx == 1) {
      selected <- rep(FALSE, p)
    } else {
      selected <- as.logical(fit$gamma[best_idx - 1, ])
    }
  }
  
  tpr <- if (p_true > 0) mean(selected[true_idx]) else NA_real_
  tnr <- if (length(inactive_idx) > 0) 1 - mean(selected[inactive_idx]) else NA_real_
  active_count <- sum(selected)
  false_positives <- sum(selected[inactive_idx])
  
  data.frame(
    TPR = tpr,
    TNR = tnr,
    active_count = active_count,
    false_positives = false_positives
  )
}

#' Evaluate Post-Hoc Posterior Estimation
#'
#' Evaluates coefficient recovery metrics (RMSE, 95% CI Coverage, Interval Length)
#' and residual variance estimation from deterministic posterior approximations.
#'
#' @param fit Model fit list from fit_approx_posterior, EP, VB, or Laplace_Approx
#' @param beta_true Numeric vector of ground-truth coefficients
#' @param alpha_true Numeric scalar ground-truth intercept (optional)
#' @param sigma2_true Numeric scalar ground-truth error variance (default 1.0)
#' @return A data.frame row with estimation metrics
#' @export
eval_posthoc_fit <- function(fit, beta_true, alpha_true = NA, sigma2_true = 1.0) {
  p <- length(beta_true)
  
  # Point estimates
  mu <- as.vector(fit$mu)
  Sigma <- fit$Sigma
  sigma2_est <- fit$sigma2
  alpha_est <- if (!is.null(fit$alpha)) fit$alpha else fit$alpha.out
  
  # Estimation RMSE
  rmse_beta <- sqrt(mean((mu - beta_true)^2))
  mae_beta <- mean(abs(mu - beta_true))
  
  # 95% Credible interval coverage & length
  se <- sqrt(pmax(1e-12, diag(Sigma)))
  ci_lower <- mu - 1.96 * se
  ci_upper <- mu + 1.96 * se
  covered <- (beta_true >= ci_lower) & (beta_true <= ci_upper)
  coverage_beta <- mean(covered)
  mean_ci_length <- mean(ci_upper - ci_lower)
  
  # Variance & intercept errors
  err_sigma2 <- abs(sigma2_est - sigma2_true)
  err_alpha <- if (!is.na(alpha_true) && !is.null(alpha_est)) abs(alpha_est - alpha_true) else NA_real_
  
  data.frame(
    rmse_beta = rmse_beta,
    mae_beta = mae_beta,
    coverage_95 = coverage_beta,
    ci_length = mean_ci_length,
    sigma2_est = sigma2_est,
    err_sigma2 = err_sigma2,
    err_alpha = err_alpha
  )
}
