################################################################################
### Out-of-Sample Prediction and Scoring Metrics
################################################################################

#' Compute Empirical Continuous Ranked Probability Score (CRPS)
#'
#' Evaluates CRPS empirically from predictive MCMC samples.
#'
#' @param preds Matrix of dimensions n x M containing M draws for n observations
#' @param y Vector of true observed response values of length n
#' @return Numeric vector of length n containing CRPS scores
#' @export
compute_empirical_crps <- function(preds, y) {
  n <- nrow(preds)
  M <- ncol(preds)
  
  crps_scores <- numeric(n)
  weights <- 2 * seq_len(M) - 1 - M
  
  for (i in seq_len(n)) {
    x_sorted <- sort(preds[i, ])
    term1 <- mean(abs(x_sorted - y[i]))
    term2 <- sum(weights * x_sorted) / (M^2)
    crps_scores[i] <- term1 - term2
  }
  
  return(crps_scores)
}

#' Compute Analytical Gaussian CRPS
#'
#' @param mu_pred Vector of predicted conditional means
#' @param sig_pred Vector of predicted conditional standard deviations
#' @param y Vector of observed response values
#' @return Numeric vector of CRPS scores
#' @export
compute_gaussian_crps <- function(mu_pred, sig_pred, y) {
  z <- (y - mu_pred) / sig_pred
  pdf_z <- stats::dnorm(z)
  cdf_z <- stats::pnorm(z)
  crps_scores <- sig_pred * (z * (2 * cdf_z - 1) + 2 * pdf_z - 1 / sqrt(pi))
  return(crps_scores)
}

#' Prediction from MCMC/Gibbs Posterior Output
#'
#' Computes out-of-sample posterior predictive draws, CRPS, MSE, RMSE, and 95\%
#' predictive intervals using the median probability model
#' (\code{colMeans(outgamma) > 0.5}) from the MCMC draws.
#'
#' @param output Output list from \code{\link{Gibbs}}, \code{\link{Gibbs_PC}}
#'   or \code{\link{skinnyGibbs}}
#' @param X Out-of-sample covariate design matrix (same columns as in training)
#' @param y Out-of-sample response vector
#' @return List with \code{predictions} (n x draws matrix) and \code{results}
#'   (summary, coverage, length, n_selected, CRPS, MSE, RMSE)
#' @export
predict_gibbs <- function(output, X, y) {
  val <- validate_predict_data(output, X, y)
  X <- val$X
  y <- val$y
  n <- length(y)
  p <- ncol(X)
  
  if (is.null(output$outgamma) || is.null(output$outbeta) || is.null(output$outsigma)) {
    stop("Output must contain 'outgamma', 'outbeta', and 'outsigma' matrices/vectors from Gibbs samplers.")
  }
  if (ncol(output$outgamma) != p) {
    stop(sprintf("Feature count mismatch: model was trained on %d features but test matrix 'X' has %d columns.",
                 ncol(output$outgamma), p))
  }
  
  gamma <- as.logical(colMeans(output$outgamma) > 0.5)
  
  predict_from_draws(output$outbeta[, gamma, drop = FALSE], output$outalpha, output$outsigma,
                     X[, gamma, drop = FALSE], y)
}

#' Posterior Predictive Summaries from Coefficient Draws
#'
#' @param post_beta Draws x k matrix of coefficients for the columns of X_sub
#' @param post_alpha Vector of intercept draws (or NULL)
#' @param post_sigma2 Vector of residual variance draws
#' @param X_sub n x k test design matrix
#' @param y Test response vector
#' @return List with predictions and results (see predict_gibbs)
#' @noRd
predict_from_draws <- function(post_beta, post_alpha, post_sigma2, X_sub, y) {
  n <- length(y)
  n_draws <- length(post_sigma2)
  n_selected <- ncol(post_beta)
  if (is.null(n_selected) || n_selected == 0) {
    post_beta <- matrix(0, nrow = n_draws, ncol = 1)
    X_sub <- matrix(0, nrow = n, ncol = 1)
    n_selected <- 0L
  }
  has_alpha <- !is.null(post_alpha) && length(post_alpha) >= n_draws
  
  mean_mat <- X_sub %*% t(post_beta)
  if (has_alpha) mean_mat <- sweep(mean_mat, 2, post_alpha[seq_len(n_draws)], "+")
  y_pred <- matrix(stats::rnorm(n * n_draws, mean = mean_mat,
                                sd = rep(sqrt(post_sigma2), each = n)),
                   nrow = n, ncol = n_draws)
  
  CRPS <- compute_empirical_crps(y_pred, y)
  
  summary_stats <- t(apply(y_pred, 1, function(x) {
    c(mean  = mean(x, na.rm = TRUE),
      q2.5  = as.numeric(stats::quantile(x, probs = 0.025, na.rm = TRUE)),
      q97.5 = as.numeric(stats::quantile(x, probs = 0.975, na.rm = TRUE)))
  }))
  
  mse <- mean((summary_stats[, "mean"] - y)^2)
  rmse <- sqrt(mse)
  covered <- (y >= summary_stats[, "q2.5"] & y <= summary_stats[, "q97.5"])
  interval_lengths <- summary_stats[, "q97.5"] - summary_stats[, "q2.5"]
  
  results <- list(
    summary = summary_stats,
    coverage = covered,
    length = interval_lengths,
    n_selected = n_selected,
    CRPS = CRPS,
    MSE = mse,
    RMSE = rmse
  )
  
  return(list(predictions = y_pred, results = results))
}

#' Prediction from Approximate Posteriors (FASS, EP, VB, LA)
#'
#' Computes out-of-sample Gaussian predictive distributions and CRPS. For a
#' \code{\link{FASS.sampler}} fit, the highest-odds model is used with its
#' post-hoc posterior (if \code{full_posterior = "Gibbs"}, predictive draws are
#' used instead of the Gaussian form). A bare posterior list from
#' \code{\link{fit_approx_posterior}} is treated as using all columns of \code{X}.
#'
#' @param output Output list from \code{FASS.sampler} (fitted with a
#'   \code{full_posterior} method) or from \code{fit_approx_posterior}
#' @param X Out-of-sample covariate design matrix (all candidate columns)
#' @param y Out-of-sample response vector
#' @return List containing summary, coverage, length, n_selected, CRPS, MSE, RMSE
#' @export
predict_approx <- function(output, X, y) {
  val <- validate_predict_data(output, X, y)
  X <- val$X
  y <- val$y
  n <- length(y)
  p <- ncol(X)
  
  if (!is.null(output$gamma)) {
    max_idx <- which.max(output$odds)
    if (max_idx == 1) {
      gamma <- rep(FALSE, p)
    } else {
      gamma <- as.logical(as.vector(output$gamma[max_idx - 1, ]))
    }
  } else {
    gamma <- rep(TRUE, p)
  }
  
  n_selected <- sum(gamma)
  X_sub <- X[, gamma, drop = FALSE]
  
  posterior <- if (!is.null(output$posterior)) output$posterior else output
  if (is.list(posterior) && !is.null(posterior$outbeta) && !is.null(posterior$outgamma)) {
    return(predict_gibbs(posterior, X, y)$results)
  }
  if (is.list(posterior) && !is.null(posterior$outsigma) && is.null(posterior$outgamma)) {
    # Post-hoc Gibbs draws from FASS.sampler(full_posterior = "Gibbs"); the null
    # model has no outbeta
    post_beta <- if (!is.null(posterior$outbeta)) posterior$outbeta else matrix(0, length(posterior$outsigma), 0)
    if (ncol(post_beta) != n_selected) {
      stop("Dimension mismatch between posterior draws and selected variables in gamma")
    }
    return(predict_from_draws(post_beta, posterior$outalpha, posterior$outsigma, X_sub, y)$results)
  }
  if (!is.list(posterior) || is.null(posterior$mu)) {
    stop("Output must contain an estimated posterior list with 'mu' and 'Sigma'.")
  }
  
  alpha_hat <- if (!is.null(posterior$alpha)) posterior$alpha else if (!is.null(posterior$alpha.out)) posterior$alpha.out else if (!is.null(output$alpha)) output$alpha else 0
  
  if (n_selected == 0) {
    y_pred_mean <- rep(alpha_hat, n)
    sigma2_noise <- if (!is.null(posterior$sigma2)) posterior$sigma2 else if (!is.null(posterior$sigma2.out)) posterior$sigma2.out else stats::var(y)
    y_pred_sd <- rep(sqrt(sigma2_noise), n)
  } else {
    mu_beta <- posterior$mu
    Sigma_beta <- posterior$Sigma
    
    if (length(mu_beta) == length(gamma)) {
      mu_beta <- mu_beta[gamma]
      Sigma_beta <- Sigma_beta[gamma, gamma, drop = FALSE]
    } else if (length(mu_beta) != sum(gamma)) {
      stop("Dimension mismatch between posterior$mu and selected variables in gamma")
    }
    
    sigma2_noise <- if (!is.null(posterior$sigma2)) posterior$sigma2 else if (!is.null(posterior$sigma2.out)) posterior$sigma2.out else stats::var(y)
    y_pred_mean <- alpha_hat + as.vector(X_sub %*% mu_beta)
    model_uncertainty <- rowSums((X_sub %*% Sigma_beta) * X_sub)
    y_pred_var <- model_uncertainty + sigma2_noise
    y_pred_sd <- sqrt(y_pred_var)
  }
  
  CRPS <- compute_gaussian_crps(y_pred_mean, y_pred_sd, y)
  lower_bound <- y_pred_mean - stats::qnorm(0.975) * y_pred_sd
  upper_bound <- y_pred_mean + stats::qnorm(0.975) * y_pred_sd
  
  model_summary <- data.frame(
    mean = as.vector(y_pred_mean), 
    q2.5 = as.vector(lower_bound),
    q97.5 = as.vector(upper_bound)
  )
  
  mse <- mean((model_summary$mean - y)^2)
  rmse <- sqrt(mse)
  covered <- (y >= model_summary$q2.5 & y <= model_summary$q97.5)
  interval_lengths <- model_summary$q97.5 - model_summary$q2.5
  
  results <- list(
    summary = model_summary,
    coverage = covered,
    length = interval_lengths,
    n_selected = n_selected,
    CRPS = CRPS,
    MSE = mse,
    RMSE = rmse
  )
  
  return(results)
}

#' @rdname predict_approx
#' @export
predict_FASS <- predict_approx
