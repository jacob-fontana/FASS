test_that("C++ calc_rss_cpp computes exact residual sums of squares", {
  set.seed(42)
  n <- 50
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  beta <- rnorm(p)
  
  rss_r <- as.numeric(sum((y - X %*% beta)^2))
  rss_cpp <- calc_rss_cpp(X, y, beta, seq_len(p))
  
  expect_equal(rss_cpp, rss_r, tolerance = 1e-12)
})

test_that("compute_gaussian_crps computes Continuous Ranked Probability Score accurately", {
  # Single observation and prediction
  y_obs <- 2.0
  pred_mean <- 1.5
  pred_sd <- 1.0
  
  z <- (y_obs - pred_mean) / pred_sd
  crps_formula <- pred_sd * (z * (2 * stats::pnorm(z) - 1) + 2 * stats::dnorm(z) - 1 / sqrt(pi))
  crps_fn <- compute_gaussian_crps(pred_mean, pred_sd, y_obs)
  
  expect_equal(crps_fn, crps_formula, tolerance = 1e-10)
  
  # Vector calculation
  y_vec <- c(1.0, 2.0, -0.5)
  means <- c(1.2, 1.8, 0.0)
  sds <- c(0.8, 1.1, 0.9)
  crps_vec <- compute_gaussian_crps(means, sds, y_vec)
  expect_equal(length(crps_vec), 3)
})

test_that("predict_gibbs and predict_approx generate valid predictive distributions", {
  set.seed(123)
  n_train <- 80
  n_test <- 20
  p <- 4
  X_all <- matrix(rnorm((n_train + n_test) * p), n_train + n_test, p)
  y_all <- X_all[, 1] * 2 + rnorm(n_train + n_test)
  
  X_train <- X_all[1:n_train, ]
  y_train <- y_all[1:n_train]
  X_test <- X_all[(n_train + 1):(n_train + n_test), ]
  y_test <- y_all[(n_train + 1):(n_train + n_test)]
  
  # 1. Prediction from Gibbs
  fit_gibbs <- Gibbs_PC(X_train, y_train, tau0_2 = 1/n_train, tau1_2 = n_train,
                        a_pi = 1, b_pi = 1, nburn = 30, niter = 50, family = "Laplace")
  pred_g <- predict_gibbs(fit_gibbs, X_test, y_test)
  expect_true(is.list(pred_g))
  expect_equal(length(pred_g$results$summary[, "mean"]), n_test)
  expect_true(is.numeric(mean(pred_g$results$CRPS)))
  expect_true(mean(pred_g$results$CRPS) > 0)
  
  # 2. Prediction from FASS / approx
  fit_fass <- FASS.sampler(y_train, X_train, niter = 50, family = "Laplace",
                           tau1 = 1, tau0 = sqrt(n_train), full_posterior = "Block-VB")
  pred_fass <- predict_approx(fit_fass, X_test, y_test)
  expect_true(is.list(pred_fass))
  expect_equal(length(pred_fass$summary$mean), n_test)
  expect_true(is.numeric(mean(pred_fass$CRPS)))
  expect_true(mean(pred_fass$CRPS) > 0)
})
