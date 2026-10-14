test_that("validate_regression_data catches invalid data inputs", {
  # Missing y or X
  expect_error(FASS.sampler(X = matrix(1, 10, 2)), "missing")
  expect_error(FASS.sampler(y = 1:10), "missing")
  
  # Non-matrix, non-numeric
  expect_error(FASS.sampler(y = "invalid", X = matrix(1, 10, 2), tau1 = 1, tau0 = 10, niter = 10), "numeric")
  expect_error(FASS.sampler(y = 1:10, X = list(1:10), tau1 = 1, tau0 = 10, niter = 10), "numeric matrix")
  
  # Dimension mismatch
  expect_error(FASS.sampler(y = 1:10, X = matrix(1, 5, 2), tau1 = 1, tau0 = 10, niter = 10), "Dimension mismatch")
  
  # Too few samples or features
  expect_error(FASS.sampler(y = 1:2, X = matrix(1, 2, 1), tau1 = 1, tau0 = 10, niter = 10), "at least 3")
  
  # NA / NaN / Inf values
  expect_error(FASS.sampler(y = c(1:9, NA), X = matrix(1, 10, 2), tau1 = 1, tau0 = 10, niter = 10), "Missing values")
  expect_error(FASS.sampler(y = c(1:9, Inf), X = matrix(1, 10, 2), tau1 = 1, tau0 = 10, niter = 10), "Infinite values")
  
  # Zero variance response
  expect_error(FASS.sampler(y = rep(5, 10), X = matrix(rnorm(20), 10, 2), tau1 = 1, tau0 = 10, niter = 10), "zero variance")
  
  # Invertibility check: n <= m
  expect_error(FASS.sampler(y = 1:5, X = matrix(rnorm(50), 5, 10), tau1 = 1, tau0 = 10, niter = 10), "exceed number of candidate features")
  
  # data.frame coercion works
  df_X <- as.data.frame(matrix(rnorm(50), 25, 2))
  y_val <- rnorm(25)
  expect_error(FASS.sampler(y = y_val, X = df_X, tau1 = 1, tau0 = 10, niter = 5, verbose = FALSE), NA)
})

test_that("validate_fass_params catches invalid prior and sampler arguments", {
  n <- 50
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  
  # Missing tau1 or tau0
  expect_error(FASS.sampler(y, X, niter = 10, tau0 = 10), "tau1")
  expect_error(FASS.sampler(y, X, niter = 10, tau1 = 1), "tau0")
  
  # Non-positive parameters
  expect_error(FASS.sampler(y, X, niter = 10, tau1 = -1, tau0 = 10), "positive")
  expect_error(FASS.sampler(y, X, niter = 10, tau1 = 1, tau0 = 0), "positive")
  expect_error(FASS.sampler(y, X, niter = 10, tau1 = 1, tau0 = 10, a_pi = -1), "positive")
  expect_error(FASS.sampler(y, X, niter = 10, tau1 = 1, tau0 = 10, b_pi = 0), "positive")
  
  # Student-t missing or invalid nu
  expect_error(FASS.sampler(y, X, niter = 10, tau1 = 1, tau0 = 10, family = "Student"), "degrees of freedom")
  expect_error(FASS.sampler(y, X, niter = 10, tau1 = 1, tau0 = 10, family = "Student", nu = -2), "degrees of freedom")
  
  # Invalid family name
  expect_error(FASS.sampler(y, X, niter = 10, tau1 = 1, tau0 = 10, family = "Cauchy"), "Invalid prior family")
  
  # Invalid full_posterior
  expect_error(FASS.sampler(y, X, niter = 10, tau1 = 1, tau0 = 10, full_posterior = "NonExistent"), "Invalid full_posterior")
  
  # Gibbs invalid nburn or niter2
  expect_error(FASS.sampler(y, X, niter = 10, tau1 = 1, tau0 = 10, full_posterior = "Gibbs", nburn = -5), "nburn")
})

test_that("validate_gibbs_params catches invalid arguments for MCMC benchmarks", {
  n <- 30
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  
  expect_error(Gibbs(X, y, tau0_2 = -1, tau1_2 = 10, a_pi = 1, b_pi = 1), "positive")
  expect_error(Gibbs(X, y, tau0_2 = 1, tau1_2 = -10, a_pi = 1, b_pi = 1), "positive")
  expect_error(Gibbs(X, y, tau0_2 = 1, tau1_2 = 10, a_pi = 1, b_pi = 1, family = "Student"), "degrees of freedom")
  expect_error(Gibbs(X, y, tau0_2 = 1, tau1_2 = 10, a_pi = 1, b_pi = 1, family = "Unknown"), "Invalid prior family")
})

test_that("validate_predict_data catches invalid prediction inputs", {
  expect_error(predict_gibbs(list(bad = 1), matrix(1, 5, 2), 1:5), "matrices/vectors")
  expect_error(predict_approx(list(bad = 1), matrix(1, 5, 2), 1:5), "mu")
  
  # Dimension mismatch between test X and posterior coefficients
  dummy_approx <- list(mu = c(1, 2, 3), Sigma = diag(3), sigma2 = 1, alpha.out = 0)
  expect_error(predict_approx(dummy_approx, matrix(1, 5, 2), 1:5), "Dimension mismatch")
})
