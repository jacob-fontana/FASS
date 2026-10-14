test_that("FASS handles null model edge case where no predictors are selected", {
  set.seed(999)
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  # Pure noise response completely independent of X
  y <- rnorm(n)
  
  # Run FASS with very strict prior penalization to encourage null model
  fit_null <- FASS.sampler(y, X, niter = 100, family = "Laplace", tau1 = 0.01, tau0 = 1000,
                           full_posterior = "Block-VB")
  
  expect_true(is.list(fit_null))
  expect_true(is.numeric(fit_null$odds))
  # Should successfully handle posterior without crashing
  expect_true(is.list(fit_null$posterior))
})

test_that("FASS handles all-active edge case where all predictors have massive signal", {
  set.seed(42)
  n <- 100
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  # Giant signals so all variables must be selected
  y <- X %*% c(100, -100, 80, -80) + rnorm(n, sd = 0.1)
  
  fit_all <- FASS.sampler(y, X, niter = 100, family = "Laplace", tau1 = 1, tau0 = sqrt(n),
                          full_posterior = "LA")
  
  top_idx <- which.max(fit_all$odds)
  selected <- which(fit_all$gamma[top_idx - 1, ] == 1)
  expect_equal(length(selected), p)
  expect_true(is.list(fit_all$posterior))
  expect_equal(length(fit_all$posterior$mu), p)
})

test_that("FASS handles single-predictor model (p = 1)", {
  set.seed(123)
  n <- 50
  p <- 1
  X <- matrix(rnorm(n * p), n, p)
  y <- X * 3.5 + rnorm(n)
  
  expect_error({
    fit_p1 <- FASS.sampler(y, X, niter = 30, family = "Laplace", tau1 = 1, tau0 = sqrt(n),
                           full_posterior = "Block-VB")
  }, NA)
})

test_that("FASS handles high-dimensional boundary where n is close to p", {
  set.seed(777)
  p <- 15
  n <- 20 # n - p = 5, n - p - 3 = 2
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 2 + rnorm(n)
  
  expect_error({
    fit_close <- FASS.sampler(y, X, niter = 40, family = "Gaussian", tau1 = 1, tau0 = n,
                              full_posterior = "LA")
  }, NA)
})

test_that("FASS handles negative coefficients and sign-switching with Laplace prior", {
  set.seed(888)
  n <- 80
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  # Explicitly strong negative coefficients
  y <- X %*% c(-5, -4, -3, 0, 0, 0) + rnorm(n)
  
  fit_neg <- FASS.sampler(y, X, niter = 100, family = "Laplace", tau1 = 1, tau0 = sqrt(n),
                          full_posterior = "LA")
  
  top_idx <- which.max(fit_neg$odds)
  selected <- which(fit_neg$gamma[top_idx - 1, ] == 1)
  expect_true(all(c(1, 2, 3) %in% selected))
  # Point estimates should be negative
  expect_true(all(fit_neg$posterior$mu[1:3] < 0))
})
