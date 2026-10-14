test_that("Golub-Welsch Gauss-Hermite quadrature integrates polynomial moments accurately", {
  # Direct access to internal quadrature function
  gh <- fass:::gauss_hermite_nodes_weights(21)
  
  # Integral of e^{-x^2} dx = sqrt(pi)
  expect_equal(sum(gh$weights), sqrt(pi), tolerance = 1e-12)
  
  # Integral of x^2 e^{-x^2} dx = sqrt(pi) / 2
  mom2 <- sum(gh$weights * (gh$nodes^2))
  expect_equal(mom2, sqrt(pi) / 2, tolerance = 1e-12)
  
  # Integral of x^4 e^{-x^2} dx = 3 * sqrt(pi) / 4
  mom4 <- sum(gh$weights * (gh$nodes^4))
  expect_equal(mom4, 3 * sqrt(pi) / 4, tolerance = 1e-12)
})

test_that("fit_approx_posterior executes all methods and prior types", {
  set.seed(42)
  n <- 80
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(1.5, -1.0, 0.5, 0.0)
  y <- X %*% beta_true + rnorm(n)
  
  XtX <- crossprod(X)
  Xy <- as.vector(crossprod(X, y))
  sigma2_hat <- 1.0
  
  priors <- list(
    list(type = "laplace", hyper = list(tau = 1)),
    list(type = "gaussian", hyper = list(tau = 1)),
    list(type = "student_t", hyper = list(tau = 1, nu = 4))
  )
  
  methods <- c("LA", "Block-VB", "MFVB-CAVI", "MMFVB", "EP")
  
  for (pr in priors) {
    for (m in methods) {
      fit <- fit_approx_posterior(
        y = y, X = X, method = m,
        prior_type = pr$type,
        hyper_params = pr$hyper,
        intercept = TRUE,
        sigma2.hat = sigma2_hat,
        XtX = XtX, X.y = Xy
      )
      
      expect_true(is.list(fit))
      expect_equal(length(fit$mu), p)
      expect_equal(dim(fit$Sigma), c(p, p))
      expect_true(is.numeric(fit$sigma2))
      expect_true(fit$sigma2 > 0)
      
      # Positive definite covariance
      eigen_vals <- eigen(fit$Sigma, symmetric = TRUE, only.values = TRUE)$values
      expect_true(all(eigen_vals > 0))
    }
  }
})

test_that("fit_approx_posterior handles intercept = FALSE properly", {
  set.seed(99)
  n <- 60
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- X %*% c(2, -1, 1) + rnorm(n)
  
  fit_noint <- fit_approx_posterior(
    y = y, X = X, method = "Block-VB",
    prior_type = "laplace",
    hyper_params = list(tau = 1),
    intercept = FALSE,
    sigma2.hat = 1.0
  )
  expect_true(is.list(fit_noint))
  expect_equal(length(fit_noint$mu), p)
})
