test_that("FASS.sampler runs across all families and selects active variables", {
  set.seed(42)
  n <- 120
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 2.5 - X[, 2] * 2.0 + rnorm(n)
  
  # 1. Laplace family
  fit_lap <- FASS.sampler(y, X, niter = 150, family = "Laplace", tau1 = 1, tau0 = sqrt(n), rescale_y = FALSE)
  expect_true(is.list(fit_lap))
  expect_equal(dim(fit_lap$gamma), c(150, p))
  top_model_lap <- which.max(fit_lap$odds)
  selected_lap <- which(fit_lap$gamma[top_model_lap - 1, ] == 1)
  expect_true(all(c(1, 2) %in% selected_lap))
  
  # 2. Gaussian family
  fit_gauss <- FASS.sampler(y, X, niter = 150, family = "Gaussian", tau1 = 1, tau0 = n, rescale_y = FALSE)
  expect_true(is.list(fit_gauss))
  top_model_gauss <- which.max(fit_gauss$odds)
  selected_gauss <- which(fit_gauss$gamma[top_model_gauss - 1, ] == 1)
  expect_true(all(c(1, 2) %in% selected_gauss))
  
  # 3. Student-t family (nu = 5)
  fit_stud <- FASS.sampler(y, X, niter = 150, family = "Student", tau1 = 1, tau0 = n, nu = 5, rescale_y = FALSE)
  expect_true(is.list(fit_stud))
  top_model_stud <- which.max(fit_stud$odds)
  selected_stud <- which(fit_stud$gamma[top_model_stud - 1, ] == 1)
  expect_true(all(c(1, 2) %in% selected_stud))
})

test_that("FASS.sampler correctly respects intercept = FALSE and rescale_y = TRUE", {
  set.seed(123)
  n <- 80
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  # Data with no intercept and large scale
  y <- (X[, 1] * 3.0 + X[, 2] * 2.0 + rnorm(n)) * 10
  
  # Intercept = FALSE
  fit_no_int <- FASS.sampler(y, X, niter = 100, family = "Laplace", tau1 = 1, tau0 = sqrt(n),
                             intercept = FALSE, full_posterior = "LA")
  expect_true(is.list(fit_no_int$posterior))
  
  # Rescale_y = TRUE
  fit_rescale <- FASS.sampler(y, X, niter = 100, family = "Laplace", tau1 = 1, tau0 = sqrt(n),
                              rescale_y = TRUE, full_posterior = "Block-VB")
  expect_true(is.list(fit_rescale$posterior))
  # Check that posterior point estimate was scaled back to original scale (approx > 5, not normalized ~ 0.5)
  expect_true(abs(fit_rescale$posterior$mu[1]) > 5)
})

test_that("FASS.sampler handles degree of freedom boundary without crashing", {
  set.seed(99)
  p <- 10
  n <- p + 4 # Extreme boundary where n - p - 3 = 1
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  
  expect_error({
    fit_bound <- FASS.sampler(y, X, niter = 50, family = "Laplace", tau1 = 1, tau0 = sqrt(n),
                              verbose = FALSE, intercept = TRUE)
  }, NA)
})

test_that("FASS.sampler handles collinear design matrices gracefully", {
  n <- 50
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  X[, 5] <- X[, 1] + X[, 2] # Exactly collinear column
  y <- rnorm(n)
  
  expect_error(
    FASS.sampler(y, X, niter = 50, family = "Laplace", tau1 = 1, tau0 = sqrt(n)),
    "collinear"
  )
})

test_that("FASS.sampler accepts features with small natural scales without false collinearity trigger", {
  set.seed(42)
  n <- 50
  p <- 5
  # Columns with natural standard deviation on the order of 1e-6
  X <- matrix(rnorm(n * p, sd = 1e-6), n, p)
  y <- rnorm(n)
  
  # Old check min(diag(R)) < 1e-5 would falsely fail here; scale-invariant check succeeds
  expect_error({
    fit <- FASS.sampler(y, X, niter = 30, family = "Laplace", tau1 = 1, tau0 = sqrt(n))
  }, NA)
})

test_that("FASS.sampler integrates all post-hoc approximation methods", {
  set.seed(111)
  n <- 100
  p <- 8
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 2 + rnorm(n)
  
  for (method in c("LA", "Block-VB", "MFVB-CAVI", "MMFVB", "EP")) {
    fit <- FASS.sampler(y, X, niter = 50, family = "Laplace", tau1 = 1, tau0 = sqrt(n),
                        full_posterior = method)
    expect_true(is.list(fit$posterior))
    expect_true(!is.null(fit$posterior$mu))
    expect_true(!is.null(fit$posterior$Sigma))
  }
  
  # Gibbs post-hoc
  fit_gibbs <- FASS.sampler(y, X, niter = 50, family = "Laplace", tau1 = 1, tau0 = sqrt(n),
                            full_posterior = "Gibbs", nburn = 50, niter2 = 50)
  expect_true(is.list(fit_gibbs$posterior))
  expect_true(!is.null(fit_gibbs$posterior$outbeta))
})
