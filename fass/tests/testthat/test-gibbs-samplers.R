test_that("Full Gibbs benchmark samplers run across families and produce valid chains", {
  set.seed(42)
  n <- 50
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 2 + rnorm(n)
  
  # Standard Gibbs
  fit_g <- Gibbs(X, y, tau0_2 = 1/n, tau1_2 = n, a_pi = 1, b_pi = 1,
                 nburn = 50, niter = 100, family = "Laplace", rescale_y = FALSE)
  expect_equal(dim(fit_g$outbeta), c(100, p))
  expect_equal(dim(fit_g$outgamma), c(100, p))
  expect_equal(length(fit_g$outsigma), 100)
  
  # Partially Collapsed Gibbs
  fit_pc <- Gibbs_PC(X, y, tau0_2 = 1/n, tau1_2 = n, a_pi = 1, b_pi = 1,
                     nburn = 50, niter = 100, family = "Gaussian", rescale_y = FALSE)
  expect_equal(dim(fit_pc$outbeta), c(100, p))
  expect_equal(dim(fit_pc$outgamma), c(100, p))
  
  # Skinny Gibbs
  fit_skinny <- skinnyGibbs(X, y, tau0_2 = 1/n, tau1_2 = n, a_pi = 1, b_pi = 1,
                           nburn = 50, niter = 100, family = "Student", nu = 4, rescale_y = FALSE)
  expect_equal(dim(fit_skinny$outbeta), c(100, p))
  expect_equal(dim(fit_skinny$outgamma), c(100, p))
})

test_that("Gibbs samplers respect rescale_y = TRUE and intercept = FALSE", {
  set.seed(77)
  n <- 40
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- (X[, 1] * 5 + rnorm(n)) * 10
  
  # Rescale y = TRUE
  fit_rescale <- Gibbs_PC(X, y, tau0_2 = 1/n, tau1_2 = n, a_pi = 1, b_pi = 1,
                          nburn = 30, niter = 50, family = "Laplace",
                          rescale_y = TRUE, intercept = FALSE)
  
  expect_equal(dim(fit_rescale$outbeta), c(50, p))
  # Point estimate should be scaled back to original data scale (mean beta_1 > 10)
  expect_true(mean(fit_rescale$outbeta[, 1]) > 10)
})
