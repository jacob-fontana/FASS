test_that("centering prevents false positives with shifted covariates and intercept", {
  set.seed(42)
  n <- 100
  p <- 15
  # Covariates with non-zero mean = 5
  X <- matrix(rnorm(n * p, mean = 5, sd = 1), n, p)
  beta_true <- c(3, -2.5, rep(0, p - 2))
  alpha_true <- 20
  y <- as.vector(alpha_true + X %*% beta_true + rnorm(n, 0, 1))

  fit <- FASS.sampler(
    y = y,
    X = X,
    niter = 200,
    family = "Laplace",
    tau1 = 1,
    tau0 = sqrt(n),
    intercept = TRUE,
    full_posterior = "LA"
  )

  # Check that posterior and intercept are present
  expect_true(!is.null(fit$posterior))
  expect_true(!is.null(fit$posterior$alpha))
  expect_true(abs(fit$posterior$alpha - alpha_true) < 5)

  # Active variables (1 and 2) should be selected in top model
  top_model <- which.max(fit$odds)
  selected <- which(fit$gamma[top_model - 1, ] == 1)
  expect_true(all(c(1, 2) %in% selected))
})

test_that("rinversegauss is numerically stable for extreme small beta and mu", {
  set.seed(123)
  # beta very small -> mu is huge (mu = 1 / (tau * beta))
  # Should never return NaN or Inf
  samples_huge_mu <- fass:::rinversegauss(100, mu = 1e8, lambda = 1)
  expect_false(any(is.na(samples_huge_mu)))
  expect_false(any(is.nan(samples_huge_mu)))
  expect_true(all(samples_huge_mu > 0))

  # mu very small
  samples_small_mu <- fass:::rinversegauss(100, mu = 1e-10, lambda = 1)
  expect_false(any(is.na(samples_small_mu)))
  expect_false(any(is.nan(samples_small_mu)))
  expect_true(all(samples_small_mu > 0))
})

test_that("approximation methods accept partial opts without crashing via fit_approx_posterior", {
  set.seed(42)
  n <- 50
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- as.vector(X %*% c(1, 2, -1) + rnorm(n))

  partial_opts <- list(max_iter = 15)

  for (m in c("LA", "Block-VB", "MFVB-CAVI", "MMFVB", "EP")) {
    expect_error(
      fit_approx_posterior(
        y = y, X = X, method = m,
        prior_type = "laplace",
        hyper_params = list(tau = 1),
        intercept = TRUE,
        opts = partial_opts
      ),
      NA
    )
  }
})

test_that("predict_approx accurately applies intercept and handles null selection", {
  set.seed(42)
  n <- 50
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  y <- as.vector(10 + X %*% c(2, 0, 0, 0) + rnorm(n, 0, 0.1))

  fit <- FASS.sampler(
    y = y, X = X, niter = 100,
    family = "Laplace", tau1 = 1, tau0 = sqrt(n),
    intercept = TRUE, full_posterior = "Block-VB"
  )

  X_new <- matrix(rnorm(5 * p), 5, p)
  y_new <- rnorm(5, mean = 10)
  pred <- predict_approx(fit, X = X_new, y = y_new)
  expect_equal(length(pred$summary$mean), 5)
  # Predictions should be in the ballpark of y (around 10)
  expect_true(all(abs(pred$summary$mean - 10) < 15))

  # Null selection edge case (model where empty model has highest odds)
  null_post <- list(
    gamma = matrix(0, 2, p),
    odds = c(10, 1), # max_idx == 1 is the null model
    posterior = list(mu = numeric(0), Sigma = matrix(0, 0, 0), alpha = 5, sigma2 = 1)
  )
  null_pred <- predict_approx(null_post, X = X_new, y = y_new)
  expect_equal(null_pred$summary$mean, rep(5, 5))
})

test_that("calc_rss_cpp bounds checking protects against invalid indices", {
  set.seed(42)
  n <- 20
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  beta <- rnorm(p)

  # Column out of bounds (1-based: col 10 does not exist; negative col)
  expect_error(calc_rss_cpp(X, y, beta, c(1, 2, 10)), NA)
  expect_error(calc_rss_cpp(X, y, beta, c(-1, 0, 2)), NA)
})

test_that("Gibbs samplers handle the intercept correctly with uncentred covariates", {
  set.seed(2026)
  n <- 400
  p <- 4
  X <- matrix(rnorm(n * p, mean = 5), n, p)
  y <- as.vector(20 + X %*% c(2, 0, 0, -1) + rnorm(n))
  for (sampler in list(Gibbs, Gibbs_PC, skinnyGibbs)) {
    fit <- sampler(X, y, tau0_2 = 1 / n, tau1_2 = 1, a_pi = 1, b_pi = 1,
                   nburn = 200, niter = 300, family = "Gaussian")
    expect_equal(mean(fit$outbeta[, 1]), 2, tolerance = 0.05)
    expect_equal(mean(fit$outalpha), 20, tolerance = 0.05)
  }
})

test_that("fit_approx_posterior returns the intercept for uncentred covariates", {
  set.seed(7)
  n <- 500
  X <- matrix(rnorm(n * 3, mean = 3), n, 3)
  y <- as.vector(-4 + X %*% c(1, -2, 0.5) + rnorm(n))
  fit <- fit_approx_posterior(y, X, method = "Block-VB", prior_type = "student_t",
                              hyper_params = list(tau = 1, nu = 5))
  ols <- coef(lm(y ~ X))
  expect_equal(fit$alpha, unname(ols[1]), tolerance = 0.05)
  expect_equal(as.vector(fit$mu), unname(ols[-1]), tolerance = 0.05)
})

test_that("predict_approx works with FASS post-hoc Gibbs draws", {
  set.seed(11)
  n <- 300
  p <- 5
  X <- matrix(rnorm(n * p, mean = 2), n, p)
  y <- as.vector(3 + 2 * X[, 1] + rnorm(n))
  fit <- FASS.sampler(y, X, niter = 100, family = "Student", nu = 5, tau1 = 1, tau0 = n,
                      full_posterior = "Gibbs", nburn = 100, niter2 = 200)
  X_new <- matrix(rnorm(20 * p, mean = 2), 20, p)
  y_new <- as.vector(3 + 2 * X_new[, 1] + rnorm(20))
  pred <- predict_approx(fit, X_new, y_new)
  expect_lt(pred$RMSE, 2)
})

test_that("full_posterior accepts NA, NULL and FALSE as selection only", {
  set.seed(3)
  X <- matrix(rnorm(200), 50, 4)
  y <- rnorm(50)
  for (fp in list(NA, NULL, FALSE)) {
    fit <- FASS.sampler(y, X, niter = 20, family = "Laplace", tau1 = 1, tau0 = sqrt(50),
                        full_posterior = fp)
    expect_true(identical(fit$posterior, NA))
  }
})
