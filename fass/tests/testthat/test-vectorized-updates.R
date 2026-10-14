test_that("Vectorized updates produce correct dimensions, terminal zero, and no warnings", {
  n <- 100
  m <- 15
  a_pi <- 1
  b_pi <- 1
  tau0 <- sqrt(n)
  tau1 <- 1
  
  set.seed(42)
  X <- matrix(rnorm(n * m), n, m)
  y <- rnorm(n)
  
  X2 <- crossprod(X)
  Q <- n * chol2inv(chol(X2))
  beta.hat <- as.vector((1 / n) * Q %*% crossprod(X, y))
  sigma.hat <- sum((y - mean(y) - X %*% beta.hat)^2) / (n - m - 3)
  
  gamma <- sample(c(0, 1), m, replace = TRUE)
  gamma.sum <- sum(gamma == 1)
  tau <- rep(tau0, m)
  tau[gamma == 1] <- tau1
  
  # 1. Laplace Update
  tau.star.lap <- tau * sign(beta.hat)
  Q.tau1.lap <- tau.star.lap %*% Q
  Q.tau.lap <- as.numeric(Q.tau1.lap %*% tau.star.lap)
  
  expect_warning({
    odds_lap <- update.laplace(n, m, a_pi, b_pi, tau0, tau1, gamma, gamma.sum,
                               beta.hat, sigma.hat, Q, Q.tau1.lap, Q.tau.lap, tau.star.lap)
  }, NA)
  expect_equal(length(odds_lap), m + 1)
  expect_equal(odds_lap[m + 1], 0)
  expect_false(any(is.nan(odds_lap)))
  expect_false(any(is.infinite(odds_lap)))
  
  # 2. Gaussian Update
  tau.star.gauss <- tau * beta.hat
  Q.tau1.gauss <- tau.star.gauss %*% Q
  Q.tau.gauss <- as.numeric(Q.tau1.gauss %*% tau.star.gauss - sum(tau * diag(Q)))
  
  expect_warning({
    odds_gauss <- update.gauss(n, m, a_pi, b_pi, tau0, tau1, gamma, gamma.sum,
                               beta.hat, sigma.hat, Q, Q.tau1.gauss, Q.tau.gauss, tau.star.gauss)
  }, NA)
  expect_equal(length(odds_gauss), m + 1)
  expect_equal(odds_gauss[m + 1], 0)
  expect_false(any(is.nan(odds_gauss)))
  expect_false(any(is.infinite(odds_gauss)))
  
  # 3. Student-t Update (testing df = 2, 5, 10 for numerical stability and zero NaNs)
  for (nu in c(2, 5, 10)) {
    tau.star.stud <- (nu + 1) * beta.hat / (beta.hat^2 + nu / tau)
    Q.tau1.stud <- tau.star.stud %*% Q
    Q.tau.stud <- as.numeric(Q.tau1.stud %*% tau.star.stud - sum(diag(Q) * (nu + 1) * (nu / tau - beta.hat^2) / ((beta.hat^2 + nu / tau)^2)))
    
    expect_warning({
      odds_stud <- update.student(n, m, a_pi, b_pi, tau0, tau1, gamma, gamma.sum,
                                  beta.hat, sigma.hat, Q, Q.tau1.stud, Q.tau.stud, tau.star.stud, nu = nu)
    }, NA)
    expect_equal(length(odds_stud), m + 1)
    expect_equal(odds_stud[m + 1], 0)
    expect_false(any(is.nan(odds_stud)))
    expect_false(any(is.infinite(odds_stud)))
  }
})

test_that("Vectorized Laplace odds match unvectorized reference calculation", {
  n <- 80
  m <- 10
  a_pi <- 1
  b_pi <- 1
  tau0 <- sqrt(n)
  tau1 <- 1
  
  set.seed(123)
  X <- matrix(rnorm(n * m), n, m)
  y <- rnorm(n)
  
  Q <- n * chol2inv(chol(crossprod(X)))
  beta.hat <- as.vector((1 / n) * Q %*% crossprod(X, y))
  sigma.hat <- sum((y - mean(y) - X %*% beta.hat)^2) / (n - m - 3)
  
  gamma <- c(1, 0, 1, 0, 0, 1, 0, 0, 0, 1)
  gamma.sum <- sum(gamma)
  tau <- rep(tau0, m)
  tau[gamma == 1] <- tau1
  
  tau.star <- tau * sign(beta.hat)
  Q.tau1 <- tau.star %*% Q
  Q.tau <- as.numeric(Q.tau1 %*% tau.star)
  
  odds_vec <- update.laplace(n, m, a_pi, b_pi, tau0, tau1, gamma, gamma.sum,
                             beta.hat, sigma.hat, Q, Q.tau1, Q.tau, tau.star)
  
  # Manual element-by-element loop calculation
  temp2 <- 1 + Q.tau * sigma.hat / (2 * n)
  odds_loop <- numeric(m)
  for (j in 1:m) {
    if (gamma[j] == 0) {
      temp1 <- log(a_pi + gamma.sum) - log(b_pi + m - gamma.sum - 1) + 
        log(tau1) - log(tau0) - (tau1 - tau0) * abs(beta.hat[j])
      t3 <- 1 + (Q.tau + 2 * (Q.tau1[j] - tau.star[j] * Q[j, j]) * (tau1 - tau0) * sign(beta.hat[j]) + 
                   (tau1^2 - tau0^2) * Q[j, j]) * sigma.hat / (2 * n)
      odds_loop[j] <- temp1 + log(max(1e-12, t3)) - log(max(1e-12, temp2))
    } else {
      temp1 <- log(a_pi + gamma.sum - 1) - log(b_pi + m - gamma.sum) + 
        log(tau1) - log(tau0) - (tau1 - tau0) * abs(beta.hat[j])
      t3 <- 1 + (Q.tau + 2 * (Q.tau1[j] - tau.star[j] * Q[j, j]) * (tau0 - tau1) * sign(beta.hat[j]) + 
                   (tau0^2 - tau1^2) * Q[j, j]) * sigma.hat / (2 * n)
      odds_loop[j] <- -temp1 - log(max(1e-12, temp2)) + log(max(1e-12, t3))
    }
  }
  
  expect_equal(odds_vec[1:m], odds_loop, tolerance = 1e-10)
})
