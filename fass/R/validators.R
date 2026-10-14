################################################################################
### Input Validation Helpers (Internal)
################################################################################

#' Validate Regression Data Matrix and Response
#'
#' @param y Numeric response vector
#' @param X Candidate design matrix or data.frame
#' @param check_invertibility Logical, whether to check that n > p
#' @return List containing coerced numeric matrix X and numeric vector y
#' @noRd
validate_regression_data <- function(y, X, check_invertibility = TRUE) {
  if (missing(y) || is.null(y)) {
    stop("Response vector 'y' is missing with no default.")
  }
  if (missing(X) || is.null(X)) {
    stop("Design matrix 'X' is missing with no default.")
  }
  
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(X)) {
    stop("'X' must be a numeric matrix or data.frame.")
  }
  if (!is.numeric(y)) {
    stop("Response 'y' must be a numeric vector.")
  }
  if (!is.numeric(X)) {
    stop("Design matrix 'X' must contain numeric values.")
  }
  
  n <- length(y)
  if (nrow(X) != n) {
    stop(sprintf("Dimension mismatch: length of 'y' (%d) does not match rows of 'X' (%d).", n, nrow(X)))
  }
  if (n < 3) {
    stop("Sample size 'n' must be at least 3.")
  }
  if (ncol(X) < 1) {
    stop("Design matrix 'X' must contain at least 1 column.")
  }
  
  if (anyNA(y) || anyNA(X)) {
    stop("Missing values (NA / NaN) detected in 'y' or 'X'. Please clean or impute before fitting.")
  }
  if (any(is.infinite(y)) || any(is.infinite(X))) {
    stop("Infinite values (Inf / -Inf) detected in 'y' or 'X'.")
  }
  
  if (stats::sd(y) <= 0) {
    stop("Response 'y' has zero variance (constant response).")
  }
  
  if (check_invertibility && nrow(X) <= ncol(X)) {
    stop(sprintf(
      "Sample size n (%d) must exceed number of candidate features m (%d) for standard FASS Gram inversion.",
      nrow(X), ncol(X)
    ))
  }
  
  return(list(y = as.vector(y), X = X))
}

#' Validate FASS Hyperparameters and Options
#'
#' @param niter Number of shotgun iterations
#' @param family Prior slab distribution family
#' @param tau1 Slab precision
#' @param tau0 Spike precision
#' @param nu Degrees of freedom for Student-t
#' @param a_pi Beta prior hyperparameter
#' @param b_pi Beta prior hyperparameter
#' @param full_posterior Post-hoc estimation method
#' @param nburn Burn-in iterations
#' @param niter2 Post-burn iterations
#' @return List of validated parameters
#' @noRd
validate_fass_params <- function(niter, family, tau1, tau0, nu, a_pi, b_pi, 
                                 full_posterior, nburn, niter2) {
  if (!is.numeric(niter) || length(niter) != 1 || niter < 1) {
    stop("'niter' must be a single positive integer >= 1.")
  }
  
  if (missing(tau1) || is.null(tau1)) {
    stop("Prior parameter 'tau1' is missing. Please specify slab precision/rate 'tau1'.")
  }
  if (missing(tau0) || is.null(tau0)) {
    stop("Prior parameter 'tau0' is missing. Please specify spike precision/rate 'tau0'.")
  }
  if (!is.numeric(tau1) || tau1 <= 0) {
    stop("Prior parameter 'tau1' must be a positive number (> 0).")
  }
  if (!is.numeric(tau0) || tau0 <= 0) {
    stop("Prior parameter 'tau0' must be a positive number (> 0).")
  }
  
  family_clean <- tolower(trimws(family[1]))
  valid_families <- c("laplace", "gaussian", "student")
  if (!family_clean %in% valid_families) {
    stop(sprintf("Invalid prior family '%s'. Must be one of: 'Laplace', 'Gaussian', 'Student'.", family[1]))
  }
  # Capitalize canonically
  family_canonical <- switch(family_clean,
                             "laplace" = "Laplace",
                             "gaussian" = "Gaussian",
                             "student" = "Student")
  
  if (family_canonical == "Student") {
    if (is.null(nu) || is.na(nu) || !is.numeric(nu) || nu <= 0) {
      stop("When family = 'Student', 'nu' (degrees of freedom) must be a positive number (> 0).")
    }
  }
  
  if (!is.numeric(a_pi) || a_pi <= 0 || !is.numeric(b_pi) || b_pi <= 0) {
    stop("Beta inclusion prior parameters 'a_pi' and 'b_pi' must be positive numbers (> 0).")
  }
  
  valid_posteriors <- c("Gibbs", "EP", "Block-VB", "MMFVB", "MFVB-CAVI", "LA")
  # NULL, NA and FALSE all mean "selection only"
  if (is.null(full_posterior) || length(full_posterior) != 1 || isFALSE(full_posterior) || is.na(full_posterior)) {
    if (length(full_posterior) > 1) {
      stop(sprintf("'full_posterior' must be a single value: NA or one of %s.",
                   paste(valid_posteriors, collapse = ", ")))
    }
    full_posterior <- NA
  } else {
    if (!full_posterior %in% valid_posteriors) {
      stop(sprintf("Invalid full_posterior '%s'. Must be one of: %s.", 
                   full_posterior, paste(valid_posteriors, collapse = ", ")))
    }
    if (identical(full_posterior, "Gibbs")) {
      if (is.null(nburn) || is.na(nburn)) nburn <- 1000
      if (is.null(niter2) || is.na(niter2)) niter2 <- 1000
      if (!is.numeric(nburn) || nburn < 0 || !is.numeric(niter2) || niter2 < 1) {
        stop("For full_posterior = 'Gibbs', 'nburn' must be >= 0 and 'niter2' must be >= 1.")
      }
    }
  }
  
  return(list(
    niter = as.integer(niter),
    family = family_canonical,
    tau1 = as.numeric(tau1),
    tau0 = as.numeric(tau0),
    nu = if (is.numeric(nu)) as.numeric(nu) else NA,
    a_pi = as.numeric(a_pi),
    b_pi = as.numeric(b_pi),
    full_posterior = full_posterior,
    nburn = nburn,
    niter2 = niter2
  ))
}

#' Validate Gibbs Sampler Parameters
#'
#' @inheritParams Gibbs
#' @noRd
validate_gibbs_params <- function(tau0_2, tau1_2, a_pi, b_pi, nburn, niter, family, nu) {
  if (missing(tau0_2) || !is.numeric(tau0_2) || tau0_2 <= 0) {
    stop("'tau0_2' (spike prior variance) must be a positive number (> 0).")
  }
  if (missing(tau1_2) || !is.numeric(tau1_2) || tau1_2 <= 0) {
    stop("'tau1_2' (slab prior variance) must be a positive number (> 0).")
  }
  if (missing(a_pi) || !is.numeric(a_pi) || a_pi <= 0 || missing(b_pi) || !is.numeric(b_pi) || b_pi <= 0) {
    stop("Beta prior parameters 'a_pi' and 'b_pi' must be positive numbers (> 0).")
  }
  if (!is.numeric(nburn) || nburn < 0) {
    stop("'nburn' (burn-in iterations) must be a non-negative integer >= 0.")
  }
  if (!is.numeric(niter) || niter < 1) {
    stop("'niter' (sampling iterations) must be a positive integer >= 1.")
  }
  
  family_clean <- tolower(trimws(family[1]))
  valid_families <- c("laplace", "gaussian", "student")
  if (!family_clean %in% valid_families) {
    stop(sprintf("Invalid prior family '%s'. Must be one of: 'Laplace', 'Gaussian', 'Student'.", family[1]))
  }
  family_canonical <- switch(family_clean,
                             "laplace" = "Laplace",
                             "gaussian" = "Gaussian",
                             "student" = "Student")
  
  if (family_canonical == "Student") {
    if (is.null(nu) || is.na(nu) || !is.numeric(nu) || nu <= 0) {
      stop("When family = 'Student', 'nu' (degrees of freedom) must be a positive number (> 0).")
    }
  }
  
  return(list(
    tau0_2 = as.numeric(tau0_2),
    tau1_2 = as.numeric(tau1_2),
    a_pi = as.numeric(a_pi),
    b_pi = as.numeric(b_pi),
    nburn = as.integer(nburn),
    niter = as.integer(niter),
    family = family_canonical,
    nu = if (is.numeric(nu)) as.numeric(nu) else NA
  ))
}

#' Validate Out-of-Sample Prediction Inputs
#'
#' @param output Fitted model output list
#' @param X Out-of-sample design matrix
#' @param y Out-of-sample response vector
#' @noRd
validate_predict_data <- function(output, X, y) {
  if (missing(output) || is.null(output) || !is.list(output)) {
    stop("Model 'output' must be a valid list from FASS.sampler, Gibbs, or fit_approx_posterior.")
  }
  if (missing(X) || is.null(X)) {
    stop("Test design matrix 'X' is missing.")
  }
  if (missing(y) || is.null(y)) {
    stop("Test response vector 'y' is missing.")
  }
  
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("Test 'X' must be a numeric matrix or data.frame.")
  }
  if (!is.numeric(y)) {
    stop("Test 'y' must be a numeric vector.")
  }
  if (length(y) != nrow(X)) {
    stop(sprintf("Length of test response 'y' (%d) does not match rows of test 'X' (%d).", length(y), nrow(X)))
  }
  if (anyNA(y) || anyNA(X)) {
    stop("Missing values (NA / NaN) detected in test 'X' or 'y'.")
  }
  
  return(list(y = as.vector(y), X = X))
}
