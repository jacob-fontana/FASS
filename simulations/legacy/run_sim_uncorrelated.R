source("EP_functions.R")
source("VB_functions.R")
source("Laplace_Approx.R")
library(mvtnorm)

set.seed(4747)

n_seq <- c(5000, 10000, 50000, 100000, 1000000)
p_seq <- c(10, 100, 1000)
models <- list(
  laplace = list(prior_type="laplace", nu=NA),
  student_2 = list(prior_type="student_t", nu=2),
  student_5 = list(prior_type="student_t", nu=5),
  student_10 = list(prior_type="student_t", nu=10)
)
rho <- 0.0

# Helper to calculate quantiles string
format_quantiles <- function(x) {
  q <- quantile(x, probs=c(0.025, 0.975))
  m <- mean(x)
  sprintf("%.4f [%.4f, %.4f]", m, q[1], q[2])
}

# Main simulation loops
for (model_name in names(models)) {
  prior_type <- models[[model_name]]$prior_type
  nu_val <- models[[model_name]]$nu
  hyper_params <- list(tau=1, nu=nu_val)
  
  for (p in p_seq) {
    cat(sprintf("\n=== Running Model: %s, p = %d ===\n", model_name, p))
    
    # Store results for this p and model across N
    if (p == 10) {
      results_p10 <- list() # Will be structured as a data frame later
    } else {
      results_p_large <- list()
    }
    
    # Data Generation (Generate once for max N, then subset)
    set.seed(4747 + p)
    beta_true <- numeric(0)
    while(length(beta_true) < p) {
      cand <- rnorm(p, 0, 10)
      beta_true <- c(beta_true, cand[abs(cand) >= 0.1])
    }
    beta_true <- beta_true[1:p]
    
    # AR(1) Covariance Matrix
    Sigma_X <- matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        Sigma_X[i,j] <- rho^abs(i-j)
      }
    }
    
    current_n_seq <- if (p == 1000) n_seq[n_seq <= 100000] else n_seq
    n_max <- max(current_n_seq)
    K_full <- rmvnorm(n_max, mean=rep(0, p), sigma=Sigma_X)
    alpha <- rnorm(1, 0, sqrt(10))
    y_full <- rnorm(n_max, alpha + K_full %*% beta_true, 1)
    
    for (n in current_n_seq) {
      cat(sprintf("  N = %d...\n", n))
      
      # Subset data for this N
      K_raw <- K_full[1:n, , drop=FALSE]
      y_raw <- y_full[1:n]
      
      # Center data
      K <- scale(K_raw, center=TRUE, scale=FALSE)
      y <- y_raw - mean(y_raw)
      
      # Precomputations
      XtX <- t(K) %*% K
      X.y <- t(K) %*% y
      sigma2.hat <- 1
      
      # 1. EP
      t0 <- Sys.time()
      ep_res <- EP(K, y, XtX, X.y, prior_type=prior_type, hyper_params=hyper_params, intercept=TRUE)
      t1 <- Sys.time()
      time_EP <- as.numeric(difftime(t1, t0, units="secs"))
      mu_EP <- ep_res$mu
      sig_EP <- sqrt(diag(ep_res$Sigma))
      
      # 2. VB Full
      t0 <- Sys.time()
      vb_res <- VB(y, K, sigma2.hat, XtX, X.y, prior_type=prior_type, hyper_params=hyper_params, intercept=TRUE)
      t1 <- Sys.time()
      time_VB <- as.numeric(difftime(t1, t0, units="secs"))
      mu_VB <- vb_res$mu
      sig_VB <- sqrt(diag(vb_res$Sigma))
      
      # 3. VB MF
      t0 <- Sys.time()
      vbmf_res <- VB_MF(y, K, sigma2.hat, XtX, X.y, prior_type=prior_type, hyper_params=hyper_params, intercept=TRUE)
      t1 <- Sys.time()
      time_VB_MF <- as.numeric(difftime(t1, t0, units="secs"))
      mu_VB_MF <- vbmf_res$mu
      sig_VB_MF <- sqrt(diag(vbmf_res$Sigma))
      
      # 4. VB Direct
      t0 <- Sys.time()
      vbd_res <- VB_Direct(y, K, sigma2.hat, XtX, X.y, prior_type=prior_type, hyper_params=hyper_params, intercept=TRUE)
      t1 <- Sys.time()
      time_VB_Direct <- as.numeric(difftime(t1, t0, units="secs"))
      mu_VB_Direct <- vbd_res$mu
      sig_VB_Direct <- sqrt(diag(vbd_res$Sigma))
      
      # 5. Laplace Approx
      t0 <- Sys.time()
      la_res <- Laplace_Approx(y, K, sigma2.hat, XtX, X.y, prior_type=prior_type, hyper_params=hyper_params, intercept=TRUE)
      t1 <- Sys.time()
      time_LA <- as.numeric(difftime(t1, t0, units="secs"))
      mu_LA <- la_res$mu
      sig_LA <- sqrt(diag(la_res$Sigma))
      
      # Process Results for this N
      methods <- c("EP", "VB", "VB_MF", "VB_Direct", "LA")
      times <- c(time_EP, time_VB, time_VB_MF, time_VB_Direct, time_LA)
      mus <- list(mu_EP, mu_VB, mu_VB_MF, mu_VB_Direct, mu_LA)
      sigs <- list(sig_EP, sig_VB, sig_VB_MF, sig_VB_Direct, sig_LA)
      
      if (p == 10) {
        # Format similar to Posterior_Summary_Table
        col_name <- paste0("n=", formatC(n, format="d", big.mark=","))
        for (m_idx in 1:length(methods)) {
          m <- methods[m_idx]
          mu_vals <- mus[[m_idx]]
          sig_vals <- sigs[[m_idx]]
          lower <- mu_vals - 1.96 * sig_vals
          upper <- mu_vals + 1.96 * sig_vals
          
          for (j in 1:p) {
            val_str <- sprintf("%.4f [%.4f, %.4f]", mu_vals[j], lower[j], upper[j])
            
            # Add to list
            row_key <- paste(j, m, sep="_")
            if (is.null(results_p10[[row_key]])) {
              results_p10[[row_key]] <- list(Parameter = paste0("beta_", j), True_Value = sprintf("%.4f", beta_true[j]), Method = m)
            }
            results_p10[[row_key]][[col_name]] <- val_str
          }
          
          # Runtime row
          rt_key <- paste("Runtime", m, sep="_")
          if (is.null(results_p10[[rt_key]])) {
            results_p10[[rt_key]] <- list(Parameter = "Runtime", True_Value = "", Method = m)
          }
          results_p10[[rt_key]][[col_name]] <- sprintf("%.4fs", times[m_idx])
        }
      } else { # p = 100 or 1000
        col_name <- paste0("n=", formatC(n, format="d", big.mark=","))
        for (m_idx in 1:length(methods)) {
          m <- methods[m_idx]
          mu_vals <- mus[[m_idx]]
          sig_vals <- sigs[[m_idx]]
          
          # RMSE across all p parameters
          rmse_val <- sqrt(mean((mu_vals - beta_true)^2))
          
          # Coverage across all p parameters
          lowers <- mu_vals - 1.96 * sig_vals
          uppers <- mu_vals + 1.96 * sig_vals
          cov_val <- mean(beta_true >= lowers & beta_true <= uppers)
          
          # Length
          len_vec <- uppers - lowers
          len_q <- format_quantiles(len_vec)
          
          # Add to list
          metrics <- c("Runtime", "Coverage", "RMSE", "Length (Mean [95% CI])")
          vals <- c(sprintf("%.4fs", times[m_idx]), sprintf("%.4f", cov_val), sprintf("%.4f", rmse_val), len_q)
          
          for (k in 1:4) {
            row_key <- paste(m, metrics[k], sep="_")
            if (is.null(results_p_large[[row_key]])) {
              results_p_large[[row_key]] <- list(Method = m, Metric = metrics[k])
            }
            results_p_large[[row_key]][[col_name]] <- vals[k]
          }
        }
      }
    }
    
    # Save CSVs
    if (p == 10) {
      df <- do.call(rbind, lapply(results_p10, function(x) as.data.frame(x, check.names=FALSE, stringsAsFactors=FALSE)))
      write.csv(df, file.path("csv", "uncorrelated", paste0("Uncorrelated_Sim_", model_name, "_p10.csv")), row.names=FALSE)
    } else {
      df <- do.call(rbind, lapply(results_p_large, function(x) as.data.frame(x, check.names=FALSE, stringsAsFactors=FALSE)))
      write.csv(df, file.path("csv", "uncorrelated", paste0("Uncorrelated_Sim_", model_name, "_p", p, ".csv")), row.names=FALSE)
    }
  }
}

cat("\nAll simulations complete.\n")
