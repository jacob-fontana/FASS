################################################################################
### Simulation 03: Post-Hoc Posterior Inference Benchmark
### Compares EP, Block-VB, MFVB-CAVI, MMFVB, and Laplace Approximation for parameter estimation
################################################################################

suppressPackageStartupMessages({
  library(fass)
})

# Source shared simulation helpers
source("simulations/R/dgp_generators.R")
source("simulations/R/sim_helpers.R")

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
# Set run_mode: "test" for quick verification, "full" for the complete paper grid.
# Can be overridden from the command line: Rscript simulations/03_posthoc_inference.R full
args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args) >= 1) args[1] else "test"
stopifnot(run_mode %in% c("test", "full"))

if (run_mode == "test") {
  n_seq <- c(500, 2000)
  p_seq <- c(10, 20)
  rho_seq <- c(0.0, 0.8) # Uncorrelated and Correlated
} else {
  n_seq <- c(5000, 10000, 50000, 100000)
  p_seq <- c(10, 100)
  rho_seq <- c(0.0, 0.8)
}

models <- list(
  list(name = "Laplace", prior_type = "laplace", nu = NULL),
  list(name = "Student_df2", prior_type = "student_t", nu = 2),
  list(name = "Student_df5", prior_type = "student_t", nu = 5),
  list(name = "Student_df10", prior_type = "student_t", nu = 10)
)

approx_methods <- c("EP", "Block-VB", "MFVB-CAVI", "MMFVB", "LA")

cat(sprintf("=== Running Post-Hoc Posterior Inference Benchmark (mode = %s) ===\n", run_mode))

results_list <- list()

for (rho in rho_seq) {
  setting_name <- if (abs(rho) < 1e-4) "uncorrelated" else "correlated_ar1"
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("Setting: %s (rho = %.2f)\n", setting_name, rho))
  cat(sprintf("======================================================\n"))
  
  for (model_spec in models) {
    cat(sprintf("\n--- Prior: %s ---\n", model_spec$name))
    
    for (p in p_seq) {
      cat(sprintf("  Dimension: p = %d\n", p))
      
      # Generate base data for max sample size
      n_max <- max(n_seq)
      base_data <- generate_posthoc_data(
        n = n_max,
        p = p,
        rho = rho,
        seed = 4747 + p + round(rho * 100)
      )
      
      for (n in n_seq) {
        if (n <= p) next
        
        y_sub <- base_data$y[1:n]
        X_sub <- base_data$X[1:n, , drop = FALSE]
        XtX_sub <- crossprod(X_sub)
        Xy_sub <- as.vector(crossprod(X_sub, y_sub))
        
        # OLS (with intercept) estimate for sigma2.hat initialization
        sigma2_hat <- sum(stats::lm.fit(cbind(1, X_sub), y_sub)$residuals^2) / (n - p - 1)
        
        hyper_params <- list(tau = 1, nu = model_spec$nu)
        
        for (method in approx_methods) {
          t0 <- Sys.time()
          fit <- tryCatch({
            fit_approx_posterior(
              y = y_sub,
              X = X_sub,
              method = method,
              prior_type = model_spec$prior_type,
              hyper_params = hyper_params,
              intercept = TRUE,
              sigma2.hat = sigma2_hat,
              XtX = XtX_sub,
              X.y = Xy_sub
            )
          }, error = function(e) {
            warning(sprintf("Method %s failed: %s", method, e$message))
            NULL
          })
          runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
          
          if (!is.null(fit)) {
            eval_res <- eval_posthoc_fit(fit, base_data$beta_true, base_data$alpha_true, base_data$sigma2_true)
            
            results_list[[length(results_list) + 1]] <- data.frame(
              setting = setting_name,
              rho = rho,
              prior = model_spec$name,
              method = method,
              n = n,
              p = p,
              runtime_sec = runtime,
              rmse_beta = eval_res$rmse_beta,
              mae_beta = eval_res$mae_beta,
              coverage_95 = eval_res$coverage_95,
              ci_length = eval_res$ci_length,
              err_sigma2 = eval_res$err_sigma2
            )
            
            cat(sprintf("    %-11s (n=%5d): Time %5.3fs | RMSE: %.4f | 95%% Cov: %.3f\n",
                        method, n, runtime, eval_res$rmse_beta, eval_res$coverage_95))
          }
        }
      }
    }
  }
}

# ------------------------------------------------------------------------------
# Save Summary
# ------------------------------------------------------------------------------
results_df <- do.call(rbind, results_list)
out_path <- "simulations/results/posthoc_results.rds"
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(results_df, out_path)

cat(sprintf("\n=== Post-Hoc Inference Benchmark Complete! Saved to '%s' ===\n", out_path))
print(head(results_df))
