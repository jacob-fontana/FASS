################################################################################
### Simulation 02: M-Closed Variable Selection Benchmark
### Evaluates FASS when the true data-generating model is strictly inside the model space (epsilon = 0)
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
# Can be overridden from the command line: Rscript simulations/02_m_closed.R full
args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args) >= 1) args[1] else "test"
stopifnot(run_mode %in% c("test", "full"))

# Each configuration is a (p, p_true) pair
if (run_mode == "test") {
  configs <- list(c(p = 50, p_true = 5), c(p = 100, p_true = 10))
  n_seq <- c(500, 2000)
  niter_fass <- 1000
} else {
  # Paper Section 5: p = 200 with 100 true signals, and p = 500 with 10 true signals
  configs <- list(c(p = 200, p_true = 100), c(p = 500, p_true = 10))
  n_seq <- c(5000, 10000, 25000, 50000, 75000, 100000, 250000, 500000, 10^6)
  niter_fass <- 10000
}

# Only recovery is evaluated here, so no post-hoc inference is needed
post_hoc_method <- NA
epsilon <- 0.00          # M-Closed setting (zero latent leakage)
seed_base <- 4747

cat(sprintf("=== Running M-Closed Simulation Benchmark (mode = %s) ===\n", run_mode))

results_list <- list()

# ------------------------------------------------------------------------------
# Grid over (p, p_true) configurations
# ------------------------------------------------------------------------------
for (cfg in configs) {
  p <- unname(cfg["p"])
  p_true <- unname(cfg["p_true"])
  
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("Configuration: p = %d, p_true = %d\n", p, p_true))
  cat(sprintf("======================================================\n"))
  
  # Generate maximum sample size dataset
  n_max <- max(n_seq)
  cat(sprintf("Generating base dataset with n_max = %d ... ", n_max))
  full_data <- generate_selection_data(
    n = n_max,
    p = p,
    p_true = p_true,
    epsilon = epsilon,
    seed = seed_base + p + p_true
  )
  cat("Done.\n")
  
  for (n in n_seq) {
    if (n <= p) next
    cat(sprintf("\n--- Sample Size: n = %d ---\n", n))
    y_sub <- full_data$y[1:n]
    X_sub <- full_data$X[1:n, , drop = FALSE]
    
    methods_to_run <- list(
      list(name = "FASS_Laplace", family = "Laplace", tau1 = 1, tau0 = sqrt(n), nu = NA),
      list(name = "FASS_Gaussian", family = "Gaussian", tau1 = 1, tau0 = n, nu = NA),
      list(name = "FASS_Student_df2", family = "Student", tau1 = 1, tau0 = n, nu = 2),
      list(name = "FASS_Student_df5", family = "Student", tau1 = 1, tau0 = n, nu = 5),
      list(name = "FASS_Student_df10", family = "Student", tau1 = 1, tau0 = n, nu = 10)
    )
    
    for (m_spec in methods_to_run) {
      cat(sprintf("  Fitting %-18s ... ", m_spec$name))
      t0 <- Sys.time()
      fit <- FASS.sampler(
        y = y_sub,
        X = X_sub,
        niter = niter_fass,
        family = m_spec$family,
        tau1 = m_spec$tau1,
        tau0 = m_spec$tau0,
        nu = m_spec$nu,
        full_posterior = post_hoc_method,
        verbose = FALSE,
        intercept = TRUE,
        rescale_y = FALSE
      )
      runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      metrics <- eval_variable_selection(fit, full_data$true_idx, p, is_gibbs = FALSE)
      
      cat(sprintf("Time: %6.2fs | TPR: %.3f | TNR: %.3f | Active: %d\n",
                  runtime, metrics$TPR, metrics$TNR, metrics$active_count))
      
      results_list[[length(results_list) + 1]] <- data.frame(
        experiment = "m_closed",
        method = m_spec$name,
        n = n,
        p = p,
        p_true = p_true,
        runtime_sec = runtime,
        TPR = metrics$TPR,
        TNR = metrics$TNR,
        active_count = metrics$active_count,
        false_positives = metrics$false_positives
      )
    }
  }
}

# ------------------------------------------------------------------------------
# Save Summary
# ------------------------------------------------------------------------------
results_df <- do.call(rbind, results_list)
out_path <- "simulations/results/m_closed_results.rds"
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(results_df, out_path)

cat(sprintf("\n=== M-Closed Benchmark Complete! Saved to '%s' ===\n", out_path))
print(head(results_df[, c("method", "p", "p_true", "n", "runtime_sec", "TPR", "TNR")]))
