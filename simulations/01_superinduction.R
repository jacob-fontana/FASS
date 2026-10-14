################################################################################
### Simulation 01: Superinduction Variable Selection Benchmark
### Compares FASS across prior families against dense latent confounding (epsilon = 0.01)
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
# Can be overridden from the command line: Rscript simulations/01_superinduction.R full
args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args) >= 1) args[1] else "test"
stopifnot(run_mode %in% c("test", "full"))

if (run_mode == "test") {
  n_seq <- c(1000, 5000)
  p <- 100
  p_true <- 10
  niter_fass <- 1000
  include_mcmc_baselines <- TRUE
  gibbs_burn <- 250
  gibbs_iters <- 500
} else {
  n_seq <- c(5000, 10000, 25000, 50000, 75000, 100000, 250000, 500000, 10^6)
  p <- 500
  p_true <- 10
  niter_fass <- 10000
  # The Gibbs baselines (Figures 2-3 of the paper) are very slow at large n:
  # set to FALSE to run the FASS part of the grid only.
  include_mcmc_baselines <- TRUE
  gibbs_burn <- 2000
  gibbs_iters <- 8000
}
mcmc_samplers <- c("Gibbs", "Gibbs_PC", "skinnyGibbs")

post_hoc_method <- "LA"  # "LA", "Block-VB", "EP", or NA (selection only)
epsilon <- 0.01          # Superinduction contamination weight
seed <- 4747

cat(sprintf("=== Running Superinduction Benchmark (mode = %s, p = %d, p_true = %d) ===\n", run_mode, p, p_true))

# ------------------------------------------------------------------------------
# 1. Generate Maximum Dataset & Subset
# ------------------------------------------------------------------------------
n_max <- max(n_seq)
cat(sprintf("Generating base dataset with n_max = %d ... ", n_max))
full_data <- generate_selection_data(n = n_max, p = p, p_true = p_true, epsilon = epsilon, seed = seed)
cat("Done.\n\n")

results_list <- list()

# ------------------------------------------------------------------------------
# 2. Iterate Over Sample Sizes
# ------------------------------------------------------------------------------
for (n in n_seq) {
  cat(sprintf("\n--- Sample Size: n = %d ---\n", n))
  y_sub <- full_data$y[1:n]
  X_sub <- full_data$X[1:n, , drop = FALSE]
  
  # Method definitions
  methods_to_run <- list(
    list(name = "FASS_Laplace", family = "Laplace", tau1 = 1, tau0 = sqrt(n), nu = NA),
    list(name = "FASS_Gaussian", family = "Gaussian", tau1 = 1, tau0 = n, nu = NA),
    list(name = "FASS_Student_df2", family = "Student", tau1 = 1, tau0 = n, nu = 2),
    list(name = "FASS_Student_df5", family = "Student", tau1 = 1, tau0 = n, nu = 5),
    list(name = "FASS_Student_df10", family = "Student", tau1 = 1, tau0 = n, nu = 10)
  )
  
  # A. Run FASS Models
  for (m_spec in methods_to_run) {
    cat(sprintf("  Fitting %-24s ... ", m_spec$name))
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
      experiment = "superinduction",
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
  
  # B. Run MCMC Benchmarks (if enabled)
  # Same priors as FASS on the variance scale: spike tau0_2 = 1/n, slab tau1_2 = 1.
  # Selection uses the median probability model.
  if (include_mcmc_baselines) {
    for (sampler_name in mcmc_samplers) {
      sampler_fun <- match.fun(sampler_name)
      for (m_spec in methods_to_run) {
        label <- sub("^FASS", sampler_name, m_spec$name)
        cat(sprintf("  Fitting %-24s ... ", label))
        t0 <- Sys.time()
        fit_mcmc <- sampler_fun(
          X = X_sub,
          y = y_sub,
          tau0_2 = 1 / n,
          tau1_2 = 1,
          a_pi = 1,
          b_pi = 1,
          nburn = gibbs_burn,
          niter = gibbs_iters,
          family = m_spec$family,
          nu = m_spec$nu,
          intercept = TRUE,
          rescale_y = FALSE
        )
        runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        metrics <- eval_variable_selection(fit_mcmc, full_data$true_idx, p, is_gibbs = TRUE)
        cat(sprintf("Time: %6.2fs | TPR: %.3f | TNR: %.3f | Active: %d\n",
                    runtime, metrics$TPR, metrics$TNR, metrics$active_count))
        
        results_list[[length(results_list) + 1]] <- data.frame(
          experiment = "superinduction",
          method = label,
          n = n,
          p = p,
          p_true = p_true,
          runtime_sec = runtime,
          TPR = metrics$TPR,
          TNR = metrics$TNR,
          active_count = metrics$active_count,
          false_positives = metrics$false_positives
        )
        rm(fit_mcmc)
      }
    }
  }
}

# ------------------------------------------------------------------------------
# 3. Save Summary Table
# ------------------------------------------------------------------------------
results_df <- do.call(rbind, results_list)
out_path <- "simulations/results/superinduction_results.rds"
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(results_df, out_path)

cat(sprintf("\n=== Superinduction Benchmark Complete! Saved to '%s' ===\n", out_path))
print(results_df[, c("method", "n", "runtime_sec", "TPR", "TNR", "active_count")])
