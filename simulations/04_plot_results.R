################################################################################
### Simulation 04: Visualize and Summarize Benchmark Results
### Generates comparative plots and tables from simulation output files
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

figures_dir <- "simulations/figures"
results_dir <- "simulations/results"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Custom clean plotting theme
theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#f0f0f0", color = NA),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14)
    )
}

# Split method labels such as "FASS_Student_df5" or "Gibbs_PC_Laplace" into the
# algorithm and the prior family
prior_pattern <- "_(Laplace|Gaussian|Student_df[0-9]+)$"
split_method <- function(df) {
  has_prior <- grepl(prior_pattern, df$method)
  df$algorithm <- ifelse(has_prior, sub(prior_pattern, "", df$method), df$method)
  df$prior <- ifelse(has_prior, sub(paste0("^.*", prior_pattern), "\\1", df$method), "n/a")
  df
}

to_long_rates <- function(df) {
  long <- reshape(
    df[, c("algorithm", "prior", "method", "n", "TPR", "TNR")],
    direction = "long",
    varying = c("TPR", "TNR"),
    v.names = "Rate",
    timevar = "Metric",
    times = c("True Positive Rate (TPR)", "True Negative Rate (TNR)")
  )
  long$Metric <- factor(long$Metric, levels = c("True Positive Rate (TPR)", "True Negative Rate (TNR)"))
  long
}

# ------------------------------------------------------------------------------
# 1. Superinduction Plots
# ------------------------------------------------------------------------------
super_path <- file.path(results_dir, "superinduction_results.rds")
if (file.exists(super_path)) {
  cat("Plotting Superinduction Results...\n")
  df_super <- split_method(readRDS(super_path))
  
  # A. Runtime Curve (Figure 2, right panel)
  p_time <- ggplot(df_super, aes(x = n, y = runtime_sec, color = prior, linetype = algorithm,
                                 group = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    labs(
      title = "Runtime Scaling under Superinduction",
      subtitle = sprintf("Feature dimension p = %d, active signals p_true = %d", df_super$p[1], df_super$p_true[1]),
      x = "Sample Size (n, log scale)",
      y = "Runtime in Seconds (log scale)",
      color = "Prior",
      linetype = "Algorithm"
    ) +
    theme_pub()
  
  ggsave(file.path(figures_dir, "superinduction_runtime.png"), p_time, width = 8, height = 5, dpi = 300)
  
  # B. TPR & TNR Curves by algorithm (Figures 2-3)
  df_acc <- to_long_rates(df_super)
  
  p_acc <- ggplot(df_acc, aes(x = n, y = Rate, color = prior, group = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_grid(Metric ~ algorithm) +
    scale_x_log10(labels = scales::comma) +
    ylim(0, 1.05) +
    labs(
      title = "Variable Selection Accuracy under Superinduction",
      x = "Sample Size (n, log scale)",
      y = "Rate",
      color = "Prior"
    ) +
    theme_pub()
  
  ggsave(file.path(figures_dir, "superinduction_accuracy.png"), p_acc,
         width = 3 + 3 * length(unique(df_acc$algorithm)), height = 6, dpi = 300)
}

# ------------------------------------------------------------------------------
# 2. M-Closed Plots
# ------------------------------------------------------------------------------
mclosed_path <- file.path(results_dir, "m_closed_results.rds")
if (file.exists(mclosed_path)) {
  cat("Plotting M-Closed Results...\n")
  df_mc <- split_method(readRDS(mclosed_path))
  df_mc$config <- sprintf("p = %d, p_true = %d", df_mc$p, df_mc$p_true)
  
  # A. Recovery rates (Figure 1)
  df_mc_acc <- to_long_rates(df_mc)
  df_mc_acc$config <- sprintf("p = %d, p_true = %d", df_mc$p, df_mc$p_true)[df_mc_acc$id]
  p_mc_acc <- ggplot(df_mc_acc, aes(x = n, y = Rate, color = prior, group = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_grid(Metric ~ config) +
    scale_x_log10(labels = scales::comma) +
    ylim(0, 1.05) +
    labs(
      title = "FASS Recovery Rates in the M-Closed Setting",
      x = "Sample Size (n, log scale)",
      y = "Rate",
      color = "Prior"
    ) +
    theme_pub()
  
  ggsave(file.path(figures_dir, "m_closed_accuracy.png"), p_mc_acc, width = 9, height = 6, dpi = 300)
  
  # B. Runtime
  p_mc_time <- ggplot(df_mc, aes(x = n, y = runtime_sec, color = prior, group = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ config, scales = "free_y") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    labs(
      title = "FASS Runtime Scaling in M-Closed Setting",
      x = "Sample Size (n, log scale)",
      y = "Runtime in Seconds (log scale)",
      color = "Prior Family"
    ) +
    theme_pub()
  
  ggsave(file.path(figures_dir, "m_closed_runtime.png"), p_mc_time, width = 9, height = 6, dpi = 300)
}

# ------------------------------------------------------------------------------
# 3. Post-Hoc Estimation Plots
# ------------------------------------------------------------------------------
posthoc_path <- file.path(results_dir, "posthoc_results.rds")
if (file.exists(posthoc_path)) {
  cat("Plotting Post-Hoc Inference Results...\n")
  df_ph <- readRDS(posthoc_path)
  
  p_rmse <- ggplot(df_ph, aes(x = n, y = rmse_beta, color = method, linetype = setting, group = interaction(method, setting))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_grid(prior ~ p, labeller = label_both) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    labs(
      title = "Post-Hoc Estimation Error (RMSE)",
      x = "Sample Size (n, log scale)",
      y = "Coefficient RMSE (log scale)",
      color = "Method",
      linetype = "Design Setting"
    ) +
    theme_pub()
  
  ggsave(file.path(figures_dir, "posthoc_rmse.png"), p_rmse, width = 10, height = 7, dpi = 300)
}

cat(sprintf("All figures successfully exported to '%s/'!\n", figures_dir))
