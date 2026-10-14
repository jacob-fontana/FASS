library(ggplot2)
library(tidyr)
library(dplyr)

suffix_map <- c(
  "Gibbs" = "Gibbs",
  "PC" = "Gibbs (PC)",
  "Skinny" = "Skinny Gibbs",
  "EP" = "FASS (EP)",
  "VB" = "FASS (VB)",
  "LA" = "FASS (LA)"
)

dist_labels <- c("Laplace", "Gaussian", "Student (df=2)", "Student (df=5)", "Student (df=10)")

load_metric <- function(metric_name, metric_file_prefix) {
  df_list <- lapply(names(suffix_map), function(s) {
    
    file_path <- file.path("Superinduction_Results", paste0("_results_", metric_file_prefix, "_", s))
    
    if (file.exists(file_path)) {
        mat <- readRDS(file_path)
    } else {
        cat("Missing:", file_path, "\n")
        return(NULL)
    }
    
    # Check dimensions or assume standard sizes
    # In Make_plots2.R, Superinduction2, 3, 4 (FASS EP, VB, LA) ran on 14 values of n.seq
    # Others ran on 8.
    if (s %in% c("EP", "VB", "LA")) {
      n_seq <- c(5000, 10000, 25000, 50000, 75000, 100000, 150000, 200000, 250000, 300000, 400000, 500000, 750000, 10^6)
    } else {
      # In some runs, they might have extended them to 14 as well! Let's check rows of mat
      n_seq_full <- c(5000, 10000, 25000, 50000, 75000, 100000, 150000, 200000, 250000, 300000, 400000, 500000, 750000, 10^6)
      if(ncol(mat) == 14) {
          n_seq <- n_seq_full
      } else {
          n_seq <- c(5000, 10000, 25000, 50000, 100000, 250000, 500000, 10^6)
      }
    }
    
    df <- as.data.frame(t(mat))
    colnames(df) <- dist_labels
    df$n <- n_seq
    df$Algorithm <- suffix_map[s]
    
    df_long <- pivot_longer(df, cols = all_of(dist_labels), names_to = "Distribution", values_to = "Value")
    df_long$Metric <- metric_name
    
    return(df_long)
  })
  return(bind_rows(df_list))
}

df_times <- load_metric("Runtime (secs)", "times")
df_TPR <- load_metric("TPR", "TPR")
df_TNR <- load_metric("TNR", "TNR")

df_all <- rbind(df_times, df_TPR, df_TNR)

df_all$Algorithm <- factor(df_all$Algorithm, levels = c("Gibbs", "Gibbs (PC)", "Skinny Gibbs", "FASS (EP)", "FASS (VB)", "FASS (LA)"))
df_all$Distribution <- factor(df_all$Distribution, levels = dist_labels)

df_all$Group <- paste(df_all$Algorithm, df_all$Distribution)

plot_metric <- function(df, metric_name, plot_title, y_label, log_y = FALSE) {
  df_sub <- df %>% filter(Metric == metric_name)
  
  p_out <- ggplot(df_sub, aes(x = n, y = Value, color = Algorithm, linetype = Distribution, shape = Algorithm, group = Group)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_x_log10(labels = scales::comma) + 
    theme_bw() +
    labs(title = plot_title,
         subtitle = "True Covariates: 10 | Total Covariates: 500",
         x = "Sample Size (n) [log scale]",
         y = y_label) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          legend.key.width = unit(1.5, "cm"))
  
  if(log_y) {
    p_out <- p_out + scale_y_log10(labels = scales::comma)
  }
  
  return(p_out)
}

p_times <- plot_metric(df_all, "Runtime (secs)", 
                       plot_title = "Algorithm Runtime vs. Sample Size", 
                       y_label = "Runtime (seconds)", log_y = TRUE)

p_TPR <- plot_metric(df_all, "TPR", 
                     plot_title = "True Positive Rate vs. Sample Size", 
                     y_label = "True Positive Rate", log_y = FALSE)

p_TNR <- plot_metric(df_all, "TNR", 
                     plot_title = "True Negative Rate vs. Sample Size", 
                     y_label = "True Negative Rate", log_y = FALSE)

ggsave("Compare_Runtime2.png", p_times, width = 12, height = 8, dpi = 300)
ggsave("Compare_TPR2.png", p_TPR, width = 12, height = 8, dpi = 300)
ggsave("Compare_TNR2.png", p_TNR, width = 12, height = 8, dpi = 300)

cat("Plots generated and saved successfully.\n")
