library(xtable)
library(dplyr)

# Function to generate standalone TeX
generate_tex <- function(xtable_obj, filename, p_val) {
  sink(filename)
  cat("\\documentclass[10pt]{article}\n")
  if (p_val == 10) {
    cat("\\usepackage[a3paper, landscape, margin=0.5in]{geometry}\n") 
  } else {
    cat("\\usepackage[a4paper, landscape, margin=1in]{geometry}\n") 
  }
  cat("\\usepackage{booktabs}\n")
  cat("\\usepackage{longtable}\n")
  cat("\\usepackage{multirow}\n")
  cat("\\begin{document}\n")
  cat("\\setlength{\\tabcolsep}{4pt}\n") 
  
  print(xtable_obj, 
        tabular.environment="longtable", 
        floating=FALSE, 
        include.rownames=FALSE, 
        size="\\small",
        sanitize.text.function = identity) 
  
  cat("\\end{document}\n")
  sink()
}

# Find all generated CSVs
files <- list.files(path="csv", pattern="^(Correlated|Uncorrelated)_Sim_.*\\.csv$", full.names=TRUE, recursive=TRUE)

for (f in files) {
  df <- read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
  model_name <- gsub("(Correlated|Uncorrelated)_Sim_(.*)_p.*\\.csv", "\\2", f)
  sim_type <- if (grepl("Uncorrelated", f)) "Uncorrelated" else "Correlated"
  p_val <- if (grepl("p1000", f)) 1000 else if (grepl("p100", f)) 100 else 10
  
  if (p_val == 10) {
    df$Method[df$Method == "VB"] <- "Block-VB"
    df$Method[df$Method == "VB_Direct"] <- "MMFVB"
    df$Method[df$Method == "VB_MF"] <- "MFVB-CAVI"
    # It is grouped by Method, then Parameter. We want to group by Parameter, then Method.
    # Split into Runtime and Parameters
    df_params <- df[df$Parameter != "Runtime", ]
    df_runtime <- df[df$Parameter == "Runtime", ]
    
    # Sort df_params by parameter number
    # Parameter looks like "beta_1", "beta_2"
    df_params$ParamNum <- as.numeric(gsub("beta_", "", df_params$Parameter))
    # Define method order
    methods_order <- c("EP", "Block-VB", "MFVB-CAVI", "MMFVB", "LA")
    df_params$MethodNum <- match(df_params$Method, methods_order)
    df_runtime$MethodNum <- match(df_runtime$Method, methods_order)
    
    df_params <- df_params[order(df_params$ParamNum, df_params$MethodNum), ]
    df_runtime <- df_runtime[order(df_runtime$MethodNum), ]
    
    # Clean up columns used for sorting
    df_params$ParamNum <- NULL
    df_params$MethodNum <- NULL
    df_runtime$MethodNum <- NULL
    
    # Apply multirow to parameters
    n_methods <- length(methods_order)
    for (i in 1:nrow(df_params)) {
      is_first <- (i %% n_methods) == 1
      if (is_first) {
        param_num <- gsub("beta_", "", df_params$Parameter[i])
        df_params$Parameter[i] <- paste0("\\multirow{5}{*}{$\\beta_{", param_num, "}$}")
        true_val <- df_params$True_Value[i]
        df_params$True_Value[i] <- paste0("\\multirow{5}{*}{", true_val, "}")
      } else {
        df_params$Parameter[i] <- ""
        df_params$True_Value[i] <- ""
      }
    }
    
    # Combine back (put runtimes at the bottom, or just append)
    df_final <- rbind(df_params, df_runtime)
    names(df_final)[names(df_final) == "True_Value"] <- "True Value"
    
    # Escape underscores in Method column
    df_final$Method <- gsub("_", "\\\\_", df_final$Method)
    
    # Make xtable
    safe_model_name <- gsub("_", "\\\\_", toupper(model_name))
    cap <- sprintf("%s Simulation (p=10) - %s Prior", sim_type, safe_model_name)
    ltx <- xtable(df_final, caption=cap, label=paste0("tab:", tolower(sim_type), "_sim_", model_name, "_p10"))
    
  } else {
    df$Method[df$Method == "VB"] <- "Block-VB"
    df$Method[df$Method == "VB_Direct"] <- "MMFVB"
    df$Method[df$Method == "VB_MF"] <- "MFVB-CAVI"
    # p = 100
    # Grouped by Method, Metric
    # Apply multirow to Method
    methods_order <- c("EP", "Block-VB", "MFVB-CAVI", "MMFVB", "LA")
    df$MethodNum <- match(df$Method, methods_order)
    
    # Define metric order
    metrics_order <- c("Runtime", "Coverage", "RMSE", "Length (Min, 25%, 75%, Max)")
    df$MetricNum <- match(df$Metric, metrics_order)
    
    df <- df[order(df$MethodNum, df$MetricNum), ]
    df$MethodNum <- NULL
    df$MetricNum <- NULL
    
    n_metrics <- length(metrics_order)
    for (i in 1:nrow(df)) {
      is_first <- (i %% n_metrics) == 1
      if (is_first) {
        method_name <- df$Method[i]
        df$Method[i] <- paste0("\\multirow{4}{*}{", gsub("_", "\\\\_", method_name), "}")
      } else {
        df$Method[i] <- ""
      }
    }
    
    # Escape percentages in Metric
    df$Metric <- gsub("%", "\\\\%", df$Metric)
    
    safe_model_name <- gsub("_", "\\\\_", toupper(model_name))
    cap <- sprintf("%s Simulation (p=%d Summary) - %s Prior", sim_type, p_val, safe_model_name)
    ltx <- xtable(df, caption=cap, label=paste0("tab:", tolower(sim_type), "_sim_", model_name, "_p", p_val))
  }
  
  tex_file <- sub("^csv/", "tables/", gsub("\\.csv", ".tex", f))
  generate_tex(ltx, tex_file, p_val)
}

cat("LaTeX tables generated successfully.\n")
