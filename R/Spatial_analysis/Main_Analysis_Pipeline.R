
# ===========================================================
# Author: [Samantha Reynoso]
# Date: [09.09.2025]
# Description: 
# Comprehensive statistical analysis of DLPFC samples using STEAM
# with nested cross-validation and robust statistical testing
# 
# Citation: Include relevant STEAM and method citations
# Data: DLPFC12.RDS contains 12 samples for analysis
# ===========================================================

# Environment Setup ----
# Set reproducible seed
set.seed(0)

# Package Management with Automatic Installation
required_packages <- c(
  # Core analysis
  "ggplot2", "scales", "dplyr",
  # Statistical analysis
  "effectsize", "multcomp", "emmeans", "broom", "PMCMRplus",
  # Optional advanced packages
  "lme4", "car", "coin", "MBESS", "rstanarm",
  # Visualization
  "patchwork"
)

# Install missing packages
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    install.packages(new_packages, dependencies = TRUE)
  }
}

install_if_missing(required_packages)

# Load required packages (suppress startup messages for cleaner output)
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(effectsize)
  library(multcomp)
  library(emmeans)
  library(broom)
  library(PMCMRplus)
  library(patchwork)
})

# Data Loading ----
if (!file.exists("DLPFC12.RDS")) {
  stop("DLPFC12.RDS file not found. Please check the file path.")
}

DLPFCs <- readRDS("DLPFC12.RDS")
message("Loaded DLPFC data with ", length(DLPFCs), " samples")

# Configuration ----
SAMPLE_IDS <- c("151507","151508","151509","151510",
                "151669","151670","151671","151672", 
                "151673","151674","151675","151676")

# STEAM parameters
STEAM_CONFIG <- list(
  n_outer_folds = 5,
  n_inner_folds = 3,
  cv_cores = 14,
  allow_parallel = FALSE,
  n_size = 5,
  model = "rf",
  metric = "Kappa",
  max_weights = 15000,
  label_column = "Labels",
  assay = "SCT"
)

# Analysis parameters
ANALYSIS_CONFIG <- list(
  confidence_level = 0.95,
  correction_method = "bonferroni",
  metrics_to_analyze = c("Accuracy", "Kappa", "ARI", "NMI", "AMI"),
  output_prefix = "DLPFC_STEAM_analysis"
)

# Helper Functions ----

# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Compute per-fold metrics from STEAM outer results
compute_outer_fold_metrics <- function(outer_res) {
  stopifnot(!is.null(outer_res), length(outer_res) > 0)
  
  fold_metrics <- lapply(seq_along(outer_res), function(i) {
    fold_result <- outer_res[[i]]
    predictions <- fold_result$preds
    
    # Clean and align predictions
    predicted <- as.character(predictions$predy)
    true_labels <- as.character(predictions$testy)
    
    # Remove missing values
    valid_indices <- !(is.na(predicted) | is.na(true_labels))
    predicted <- predicted[valid_indices]
    true_labels <- true_labels[valid_indices]
    
    # Ensure common factor levels
    all_levels <- sort(unique(c(predicted, true_labels)))
    pred_factor <- factor(predicted, levels = all_levels)
    true_factor <- factor(true_labels, levels = all_levels)
    
    # Calculate metrics using caret
    confusion_matrix <- caret::confusionMatrix(pred_factor, true_factor)
    
    accuracy <- unname(confusion_matrix$overall["Accuracy"])
    kappa <- unname(confusion_matrix$overall["Kappa"])
    
    # Calculate F1 Score (macro average)
    class_metrics <- confusion_matrix$byClass
    if (is.matrix(class_metrics)) {
      f1_scores <- class_metrics[, "F1"]
      macro_f1 <- mean(f1_scores, na.rm = TRUE)
    } else if (is.vector(class_metrics) && "F1" %in% names(class_metrics)) {
      macro_f1 <- class_metrics["F1"]
    } else {
      # Fallback calculation
      precision <- class_metrics[, "Precision"]
      recall <- class_metrics[, "Recall"]
      f1_scores <- 2 * (precision * recall) / (precision + recall)
      macro_f1 <- mean(f1_scores, na.rm = TRUE)
    }
    
    # External clustering indices
    ari <- mclust::adjustedRandIndex(true_labels, predicted)
    nmi <- aricode::NMI(true_labels, predicted)
    ami <- aricode::AMI(true_labels, predicted)
    
    data.frame(
      Outer_Fold = i,
      Accuracy = as.numeric(accuracy),
      Kappa = as.numeric(kappa),
      Macro_F1 = as.numeric(macro_f1),
      ARI = as.numeric(ari),
      NMI = as.numeric(nmi),
      AMI = as.numeric(ami),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, fold_metrics)
}

# Enhanced statistical summary with confidence intervals and effect sizes
summarize_fold_metrics <- function(fold_df, confidence_level = 0.95) {
  stopifnot(is.data.frame(fold_df), "Outer_Fold" %in% names(fold_df))
  
  metric_columns <- setdiff(names(fold_df), c("Outer_Fold", "Sample"))
  
  summary_stats <- lapply(metric_columns, function(metric) {
    values <- fold_df[[metric]]
    values <- values[!is.na(values)]
    n <- length(values)
    
    if (n == 0) {
      return(create_empty_summary_row(metric, n))
    }
    
    if (!is.numeric(values)) {
      warning("Non-numeric data in column: ", metric)
      return(create_empty_summary_row(metric, n))
    }
    
    # Handle constant values
    if (sd(values) == 0 || var(values) < .Machine$double.eps) {
      warning("Constant values in column: ", metric)
      return(data.frame(
        Metric = metric,
        Mean = mean(values),
        SD = 0,
        SE = 0,
        CI_Lower = mean(values),
        CI_Upper = mean(values),
        Cohen_d = NA,
        Cohen_d_CI_Lower = NA,
        Cohen_d_CI_Upper = NA,
        n_folds = n,
        stringsAsFactors = FALSE
      ))
    }
    
    # Compute statistics
    t_test_result <- tryCatch({
      t.test(values, conf.level = confidence_level)
    }, error = function(e) {
      warning("t.test failed for ", metric, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(t_test_result)) {
      return(create_empty_summary_row(metric, n))
    }
    
    tidy_result <- broom::tidy(t_test_result)
    
    # Effect size calculation
    cohens_d_result <- tryCatch({
      if (sd(values) > 0 && is.finite(mean(values))) {
        effectsize::cohens_d(values, mu = 0)
      } else {
        list(Cohens_d = NA, CI_low = NA, CI_high = NA)
      }
    }, error = function(e) {
      # Fallback calculation
      if (sd(values) > 0) {
        list(Cohens_d = mean(values) / sd(values), CI_low = NA, CI_high = NA)
      } else {
        list(Cohens_d = NA, CI_low = NA, CI_high = NA)
      }
    })
    
    # Extract Cohen's d values
    cohen_d_val <- if (is.list(cohens_d_result)) cohens_d_result$Cohens_d else cohens_d_result[1]
    cohen_d_ci_low <- if (is.list(cohens_d_result)) cohens_d_result$CI_low else NA
    cohen_d_ci_high <- if (is.list(cohens_d_result)) cohens_d_result$CI_high else NA
    
    # Standard error calculation
    se_val <- if (is.finite(tidy_result$statistic) && tidy_result$statistic != 0) {
      abs(tidy_result$estimate / tidy_result$statistic)
    } else {
      sd(values) / sqrt(n)
    }
    
    data.frame(
      Metric = metric,
      Mean = tidy_result$estimate,
      SD = sd(values),
      SE = se_val,
      CI_Lower = tidy_result$conf.low,
      CI_Upper = tidy_result$conf.high,
      Cohen_d = cohen_d_val,
      Cohen_d_CI_Lower = cohen_d_ci_low,
      Cohen_d_CI_Upper = cohen_d_ci_high,
      n_folds = n,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, summary_stats)
}

# Helper function for empty summary rows
create_empty_summary_row <- function(metric, n) {
  data.frame(
    Metric = metric, Mean = NA, SD = NA, SE = NA,
    CI_Lower = NA, CI_Upper = NA, Cohen_d = NA,
    Cohen_d_CI_Lower = NA, Cohen_d_CI_Upper = NA,
    n_folds = n, stringsAsFactors = FALSE
  )
}

# Extract pooled metrics from STEAM results
extract_pooled_metrics <- function(steam_result) {
  # Initialize result with expected structure
  result <- data.frame(
    Accuracy = NA_real_, Kappa = NA_real_, ARI = NA_real_,
    NMI = NA_real_, AMI = NA_real_, Macro_F1 = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # Try different paths to find metrics
  metrics_source <- NULL
  if (!is.null(steam_result[["nested"]][["metrics"]])) {
    metrics_source <- steam_result[["nested"]][["metrics"]]
  } else if (!is.null(steam_result[["train"]][["ncv"]][["summary"]])) {
    metrics_source <- steam_result[["train"]][["ncv"]][["summary"]]
  }
  
  if (is.null(metrics_source)) {
    warning("No pooled metrics found in STEAM results")
    return(result)
  }
  
  # Extract metrics with flexible naming
  if (is.list(metrics_source)) {
    result$Accuracy <- as.numeric(metrics_source$Accuracy %||% metrics_source$accuracy %||% NA_real_)
    result$Kappa <- as.numeric(metrics_source$Kappa %||% metrics_source$kappa %||% NA_real_)
    result$ARI <- as.numeric(metrics_source$ARI %||% metrics_source$ari %||% NA_real_)
    result$NMI <- as.numeric(metrics_source$NMI %||% metrics_source$nmi %||% NA_real_)
    result$AMI <- as.numeric(metrics_source$AMI %||% metrics_source$ami %||% NA_real_)
    result$Macro_F1 <- as.numeric(metrics_source$Macro_F1 %||% metrics_source$f1 %||% NA_real_)
  }
  
  return(result)
}

# Pairwise sample comparisons using established statistical methods
perform_pairwise_comparisons <- function(all_fold_data, metric = "Accuracy", 
                                        correction_method = "bonferroni") {
  
  # Prepare long-format data
  long_data <- data.frame()
  for (sample_id in names(all_fold_data)) {
    if (metric %in% names(all_fold_data[[sample_id]])) {
      sample_df <- data.frame(
        Sample = sample_id,
        Metric_Value = all_fold_data[[sample_id]][[metric]],
        Fold = all_fold_data[[sample_id]]$Outer_Fold,
        stringsAsFactors = FALSE
      )
      long_data <- rbind(long_data, sample_df)
    }
  }
  
  if (nrow(long_data) == 0 || length(unique(long_data$Sample)) < 2) {
    warning("Insufficient data for pairwise comparisons")
    return(NULL)
  }
  
  # Primary method: emmeans for pairwise comparisons
  tryCatch({
    model <- lm(Metric_Value ~ Sample, data = long_data)
    emm <- emmeans::emmeans(model, ~ Sample)
    pairwise_result <- pairs(emm, adjust = correction_method)
    pairwise_tidy <- broom::tidy(pairwise_result)
    
    # Add effect sizes
    effect_sizes <- calculate_effect_sizes(long_data)
    final_results <- merge(pairwise_tidy, effect_sizes, by = "contrast", all.x = TRUE)
    
    # Add interpretation
    final_results$Interpretation <- effectsize::interpret_cohens_d(abs(final_results$Effect_Size))
    final_results$Metric <- metric
    final_results$Correction_Method <- correction_method
    final_results$Significant <- final_results$adj.p.value < 0.05
    
    return(final_results)
    
  }, error = function(e) {
    warning("emmeans approach failed: ", e$message)
    return(NULL)
  })
}

# Helper function to calculate effect sizes
calculate_effect_sizes <- function(long_data) {
  effect_sizes <- data.frame()
  samples <- unique(long_data$Sample)
  
  for (i in 1:(length(samples)-1)) {
    for (j in (i+1):length(samples)) {
      s1 <- samples[i]
      s2 <- samples[j]
      
      x1 <- long_data$Metric_Value[long_data$Sample == s1]
      x2 <- long_data$Metric_Value[long_data$Sample == s2]
      
      cohens_d_result <- tryCatch({
        effectsize::cohens_d(x1, x2)
      }, error = function(e) {
        # Fallback calculation
        pooled_sd <- sqrt(((length(x1)-1)*var(x1) + (length(x2)-1)*var(x2))/(length(x1)+length(x2)-2))
        list(Cohens_d = (mean(x1) - mean(x2)) / pooled_sd, CI_low = NA, CI_high = NA)
      })
      
      effect_size_row <- data.frame(
        contrast = paste(s1, "-", s2),
        Effect_Size = if(is.list(cohens_d_result)) cohens_d_result$Cohens_d else cohens_d_result,
        Effect_Size_CI_Low = if(is.list(cohens_d_result)) cohens_d_result$CI_low else NA,
        Effect_Size_CI_High = if(is.list(cohens_d_result)) cohens_d_result$CI_high else NA,
        stringsAsFactors = FALSE
      )
      effect_sizes <- rbind(effect_sizes, effect_size_row)
    }
  }
  
  return(effect_sizes)
}

# Visualization Functions ----

# Create metric plots
create_plots <- function(per_fold_data, summary_data, pooled_data = NULL) {
  
  # Ensure proper factor levels
  per_fold_data$Sample <- factor(per_fold_data$Sample, levels = unique(per_fold_data$Sample))
  summary_data$Sample <- factor(summary_data$Sample, levels = unique(per_fold_data$Sample))
  
  # Color palette
  n_samples <- length(levels(per_fold_data$Sample))
  custom_colors <- scales::hue_pal()(n_samples)
  names(custom_colors) <- levels(per_fold_data$Sample)
  
  # Metric labels
  metric_labels <- c(
    "Accuracy" = "Accuracy",
    "Kappa" = "Kappa",
    "ARI" = "Adjusted Rand Index", 
    "Macro_F1" = "F1 Score"
  )
  
  plots <- list()
  
  # Create plots for each metric
  metrics_to_plot <- c("Accuracy", "Kappa", "ARI", "Macro_F1")
  
  for (metric in metrics_to_plot) {
    
    # Base plot theme
    base_theme <- theme_minimal() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        axis.text.y = element_text(size = 12, color = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.6),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
      )
    
    if (metric == "Macro_F1" && !is.null(pooled_data)) {
      # Special handling for F1 Score
      plots[[metric]] <- ggplot(pooled_data, aes(x = Sample, y = Macro_F1, fill = Sample)) +
        geom_col(alpha = 0.7, width = 0.6, color = "black", linewidth = 0.3) +
        geom_text(aes(label = sprintf("%.3f", Macro_F1)), 
                  vjust = -0.5, size = 3.5, fontface = "bold") +
        scale_fill_manual(values = custom_colors) +
        base_theme +
        coord_cartesian(ylim = c(0, max(pooled_data$Macro_F1, na.rm = TRUE) * 1.1)) +
        labs(y = metric_labels[metric],
             title = paste(metric_labels[metric], "Performance"),
             subtitle = "Values from pooled STEAM results")
    } else {
      # Standard plots with confidence intervals
      plots[[metric]] <- ggplot() +
        # Boxplots
        geom_boxplot(data = per_fold_data, 
                     aes(x = Sample, y = .data[[metric]], fill = Sample), 
                     outlier.shape = NA, width = 0.6, alpha = 0.7, 
                     linewidth = 0.4, color = "black") +
        # Individual points
        geom_jitter(data = per_fold_data,
                    aes(x = Sample, y = .data[[metric]]), 
                    color = "black", width = 0.15, alpha = 0.6, size = 2) +
        # Confidence intervals
        geom_errorbar(data = summary_data[summary_data$Metric == metric, ],
                      aes(x = Sample, ymin = CI_Lower, ymax = CI_Upper),
                      width = 0.3, linewidth = 0.7, color = "red", alpha = 0.8) +
        # Means
        geom_point(data = summary_data[summary_data$Metric == metric, ],
                   aes(x = Sample, y = Mean),
                   color = "red", size = 4, shape = 18) +
        scale_fill_manual(values = custom_colors) +
        base_theme +
        coord_cartesian(ylim = c(0, 1)) +
        labs(y = metric_labels[metric] %||% metric,
             title = paste(metric_labels[metric], "Performance"),
             subtitle = "Red diamonds = means, Red bars = 95% CI")
    }
  }
  
  return(plots)
}

# Save plots with consistent formatting
save_plots <- function(plots, prefix = "DLPFC_metrics", 
                                  width = 10, height = 8, dpi = 300) {
  
  if (!dir.exists("figures")) {
    dir.create("figures")
  }
  
  for (metric_name in names(plots)) {
    filename <- file.path("figures", paste0(prefix, "_", metric_name, ".png"))
    ggsave(filename, plots[[metric_name]], 
           width = width, height = height, dpi = dpi, bg = "white")
    message("Saved: ", filename)
  }
  
  # Combined plot
  if (length(plots) >= 4) {
    combined_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] +
      plot_layout(ncol = 2, nrow = 2)
    
    combined_filename <- file.path("figures", paste0(prefix, "_combined.png"))
    ggsave(combined_filename, combined_plot, 
           width = 16, height = 12, dpi = dpi, bg = "white")
    message("Saved: ", combined_filename)
  }
}

# Main Analysis Pipeline ----

main_analysis <- function() {
  
  message("=== Starting DLPFC STEAM Analysis ===")
  
  # Initialize storage
  all_per_fold <- list()
  all_summary <- list()
  all_pooled <- list()
  all_steam_results <- vector("list", length(SAMPLE_IDS))
  names(all_steam_results) <- SAMPLE_IDS
  
  # Process each sample
  for (sample_id in SAMPLE_IDS) {
    message("Processing sample ", sample_id, " (", which(SAMPLE_IDS == sample_id), "/", length(SAMPLE_IDS), ")")
    
    # Check if sample exists
    if (!sample_id %in% names(DLPFCs)) {
      warning("Sample ", sample_id, " not found in DLPFC data")
      next
    }
    
    dlpfc_sample <- DLPFCs[[sample_id]]
    
    # Load STEAM
    steam_obj <- LoadSTEAM(
      Seurat.obj = dlpfc_sample,
      label.column = STEAM_CONFIG$label_column,
      assay = STEAM_CONFIG$assay
    )
    
    # Run nested CV
    steam_result <- RunSTEAM(
      steam_obj, 
      mode = "nested",
      n_outer_folds = STEAM_CONFIG$n_outer_folds,
      n_inner_folds = STEAM_CONFIG$n_inner_folds,
      cv.cores = STEAM_CONFIG$cv_cores,
      allowParallel = STEAM_CONFIG$allow_parallel,
      n.size = STEAM_CONFIG$n_size,
      model = STEAM_CONFIG$model,
      metric = STEAM_CONFIG$metric,
      maxnweights = STEAM_CONFIG$max_weights
    )
    
    all_steam_results[[sample_id]] <- steam_result
    
    # Extract per-fold metrics
    if (!is.null(steam_result$train$ncv$outer_result)) {
      fold_metrics <- compute_outer_fold_metrics(steam_result$train$ncv$outer_result)
      fold_metrics$Sample <- sample_id
      all_per_fold[[sample_id]] <- fold_metrics
      
      # Statistical summary
      fold_summary <- summarize_fold_metrics(fold_metrics, ANALYSIS_CONFIG$confidence_level)
      fold_summary$Sample <- sample_id
      all_summary[[sample_id]] <- fold_summary
    } else {
      warning("No outer fold results found for sample ", sample_id)
    }
    
    # Extract pooled metrics
    pooled_metrics <- extract_pooled_metrics(steam_result)
    pooled_metrics$Sample <- sample_id
    all_pooled[[sample_id]] <- pooled_metrics
    
    message("Completed sample ", sample_id)
  }
  
  # Combine results
  message("Combining results across samples...")
  
  per_fold_table <- do.call(rbind, all_per_fold)
  summary_table <- do.call(rbind, all_summary)
  pooled_table <- do.call(rbind, all_pooled)
  
  # Add model labels
  per_fold_table$Model <- "STEAM_RF"
  summary_table$Model <- "STEAM_RF"
  pooled_table$Model <- "STEAM_RF"
  
  # Statistical Analysis
  message("Performing statistical analysis...")
  
  # Pairwise comparisons
  pairwise_results <- list()
  for (metric in ANALYSIS_CONFIG$metrics_to_analyze[1:2]) {  # Focus on key metrics
    pairwise_results[[metric]] <- perform_pairwise_comparisons(
      all_per_fold, metric, ANALYSIS_CONFIG$correction_method
    )
  }
  
  # Create visualizations
  message("Creating plots...")
  
  summary_plot_data <- summary_table[!is.na(summary_table$Sample), ]
  metric_plots <- create_plots(per_fold_table, summary_plot_data, pooled_table)
  
  # Save plots
  save_plots(metric_plots, ANALYSIS_CONFIG$output_prefix)
  
  # Save results
  message("Saving analysis results...")
  
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  results_list <- list(
    per_fold_metrics = per_fold_table,
    summary_statistics = summary_table,
    pooled_metrics = pooled_table,
    pairwise_comparisons = pairwise_results,
    configuration = list(
      steam_config = STEAM_CONFIG,
      analysis_config = ANALYSIS_CONFIG,
      sample_ids = SAMPLE_IDS
    ),
    session_info = sessionInfo()
  )
  
  saveRDS(results_list, file.path("results", paste0(ANALYSIS_CONFIG$output_prefix, "_complete.rds")))
  
  # Summary report
  message("=== Analysis Complete ===")
  message("Samples processed: ", length(SAMPLE_IDS))
  message("Total CV folds: ", nrow(per_fold_table))
  message("Metrics analyzed: ", paste(ANALYSIS_CONFIG$metrics_to_analyze, collapse = ", "))
  message("Results saved to: results/", ANALYSIS_CONFIG$output_prefix, "_complete.rds")
  message("Plots saved to: figures/")
  
  return(results_list)
}

# Execute Analysis ----

# Usage Examples and Options
cat("=== STEAM Multi-Dataset Analysis Framework ===\n")
cat("Available execution modes:\n")
cat("1. Single dataset analysis (default)\n")
cat("2. Multi-dataset comparison\n")
cat("3. Dataset validation and listing\n\n")

# Check available datasets
list_datasets()

# Execution Control Variables
# Set these variables to control analysis execution:

# For single dataset analysis (default)
RUN_SINGLE_ANALYSIS <- TRUE
SINGLE_DATASET_NAME <- CURRENT_DATASET  # Change this to analyze different datasets

# For multi-dataset comparison
RUN_MULTI_ANALYSIS <- FALSE  # Set to TRUE to compare multiple datasets

# For validation only (no analysis)
VALIDATION_ONLY <- FALSE  # Set to TRUE to just validate datasets

# Main Execution Logic ----

if (VALIDATION_ONLY) {
  # Just validate datasets without running analysis
  cat("\n=== VALIDATION MODE ===\n")
  
  for (dataset_name in names(DATASET_CONFIGS)) {
    cat("\nValidating", dataset_name, "...\n")
    validate_dataset(dataset_name)
  }
  
} else if (RUN_MULTI_ANALYSIS) {
  # Multi-dataset comparison mode
  cat("\n=== MULTI-DATASET ANALYSIS MODE ===\n")
  
  # Run multi-dataset analysis
  multi_results <- main_analysis(multi_dataset_mode = TRUE)
  
  # Display summary
  cat("\nMulti-dataset analysis completed!\n")
  cat("Access results using 'multi_results' object\n")
  cat("Individual dataset results: multi_results$individual_results\n")
  cat("Combined metrics: multi_results$per_fold_metrics\n")
  
} else if (RUN_SINGLE_ANALYSIS) {
  # Single dataset analysis mode (default)
  cat("\n=== SINGLE DATASET ANALYSIS MODE ===\n")
  cat("Analyzing dataset:", SINGLE_DATASET_NAME, "\n")
  
  # Validate the selected dataset first
  if (!validate_dataset(SINGLE_DATASET_NAME)) {
    stop("Dataset validation failed for: ", SINGLE_DATASET_NAME)
  }
  
  # Run single dataset analysis
  single_results <- main_analysis(dataset_name = SINGLE_DATASET_NAME, multi_dataset_mode = FALSE)
  
  # Display summary
  cat("\nSingle dataset analysis completed!\n")
  cat("Access results using 'single_results' object\n")
  cat("Per-fold metrics: single_results$per_fold_metrics\n")
  cat("Summary statistics: single_results$summary_statistics\n")
  cat("Pooled metrics: single_results$pooled_metrics\n")
  
} else {
  cat("\n=== NO ANALYSIS SELECTED ===\n")
  cat("Set one of the execution flags to TRUE:\n")
  cat("- RUN_SINGLE_ANALYSIS: Analyze single dataset\n")
  cat("- RUN_MULTI_ANALYSIS: Compare multiple datasets\n")
  cat("- VALIDATION_ONLY: Just validate datasets\n")
}

# Additional Usage Examples ----

cat("\n=== USAGE EXAMPLES ===\n")
cat("# To analyze a specific dataset:\n")
cat("# results <- main_analysis('DLPFC')\n")
cat("# results <- main_analysis('HIPPOCAMPUS')\n\n")

cat("# To compare multiple datasets:\n")
cat("# multi_results <- main_analysis(multi_dataset_mode = TRUE)\n\n")

cat("# To validate a specific dataset:\n")
cat("# validate_dataset('DLPFC')\n\n")

cat("# To list all available datasets:\n")
cat("# list_datasets()\n\n")

cat("# To add a new dataset, modify DATASET_CONFIGS:\n")
cat("# DATASET_CONFIGS$NEW_DATASET <- list(\n")
cat("#   name = 'NEW_DATASET',\n")
cat("#   description = 'Description here',\n")
cat("#   file_path = '/path/to/data.RDS',\n")
cat("#   sample_ids = c('sample1', 'sample2'),\n")
cat("#   label_column = 'cluster_column',\n")
cat("#   assay = 'SCT',\n")
cat("#   type = 'spatial_transcriptomics'\n")
cat("# )\n\n")

# Output Summary ----

if (exists("single_results") || exists("multi_results")) {
  cat("=== OUTPUT SUMMARY ===\n")
  cat("Generated files:\n")
  
  if (dir.exists("results")) {
    result_files <- list.files("results", pattern = "\\.rds$", full.names = TRUE)
    for (file in result_files) {
      cat("- Results:", file, "\n")
    }
  }
  
  if (dir.exists("figures")) {
    figure_files <- list.files("figures", pattern = "\\.png$", full.names = TRUE)
    for (file in head(figure_files, 10)) {  # Show first 10 files
      cat("- Figure:", file, "\n")
    }
    if (length(figure_files) > 10) {
      cat("- ... and", length(figure_files) - 10, "more figures\n")
    }
  }

# Session Information for Reproducibility ----
cat("\n=== SESSION INFORMATION ===\n")
print(sessionInfo())

}
