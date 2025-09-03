#!/usr/bin/env Rscript

#' Plot Accuracy Before and After Iterative Analysis
#' 
#' Creates a box plot showing accuracy improvement from iterative spatial analysis
#' 
#' @param STEAM.obj STEAM object with iterative analysis results
#' @return ggplot object
plotAccuracyBeforeAfter <- function(STEAM.obj) {
  
  library(ggplot2)
  library(dplyr)
  
  if (is.null(STEAM.obj$iterative_anchor_analysis)) {
    stop("No iterative analysis found. Run STEAM_anchor with iterative=TRUE first.")
  }
  
  # Get number of folds
  n_folds <- length(STEAM.obj$nested$ncv$outer_result)
  
  # Calculate baseline accuracy (before any corrections)
  baseline_accuracy <- numeric(n_folds)
  for (f in 1:n_folds) {
    fold_result <- STEAM.obj$nested$ncv$outer_result[[f]]
    if (!is.null(fold_result$preds)) {
      preds_df <- fold_result$preds
      baseline_accuracy[f] <- mean(preds_df$testy == preds_df$predy, na.rm = TRUE)
    }
  }
  
  # Calculate final accuracy (after all corrections)
  final_accuracy <- numeric(n_folds)
  
  # Get all corrections from all iterations
  all_corrections <- NULL
  if (!is.null(STEAM.obj$spatial_anchor_analysis$corrections)) {
    all_corrections <- STEAM.obj$spatial_anchor_analysis$corrections
  }
  
  for (f in 1:n_folds) {
    fold_result <- STEAM.obj$nested$ncv$outer_result[[f]]
    if (!is.null(fold_result$preds)) {
      preds_df <- fold_result$preds
      
      # Start with original predictions
      current_preds <- setNames(preds_df$predy, rownames(preds_df))
      true_labels <- setNames(preds_df$testy, rownames(preds_df))
      
      # Apply all corrections
      if (!is.null(all_corrections)) {
        fold_corrections <- all_corrections[all_corrections$fold == f, ]
        if (nrow(fold_corrections) > 0) {
          current_preds[fold_corrections$cell_id] <- fold_corrections$suggested_correction
        }
      }
      
      # Calculate final accuracy
      final_accuracy[f] <- mean(current_preds == true_labels, na.rm = TRUE)
    }
  }
  
  # Create data frame for plotting
  plot_data <- data.frame(
    Fold = rep(1:n_folds, 2),
    Accuracy = c(baseline_accuracy, final_accuracy),
    Stage = factor(rep(c("Before Spatial Analysis", "After Iterative Analysis"), each = n_folds),
                   levels = c("Before Spatial Analysis", "After Iterative Analysis"))
  )
  
  # Calculate improvement
  improvement <- final_accuracy - baseline_accuracy
  mean_improvement <- mean(improvement)
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = Stage, y = Accuracy)) +
    geom_boxplot(aes(fill = Stage), alpha = 0.7, outlier.shape = NA) +
    geom_point(aes(color = factor(Fold)), size = 3, position = position_jitter(width = 0.2)) +
    scale_fill_manual(values = c("Before Spatial Analysis" = "#E74C3C", 
                                 "After Iterative Analysis" = "#27AE60")) +
    scale_color_discrete(name = "Fold") +
    labs(
      title = "Accuracy Improvement from Iterative Spatial Analysis",
      subtitle = sprintf("Mean improvement: %.4f (%.2f%%)", mean_improvement, mean_improvement * 100),
      x = "Analysis Stage",
      y = "Accuracy",
      caption = sprintf("Based on %d-fold cross-validation", n_folds)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray50"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(title = "Stage"))
  
  # Add improvement statistics as text
  stats_text <- sprintf(
    "Improvements by fold:\n%s\nMean: %.4f\nMedian: %.4f\nRange: %.4f to %.4f",
    paste(sprintf("Fold %d: %+.4f", 1:n_folds, improvement), collapse = "\n"),
    mean(improvement),
    median(improvement),
    min(improvement),
    max(improvement)
  )
  
  # Print summary statistics
  cat("=== Accuracy Improvement Summary ===\n")
  cat("Baseline accuracy (before spatial analysis):\n")
  for (f in 1:n_folds) {
    cat(sprintf("  Fold %d: %.4f\n", f, baseline_accuracy[f]))
  }
  cat(sprintf("  Mean: %.4f\n\n", mean(baseline_accuracy)))
  
  cat("Final accuracy (after iterative analysis):\n")
  for (f in 1:n_folds) {
    cat(sprintf("  Fold %d: %.4f\n", f, final_accuracy[f]))
  }
  cat(sprintf("  Mean: %.4f\n\n", mean(final_accuracy)))
  
  cat("Improvement by fold:\n")
  for (f in 1:n_folds) {
    cat(sprintf("  Fold %d: %+.4f (%+.2f%%)\n", f, improvement[f], improvement[f] * 100))
  }
  cat(sprintf("  Mean improvement: %+.4f (%+.2f%%)\n", mean_improvement, mean_improvement * 100))
  
  return(p)
}

#' Plot Accuracy Progression Across Iterations
#' 
#' Shows how accuracy improves at each iteration
#' 
#' @param STEAM.obj STEAM object with iterative analysis results
#' @return ggplot object
plotAccuracyProgression <- function(STEAM.obj) {
  
  library(ggplot2)
  library(dplyr)
  
  if (is.null(STEAM.obj$iterative_anchor_analysis)) {
    stop("No iterative analysis found. Run STEAM_anchor with iterative=TRUE first.")
  }
  
  iter_analysis <- STEAM.obj$iterative_anchor_analysis
  n_iterations <- iter_analysis$iterations_completed
  n_folds <- length(STEAM.obj$nested$ncv$outer_result)
  
  # Create data frame to store accuracy at each iteration
  accuracy_data <- data.frame()
  
  for (iter in 1:n_iterations) {
    if (!is.null(iter_analysis$iteration_results[[iter]]$improvement_metrics)) {
      fold_accuracies <- iter_analysis$iteration_results[[iter]]$improvement_metrics$updated_accuracy_by_fold
      
      iter_data <- data.frame(
        Iteration = iter,
        Fold = 1:n_folds,
        Accuracy = fold_accuracies
      )
      accuracy_data <- rbind(accuracy_data, iter_data)
    }
  }
  
  if (nrow(accuracy_data) == 0) {
    stop("No accuracy data found in iteration results")
  }
  
  # Create the progression plot
  p <- ggplot(accuracy_data, aes(x = Iteration, y = Accuracy)) +
    geom_boxplot(aes(group = Iteration), alpha = 0.7, fill = "#3498DB") +
    geom_point(aes(color = factor(Fold)), size = 2) +
    scale_x_continuous(breaks = 1:n_iterations) +
    scale_color_discrete(name = "Fold") +
    labs(
      title = "Accuracy Progression Across Iterations",
      x = "Iteration",
      y = "Accuracy",
      caption = sprintf("Based on %d iterations, %d folds", n_iterations, n_folds)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray50"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}
