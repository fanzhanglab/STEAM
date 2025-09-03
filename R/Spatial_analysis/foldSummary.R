# Helper operator for null coalescing
`%||%` <- function(lhs, rhs) {
  if (!is.null(lhs) && length(lhs) > 0) lhs else rhs
}

# Function to show detailed fold-by-fold, iteration-by-iteration summary
showIterativeSummary <- function(STEAM.obj) {
  
  # Check if iterative analysis exists
  if (is.null(STEAM.obj$iterative_anchor_analysis)) {
    cat("No iterative analysis found. Run STEAM_anchor with iterative=TRUE first.\n")
    cat("Example: STEAM.obj <- STEAM_anchor(STEAM.obj, iterative=TRUE, max_iterations=3)\n")
    return(invisible(NULL))
  }
  
  iter_analysis <- STEAM.obj$iterative_anchor_analysis
  cat("=== Iterative Analysis: Fold-by-Fold Summary ===\n\n")
  
  cat("Method:", iter_analysis$method, "\n")
  cat("Iterations completed:", iter_analysis$iterations_completed, "\n")
  cat("Total corrections:", iter_analysis$summary$total_corrections, "\n\n")
  
  # Get number of folds
  n_folds <- length(STEAM.obj$nested$ncv$outer_result)
  
  if (iter_analysis$method == "iterative_with_ground_truth") {
    
    # Create summary table for ground truth method
    cat("=== Per-Fold, Per-Iteration Breakdown (Ground Truth Method) ===\n")
    
    # Initialize summary matrix
    summary_matrix <- matrix(0, nrow = iter_analysis$iterations_completed, ncol = n_folds + 1)
    colnames(summary_matrix) <- c(paste("Fold", 1:n_folds), "Total")
    rownames(summary_matrix) <- paste("Iteration", 1:iter_analysis$iterations_completed)
    
    # Fill in the corrections for each iteration and fold
    for (iter in seq_along(iter_analysis$iteration_results)) {
      if (!is.null(iter_analysis$iteration_results[[iter]]$corrections)) {
        corrections <- iter_analysis$iteration_results[[iter]]$corrections
        
        # Count corrections by fold
        fold_counts <- table(factor(corrections$fold, levels = 1:n_folds))
        summary_matrix[iter, 1:n_folds] <- as.numeric(fold_counts)
        summary_matrix[iter, n_folds + 1] <- sum(fold_counts)
      }
    }
    
    print(summary_matrix)
    
    # Calculate baseline accuracy (RunSTEAM results)
    n_folds <- length(STEAM.obj$nested$ncv$outer_result)
    baseline_accuracy <- numeric(n_folds)
    
    for (f in 1:n_folds) {
      fold_result <- STEAM.obj$nested$ncv$outer_result[[f]]
      if (!is.null(fold_result$preds)) {
        preds_df <- fold_result$preds
        baseline_accuracy[f] <- mean(preds_df$testy == preds_df$predy, na.rm = TRUE)
      }
    }
    
    # Get all corrections
    all_corrections <- STEAM.obj$spatial_anchor_analysis$corrections
    
    # Calculate corrected accuracy (after all corrections)
    corrected_accuracy <- numeric(n_folds)
    
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
        
        # Calculate corrected accuracy
        corrected_accuracy[f] <- mean(current_preds == true_labels, na.rm = TRUE)
      }
    }
    
    # Show baseline vs corrected accuracy
    cat("\n=== BASELINE VS CORRECTED ACCURACY ===\n")
    cat("Baseline accuracy (RunSTEAM):\n")
    for (f in 1:n_folds) {
      cat(sprintf("  Fold %d: %.4f\n", f, baseline_accuracy[f]))
    }
    cat(sprintf("  Mean: %.4f\n\n", mean(baseline_accuracy)))
    
    cat("Corrected accuracy (after all corrections):\n")
    for (f in 1:n_folds) {
      cat(sprintf("  Fold %d: %.4f\n", f, corrected_accuracy[f]))
    }
    cat(sprintf("  Mean: %.4f\n\n", mean(corrected_accuracy)))
    
    # Calculate and show improvement
    improvement <- corrected_accuracy - baseline_accuracy
    total_improvement <- mean(corrected_accuracy) - mean(baseline_accuracy)
    
    cat("ACCURACY IMPROVEMENT:\n")
    for (f in 1:n_folds) {
      cat(sprintf("  Fold %d: %+.4f (%+.2f%%)\n", f, improvement[f], improvement[f] * 100))
    }
    cat(sprintf("  TOTAL: %+.4f (%+.2f%%)\n", total_improvement, total_improvement * 100))
    
  } else {
    # Alternative validation method
    cat("=== Per-Fold, Per-Iteration Breakdown (Alternative Validation) ===\n")
    
    # Initialize summary matrix
    summary_matrix <- matrix(0, nrow = iter_analysis$iterations_completed, ncol = n_folds + 1)
    colnames(summary_matrix) <- c(paste("Fold", 1:n_folds), "Total")
    rownames(summary_matrix) <- paste("Iteration", 1:iter_analysis$iterations_completed)
    
    # Fill in the corrections for each iteration
    for (iter in seq_along(iter_analysis$correction_history)) {
      if (!is.null(iter_analysis$correction_history[[iter]])) {
        corrections <- iter_analysis$correction_history[[iter]]
        
        # Count corrections by fold
        fold_counts <- table(factor(corrections$fold, levels = 1:n_folds))
        summary_matrix[iter, 1:n_folds] <- as.numeric(fold_counts)
        summary_matrix[iter, n_folds + 1] <- sum(fold_counts)
      }
    }
    
    print(summary_matrix)
    
    # Show metrics by iteration if available
    if (!is.null(iter_analysis$iteration_metrics)) {
      cat("\n=== Validation Metrics by Iteration ===\n")
      for (iter in seq_along(iter_analysis$iteration_metrics)) {
        metrics <- iter_analysis$iteration_metrics[[iter]]
        cat(sprintf("Iteration %d: spatial_coherence=%.3f, expression_coherence=%.3f\n",
                   iter, 
                   metrics$spatial_coherence %||% 0,
                   metrics$expression_coherence %||% 0))
      }
    }
  }
  
  # Show cumulative summary
  cat("\n=== Cumulative Summary ===\n")
  cumulative <- colSums(summary_matrix[, 1:n_folds, drop = FALSE])
  total_cumulative <- sum(cumulative)
  
  for (f in 1:n_folds) {
    cat(sprintf("Fold %d: %d total corrections\n", f, cumulative[f]))
  }
  cat(sprintf("Overall: %d total corrections\n", total_cumulative))
  
  cat("\n")
}

# Function to get detailed corrections by iteration and fold
getDetailedCorrections <- function(STEAM.obj, iteration = NULL, fold = NULL) {
  
  if (is.null(STEAM.obj$iterative_anchor_analysis)) {
    stop("No iterative analysis found. Run STEAM_anchor with iterative=TRUE first.")
  }
  
  iter_analysis <- STEAM.obj$iterative_anchor_analysis
  
  if (iter_analysis$method == "iterative_with_ground_truth") {
    all_corrections <- data.frame()
    
    iter_range <- if (is.null(iteration)) {
      seq_along(iter_analysis$iteration_results)
    } else {
      iteration
    }
    
    for (i in iter_range) {
      if (!is.null(iter_analysis$iteration_results[[i]]$corrections)) {
        iter_corr <- iter_analysis$iteration_results[[i]]$corrections
        iter_corr$iteration <- i
        all_corrections <- rbind(all_corrections, iter_corr)
      }
    }
  } else {
    all_corrections <- data.frame()
    
    iter_range <- if (is.null(iteration)) {
      seq_along(iter_analysis$correction_history)
    } else {
      iteration
    }
    
    for (i in iter_range) {
      if (!is.null(iter_analysis$correction_history[[i]])) {
        iter_corr <- iter_analysis$correction_history[[i]]
        iter_corr$iteration <- i
        all_corrections <- rbind(all_corrections, iter_corr)
      }
    }
  }
  
  # Filter by fold if specified
  if (!is.null(fold)) {
    all_corrections <- all_corrections[all_corrections$fold == fold, ]
  }
  
  return(all_corrections)
}
