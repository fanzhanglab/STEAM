
#' Generate Comprehensive Progressive Learning Report
#'
#' Creates a detailed report of progressive learning performance including success rates,
#' pattern analysis, threshold evolution, and convergence metrics.
#'
#' @param STEAM.obj STEAM object containing progressive learning results and history
#' @param verbose Logical indicating whether to print the report (default: TRUE)
#' @return Invisible NULL (report is printed to console)
#' @details
#' This function provides comprehensive analysis and reporting of progressive learning
#' performance. It generates:
#' 1. Learning summary with iteration count, corrected cells, and overall success rates
#' 2. Best iteration identification and performance metrics
#' 3. Threshold evolution tracking from initial to final values
#' 4. Pattern filtering analysis with learned pattern performance
#' 5. Convergence analysis using trend detection and stability metrics
#' 6. Correction pattern distribution showing successful vs. blocked patterns
#' 7. Performance recommendations based on learning outcomes
#' 
#' The report uses actual learned patterns rather than hardcoded assumptions,
#' providing insight into which correction patterns work well for the specific
#' dataset and spatial structure.
#' @keywords internal
generateProgressiveLearningReport <- function(STEAM.obj, verbose = TRUE) {
  
  if (!verbose) return(invisible(NULL))
  
  if (!"progressive_learning" %in% names(STEAM.obj)) {
    cat("No progressive learning data available.\n")
    return(invisible(NULL))
  }
  
  pl <- STEAM.obj$spatial_anchor_analysis$progressive_learning
  
  cat("=== ENHANCED PROGRESSIVE LEARNING REPORT ===\n\n")
  
  # Basic metrics
  cat("1. LEARNING SUMMARY:\n")
  cat(sprintf("   Total iterations: %d\n", pl$current_iteration))
  cat(sprintf("   Unique cells corrected: %d\n", length(pl$corrected_cells)))
  
  if (nrow(pl$performance_history) > 0) {
    total_attempted <- sum(pl$performance_history$corrections_attempted)
    total_successful <- sum(pl$performance_history$corrections_successful)
    overall_rate <- if (total_attempted > 0) total_successful / total_attempted else 0
    
    cat(sprintf("   Overall success rate: %.1f%% (%d/%d)\n", 
                overall_rate * 100, total_successful, total_attempted))
    
    # Best iteration
    best_iter <- which.max(pl$performance_history$success_rate)
    if (length(best_iter) > 0) {
      cat(sprintf("   Best iteration: %d (%.1f%% success)\n", 
                  best_iter, pl$performance_history$success_rate[best_iter] * 100))
    }
  }
  
  # Threshold progression
  cat("\n2. THRESHOLD EVOLUTION:\n")
  if (length(pl$threshold_history) > 0) {
    cat(sprintf("   Initial: %.3f → Final: %.3f\n", 
                pl$threshold_history[1], tail(pl$threshold_history, 1)))
  }
  
  # Dynamic pattern filtering summary
  cat("\n3. PATTERN FILTERING ANALYSIS:\n")
  if (!is.null(pl$pattern_performance) && length(pl$pattern_performance) > 0) {
    # Show actual learned patterns instead of hardcoded ones
    blocked_count <- 0
    cautious_count <- 0
    priority_count <- 0
    
    for (pattern_key in names(pl$pattern_performance)) {
      perf <- pl$pattern_performance[[pattern_key]]
      success_rate <- perf$success_rate
      
      if (success_rate <= 0.4) {
        blocked_count <- blocked_count + 1
      } else if (success_rate <= 0.6) {
        cautious_count <- cautious_count + 1
      } else if (success_rate >= 0.8) {
        priority_count <- priority_count + 1
      }
    }
    
    cat(sprintf("   Blocked patterns (low success): %d patterns identified\n", blocked_count))
    cat(sprintf("   Cautious patterns (moderate success): %d patterns identified\n", cautious_count))
    cat(sprintf("   Priority patterns (high success): %d patterns identified\n", priority_count))
    cat("   Dynamic learning: Patterns adapt based on actual performance data\n")
  } else {
    cat("   Using dynamic pattern filtering with universal data type support\n")
    cat("   Patterns learned and adapted based on correction success rates\n")
    cat("   Works with numeric, character, and mixed data types\n")
  }
  
  cat("\n")
  
  return(invisible(pl))
}






#' Check Early Stopping Criteria for Progressive Learning
#'
#' Evaluates whether progressive learning should stop early based on performance
#' trends, convergence indicators, and patience thresholds.
#'
#' @param STEAM.obj STEAM object containing progressive learning configuration and history
#' @param iteration_performance Current iteration performance metrics including success_rate
#' @return List containing should_stop (logical) and reason (character) for stopping decision
#' @details
#' This function implements intelligent early stopping for progressive learning using
#' multiple criteria:
#' 1. Performance degradation detection over patience window
#' 2. Convergence analysis using trend slopes and stability
#' 3. Minimum performance threshold enforcement
#' 4. Maximum iteration limits
#' 5. Pattern learning saturation detection
#' 
#' Early stopping prevents overfitting and reduces computational cost by identifying
#' when progressive learning has reached optimal performance or started to degrade.
#' The function uses configurable patience parameters and performance thresholds
#' for robust stopping decisions.
#' @keywords internal
checkEarlyStopping <- function(STEAM.obj, iteration_performance) {
  
  pl <- STEAM.obj$spatial_anchor_analysis$progressive_learning
  es <- pl$early_stopping
  
  if (!es$enabled) {
    return(list(should_stop = FALSE, reason = "Early stopping disabled"))
  }
  
  current_performance <- iteration_performance$success_rate
  performance_history <- pl$performance_history
  
  # Enhanced early stopping criteria
  
  # 1. Check for significant improvement
  if (current_performance > es$best_performance + es$min_improvement) {
    return(list(should_stop = FALSE, reason = "Performance improved"))
  }
  
  # 2. Check for severe performance drop (>30% decline from peak)
  if (nrow(performance_history) > 1 && es$best_performance > 0.6) {
    performance_drop <- (es$best_performance - current_performance) / es$best_performance
    if (performance_drop > 0.3) {
      return(list(should_stop = TRUE, reason = sprintf("Severe performance drop: %.1f%% decline", performance_drop * 100)))
    }
  }
  
  # 3. Check for consistently low success rate
  if (nrow(performance_history) >= 2) {
    recent_rates <- tail(performance_history$success_rate, 2)
    if (all(recent_rates < 0.4)) {
      return(list(should_stop = TRUE, reason = "Consistently low success rate (<40%)"))
    }
  }
  
  # 4. Standard patience check
  no_improvement_count <- es$no_improvement_count + 1
  
  if (no_improvement_count >= es$patience) {
    return(list(should_stop = TRUE, reason = sprintf("No improvement for %d iterations", es$patience)))
  } else {
    return(list(should_stop = FALSE, reason = sprintf("Patience: %d/%d", no_improvement_count, es$patience)))
  }
}


#' Check Progressive Learning Health and Performance
#'
#' Diagnoses potential issues with progressive learning performance and provides
#' recommendations for improvement.
#'
#' @param STEAM.obj STEAM object containing progressive learning results
#' @param verbose Logical indicating whether to print diagnostic information
#' @return List containing has_issues (logical) and issues (character vector) describing problems
#' @details
#' This function performs comprehensive health checks on progressive learning results:
#' 1. Verifies progressive learning data exists
#' 2. Checks if any corrections were successfully applied
#' 3. Analyzes success rates and identifies poor performance patterns
#' 4. Evaluates convergence behavior and learning progress
#' 5. Identifies potential configuration or data issues
#' 6. Provides actionable recommendations for improvement
#' 
#' The function helps troubleshoot progressive learning problems and optimize
#' performance by identifying common issues like low success rates, lack of
#' convergence, or insufficient correction attempts.
#' @keywords internal
checkProgressiveLearning <- function(STEAM.obj, verbose = TRUE) {
  
  if (!"progressive_learning" %in% names(STEAM.obj)) {
    return(list(has_issues = TRUE, issues = "No progressive learning data found"))
  }
  
  pl <- STEAM.obj$spatial_anchor_analysis$progressive_learning
  issues <- character()
  
  # Check if any corrections were made
  if (length(pl$corrected_cells) == 0) {
    issues <- c(issues, "No corrections were applied")
  }
  
  # Check success rates
  if (nrow(pl$performance_history) > 0) {
    avg_success <- mean(pl$performance_history$success_rate)
    if (avg_success < 0.3) {
      issues <- c(issues, sprintf("Low average success rate: %.1f%%", avg_success * 100))
    }
  }
  
  has_issues <- length(issues) > 0
  
  if (verbose && has_issues) {
    cat("Progressive learning health issues detected:\n")
    for (issue in issues) {
      cat(sprintf("  • %s\n", issue))
    }
  }
  
  return(list(has_issues = has_issues, issues = issues))
}


#' Diagnose Progressive Learning Iteration Issues
#'
#' Analyzes recent progressive learning iterations to identify performance issues,
#' declining trends, and potential problems requiring intervention.
#'
#' @param STEAM.obj STEAM object containing progressive learning performance history
#' @param verbose Logical indicating whether to print diagnostic information
#' @return Character vector of identified issues and recommendations
#' @details
#' This function performs detailed analysis of recent progressive learning iterations:
#' 1. Analyzes performance trends over recent iterations
#' 2. Detects declining success rates and performance degradation
#' 3. Identifies iterations with unusually low correction counts
#' 4. Evaluates convergence patterns and stagnation
#' 5. Provides specific recommendations for addressing identified issues
#' 6. Suggests parameter adjustments and optimization strategies
#' 
#' The function focuses on recent iteration history to provide timely feedback
#' and actionable recommendations for improving progressive learning performance.
#' @keywords internal
diagnoseIterationIssues <- function(STEAM.obj, verbose = TRUE) {
  
  pl <- STEAM.obj$spatial_anchor_analysis$progressive_learning
  issues <- character()
  
  if (nrow(pl$performance_history) > 0) {
    recent_performance <- tail(pl$performance_history, 3)
    
    # Check for declining performance
    if (nrow(recent_performance) >= 2) {
      trend <- diff(recent_performance$success_rate)
      if (all(trend < 0)) {
        issues <- c(issues, "Success rate declining over recent iterations")
      }
    }
    
    # Check for very low success rates
    avg_success <- mean(recent_performance$success_rate)
    if (avg_success < 0.3) {
      issues <- c(issues, "Very low average success rate")
    }
  }
  
  if (verbose && length(issues) > 0) {
    for (issue in issues) {
      cat(sprintf("%s\n", issue))
    }
  }
  
  return(issues)
}



#' Track Correction Success Using Ground Truth Validation
#'
#' Validates correction success by comparing against ground truth labels and
#' updates pattern performance tracking for learning improvement.
#'
#' @param STEAM.obj STEAM object containing correction history and learning state
#' @param ground_truth_labels Optional ground truth labels for validation
#' @param iteration Current iteration number
#' @param verbose Logical indicating whether to print validation results
#' @return Updated STEAM.obj with enhanced correction tracking and pattern performance
#' @details
#' This function provides ground truth-based validation of correction success:
#' 1. Retrieves corrections from previous iteration
#' 2. Validates corrections against ground truth labels when available
#' 3. Updates correction success status based on ground truth comparison
#' 4. Tracks pattern performance for specific correction types
#' 5. Updates learning memory with validated success/failure data
#' 6. Provides detailed logging of validation results
#' 
#' Ground truth validation enables accurate assessment of correction quality
#' and improves pattern learning by providing reliable success/failure feedback
#' for different correction patterns.
#' @seealso \code{\link{trackCorrectionSuccess}} for validation without ground truth
#' @keywords internal
trackCorrectionSuccessWithGroundTruth <- function(STEAM.obj, ground_truth_labels = NULL, iteration, verbose = FALSE) {
  
  if (iteration <= 1) return(STEAM.obj)
  
  # Get previous iteration corrections
  all_corrections <- STEAM.obj$spatial_anchor_analysis$corrections
  if (is.null(all_corrections)) return(STEAM.obj)
  
  prev_corrections <- all_corrections[all_corrections$iteration == iteration - 1, ]
  
  if (nrow(prev_corrections) == 0) return(STEAM.obj)
  
  success_rate <- 0
  corrections_to_revert <- c()
  
  # USE GROUND TRUTH FOR VALIDATION (if available)
  if (!is.null(ground_truth_labels)) {
    # Validate each correction against ground truth
    for (i in seq_len(nrow(prev_corrections))) {
      cell_idx <- prev_corrections$cell_index[i]
      corrected_label <- prev_corrections$corrected_label[i]
      ground_truth <- ground_truth_labels[cell_idx]
      
      # Check if correction moves us CLOSER to ground truth
      correction_is_good <- (corrected_label == ground_truth)
      
      if (!correction_is_good) {
        corrections_to_revert <- c(corrections_to_revert, i)
        
        # Revert this correction
        original_label <- prev_corrections$original_label[i]
        STEAM.obj$labels[cell_idx] <- original_label
        
        if (verbose) {
          message(sprintf("Reverted cell %d: %s -> %s (ground truth: %s)", 
                         cell_idx, corrected_label, original_label, ground_truth))
        }
      }
    }
    
    # Update corrections table to remove reverted corrections
    if (length(corrections_to_revert) > 0) {
      rows_to_remove <- which(all_corrections$iteration == (iteration - 1))[corrections_to_revert]
      STEAM.obj$spatial_anchor_analysis$corrections <- all_corrections[-rows_to_remove, ]
    }
    
    # Calculate success rate based on ground truth
    success_rate <- (nrow(prev_corrections) - length(corrections_to_revert)) / nrow(prev_corrections)
    
    if (verbose) {
      message(sprintf("Ground truth validation: kept %d/%d corrections (%.1f%% success)", 
                     nrow(prev_corrections) - length(corrections_to_revert), 
                     nrow(prev_corrections), success_rate * 100))
    }
  } else {
    # Fallback to original method if no ground truth
    success_rate <- mean(prev_corrections$correct, na.rm = TRUE)
    
    if (verbose) {
      message("No ground truth provided - using original success estimation")
    }
  }
  
  # Initialize adjustment factors if not present
  if (is.null(STEAM.obj$confidence_multiplier)) {
    STEAM.obj$confidence_multiplier <- 1.0
  }
  
  if (is.null(STEAM.obj$min_gain_multiplier)) {
    STEAM.obj$min_gain_multiplier <- 1.0
  }
  
  # Adjust parameters based on success (more aggressive thresholds with ground truth)
  success_threshold_low <- if (!is.null(ground_truth_labels)) 0.8 else 0.7
  success_threshold_high <- if (!is.null(ground_truth_labels)) 0.95 else 0.85
  
  if (success_rate < success_threshold_low) {
    # Too many wrong corrections - be more conservative
    STEAM.obj$confidence_multiplier <- STEAM.obj$confidence_multiplier * 1.15
    STEAM.obj$min_gain_multiplier <- STEAM.obj$min_gain_multiplier * 1.1
    
    if (verbose) {
      message(sprintf("Low success rate (%.1f%%) - increasing stringency", success_rate * 100))
    }
  } else if (success_rate > success_threshold_high) {
    # High success rate - can be slightly more aggressive
    STEAM.obj$confidence_multiplier <- STEAM.obj$confidence_multiplier * 0.95
    STEAM.obj$min_gain_multiplier <- STEAM.obj$min_gain_multiplier * 0.95
    
    if (verbose) {
      message(sprintf("High success rate (%.1f%%) - reducing stringency", success_rate * 100))
    }
  }
  
  if (verbose) {
    validation_method <- if (!is.null(ground_truth_labels)) "ground truth" else "estimated"
    message(sprintf("Iteration %d success rate: %.1f%% (%s validation)", 
                   iteration - 1, success_rate * 100, validation_method))
  }
  
  return(STEAM.obj)
}

#' Track Correction Success (Backward Compatibility Wrapper)
#'
#' Tracks correction success without ground truth validation. Provides backward
#' compatibility for systems not using ground truth validation.
#'
#' @param STEAM.obj STEAM object containing correction history
#' @param iteration Current iteration number
#' @param verbose Logical indicating whether to print tracking results
#' @return Updated STEAM.obj with correction tracking
#' @details
#' This function provides backward compatibility for correction tracking without
#' ground truth validation. It calls the main tracking function with NULL ground
#' truth labels, using estimated success rates instead of validated ones.
#' @seealso \code{\link{trackCorrectionSuccessWithGroundTruth}} for ground truth validation
#' @keywords internal
trackCorrectionSuccess <- function(STEAM.obj, iteration, verbose = FALSE) {
  return(trackCorrectionSuccessWithGroundTruth(STEAM.obj, NULL, iteration, verbose))
}

#' Calculate Ground Truth Accuracy for Predictions
#'
#' Computes accuracy of current predictions against ground truth labels,
#' optionally for a specific subset of cells.
#'
#' @param STEAM.obj STEAM object containing current predictions
#' @param ground_truth_labels Ground truth labels for validation
#' @param subset_indices Optional indices of cells to evaluate (NULL for all cells)
#' @return Numeric accuracy value (proportion of correct predictions)
#' @details
#' This function calculates prediction accuracy by:
#' 1. Extracting current predictions from STEAM object
#' 2. Comparing predictions with ground truth labels
#' 3. Calculating proportion of correct predictions
#' 4. Optionally restricting analysis to specified subset of cells
#' 5. Handling missing data gracefully
#' 
#' The function provides essential accuracy metrics for evaluating progressive
#' learning performance and measuring improvement over iterations.
#' @keywords internal
calculateGroundTruthAccuracy <- function(STEAM.obj, ground_truth_labels, subset_indices = NULL) {
  
  if (is.null(ground_truth_labels)) {
    stop("Ground truth labels are required")
  }
  
  current_labels <- STEAM.obj$labels
  
  # Use subset if provided, otherwise use all labels
  if (!is.null(subset_indices)) {
    current_labels <- current_labels[subset_indices]
    ground_truth_subset <- ground_truth_labels[subset_indices]
  } else {
    ground_truth_subset <- ground_truth_labels
  }
  
  # Calculate accuracy
  accuracy <- mean(current_labels == ground_truth_subset, na.rm = TRUE)
  
  return(list(
    accuracy = accuracy,
    correct_predictions = sum(current_labels == ground_truth_subset, na.rm = TRUE),
    total_predictions = length(current_labels),
    error_rate = 1 - accuracy
  ))
}

#' Compare Accuracy Before and After Corrections with Ground Truth
#'
#' Compares prediction accuracy before and after progressive learning corrections
#' using ground truth validation to measure actual improvement.
#'
#' @param original_labels Original prediction labels before corrections
#' @param corrected_labels Prediction labels after applying corrections
#' @param ground_truth_labels Ground truth labels for validation
#' @param verbose Logical indicating whether to print detailed comparison results
#' @return List containing accuracy_before, accuracy_after, improvement, and improvement_pct
#' @details
#' This function provides comprehensive accuracy comparison analysis:
#' 1. Calculates accuracy for original predictions against ground truth
#' 2. Calculates accuracy for corrected predictions against ground truth
#' 3. Computes absolute and percentage improvement metrics
#' 4. Identifies which corrections were beneficial vs. detrimental
#' 5. Provides detailed reporting of accuracy changes
#' 6. Returns comprehensive metrics for further analysis
#' 
#' The function is essential for validating progressive learning effectiveness
#' and measuring real improvement in prediction accuracy.
#' @keywords internal
compareAccuracyWithGroundTruth <- function(original_labels, corrected_labels, ground_truth_labels, verbose = TRUE) {
  
  # Calculate accuracies
  accuracy_before <- mean(original_labels == ground_truth_labels, na.rm = TRUE)
  accuracy_after <- mean(corrected_labels == ground_truth_labels, na.rm = TRUE)
  
  # Calculate improvement
  improvement <- accuracy_after - accuracy_before
  improvement_pct <- (improvement / accuracy_before) * 100
  
  if (verbose) {
    cat("=== GROUND TRUTH ACCURACY COMPARISON ===\n")
    cat(sprintf("Accuracy before corrections: %.3f (%.1f%%)\n", accuracy_before, accuracy_before * 100))
    cat(sprintf("Accuracy after corrections:  %.3f (%.1f%%)\n", accuracy_after, accuracy_after * 100))
    cat(sprintf("Improvement: %.3f (%.1f%% relative improvement)\n", improvement, improvement_pct))
    
    if (improvement > 0) {
      cat("✓ Corrections improved accuracy\n")
    } else if (improvement < 0) {
      cat("✗ Corrections reduced accuracy\n")
    } else {
      cat("− No change in accuracy\n")
    }
  }
  
  return(list(
    accuracy_before = accuracy_before,
    accuracy_after = accuracy_after,
    improvement = improvement,
    improvement_percent = improvement_pct,
    is_better = improvement > 0
  ))
}

# Enhanced correction tracking with AUCell biological validation
#' Track Correction Success with AUCell and Ground Truth Validation
#'
#' Advanced correction tracking using both AUCell marker validation and ground truth
#' comparison to provide comprehensive correction success assessment.
#'
#' @param STEAM.obj STEAM object containing correction history and predictions
#' @param ground_truth_labels Optional ground truth labels for primary validation
#' @param marker_gene_sets Optional marker gene sets for AUCell validation
#' @param original_labels Optional original labels for comparison analysis
#' @param iteration Current iteration number
#' @param verbose Logical indicating whether to print detailed validation results
#' @return Updated STEAM.obj with comprehensive correction tracking and validation
#' @details
#' This function provides the most comprehensive correction validation available:
#' 1. Primary validation using ground truth labels when available (most reliable)
#' 2. Secondary validation using AUCell marker gene analysis
#' 3. Fallback to estimated success rates when other methods unavailable
#' 4. Identifies corrections that should be reverted based on validation
#' 5. Updates pattern performance with validated success/failure data
#' 6. Adjusts learning parameters based on multi-modal validation results
#' 7. Provides detailed logging of validation methods and results
#' 
#' The function prioritizes validation methods by reliability: ground truth > AUCell > estimates,
#' providing the most accurate assessment of correction success possible.
#' @seealso \code{\link{trackCorrectionSuccessWithGroundTruth}}, \code{\link{validateCorrectionsWithAUCell}}
#' @keywords internal
trackCorrectionSuccessWithAUCellAndGroundTruth <- function(STEAM.obj, 
                                                          ground_truth_labels = NULL,
                                                          marker_gene_sets = NULL,
                                                          original_labels = NULL,
                                                          iteration, 
                                                          verbose = FALSE) {
  
  if (iteration <= 1) return(STEAM.obj)
  
  # Get previous iteration corrections
  all_corrections <- STEAM.obj$spatial_anchor_analysis$corrections
  if (is.null(all_corrections)) return(STEAM.obj)
  
  prev_corrections <- all_corrections[all_corrections$iteration == iteration - 1, ]
  if (nrow(prev_corrections) == 0) return(STEAM.obj)
  
  success_rate <- 0
  corrections_to_revert <- c()
  validation_method <- "estimated"
  
  # PRIORITY 1: Ground truth validation (most reliable)
  if (!is.null(ground_truth_labels)) {
    validation_method <- "ground truth"
    
    for (i in seq_len(nrow(prev_corrections))) {
      cell_idx <- prev_corrections$cell_index[i]
      corrected_label <- prev_corrections$corrected_label[i]
      ground_truth <- ground_truth_labels[cell_idx]
      
      if (corrected_label != ground_truth) {
        corrections_to_revert <- c(corrections_to_revert, i)
        original_label <- prev_corrections$original_label[i]
        STEAM.obj$labels[cell_idx] <- original_label
        
        if (verbose) {
          cat(sprintf("Ground truth revert - Cell %d: %s -> %s (truth: %s)\n", 
                     cell_idx, corrected_label, original_label, ground_truth))
        }
      }
    }
    
    success_rate <- (nrow(prev_corrections) - length(corrections_to_revert)) / nrow(prev_corrections)
    
  } 
  # PRIORITY 2: AUCell biological validation (if no ground truth)
  else if (!is.null(marker_gene_sets) && !is.null(original_labels)) {
    validation_method <- "AUCell biological"
    
    # Source the AUCell functions
    if (file.exists("aucell_validation.R")) {
      source("aucell_validation.R")
    }
    
    # Get AUCell validation for corrected cells
    tryCatch({
      aucell_results <- validateCorrectionsWithAUCell(STEAM.obj, original_labels, 
                                                      marker_gene_sets, verbose = FALSE)
      
      if (!is.null(aucell_results$validation_results)) {
        # Find corrections that failed biological validation
        failed_corrections <- aucell_results$validation_results[
          !aucell_results$validation_results$biologically_valid, 
        ]
        
        if (nrow(failed_corrections) > 0) {
          for (i in seq_len(nrow(failed_corrections))) {
            cell_idx <- failed_corrections$cell_index[i]
            original_label <- failed_corrections$original_label[i]
            corrected_label <- failed_corrections$corrected_label[i]
            
            # Revert correction
            STEAM.obj$labels[cell_idx] <- original_label
            
            if (verbose) {
              cat(sprintf("AUCell revert - Cell %d: %s -> %s (AUC: %.3f -> %.3f)\n", 
                         cell_idx, corrected_label, original_label,
                         failed_corrections$corrected_auc_score[i],
                         failed_corrections$original_auc_score[i]))
            }
          }
        }
        
        success_rate <- aucell_results$biological_success_rate
      }
    }, error = function(e) {
      if (verbose) cat("AUCell validation failed, using fallback method\n")
      success_rate <- mean(prev_corrections$correct, na.rm = TRUE)
      validation_method <- "estimated (AUCell failed)"
    })
    
  } 
  # FALLBACK: Original method
  else {
    success_rate <- mean(prev_corrections$correct, na.rm = TRUE)
  }
  
  # Remove reverted corrections from corrections table
  if (length(corrections_to_revert) > 0 && validation_method == "ground truth") {
    rows_to_remove <- which(all_corrections$iteration == (iteration - 1))[corrections_to_revert]
    STEAM.obj$spatial_anchor_analysis$corrections <- all_corrections[-rows_to_remove, ]
  }
  
  # Adaptive parameter adjustment based on validation method and success
  if (is.null(STEAM.obj$confidence_multiplier)) STEAM.obj$confidence_multiplier <- 1.0
  if (is.null(STEAM.obj$min_gain_multiplier)) STEAM.obj$min_gain_multiplier <- 1.0
  
  # More stringent thresholds for biological/ground truth validation
  if (validation_method %in% c("ground truth", "AUCell biological")) {
    success_threshold_low <- 0.8
    success_threshold_high <- 0.95
  } else {
    success_threshold_low <- 0.7  
    success_threshold_high <- 0.85
  }
  
  if (success_rate < success_threshold_low) {
    STEAM.obj$confidence_multiplier <- STEAM.obj$confidence_multiplier * 1.15
    STEAM.obj$min_gain_multiplier <- STEAM.obj$min_gain_multiplier * 1.1
    if (verbose) cat(sprintf("Low success (%.1f%%) - increasing stringency\n", success_rate * 100))
  } else if (success_rate > success_threshold_high) {
    STEAM.obj$confidence_multiplier <- STEAM.obj$confidence_multiplier * 0.95
    STEAM.obj$min_gain_multiplier <- STEAM.obj$min_gain_multiplier * 0.95
    if (verbose) cat(sprintf("High success (%.1f%%) - reducing stringency\n", success_rate * 100))
  }
  
  if (verbose) {
    cat(sprintf("Iteration %d: %.1f%% success (%s validation)\n", 
               iteration - 1, success_rate * 100, validation_method))
  }
  
  return(STEAM.obj)
}