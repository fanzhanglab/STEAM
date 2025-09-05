
#' Execute Progressive Learning for Single Fold Dataset
#' 
#' Runs progressive learning iterations on a single fold or entire dataset until convergence
#' criteria are met. This function handles the iterative correction process with adaptive
#' thresholding and early stopping mechanisms for optimal performance.
#' 
#' @param STEAM.obj STEAM object containing spatial data and analysis results
#' @param correction_function Function to execute corrections, typically runProgressiveIteration
#' @param verbose Logical indicating whether to print detailed progress information
#' 
#' @return Updated STEAM object with progressive learning results integrated into spatial_anchor_analysis structure
#' 
#' @details
#' This function implements the core progressive learning loop for single-fold processing:
#' - Executes iterative corrections using executeProgressiveIteration
#' - Monitors convergence through adaptive thresholds and early stopping
#' - Integrates final results into spatial_anchor_analysis structure
#' - Handles both single-dataset and fold-specific processing scenarios
#' The function continues learning until convergence criteria or maximum iterations are reached.

runProgressiveLearningSingleFold <- function(STEAM.obj, correction_function, verbose = TRUE) {
  
  # Use the original single-dataset logic
  continue_learning <- TRUE
  
  while (continue_learning) {
    iteration_result <- executeProgressiveIteration(STEAM.obj, correction_function, verbose)
    
    continue_learning <- iteration_result$continue
    if ("STEAM.obj" %in% names(iteration_result)) {
      STEAM.obj <- iteration_result$STEAM.obj
    }
    
    if (!continue_learning) {
      if (verbose) {
        cat(sprintf("Progressive learning for this fold completed: %s\n", iteration_result$reason))
      }
      break
    }
  }
  
  # Integrate results into spatial_anchor_analysis structure
  STEAM.obj <- integrateProgressiveResults(STEAM.obj)
  
  return(STEAM.obj)
}





#' Execute Single Progressive Learning Iteration
#' 
#' Performs one iteration of progressive learning by identifying uncertain cells,
#' generating corrections, and applying them with pattern learning integration.
#' This function implements the core correction logic with adaptive thresholding.
#' 
#' @param STEAM.obj STEAM object containing spatial data and analysis results
#' @param expression_matrix Expression data matrix for correction generation
#' @param threshold Confidence threshold for identifying uncertain cells
#' @param max_corrections Maximum number of corrections to attempt per iteration
#' @param iteration Current iteration number for tracking progress
#' @param progressive_state State object tracking corrected cells and patterns
#' @param k_grid Grid of k values for nearest neighbor analysis
#' @param scoring_weights Weights for scoring correction quality
#' @param verbose Logical indicating whether to print detailed progress information
#' 
#' @return List containing updated STEAM object, progressive state, performance metrics, and correction data
#' 
#' @details
#' This function implements a single iteration of the progressive learning algorithm:
#' - Extracts predictions and calculates confidence scores for all cells
#' - Identifies uncertain cells below the confidence threshold (excluding already corrected)
#' - Limits corrections to maximum allowed per iteration, prioritizing lowest confidence
#' - Generates progressive corrections using spatial and expression data
#' - Applies corrections and tracks success rates with pattern learning
#' - Updates correction memory and analyzes correction patterns
#' Returns comprehensive metrics including success rates, confidence scores, and pattern analysis.

runProgressiveIteration <- function(
    STEAM.obj, expression_matrix, threshold, max_corrections, 
    iteration, progressive_state, k_grid, scoring_weights, verbose = TRUE
) {
  
  if (verbose) {
    cat(sprintf("=== Progressive STEAM Corrections (Threshold: %.3f) ===\n", threshold))
  }
  
  # Extract predictions and identify uncertain cells
  all_predictions <- extractPreds(STEAM.obj)
  
  if (verbose) {
    cat(sprintf("Using expression data from STEAM.obj$count_exp\n"))
    pred_source <- attr(all_predictions, "source_field")
    if (is.null(pred_source)) pred_source <- "unknown"
    cat(sprintf("Predictions extracted from: %s\n", pred_source))
    cat(sprintf("Total cells: %d\n", length(all_predictions)))
    cat(sprintf("Prediction range: %d - %d\n", min(all_predictions), max(all_predictions)))
    cat(sprintf("Already corrected: %d\n", length(progressive_state$corrected_cells)))
  }
  
  # Calculate confidence scores
  confidence_scores <- calculateConfidenceScores(STEAM.obj, all_predictions)
  
  # Identify uncertain cells (excluding already corrected ones)
  uncertain_mask <- confidence_scores < threshold & 
    !(names(all_predictions) %in% progressive_state$corrected_cells)
  uncertain_cells <- names(all_predictions)[uncertain_mask]
  
  if (verbose) {
    cat(sprintf("Identified %d uncertain cells (confidence < %.3f)\n", 
                length(uncertain_cells), threshold))
  }
  
  # Limit to max corrections
  if (length(uncertain_cells) > max_corrections) {
    # Select cells with lowest confidence first
    uncertain_confidence <- confidence_scores[uncertain_cells]
    uncertain_cells <- uncertain_cells[order(uncertain_confidence)][1:max_corrections]
    
    if (verbose) {
      cat(sprintf("Limited to %d corrections this iteration\n", max_corrections))
    }
  }
  
  if (length(uncertain_cells) == 0) {
    if (verbose) {
      cat("No uncertain cells found for correction.\n")
    }
    return(list(
      STEAM.obj = STEAM.obj,
      progressive_state = progressive_state,
      corrections_attempted = 0,
      corrections_successful = 0,
      success_rate = 1.0,
      avg_confidence = mean(confidence_scores),
      improvement = 0,
      correction_patterns = list()
    ))
  }
  
  if (verbose) {
    cat(sprintf("Evaluating %d uncertain cells\n", length(uncertain_cells)))
  }
  
  # Generate corrections using progressive approach
  corrections_result <- generateProgressiveCorrections(
    STEAM.obj = STEAM.obj,
    expression_matrix = expression_matrix,
    uncertain_cells = uncertain_cells,
    k_grid = k_grid,
    scoring_weights = scoring_weights,
    verbose = verbose
  )
  
  # Apply corrections and track success
  application_result <- applyProgressiveCorrections(
    STEAM.obj = STEAM.obj,
    corrections = corrections_result$corrections,
    verbose = verbose
  )
  
  # Update pattern learning memory with correction results
  if (nrow(corrections_result$corrections) > 0) {
    application_result$STEAM.obj <- updatePatternPerformance(
      STEAM.obj = application_result$STEAM.obj,
      corrections = corrections_result$corrections,
      verbose = verbose
    )
  }
  
  # Update corrected cells list
  progressive_state$corrected_cells <- unique(c(
    progressive_state$corrected_cells,
    application_result$successful_cells
  ))
  
  # Track correction patterns
  if (length(corrections_result$corrections) > 0) {
    progressive_state$correction_patterns[[iteration]] <- analyzeCorrections(corrections_result$corrections)
  }
  
  # Enhanced tracking output
  if (verbose) {
    cat(sprintf("Corrections attempted: %d\n", nrow(corrections_result$corrections)))
    cat(sprintf("Corrections successful: %d\n", application_result$successful_count))
    cat(sprintf("Success rate: %.1f%%\n", application_result$success_rate * 100))
    
    cat("\n--- Enhanced Correction Tracking ---\n")
    cat(sprintf("Success Rate: %.1f%% (%d/%d)\n", 
                application_result$success_rate * 100,
                application_result$successful_count,
                nrow(corrections_result$corrections)))
    cat(sprintf("Average Confidence: %.3f\n", corrections_result$avg_confidence))
    cat(sprintf("Correction Types: %d different patterns\n", 
                length(unique(corrections_result$correction_types))))
    cat("Top Successful Patterns:\n")
  }
  
  return(list(
    STEAM.obj = application_result$STEAM.obj,
    progressive_state = progressive_state,
    corrections_attempted = nrow(corrections_result$corrections),
    corrections_successful = application_result$successful_count,
    success_rate = application_result$success_rate,
    avg_confidence = corrections_result$avg_confidence,
    improvement = application_result$improvement,
    correction_patterns = corrections_result$correction_types,
    corrections_data = corrections_result$corrections  # Store actual corrections dataframe
  ))
}


#' Execute Complete Progressive Learning Process
#' 
#' Main entry point for progressive learning that automatically detects and handles
#' both single-dataset and multi-fold scenarios. Orchestrates the entire progressive
#' learning workflow with comprehensive monitoring and reporting.
#' 
#' @param STEAM.obj STEAM object containing spatial data and analysis results
#' @param correction_function Function to execute corrections for each iteration
#' @param verbose Logical indicating whether to print detailed progress information
#' 
#' @return Updated STEAM object with complete progressive learning results and analytics
#' 
#' @details
#' This function serves as the main orchestrator for progressive learning:
#' - Automatically detects presence of cross-validation fold structure
#' - Routes to appropriate processing function (single-fold or multi-fold)
#' - Ensures proper initialization of progressive learning framework
#' - Executes iterative learning with convergence monitoring
#' - Provides comprehensive final reporting including success rates and coverage
#' - Performs health checks and suggests improvements for poor performance
#' - Integrates all results into the spatial_anchor_analysis structure
#' For multi-fold scenarios, delegates to runProgressiveLearningMultiFold for parallel fold processing.

runProgressiveLearning <- function(STEAM.obj, correction_function, verbose = TRUE) {
  
  if (verbose) cat("=== Starting Progressive Learning Process ===\n")
  
  # Check if we have multiple folds to process
  has_folds <- !is.null(STEAM.obj$nested) && 
    !is.null(STEAM.obj$nested$ncv) && 
    !is.null(STEAM.obj$nested$ncv$outer_result)
  
  if (has_folds) {
    n_folds <- length(STEAM.obj$nested$ncv$outer_result)
    if (verbose) cat(sprintf("Found %d folds for processing\n", n_folds))
    
    # Process each fold separately
    return(runProgressiveLearningMultiFold(STEAM.obj, correction_function, verbose))
  } else {
    if (verbose) cat("No fold structure found - processing entire dataset\n")
  }
  
  # Initialize if not already done (preserve existing settings)
  if (is.null(STEAM.obj$spatial_anchor_analysis$progressive_learning)) {
    # Only initialize if truly missing - this should not happen if STEAM_anchor called properly
    warning("Progressive learning not initialized properly - using default settings")
    STEAM.obj <- initializeProgressiveLearning(STEAM.obj)
  }
  
  continue_learning <- TRUE
  
  while (continue_learning) {
    iteration_result <- executeProgressiveIteration(STEAM.obj, correction_function, verbose)
    
    continue_learning <- iteration_result$continue
    if ("STEAM.obj" %in% names(iteration_result)) {
      STEAM.obj <- iteration_result$STEAM.obj
    }
    
    if (!continue_learning) {
      if (verbose) {
        cat(sprintf("\n=== Progressive Learning Completed ===\n"))
        cat(sprintf("Reason: %s\n", iteration_result$reason))
        cat(sprintf("Total iterations: %d\n", STEAM.obj$spatial_anchor_analysis$progressive_learning$current_iteration))
        cat(sprintf("Total unique cells corrected: %d\n", 
                    length(STEAM.obj$spatial_anchor_analysis$progressive_learning$corrected_cells)))
        
        # Enhanced Final Report
        cat("\n--- Final Analysis ---\n")
        generateProgressiveLearningReport(STEAM.obj, verbose = TRUE)
        
        # Final health check
        final_health <- checkProgressiveLearning(STEAM.obj, verbose = FALSE)
        # Validate final_health structure and has_issues field
        if (!is.null(final_health) && is.list(final_health) && 
            "has_issues" %in% names(final_health) && 
            !is.null(final_health$has_issues) && 
            length(final_health$has_issues) > 0 && 
            !is.na(final_health$has_issues) && 
            final_health$has_issues) {
          cat("\nFinal Issues Detected:\n")
          if ("issues" %in% names(final_health) && !is.null(final_health$issues)) {
            for (issue in final_health$issues) {
              cat(sprintf("   • %s\n", issue))
            }
          }
          cat("\n--- Recommendations for Future Runs ---\n")
          suggestIterationFixes(STEAM.obj, verbose = TRUE)
        } else {
          cat("\nProgressive learning completed successfully with no major issues!\n")
        }
      }
      break
    }
  }
  
  return(STEAM.obj)
}


#' Execute Progressive Learning Across Multiple Cross-Validation Folds
#' 
#' Processes progressive learning across all cross-validation folds in parallel,
#' aggregating results and providing comprehensive multi-fold performance analysis.
#' This function handles the complexity of fold-specific processing and result integration.
#' 
#' @param STEAM.obj STEAM object containing nested cross-validation structure with multiple folds
#' @param correction_function Function to execute corrections, typically runProgressiveIteration
#' @param verbose Logical indicating whether to print detailed progress information
#' 
#' @return Updated STEAM object with aggregated multi-fold progressive learning results
#' 
#' @details
#' This function orchestrates progressive learning across multiple cross-validation folds:
#' - Extracts fold structure from nested cross-validation results
#' - Processes each fold independently using runProgressiveLearningSingleFold
#' - Tracks fold-specific corrections and assigns proper fold identifiers
#' - Aggregates performance metrics across all folds
#' - Provides detailed fold-by-fold performance breakdown
#' - Handles both 'correct' and 'is_correct' column name variations for compatibility
#' - Integrates all results into unified spatial_anchor_analysis structure
#' Essential for multi-fold cross-validation scenarios where corrections need fold-specific tracking.

runProgressiveLearningMultiFold <- function(STEAM.obj, correction_function, verbose = TRUE) {
  
  outer_result <- STEAM.obj$nested$ncv$outer_result
  n_folds <- length(outer_result)
  
  if (verbose) cat(sprintf("Processing %d folds with progressive learning\n\n", n_folds))
  
  # Ensure progressive learning structure exists (preserve existing settings)
  if (is.null(STEAM.obj$spatial_anchor_analysis$progressive_learning)) {
    # Only initialize if truly missing - this should not happen if STEAM_anchor called properly
    warning("Progressive learning not initialized properly - using default settings")
    STEAM.obj <- initializeProgressiveLearning(STEAM.obj)
  }
  
  all_corrections <- data.frame()
  fold_summaries <- list()
  
  # Process each fold
  for (fold_idx in seq_len(n_folds)) {
    if (verbose) cat(sprintf("Processing Fold %d/%d\n", fold_idx, n_folds))
    cat(sprintf("=====================================\n\n"))
    
    fold_data <- outer_result[[fold_idx]]
    
    # Create fold-specific STEAM object
    fold_STEAM <- prepareFoldForProgressive(STEAM.obj, fold_data, fold_idx)
    
    # Run progressive learning on this fold
    fold_result <- runProgressiveLearningSingleFold(fold_STEAM, correction_function, verbose)
    
    # Extract corrections and assign fold number
    if (!is.null(fold_result$spatial_anchor_analysis$corrections)) {
      fold_corrections <- fold_result$spatial_anchor_analysis$corrections
      # Only assign fold if there are corrections
      if (nrow(fold_corrections) > 0) {
        fold_corrections$fold <- fold_idx
        all_corrections <- rbind(all_corrections, fold_corrections)
      }
    }
    
    # Store fold summary
    fold_summaries[[paste0("fold_", fold_idx)]] <- list(
      fold = fold_idx,
      corrections = if(exists("fold_corrections")) nrow(fold_corrections) else 0,
      success_rate = if(exists("fold_corrections") && nrow(fold_corrections) > 0) {
        # Handle both column name variations
        if ("correct" %in% colnames(fold_corrections)) {
          mean(fold_corrections$correct, na.rm = TRUE)
        } else if ("is_correct" %in% colnames(fold_corrections)) {
          mean(fold_corrections$is_correct, na.rm = TRUE)
        } else {
          1.0  # Default assumption
        }
      } else 0,
      performance = fold_result$progressive_learning
    )
    
    if (verbose) cat(sprintf("Fold %d completed\n\n", fold_idx))
  }
  
  # Aggregate results
  STEAM.obj <- aggregateMultiFoldResults(STEAM.obj, all_corrections, fold_summaries, verbose)
  
  if (verbose) {
    cat("=== Multi-Fold Progressive Learning Completed ===\n")
    cat(sprintf("Total corrections across all folds: %d\n", nrow(all_corrections)))
    if (nrow(all_corrections) > 0) {
      # Handle both column name variations
      success_col <- if ("correct" %in% colnames(all_corrections)) {
        all_corrections$correct
      } else if ("is_correct" %in% colnames(all_corrections)) {
        all_corrections$is_correct
      } else {
        rep(TRUE, nrow(all_corrections))  # Default assumption
      }
      overall_success <- mean(success_col, na.rm = TRUE)
      cat(sprintf("Overall success rate: %.1f%%\n", overall_success * 100))
      
      # Fold breakdown
      # Handle both column name variations for aggregation
      success_col <- if ("correct" %in% colnames(all_corrections)) {
        all_corrections$correct
      } else if ("is_correct" %in% colnames(all_corrections)) {
        all_corrections$is_correct
      } else {
        rep(TRUE, nrow(all_corrections))
      }
      fold_breakdown <- aggregate(success_col, 
                                  by = list(fold = all_corrections$fold), 
                                  FUN = function(x) c(count = length(x), success_rate = mean(x, na.rm = TRUE)))
      cat("\nFold-by-fold results:\n")
      for (i in seq_len(nrow(fold_breakdown))) {
        results <- fold_breakdown$x[i,]
        cat(sprintf("  Fold %d: %d corrections, %.1f%% success\n", 
                    fold_breakdown$fold[i], results[1], results[2] * 100))
      }
    }
  }
  
  # Aggregate and store all corrections in the main STEAM.obj
  STEAM.obj <- aggregateMultiFoldResults(STEAM.obj, all_corrections, fold_summaries, verbose = verbose)
  
  return(STEAM.obj)
}




#' Execute Single Progressive Learning Iteration with Comprehensive Monitoring
#' 
#' Performs one complete iteration of progressive learning including correction execution,
#' performance tracking, health monitoring, and adaptive threshold adjustment.
#' This function provides the core iteration logic with extensive analytics.
#' 
#' @param STEAM.obj STEAM object containing progressive learning framework and spatial data
#' @param correction_function Function to execute corrections for this iteration
#' @param verbose Logical indicating whether to print detailed progress information
#' 
#' @return List containing continuation flag, updated STEAM object, and iteration performance metrics
#' 
#' @details
#' This function implements comprehensive single-iteration processing:
#' - Validates progressive learning framework initialization
#' - Executes corrections using the provided correction function
#' - Records detailed performance metrics including timing and success rates
#' - Calculates improvement over previous iterations with trend analysis
#' - Updates learning memory and correction tracking systems
#' - Performs health checks and diagnostic analysis for poor performance
#' - Implements early stopping mechanisms based on convergence criteria
#' - Calculates adaptive thresholds for subsequent iterations
#' - Provides extensive progress analytics including coverage rates and convergence indicators
#' Returns comprehensive iteration results for continuation decision making.

executeProgressiveIteration <- function(STEAM.obj, correction_function, verbose = TRUE) {
  
  if (is.null(STEAM.obj$spatial_anchor_analysis$progressive_learning)) {
    stop("Progressive learning framework not initialized. Run initializeProgressiveLearning() first.")
  }
  
  pl <- STEAM.obj$spatial_anchor_analysis$progressive_learning
  iteration <- pl$current_iteration + 1
  
  if (verbose) {
    cat(sprintf("\n=== Progressive Learning Iteration %d ===\n", iteration))
    cat(sprintf("Current threshold: %.3f\n", pl$confidence_threshold))
  }
  
  # Check if we've reached max iterations
  if (iteration > pl$max_iterations) {
    cat("Maximum iterations reached. Stopping.\n")
    return(list(continue = FALSE, reason = "Max iterations reached"))
  }
  
  # Execute corrections with current threshold
  start_time <- Sys.time()
  correction_results <- correction_function(STEAM.obj, pl$confidence_threshold)
  end_time <- Sys.time()
  
  # Record iteration results
  n_attempted <- nrow(correction_results)
  n_successful <- if (n_attempted > 0) sum(correction_results$correct, na.rm = TRUE) else 0
  success_rate <- if (n_attempted > 0) mean(correction_results$correct, na.rm = TRUE) else 0
  
  iteration_performance <- data.frame(
    iteration = iteration,
    corrections_attempted = n_attempted,
    corrections_successful = n_successful,
    success_rate = success_rate,
    threshold_used = pl$confidence_threshold,
    improvement = 0,  # Will be calculated after comparison
    execution_time = as.numeric(difftime(end_time, start_time, units = "secs")),
    stringsAsFactors = FALSE
  )
  
  # Store actual corrections data for integration later
  iteration_performance$corrections_data <- list(correction_results)
  
  # Calculate improvement
  if (nrow(pl$performance_history) > 0) {
    last_success_rate <- tail(pl$performance_history$success_rate, 1)
    # Handle potential NaN values
    if (is.na(last_success_rate) || is.nan(last_success_rate)) {
      last_success_rate <- 0
    }
    iteration_performance$improvement <- iteration_performance$success_rate - last_success_rate
  }
  
  if (verbose) {
    cat(sprintf("Corrections attempted: %d\n", iteration_performance$corrections_attempted))
    cat(sprintf("Corrections successful: %d\n", iteration_performance$corrections_successful))
    cat(sprintf("Success rate: %.1f%%\n", iteration_performance$success_rate * 100))
    if (iteration_performance$improvement != 0) {
      cat(sprintf("Improvement: %+.1f%%\n", iteration_performance$improvement * 100))
    }
    
    # Enhanced Progress Analytics
    if (iteration > 1) {
      # Calculate trend over last few iterations
      recent_performance <- tail(STEAM.obj$spatial_anchor_analysis$progressive_learning$performance_history, min(3, iteration-1))
      if (nrow(recent_performance) > 1) {
        trend <- lm(success_rate ~ iteration, data = recent_performance)
        trend_slope <- coef(trend)[2]
        if (!is.na(trend_slope)) {
          if (trend_slope > 0.01) {
            cat("Trend: Improving\n")
          } else if (trend_slope < -0.01) {
            cat("Trend: Declining\n")  
          } else {
            cat("Trend: Stable\n")
          }
        }
      }
      
      # Show convergence indicators
      if (iteration > 2) {
        total_corrected <- length(STEAM.obj$spatial_anchor_analysis$progressive_learning$corrected_cells)
        coverage_rate <- total_corrected / nrow(STEAM.obj$spatial)
        cat(sprintf("Coverage: %.1f%% of cells processed\n", coverage_rate * 100))
      }
    }
  }
  
  # Update learning components
  STEAM.obj$spatial_anchor_analysis$progressive_learning$current_iteration <- iteration
  STEAM.obj$spatial_anchor_analysis$progressive_learning$performance_history <- rbind(
    pl$performance_history, iteration_performance
  )
  
  # Update learning memory
  STEAM.obj <- updateLearningMemory(STEAM.obj, correction_results)
  
  # Enhanced Success Tracking
  STEAM.obj <- trackCorrectionSuccess(STEAM.obj, iteration, verbose = verbose)
  
  # Store iteration details
  STEAM.obj$spatial_anchor_analysis$progressive_learning$iteration_history[[paste0("iteration_", iteration)]] <- list(
    performance = iteration_performance,
    corrections = correction_results,
    timestamp = Sys.time()
  )
  
  # Update corrected cells list
  new_corrected_cells <- unique(correction_results$cell_id)
  STEAM.obj$spatial_anchor_analysis$progressive_learning$corrected_cells <- unique(c(
    pl$corrected_cells, new_corrected_cells
  ))
  
  # Enhanced Analytics: Check for issues and suggest fixes if performance is poor
  if (iteration_performance$success_rate < 0.5 && verbose) {
    cat("\n--- Performance Analysis ---\n")
    diagnoseIterationIssues(STEAM.obj, verbose = TRUE)
    cat("\n--- Suggested Improvements ---\n")
    suggestIterationFixes(STEAM.obj, verbose = TRUE)
  }
  
  # Progressive Learning Health Check
  if (iteration > 2 && verbose) {
    health_check <- checkProgressiveLearning(STEAM.obj, verbose = FALSE)
    # Validate health_check structure and has_issues field
    if (!is.null(health_check) && is.list(health_check) && 
        "has_issues" %in% names(health_check) && 
        !is.null(health_check$has_issues) && 
        length(health_check$has_issues) > 0 && 
        !is.na(health_check$has_issues) && 
        health_check$has_issues) {
      cat(sprintf("\nProgressive Learning Health Check: %s\n", 
                  paste(health_check$issues, collapse = "; ")))
    }
  }
  
  # Advanced Progressive Learning Metrics
  if (!is.null(iteration) && length(iteration) > 0 && !is.na(iteration) && iteration > 1 && verbose) {
    metrics <- calculateProgressiveMetrics(STEAM.obj, iteration, verbose = TRUE)
    # Store metrics for future analysis
    STEAM.obj$spatial_anchor_analysis$progressive_learning$convergence_metrics <- metrics
  }
  
  # Check early stopping
  early_stop_result <- checkEarlyStopping(STEAM.obj, iteration_performance)
  
  if (verbose) cat(sprintf("Early stopping check: %s\n", early_stop_result$reason))
  
  if (early_stop_result$should_stop) {
    if (verbose) cat("Early stopping triggered.\n")
    return(list(continue = FALSE, reason = early_stop_result$reason, STEAM.obj = STEAM.obj))
  }
  
  # Calculate adaptive threshold for next iteration
  adaptive_conv_threshold <- calculateAdaptiveConvergenceThreshold(STEAM.obj)
  new_threshold <- calculateAdaptiveThreshold(STEAM.obj, iteration_performance)
  STEAM.obj$spatial_anchor_analysis$progressive_learning$confidence_threshold <- new_threshold
  STEAM.obj$spatial_anchor_analysis$progressive_learning$threshold_history <- c(
    pl$threshold_history, new_threshold
  )
  
  if (verbose) {
    cat(sprintf("Adaptive Convergence Threshold: %.4f\n", adaptive_conv_threshold))
    cat(sprintf("Next iteration threshold: %.3f\n", new_threshold))
  }
  
  return(list(continue = TRUE, STEAM.obj = STEAM.obj, performance = iteration_performance))
}



#' Prepare Fold-Specific STEAM Object for Progressive Learning
#' 
#' Creates a fold-specific STEAM object from cross-validation data, setting up
#' fold-specific predictions and initializing progressive learning framework
#' with inherited parent settings for consistent processing.
#' 
#' @param STEAM.obj Main STEAM object containing overall spatial data and analysis results
#' @param fold_data Cross-validation fold data containing test indices and fold-specific predictions
#' @param fold_idx Index number of the current fold being processed
#' 
#' @return Fold-specific STEAM object configured for progressive learning with proper fold context
#' 
#' @details
#' This function prepares individual folds for progressive learning processing:
#' - Creates a copy of the main STEAM object for fold-specific modifications
#' - Extracts fold-specific information including test cell indices
#' - Sets up fold-specific predictions from cross-validation results
#' - Stores fold metadata including fold number and cell mappings
#' - Initializes progressive learning framework with parent object settings
#' - Handles fallback scenarios when fold data structure varies
#' - Ensures consistency in progressive learning parameters across folds
#' Essential for maintaining fold isolation while preserving global learning settings.

prepareFoldForProgressive <- function(STEAM.obj, fold_data, fold_idx) {
  
  # Create a copy of the main STEAM object for this fold
  fold_STEAM <- STEAM.obj
  
  # Extract fold-specific information
  if ("test_index" %in% names(fold_data)) {
    test_cells <- fold_data$test_index
  } else {
    # If no test_index, assume all cells (single fold scenario)
    test_cells <- seq_len(nrow(STEAM.obj$spatial))
  }
  
  # Set up fold-specific predictions if available
  if ("preds" %in% names(fold_data)) {
    # Use fold-specific predictions
    fold_predictions <- fold_data$preds
  } else if ("predictions" %in% names(fold_data)) {
    fold_predictions <- fold_data$predictions
  } else {
    # Fallback to main predictions
    fold_predictions <- STEAM.obj$predictions
  }
  
  # Store fold information
  fold_STEAM$current_fold <- fold_idx
  fold_STEAM$fold_cells <- test_cells
  fold_STEAM$fold_predictions <- fold_predictions
  
  # Initialize progressive learning for this fold with parent settings
  if (is.null(fold_STEAM$spatial_anchor_analysis$progressive_learning)) {
    # Copy settings from parent STEAM.obj if available
    if (!is.null(STEAM.obj$spatial_anchor_analysis$progressive_learning)) {
      parent_pl <- STEAM.obj$spatial_anchor_analysis$progressive_learning
      fold_STEAM <- initializeProgressiveLearning(
        fold_STEAM,
        initial_threshold = parent_pl$initial_threshold,
        min_threshold = parent_pl$min_threshold,
        max_iterations = parent_pl$max_iterations
      )
    } else {
      fold_STEAM <- initializeProgressiveLearning(fold_STEAM)
    }
  }
  
  return(fold_STEAM)
}



#' Identify Remaining Misclassified Cells in Cross-Validation Fold
#' 
#' Identifies cells in a specific cross-validation fold that are misclassified
#' and have not yet been corrected, providing targets for progressive learning
#' iterations with explicit ground truth comparison.
#' 
#' @param STEAM.obj STEAM object containing nested cross-validation results
#' @param fold_idx Index of the cross-validation fold to analyze
#' @param already_corrected Vector of cell names that have already been corrected
#' @param ground_truth Named vector of true cell type labels for comparison
#' 
#' @return Character vector of cell names that are misclassified and available for correction
#' 
#' @details
#' This function performs fold-specific misclassification analysis:
#' - Extracts predictions from the specified cross-validation fold
#' - Matches ground truth labels to fold-specific cells
#' - Identifies discrepancies between predictions and ground truth
#' - Filters out cells that have already been corrected in previous iterations
#' - Returns remaining misclassified cells as targets for progressive learning
#' - Handles missing data and ensures robust cell name matching
#' Essential for maintaining accurate correction targeting in multi-fold progressive learning.

remainingCells <- function(STEAM.obj, fold_idx, already_corrected, ground_truth) {
  
  outer_result <- STEAM.obj$nested$ncv$outer_result[[fold_idx]]
  if (is.null(outer_result$preds)) return(character())
  
  preds_df <- outer_result$preds
  
  # Get cell names for this fold
  fold_cells <- rownames(preds_df)
  
  # Use explicit ground truth instead of testy column
  if (is.null(ground_truth)) {
    stop("Ground truth must be provided for misclassification identification")
  }
  
  # Match ground truth to fold cells
  fold_ground_truth <- ground_truth[fold_cells]
  fold_predictions <- preds_df$predy
  
  # Find misclassified cells in this fold
  misclassified_mask <- fold_ground_truth != fold_predictions
  misclassified_cells <- fold_cells[misclassified_mask & !is.na(misclassified_mask)]
  
  # Remove already corrected cells
  remaining_misclassified <- setdiff(misclassified_cells, already_corrected)
  
  return(remaining_misclassified)
}

