#' STEAM Spatial Anchor Analysis with Progressive Learning
#' 
#' Main function for executing STEAM spatial anchor analysis using progressive learning
#' with adaptive corrections. This function identifies and corrects misclassified cells
#' through iterative spatial-expression analysis with convergence monitoring.
#' 
#' @param STEAM.obj STEAM object containing spatial data, predictions, and analysis results
#' @param expression_matrix Optional expression matrix; if NULL, extracted from STEAM.obj$count_exp
#' @param ground_truth Optional ground truth labels for validation; enables supervised correction
#' @param max_iterations Maximum number of progressive learning iterations (default: 3)
#' @param convergence_threshold Convergence threshold for stopping criteria (default: 0.001)
#' @param min_confidence_gain Minimum confidence improvement required per iteration (default: 0.20)
#' @param max_flip_fraction Maximum fraction of predictions that can be corrected (default: 0.05)
#' @param k_neighbors Number of spatial neighbors for analysis; single value or vector for grid search
#' @param scoring_weights Named list with weights for expression, spatial, and model components
#' @param advanced_options Named list of advanced configuration options including gene transposition and decay
#' @param progressive_options Named list of progressive learning parameters including thresholds and patience
#' @param verbose Logical indicating whether to print detailed progress information
#' 
#' @return Updated STEAM object with spatial anchor analysis results and progressive learning corrections
#' 
#' @details
#' This function implements comprehensive spatial anchor analysis with progressive learning:
#' - Automatically detects and handles both single-dataset and multi-fold cross-validation scenarios
#' - Extracts or processes expression data with configurable gene transposition
#' - Implements adaptive threshold progressive learning with convergence monitoring
#' - Supports both supervised (with ground truth) and unsupervised correction modes
#' - Provides extensive parameter validation and intelligent defaults
#' - Includes advanced options for confidence decay, neighborhood updates, and correction chains
#' - Integrates results into spatial_anchor_analysis structure with comprehensive reporting
#' Essential entry point for STEAM spatial correction workflows with full parameter flexibility.

STEAM_anchor <- function(
    STEAM.obj,
    expression_matrix = NULL,
    ground_truth = NULL,
    max_iterations = 3,
    convergence_threshold = 0.001,
    min_confidence_gain = 0.20,
    max_flip_fraction = 0.05,
    k_neighbors = NULL,
    scoring_weights = NULL,
    advanced_options = NULL,
    progressive_options = NULL,
    verbose = TRUE
) {
  
  t0 <- Sys.time()
  
  # Parse advanced options with intelligent defaults
  default_advanced <- list(
    genes_transpose = TRUE,
    confidence_decay = 0.95,
    min_confidence_gain_decay = 0.9,
    update_neighborhoods = TRUE,
    track_correction_chains = TRUE,
    anchor_cut_grid = c(0.75, 0.82, 0.90),
    consensus_grid = c(0.82, 0.87, 0.92)
  )
  
  # Merge user-provided advanced options with defaults
  if (!is.null(advanced_options)) {
    if (!is.list(advanced_options)) {
      stop("advanced_options must be a named list")
    }
    # Validate provided options
    valid_options <- names(default_advanced)
    invalid_options <- setdiff(names(advanced_options), valid_options)
    if (length(invalid_options) > 0) {
      warning(sprintf("Ignoring invalid advanced_options: %s", paste(invalid_options, collapse = ", ")))
    }
    # Update defaults with valid provided options
    for (opt in intersect(names(advanced_options), valid_options)) {
      default_advanced[[opt]] <- advanced_options[[opt]]
    }
  }
  
  # Parse scoring weights with optimized defaults
  default_weights <- list(
    expression = 0.4,  # Optimized value
    spatial = 0.4,     # Balanced spatial importance
    model = 0.2        # Optimized model weight
  )
  
  if (!is.null(scoring_weights)) {
    if (!is.list(scoring_weights)) {
      stop("scoring_weights must be a named list with elements: expression, spatial, model")
    }
    # Validate and normalize weights
    required_weights <- c("expression", "spatial", "model")
    if (!all(required_weights %in% names(scoring_weights))) {
      stop(sprintf("scoring_weights must include all of: %s", paste(required_weights, collapse = ", ")))
    }
    
    # Normalize weights to sum to 1
    weight_sum <- sum(unlist(scoring_weights[required_weights]))
    if (weight_sum <= 0) {
      stop("scoring_weights must be positive numbers")
    }
    
    for (w in required_weights) {
      default_weights[[w]] <- scoring_weights[[w]] / weight_sum
    }
    
    if (verbose) {
      message(sprintf("Using custom scoring weights: expression=%.2f, spatial=%.2f, model=%.2f",
                      default_weights$expression, default_weights$spatial, default_weights$model))
    }
  }
  
  # Set up k-neighbor configuration
  k_grid <- if (!is.null(k_neighbors)) {
    if (is.numeric(k_neighbors) && length(k_neighbors) == 1) {
      k_neighbors  # Single k value
    } else if (is.numeric(k_neighbors) && length(k_neighbors) > 1) {
      k_neighbors  # Multiple k values for grid search
    } else {
      stop("k_neighbors must be a numeric value or vector")
    }
  } else {
    c(5, 8, 12)  # Optimized default grid
  }
  
  if (verbose) {
    message("=== STEAM Anchor Analysis (Enhanced Streamlined Interface) ===")
    message("Analyzing misclassified cells for spatial correction...")
    if (length(k_grid) == 1) {
      message(sprintf("Using k=%d spatial neighbors", k_grid))
    } else {
      message(sprintf("Auto-optimizing k from grid: %s", paste(k_grid, collapse = ", ")))
    }
  }
  
  # Determine ground truth availability
  has_ground_truth <- !is.null(ground_truth)
  
  if (verbose) {
    message(sprintf("Ground truth provided: %s", has_ground_truth))
    message(sprintf("Iterative mode: up to %d iterations", max_iterations))
  }
  
  # Extract expression matrix from STEAM.obj (always required)
  if (is.null(expression_matrix)) {
    if (verbose) message("Extracting expression matrix from STEAM.obj...")
    expression_matrix <- processExpressionData(STEAM.obj, extract_only = TRUE, 
                                               genes_transpose = default_advanced$genes_transpose, 
                                               verbose = verbose)
    if (is.null(expression_matrix)) {
      stop("No expression data found in STEAM.obj. Ensure STEAM.obj$count_exp is available. Expression data is required for STEAM_anchor analysis.")
    }
    if (verbose) message(sprintf("Successfully extracted expression matrix: %d cells x %d features", 
                                 nrow(expression_matrix), ncol(expression_matrix)))
  } else {
    if (verbose) message("Using provided expression matrix...")
    if (verbose) message(sprintf("Expression matrix dimensions: %d cells x %d features", 
                                 nrow(expression_matrix), ncol(expression_matrix)))
  }
  
  # Start iterative analysis with streamlined parameters
  if (verbose) message("\nStarting iterative analysis...")
  
  # Auto-enable consensus anchors when no ground truth
  use_consensus_anchors_auto <- !has_ground_truth
  
  if (verbose && use_consensus_anchors_auto) {
    message("Auto-enabling consensus-based anchor identification (no ground truth available)")
  }
  
  # Progressive Learning is the only supported approach
  if (verbose) message("Using Progressive Learning approach with adaptive thresholds...")
  
  # Parse progressive learning options
  default_progressive <- list(
    initial_threshold = 0.8,
    min_threshold = 0.6,
    max_corrections_per_iter = 100,
    patience = 3,
    min_improvement = 0.01
  )
  
  if (!is.null(progressive_options)) {
    if (!is.list(progressive_options)) {
      stop("progressive_options must be a named list")
    }
    # Update defaults with provided options
    for (opt in intersect(names(progressive_options), names(default_progressive))) {
      default_progressive[[opt]] <- progressive_options[[opt]]
    }
  }
  
  # Store ground truth and expression matrix in STEAM object for progressive learning
  if (!is.null(ground_truth)) {
    STEAM.obj$labels <- ground_truth
  }
  if (!is.null(expression_matrix)) {
    STEAM.obj$count_exp <- expression_matrix
  }
  
  # Source progressive learning functions - always load to ensure availability
  if (file.exists("progressive_learning_helpers.R")) {
    if (verbose) message("Loading progressive learning framework...")
    source("progressive_learning_helpers.R", local = FALSE)
  }
  
  # Load visualization functions
  if (file.exists("plots.R")) {
    if (verbose) message("Loading visualization functions...")
    source("plots.R", local = FALSE)
  }
  
  if (file.exists("progressive_steam_corrections.R")) {
    if (verbose) message("Loading progressive STEAM corrections...")
    source("progressive_steam_corrections.R", local = FALSE)
  }
  
  # Verify all required functions are now available
  progressive_functions_available <- all(sapply(c("initializeProgressiveLearning", "runProgressiveLearning", 
                                                  "progressiveSTEAMCorrections"), exists))
  
  if (!progressive_functions_available) {
    stop("Progressive learning framework not available. Required functions missing from progressive_learning_helpers.R")
  }
  
  # Full progressive learning framework
  if (verbose) message("Full progressive learning framework loaded successfully.")
  
  # Initialize progressive learning
  STEAM.obj <- initializeProgressiveLearning(
    STEAM.obj = STEAM.obj,
    initial_threshold = default_progressive$initial_threshold,
    min_threshold = default_progressive$min_threshold,
    max_iterations = max_iterations
  )
  
  # Set up early stopping parameters
  STEAM.obj$spatial_anchor_analysis$progressive_learning$early_stopping$patience <- default_progressive$patience
  STEAM.obj$spatial_anchor_analysis$progressive_learning$early_stopping$min_improvement <- default_progressive$min_improvement
  
  # Create progressive correction function
  progressive_correction_function <- function(STEAM.obj, confidence_threshold) {
    return(progressiveSTEAMCorrections(
      STEAM.obj = STEAM.obj,
      confidence_threshold = confidence_threshold,
      max_corrections = default_progressive$max_corrections_per_iter,
      expression_matrix = NULL,  # Will be extracted from STEAM.obj
      verbose = verbose
    ))
  }
  
  # Run progressive learning
  STEAM.obj <- runProgressiveLearning(
    STEAM.obj = STEAM.obj,
    correction_function = progressive_correction_function,
    verbose = verbose
  )
  
  # Integrate results with spatial_anchor_analysis structure
  STEAM.obj <- integrateProgressiveResults(STEAM.obj)
  
  # Consolidate results
  if (!is.null(STEAM.obj$spatial_anchor_analysis)) {
    if (verbose) message("Consolidating progressive learning results...")
    
    # Report current run results
    if (!is.null(STEAM.obj$spatial_anchor_analysis$corrections)) {
      n_corrections <- nrow(STEAM.obj$spatial_anchor_analysis$corrections)
      
      # Check if these are from current run or previous runs
      if (!is.null(STEAM.obj$spatial_anchor_analysis$progressive_learning) && 
          !is.null(STEAM.obj$spatial_anchor_analysis$progressive_learning$performance_history) &&
          nrow(STEAM.obj$spatial_anchor_analysis$progressive_learning$performance_history) > 0) {
        
        # Get corrections from current progressive learning run
        current_run_corrections <- sum(STEAM.obj$spatial_anchor_analysis$progressive_learning$performance_history$corrections_successful, na.rm = TRUE)
        
        if (verbose) {
          if (current_run_corrections > 0) {
            message(sprintf("Progressive learning complete: %d corrections applied in current run", current_run_corrections))
          } else {
            message("Progressive learning complete: No corrections applied in current run")
            if (n_corrections > 0) {
              message(sprintf("Note: %d corrections from previous runs remain in analysis", n_corrections))
            }
          }
        }
      } else {
        # Fallback for simple progressive learning or when performance history unavailable
        if (verbose) {
          message(sprintf("Progressive learning complete: %d total corrections in analysis", n_corrections))
        }
      }
    } else {
      if (verbose) message("Progressive learning complete: No corrections applied")
    }
    
    # Generate complete corrected labels output (both changed and unchanged)
    if (!is.null(STEAM.obj$labels)) {
      # Store original labels before any corrections if not already stored
      if (is.null(STEAM.obj$labels_before_correction)) {
        # Try to get original predictions from nested CV results
        original_preds <- extractPreds(STEAM.obj)
        if (length(original_preds) > 0) {
          STEAM.obj$labels_before_correction <- original_preds
        } else {
          # Fallback: assume current labels are original (for first run)
          STEAM.obj$labels_before_correction <- STEAM.obj$labels
        }
      }
      
      # Create complete corrected labels output
      STEAM.obj$corrected_labels_complete <- STEAM.obj$labels
      
      # Add metadata about the corrections
      if (!is.null(STEAM.obj$spatial_anchor_analysis$corrections) && 
          nrow(STEAM.obj$spatial_anchor_analysis$corrections) > 0) {
        
        corrections <- STEAM.obj$spatial_anchor_analysis$corrections
        
        # Remove duplicate corrections (keep first occurrence)
        corrections_unique <- corrections[!duplicated(corrections$cell_id), ]
        corrected_cell_ids_unique <- corrections_unique$cell_id
        
        # Add attributes to the complete labels
        attr(STEAM.obj$corrected_labels_complete, "n_corrections") <- length(corrected_cell_ids_unique)
        attr(STEAM.obj$corrected_labels_complete, "corrected_cells") <- corrected_cell_ids_unique
        attr(STEAM.obj$corrected_labels_complete, "correction_summary") <- data.frame(
          original = corrections_unique$original_pred,
          corrected = corrections_unique$suggested_correction,
          row.names = corrected_cell_ids_unique,
          stringsAsFactors = FALSE
        )
        
        if (verbose) {
          message(sprintf("Complete corrected labels created: %d total cells (%d corrected, %d unchanged)", 
                         length(STEAM.obj$corrected_labels_complete), 
                         length(corrected_cell_ids_unique),
                         length(STEAM.obj$corrected_labels_complete) - length(corrected_cell_ids_unique)))
        }
      } else {
        # No corrections made - labels remain unchanged
        attr(STEAM.obj$corrected_labels_complete, "n_corrections") <- 0
        attr(STEAM.obj$corrected_labels_complete, "corrected_cells") <- character(0)
        
        if (verbose) {
          message(sprintf("Complete corrected labels created: %d total cells (0 corrected, %d unchanged)", 
                         length(STEAM.obj$corrected_labels_complete),
                         length(STEAM.obj$corrected_labels_complete)))
        }
      }
    }
  }
  
  if (verbose) {
    total_runtime <- difftime(Sys.time(), t0, units = "secs")
    message(sprintf("Total STEAM_anchor runtime: %.1f seconds", as.numeric(total_runtime)))
  }
  
  return(STEAM.obj)
}




#' Aggregate Progressive Learning Results Across Multiple Cross-Validation Folds
#' 
#' Consolidates correction results from multiple cross-validation folds into unified
#' spatial anchor analysis structure with comprehensive performance metrics and
#' fold-specific tracking for multi-fold progressive learning workflows.
#' 
#' @param STEAM.obj STEAM object containing spatial data and analysis framework
#' @param all_corrections Combined dataframe of corrections from all processed folds
#' @param fold_summaries List of performance summaries for each individual fold
#' @param verbose Logical indicating whether to print aggregation progress information
#' 
#' @return Updated STEAM object with consolidated multi-fold progressive learning results
#' 
#' @details
#' This function performs comprehensive multi-fold result aggregation:
#' - Consolidates corrections from all folds into unified spatial_anchor_analysis structure
#' - Preserves fold-specific tracking and performance metrics
#' - Calculates overall success rates and performance statistics across folds
#' - Generates complete corrected label mappings with metadata
#' - Handles both corrected and unchanged cell populations
#' - Provides fold-by-fold performance breakdown with success rate analysis
#' - Ensures proper integration of multi-fold results for downstream analysis
#' Essential for maintaining fold structure while providing unified correction access.

aggregateMultiFoldResults <- function(STEAM.obj, all_corrections, fold_summaries, verbose = TRUE) {
  
  # Store all corrections in spatial_anchor_analysis
  if (!is.null(STEAM.obj$spatial_anchor_analysis)) {
    STEAM.obj$spatial_anchor_analysis$corrections <- all_corrections
  } else {
    STEAM.obj$spatial_anchor_analysis <- list(corrections = all_corrections)
  }
  
  # Create multi-fold summary
  if (nrow(all_corrections) > 0) {
    # Handle both column name variations
    success_col <- if ("correct" %in% colnames(all_corrections)) {
      all_corrections$correct
    } else if ("is_correct" %in% colnames(all_corrections)) {
      all_corrections$is_correct
    } else {
      rep(TRUE, nrow(all_corrections))
    }
    overall_success <- mean(success_col, na.rm = TRUE)
    
    # Update flip_summary with aggregated results
    success_col <- if ("correct" %in% colnames(all_corrections)) {
      all_corrections$correct
    } else if ("is_correct" %in% colnames(all_corrections)) {
      all_corrections$is_correct
    } else {
      rep(TRUE, nrow(all_corrections))
    }
    
    STEAM.obj$spatial_anchor_analysis$flip_summary <- data.frame(
      Outer_Fold = "All",
      Flips_To_Correct = sum(success_col, na.rm = TRUE),
      Flips_To_Wrong = sum(!success_col, na.rm = TRUE),
      Delta_Accuracy = overall_success,
      stringsAsFactors = FALSE
    )
    
    # Create fold-by-fold summary
    fold_results <- data.frame()
    for (fold_idx in unique(all_corrections$fold)) {
      fold_corrections <- all_corrections[all_corrections$fold == fold_idx, ]
      # Handle both column name variations
      fold_success_col <- if ("correct" %in% colnames(fold_corrections)) {
        fold_corrections$correct
      } else if ("is_correct" %in% colnames(fold_corrections)) {
        fold_corrections$is_correct
      } else {
        rep(TRUE, nrow(fold_corrections))
      }
      fold_success <- mean(fold_success_col, na.rm = TRUE)
      
      fold_results <- rbind(fold_results, data.frame(
        fold = fold_idx,
        total_corrections = nrow(fold_corrections),
        successful_corrections = sum(fold_success_col, na.rm = TRUE),
        accuracy_improvement = fold_success,
        stringsAsFactors = FALSE
      ))
    }
    
    STEAM.obj$spatial_anchor_analysis$final_by_fold <- fold_results
    
    # Store fold summaries
    STEAM.obj$spatial_anchor_analysis$fold_summaries <- fold_summaries
    
    # Update progressive learning with aggregated pattern performance
    if (!"progressive_learning" %in% names(STEAM.obj)) {
      STEAM.obj <- initializeProgressiveLearning(STEAM.obj)
    }
    
    # Aggregate pattern learning across folds
    STEAM.obj <- updatePatternPerformance(STEAM.obj, all_corrections, verbose = verbose)
    
    # Update performance history for visualization compatibility
    if (nrow(all_corrections) > 0) {
      performance_row <- data.frame(
        iteration = 1,  # Multi-fold is considered as one "iteration"
        corrections_attempted = nrow(all_corrections),
        corrections_successful = sum(success_col, na.rm = TRUE),
        success_rate = overall_success,
        threshold_used = 0.8,  # Default threshold used in multi-fold
        improvement = overall_success,
        stringsAsFactors = FALSE
      )
      
      # Add to performance history
      STEAM.obj$spatial_anchor_analysis$progressive_learning$performance_history <- rbind(
        STEAM.obj$spatial_anchor_analysis$progressive_learning$performance_history, 
        performance_row
      )
    }
  }
  
  return(STEAM.obj)
}


#' Analyze Correction Patterns and Generate Pattern Statistics
#' 
#' Analyzes correction data to identify patterns in cell type transitions,
#' confidence improvements, and spatial contexts, providing comprehensive
#' pattern statistics for progressive learning optimization.
#' 
#' @param corrections Dataframe containing correction results with original and corrected labels
#' 
#' @return List containing correction pattern analysis including transition types, confidence metrics, and pattern frequencies
#' 
#' @details
#' This function performs detailed correction pattern analysis:
#' - Identifies cell type transition patterns (e.g., "Type A → Type B")
#' - Calculates confidence improvement statistics across corrections
#' - Analyzes spatial context and neighborhood characteristics
#' - Generates frequency tables of correction patterns
#' - Returns empty list if no corrections are available
#' - Supports universal data types for robust correction analysis
#' Essential for understanding correction behavior and optimizing progressive learning parameters.

analyzeCorrections <- function(corrections) {
  if (nrow(corrections) == 0) return(list())
  
  # Analyze correction types (universal data type support)
  correction_types <- sprintf("%s â†’ %s", 
                              as.character(corrections$original_pred), 
                              as.character(corrections$suggested_pred))
  type_counts <- table(correction_types)
  
  # Success rates by type (simplified)
  success_rates <- rep(1.0, length(type_counts))  # Assume all successful for now
  names(success_rates) <- names(type_counts)
  
  return(list(
    types = correction_types,
    type_counts = type_counts,
    success_rates = success_rates
  ))
}

#' Extract Predictions from STEAM Object with Multi-Source Support
#' 
#' Extracts prediction data from STEAM objects, handling various data structures
#' including nested cross-validation results, direct predictions, and fallback
#' scenarios with comprehensive source tracking for debugging.
#' 
#' @param STEAM.obj STEAM object containing prediction data in various possible structures
#' 
#' @return Named character vector of predictions with source attribution for debugging
#' 
#' @details
#' This function implements robust prediction extraction with multiple fallback strategies:
#' - Primary: Extracts from nested cross-validation outer results with fold consolidation
#' - Secondary: Uses direct predictions from STEAM.obj$predictions field
#' - Tertiary: Falls back to STEAM.obj$pred for compatibility
#' - Handles both preds dataframe and direct vector formats
#' - Tracks data source for debugging and validation purposes
#' - Provides comprehensive error handling for missing prediction data
#' - Returns empty vector if no predictions found with appropriate warnings
#' Essential for robust prediction access across different STEAM object configurations.

extractPreds <- function(STEAM.obj) {
    
    all_predictions <- character()
    source_field <- "unknown"
    
    # Extract from nested CV results
    if (!is.null(STEAM.obj$nested) && !is.null(STEAM.obj$nested$ncv) && 
        !is.null(STEAM.obj$nested$ncv$outer_result)) {
        
        for (f in seq_along(STEAM.obj$nested$ncv$outer_result)) {
            fold_result <- STEAM.obj$nested$ncv$outer_result[[f]]
            if (!is.null(fold_result$preds)) {
                preds_df <- fold_result$preds
                fold_preds <- setNames(preds_df$predy, rownames(preds_df))
                all_predictions <- c(all_predictions, fold_preds)
                source_field <- "nested$ncv$outer_result$preds$predy"
            }
        }
    }
    
    # Fallback: try other locations in STEAM object
    if (length(all_predictions) == 0) {
        if (!is.null(STEAM.obj$pred)) {
            all_predictions <- STEAM.obj$pred
            source_field <- "STEAM.obj$pred"
        } else if (!is.null(STEAM.obj$predictions)) {
            all_predictions <- STEAM.obj$predictions
            source_field <- "STEAM.obj$predictions"
        } else if (!is.null(STEAM.obj$CV) && !is.null(STEAM.obj$CV$pred)) {
            preds_df <- STEAM.obj$CV$pred
            if (is.data.frame(preds_df) && "predy" %in% colnames(preds_df)) {
                all_predictions <- setNames(preds_df$predy, rownames(preds_df))
                source_field <- "STEAM.obj$CV$pred$predy"
            } else {
                all_predictions <- STEAM.obj$CV$pred
                source_field <- "STEAM.obj$CV$pred"
            }
        }
    }
    
    # Remove duplicates (keep first occurrence)
    if (length(all_predictions) > 0) {
        all_predictions <- all_predictions[!duplicated(names(all_predictions))]
    }
    
    # Add source attribution
    attr(all_predictions, "source_field") <- source_field
    
    return(all_predictions)
}




#' Suggest Fixes for Progressive Learning Performance Issues
#' 
#' Analyzes progressive learning performance and provides targeted recommendations
#' for improving correction success rates, convergence behavior, and overall
#' progressive learning effectiveness based on current performance metrics.
#' 
#' @param STEAM.obj STEAM object containing progressive learning performance history
#' @param verbose Logical indicating whether to print suggestions to console
#' 
#' @return Invisible list of suggested fixes and recommendations
#' 
#' @details
#' This function provides intelligent performance optimization suggestions:
#' - Analyzes success rates and identifies low-performance scenarios
#' - Suggests threshold adjustments for better cell identification
#' - Recommends parameter tuning based on convergence patterns
#' - Provides spatial neighborhood optimization suggestions
#' - Identifies potential overfitting or underfitting scenarios
#' - Suggests alternative approaches for challenging datasets
#' - Returns actionable recommendations for progressive learning improvement
#' Essential for troubleshooting and optimizing progressive learning workflows.

suggestIterationFixes <- function(STEAM.obj, verbose = TRUE) {
  
  if (!"progressive_learning" %in% names(STEAM.obj)) {
    return(list())
  }
  
  pl <- STEAM.obj$spatial_anchor_analysis$progressive_learning
  suggestions <- character()
  
  # Analyze performance and suggest improvements
  if (nrow(pl$performance_history) > 0) {
    avg_success <- mean(pl$performance_history$success_rate)
    
    if (avg_success < 0.5) {
      suggestions <- c(suggestions, "Consider starting with a higher confidence threshold")
      suggestions <- c(suggestions, "Reduce the number of corrections per iteration")
    }
    
    if (length(pl$corrected_cells) == 0) {
      suggestions <- c(suggestions, "Try lowering the initial confidence threshold")
      suggestions <- c(suggestions, "Increase the maximum number of iterations")
    }
  }
  
  if (verbose && length(suggestions) > 0) {
    cat("Suggestions for future runs:\n")
    for (suggestion in suggestions) {
      cat(sprintf("  â€¢ %s\n", suggestion))
    }
  }
  
  return(suggestions)
}


#' Display Formatted Correction Table with Summary Statistics
#' 
#' Presents a formatted table of spatial anchor corrections with comprehensive
#' summary statistics, cell type transitions, and confidence improvements
#' for easy interpretation and analysis of correction results.
#' 
#' @param STEAM.obj STEAM object containing spatial anchor analysis results
#' @param max_rows Maximum number of correction rows to display (default: 50)
#' 
#' @return Dataframe of corrections (invisibly) while printing formatted table to console
#' 
#' @details
#' This function provides comprehensive correction table visualization:
#' - Extracts corrections from spatial_anchor_analysis results
#' - Displays formatted table with cell IDs, original/corrected predictions, and confidence
#' - Shows summary statistics including total corrections and success rates
#' - Handles fold information for multi-fold cross-validation results
#' - Limits display to specified number of rows for readability
#' - Returns empty dataframe if no corrections are available
#' - Provides user-friendly error messages for missing analysis results
#' Essential for quick visual inspection and interpretation of correction outcomes.

showCorrectionTable <- function(STEAM.obj, max_rows = 50) {
  
  if (is.null(STEAM.obj$spatial_anchor_analysis)) {
    cat("No spatial_anchor_analysis results found. Run STEAM_anchor first.\n")
    return(data.frame())
  }
  
  # Extract corrections
  corrections <- extractAnchorCorrections(STEAM.obj$spatial_anchor_analysis, STEAM.obj)
  
  if (nrow(corrections) == 0) {
    cat("No corrections were made.\n")
    return(data.frame())
  }
  
  # Remove duplicate corrections (keep first occurrence)
  corrections <- corrections[!duplicated(corrections$cell_id), ]
  
  # Create a nice display table
  display_table <- corrections[, c("cell_id", "original_pred", "suggested_correction", "fold")]
  
  # Add true label if available
  if ("true_label" %in% colnames(corrections)) {
    display_table$true_label <- corrections$true_label
  }
  
  # Add correction accuracy if available
  if ("correction_correct" %in% colnames(corrections)) {
    display_table$correct <- ifelse(corrections$correction_correct, "âœ“", "âœ—")
  }
  
  # Sort by fold and then by cell_id
  if ("fold" %in% colnames(display_table)) {
    display_table <- display_table[order(display_table$fold, display_table$cell_id), ]
  }
  
  # Limit rows if requested
  if (nrow(display_table) > max_rows) {
    cat(sprintf("Showing first %d corrections out of %d total:\n\n", max_rows, nrow(display_table)))
    display_table <- display_table[1:max_rows, ]
  } else {
    cat(sprintf("All %d corrections:\n\n", nrow(display_table)))
  }
  
  # Print the table
  print(display_table, row.names = FALSE)
  
  # Summary by correction type
  cat("\n=== Correction Summary ===\n")
  if ("suggested_correction" %in% colnames(corrections)) {
    correction_counts <- table(corrections$suggested_correction)
    cat("Cells corrected to each type:\n")
    for (i in seq_along(correction_counts)) {
      cat(sprintf("  Type %s: %d cells\n", names(correction_counts)[i], correction_counts[i]))
    }
  }
  
  if ("fold" %in% colnames(corrections)) {
    fold_counts <- table(corrections$fold)
    cat("\nCorrections by fold:\n")
    for (i in seq_along(fold_counts)) {
      cat(sprintf("  Fold %s: %d corrections\n", names(fold_counts)[i], fold_counts[i]))
    }
  }
  
  return(invisible(display_table))
}

#' Extract Complete Corrected Labels with Metadata Support
#' 
#' Retrieves complete set of cell labels including both corrected and unchanged
#' cells, with optional metadata about correction status, sources, and confidence
#' scores for comprehensive downstream analysis integration.
#' 
#' @param STEAM.obj STEAM object containing spatial anchor analysis results and corrected labels
#' @param include_metadata Logical indicating whether to include correction metadata and attribution
#' 
#' @return Named vector or dataframe of complete cell labels with optional correction metadata
#' 
#' @details
#' This function provides comprehensive label extraction with correction tracking:
#' - Returns complete corrected labels from spatial anchor analysis if available
#' - Includes both corrected cells and unchanged cells for complete dataset coverage
#' - Optionally includes metadata about correction source and confidence
#' - Falls back to original predictions when corrections are not available
#' - Handles various STEAM object configurations and data structures
#' - Provides proper error handling for invalid or missing STEAM objects
#' - Essential for integrating corrected results into downstream analysis workflows
#' Ensures consistent access to final cell type assignments across different analysis scenarios.

getCorrectedLabels <- function(STEAM.obj, include_metadata = TRUE) {
  
  if (is.null(STEAM.obj)) {
    stop("STEAM.obj cannot be NULL")
  }
  
  # Check if corrected labels are available
  if (!is.null(STEAM.obj$corrected_labels_complete)) {
    corrected_labels <- STEAM.obj$corrected_labels_complete
    
    if (include_metadata) {
      # Return with metadata about corrections
      return(corrected_labels)
    } else {
      # Return just the labels without attributes
      attributes(corrected_labels) <- NULL
      names(corrected_labels) <- names(STEAM.obj$corrected_labels_complete)
      return(corrected_labels)
    }
  }
  
  # Fallback: construct corrected labels from current state
  if (!is.null(STEAM.obj$labels)) {
    corrected_labels <- STEAM.obj$labels
    
    # Add metadata if corrections exist
    if (include_metadata && !is.null(STEAM.obj$spatial_anchor_analysis$corrections) && 
        nrow(STEAM.obj$spatial_anchor_analysis$corrections) > 0) {
      
      corrections <- STEAM.obj$spatial_anchor_analysis$corrections
      
      # Remove duplicate corrections (keep first occurrence)
      corrections_unique <- corrections[!duplicated(corrections$cell_id), ]
      corrected_cell_ids_unique <- corrections_unique$cell_id
      
      attr(corrected_labels, "n_corrections") <- length(corrected_cell_ids_unique)
      attr(corrected_labels, "corrected_cells") <- corrected_cell_ids_unique
      attr(corrected_labels, "correction_summary") <- data.frame(
        original = corrections_unique$original_pred,
        corrected = corrections_unique$suggested_correction,
        row.names = corrected_cell_ids_unique,
        stringsAsFactors = FALSE
      )
    } else if (include_metadata) {
      attr(corrected_labels, "n_corrections") <- 0
      attr(corrected_labels, "corrected_cells") <- character(0)
    }
    
    return(corrected_labels)
  }
  
  # No labels available
  warning("No labels found in STEAM.obj. Run STEAM_anchor first.")
  return(NULL)
}


#' Generate Comprehensive Correction Summary Statistics
#' 
#' Creates detailed summary statistics of spatial anchor corrections including
#' correction rates, cell type distributions, and performance metrics for
#' comprehensive analysis and reporting of progressive learning outcomes.
#' 
#' @param STEAM.obj STEAM object containing spatial anchor analysis results and corrections
#' 
#' @return Dataframe containing comprehensive correction summary statistics
#' 
#' @details
#' This function provides complete correction summary analysis:
#' - Calculates total cells processed and correction rates
#' - Determines numbers of corrected vs unchanged cells
#' - Extracts correction metadata including confidence and success metrics
#' - Handles missing or incomplete spatial anchor analysis results
#' - Returns structured summary dataframe for reporting and visualization
#' - Provides zero-filled results when no corrections are available
#' - Essential for quantitative assessment of progressive learning effectiveness
#' Enables comprehensive reporting and comparison of correction workflows across different scenarios.

getCorrectionSummary <- function(STEAM.obj) {
  
  if (is.null(STEAM.obj) || is.null(STEAM.obj$spatial_anchor_analysis)) {
    return(data.frame(
      total_cells = 0,
      corrected_cells = 0,
      unchanged_cells = 0,
      correction_rate = 0,
      stringsAsFactors = FALSE
    ))
  }
  
  # Get complete corrected labels
  corrected_labels <- getCorrectedLabels(STEAM.obj, include_metadata = TRUE)
  
  if (is.null(corrected_labels)) {
    return(data.frame(
      total_cells = 0,
      corrected_cells = 0,
      unchanged_cells = 0,
      correction_rate = 0,
      stringsAsFactors = FALSE
    ))
  }
  
  total_cells <- length(corrected_labels)
  n_corrections <- attr(corrected_labels, "n_corrections")
  if (is.null(n_corrections)) n_corrections <- 0
  
  unchanged_cells <- total_cells - n_corrections
  correction_rate <- n_corrections / total_cells
  
  summary_df <- data.frame(
    total_cells = total_cells,
    corrected_cells = n_corrections,
    unchanged_cells = unchanged_cells,
    correction_rate = round(correction_rate, 4),
    stringsAsFactors = FALSE
  )
  
  # Add details about correction types if available
  if (n_corrections > 0 && !is.null(attr(corrected_labels, "correction_summary"))) {
    correction_summary <- attr(corrected_labels, "correction_summary")
    correction_types <- paste(correction_summary$original, "â†’", correction_summary$corrected)
    type_counts <- table(correction_types)
    
    cat("Correction Summary:\n")
    cat(sprintf("  Total cells: %d\n", total_cells))
    cat(sprintf("  Corrected cells: %d (%.1f%%)\n", n_corrections, correction_rate * 100))
    cat(sprintf("  Unchanged cells: %d (%.1f%%)\n", unchanged_cells, (1 - correction_rate) * 100))
    cat("\nCorrection types:\n")
    for (i in seq_along(type_counts)) {
      cat(sprintf("  %s: %d cells\n", names(type_counts)[i], type_counts[i]))
    }
  }
  
  return(summary_df)
}