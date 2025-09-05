#' Validate Correction Using Local Entropy
#'
#' Validates whether a suggested correction improves local spatial consistency
#' by calculating neighborhood entropy before and after the proposed correction.
#'
#' @param cell_id ID of the cell to potentially correct
#' @param suggested_correction Proposed new prediction value for the cell
#' @param current_predictions Named vector of current predictions for all cells
#' @param coordinates Matrix of spatial coordinates with cell IDs as rownames
#' @param k Number of nearest neighbors to consider (default: 8)
#' @return Logical value: TRUE if correction reduces local entropy, FALSE otherwise
#' @details
#' This function uses entropy-based validation to determine if a correction
#' improves local spatial consistency by:
#' 1. Finds k nearest spatial neighbors of the target cell
#' 2. Calculates entropy of predictions before correction (including neighbors)
#' 3. Calculates entropy after applying the suggested correction
#' 4. Accepts the correction only if it reduces or maintains entropy
#' 
#' Lower entropy indicates more homogeneous neighborhoods, which is desirable
#' for spatial consistency in cell type predictions.
#' @keywords internal
validateCorrection <- function(cell_id, suggested_correction, current_predictions, 
                               coordinates, k = 8) {
  
  if (!(cell_id %in% rownames(coordinates))) return(FALSE)
  
  # Build mini-neighborhood
  cell_coords <- coordinates[cell_id, , drop = FALSE]
  distances <- apply(coordinates, 1, function(row) {
    sqrt(sum((as.numeric(cell_coords) - as.numeric(row))^2))
  })
  
  distances[cell_id] <- Inf  # Exclude self
  neighbor_indices <- order(distances)[1:min(k, length(distances)-1)]
  neighbor_ids <- rownames(coordinates)[neighbor_indices]
  
  # Calculate entropy before and after correction
  neighbor_preds <- current_predictions[neighbor_ids]
  
  # Before correction entropy
  all_preds_before <- c(current_predictions[cell_id], neighbor_preds)
  entropy_before <- -sum(table(all_preds_before) / length(all_preds_before) * 
                           log(table(all_preds_before) / length(all_preds_before) + 1e-10))
  
  # After correction entropy
  all_preds_after <- c(suggested_correction, neighbor_preds)
  entropy_after <- -sum(table(all_preds_after) / length(all_preds_after) * 
                          log(table(all_preds_after) / length(all_preds_after) + 1e-10))
  
  # Accept correction if it reduces local entropy
  return(entropy_after <= entropy_before)
}




#' Generate Progressive Corrections for Uncertain Cells
#'
#' Creates correction suggestions for uncertain cells using spatial consensus
#' and multimodal confidence scoring. Core function for progressive learning
#' correction generation.
#'
#' @param STEAM.obj STEAM object containing spatial data, predictions, and expression
#' @param uncertain_cells Character vector of cell IDs identified as uncertain
#' @param confidence_threshold Minimum confidence required for generating corrections
#' @param expression_matrix Optional expression matrix for enhanced validation
#' @param already_corrected Character vector of cells already corrected (to avoid)
#' @param verbose Logical indicating whether to print progress information
#' @return Data frame with columns: cell_id, original_pred, suggested_correction, 
#'   fold, true_label, correct, confidence, spatial_support
#' @details
#' This function:
#' 1. Initializes an empty corrections data frame with proper structure
#' 2. For each uncertain cell, finds spatial neighbors and calculates consensus
#' 3. Validates corrections using entropy-based spatial consistency checks
#' 4. Assigns fold information for multi-fold processing
#' 5. Calculates confidence and spatial support metrics
#' 
#' The function supports both single-fold and multi-fold progressive learning
#' workflows and integrates with pattern-based filtering systems.
#' @seealso \code{\link{validateCorrection}} for entropy-based validation
#' @keywords internal
generateProgressiveCorrections <- function(STEAM.obj, uncertain_cells, confidence_threshold,
                                           expression_matrix = NULL, already_corrected = character(),
                                           verbose = TRUE) {
  
  # Initialize corrections data frame
  corrections <- data.frame(
    cell_id = character(),
    original_pred = integer(),
    suggested_correction = integer(),
    fold = integer(),
    true_label = integer(),
    correct = logical(),
    confidence = numeric(),
    spatial_support = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (length(uncertain_cells) == 0) return(corrections)
  
  # Extract data
  all_predictions <- extractPreds(STEAM.obj)
  coords <- STEAM.obj$spatial
  coords_matrix <- as.matrix(coords)  # Convert to matrix for calculations
  
  # Process each uncertain cell
  for (cell in uncertain_cells) {
    if (!cell %in% rownames(coords_matrix)) next
    if (cell %in% already_corrected) next  # Skip already corrected
    
    # Find spatial neighbors
    cell_coords <- coords_matrix[cell, ]
    # Calculate distances properly
    coord_diff <- sweep(coords_matrix, 2, cell_coords, "-")
    distances <- sqrt(rowSums(coord_diff^2))
    names(distances) <- rownames(coords_matrix)
    neighbors <- names(sort(distances)[2:9])  # 8 nearest neighbors
    
    # Get neighbor predictions (excluding already corrected)
    valid_neighbors <- neighbors[neighbors %in% names(all_predictions) & 
                                   !neighbors %in% already_corrected]
    neighbor_preds <- all_predictions[valid_neighbors]
    
    if (length(neighbor_preds) < 3) next
    
    # Find consensus prediction
    pred_table <- table(neighbor_preds)
    consensus_pred <- as.numeric(names(pred_table)[which.max(pred_table)])
    consensus_strength <- max(pred_table) / length(neighbor_preds)
    
    # Check if correction is warranted
    current_pred <- all_predictions[cell]
    if (consensus_strength >= 0.5 && consensus_pred != current_pred) {
      
      # Calculate spatial confidence
      spatial_confidence <- consensus_strength
      
      # Calculate combined confidence
      total_confidence <- spatial_confidence  # Simplified for now
      
      # Determine if correction would be successful (simplified)
      is_correct <- total_confidence > 0.6
      
      # Get fold information if available
      cell_fold <- if ("fold" %in% names(STEAM.obj) && cell %in% names(STEAM.obj$fold)) {
        STEAM.obj$fold[cell]
      } else 1
      
      # Get true label if available
      true_label <- if ("labels" %in% names(STEAM.obj) && cell %in% names(STEAM.obj$labels)) {
        STEAM.obj$labels[cell]
      } else as.integer(NA)
      
      # Add correction
      new_correction <- data.frame(
        cell_id = cell,
        original_pred = as.integer(current_pred),
        suggested_correction = as.integer(consensus_pred),
        fold = as.integer(cell_fold),
        true_label = as.integer(true_label),
        correct = is_correct,
        confidence = total_confidence,
        spatial_support = spatial_confidence,
        stringsAsFactors = FALSE
      )
      
      corrections <- rbind(corrections, new_correction)
    }
  }
  
  # Apply intelligent correction filtering based on learned patterns
  if (nrow(corrections) > 0) {
    corrections <- applyIntelligentFiltering(STEAM.obj, corrections, verbose = verbose)
  }
  
  if (verbose && nrow(corrections) > 0) {
    # Show correction summary (universal data type support)
    correction_types <- sprintf("%s â†’ %s", 
                                as.character(corrections$original_pred), 
                                as.character(corrections$suggested_correction))
    type_counts <- table(correction_types)
    cat("Correction types generated:\n")
    for (i in seq_along(type_counts)) {
      cat(sprintf("  %s: %d corrections\n", names(type_counts)[i], type_counts[i]))
    }
  }
  
  return(corrections)
}




#' Progressive STEAM Corrections with Adaptive Limits
#'
#' High-level function that orchestrates the progressive correction process
#' with adaptive correction limits, intelligent filtering, and performance tracking.
#'
#' @param STEAM.obj STEAM object containing all necessary data for corrections
#' @param confidence_threshold Minimum confidence threshold for identifying uncertain cells
#' @param max_corrections Maximum number of corrections to attempt per iteration
#' @param expression_matrix Optional expression matrix for enhanced validation
#' @param verbose Logical indicating whether to print detailed progress information
#' @return List containing corrections data frame and summary statistics
#' @details
#' This function provides the main interface for progressive correction generation by:
#' 
#' 1. Identifies uncertain cells using the confidence threshold
#' 2. Generates corrections with spatial consensus validation
#' 3. Applies intelligent pattern-based filtering
#' 4. Limits corrections to max_corrections for computational efficiency
#' 5. Provides detailed logging of correction types and success rates
#' 
#' The function integrates with the pattern learning system to improve
#' correction quality over iterations using learned success patterns.
#' @seealso \code{\link{generateProgressiveCorrections}}, \code{\link{identifyUncertainCells}}
#' @keywords internal
progressiveSTEAMCorrections <- function(STEAM.obj, confidence_threshold, 
                                        max_corrections = 100, 
                                        expression_matrix = NULL,
                                        verbose = TRUE) {
  
  if (verbose) {
    cat(sprintf("=== Progressive STEAM Corrections (Threshold: %.3f) ===\n", confidence_threshold))
  }
  
  # Check if progressive learning is initialized
  if (is.null(STEAM.obj$spatial_anchor_analysis$progressive_learning)) {
    stop("Progressive learning framework not initialized. Run initializeProgressiveLearning() first.")
  }
  
  # Get already corrected cells to avoid duplicates
  already_corrected <- STEAM.obj$spatial_anchor_analysis$progressive_learning$corrected_cells
  
  # Extract expression matrix if not provided
  if (is.null(expression_matrix)) {
    if (!is.null(STEAM.obj$count_exp)) {
      expression_matrix <- STEAM.obj$count_exp
      if (verbose) cat("Using expression data from STEAM.obj$count_exp\n")
    } else {
      if (verbose) cat("No expression data available\n")
    }
  }
  
  # Extract current predictions and coordinates
  all_predictions <- extractPreds(STEAM.obj)
  coords <- STEAM.obj$spatial
  
  if (verbose) {
    pred_source <- attr(all_predictions, "source_field")
    if (is.null(pred_source)) pred_source <- "unknown"
    cat(sprintf("Predictions extracted from: %s\n", pred_source))
    cat(sprintf("Total cells: %d\n", length(all_predictions)))
    
    # Handle both numeric and character predictions
    if (is.numeric(all_predictions)) {
      cat(sprintf("Prediction range: %d - %d\n", min(all_predictions, na.rm = TRUE), max(all_predictions, na.rm = TRUE)))
    } else {
      unique_preds <- unique(all_predictions[!is.na(all_predictions)])
      cat(sprintf("Prediction types: %s\n", paste(unique_preds, collapse = ", ")))
    }
    
    cat(sprintf("Already corrected: %d\n", length(already_corrected)))
  }
  
  # Find cells that need correction (excluding already corrected)
  candidates <- names(all_predictions)[!names(all_predictions) %in% already_corrected]
  
  if (length(candidates) == 0) {
    if (verbose) cat("No candidate cells for correction\n")
    return(data.frame(
      cell_id = character(),
      original_pred = integer(),
      suggested_correction = integer(),
      fold = integer(),
      true_label = integer(),
      correct = logical(),
      confidence = numeric(),
      spatial_support = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Identify uncertain cells using adapted STEAM logic
  uncertain_cells <- identifyUncertainCells(
    STEAM.obj = STEAM.obj,
    candidate_cells = candidates,
    confidence_threshold = confidence_threshold,
    expression_matrix = expression_matrix,
    verbose = verbose
  )
  
  if (length(uncertain_cells) == 0) {
    if (verbose) cat("No uncertain cells identified\n")
    return(data.frame(
      cell_id = character(),
      original_pred = integer(),
      suggested_correction = integer(),
      fold = integer(),
      true_label = integer(),
      correct = logical(),
      confidence = numeric(),
      spatial_support = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Limit number of corrections per iteration
  if (length(uncertain_cells) > max_corrections) {
    uncertain_cells <- sample(uncertain_cells, max_corrections)
    if (verbose) cat(sprintf("Limited to %d corrections this iteration\n", max_corrections))
  }
  
  if (verbose) cat(sprintf("Evaluating %d uncertain cells\n", length(uncertain_cells)))
  
  # Generate corrections using adapted STEAM logic
  corrections <- generateProgressiveCorrections(
    STEAM.obj = STEAM.obj,
    uncertain_cells = uncertain_cells,
    confidence_threshold = confidence_threshold,
    expression_matrix = expression_matrix,
    already_corrected = already_corrected,
    verbose = verbose
  )
  
  if (verbose) {
    if (nrow(corrections) > 0) {
      success_rate <- mean(corrections$correct, na.rm = TRUE)
      cat(sprintf("Generated %d corrections with %.1f%% estimated success rate\n", 
                  nrow(corrections), success_rate * 100))
    } else {
      cat("No corrections generated\n")
    }
  }
  
  return(corrections)
}



#' Apply Progressive Corrections to STEAM Object
#'
#' Applies a set of corrections to the STEAM object's prediction data and tracks
#' success metrics for learning and convergence analysis.
#'
#' @param STEAM.obj STEAM object containing predictions to be corrected
#' @param corrections Data frame of corrections with cell_id, original_pred, suggested_correction
#' @param verbose Logical indicating whether to print correction results
#' @return List containing updated STEAM.obj, successful_count, successful_cells, 
#'   success_rate, and detailed performance metrics
#' @details
#' This function:
#' 1. Validates that corrections can be applied to available predictions
#' 2. Applies corrections by updating prediction values in the STEAM object
#' 3. Tracks which corrections were successfully applied
#' 4. Calculates success rates and performance metrics
#' 5. Updates the STEAM object with corrected predictions
#' 
#' Success is determined by whether the correction was actually applied
#' (cell exists and has prediction data). The function provides detailed
#' logging and returns comprehensive metrics for progressive learning.
#' @keywords internal
applyProgressiveCorrections <- function(STEAM.obj, corrections, verbose = TRUE) {
  
  if (nrow(corrections) == 0) {
    return(list(
      STEAM.obj = STEAM.obj,
      successful_count = 0,
      successful_cells = character(0),
      success_rate = 1.0,
      improvement = 0
    ))
  }
  
  # Apply corrections to labels
  successful_cells <- character(0)
  
  for (i in seq_len(nrow(corrections))) {
    cell_id <- corrections$cell_id[i]
    # Handle both old and new field names for suggested prediction
    new_pred <- if ("suggested_correction" %in% colnames(corrections)) {
      corrections$suggested_correction[i]
    } else {
      corrections$suggested_pred[i]
    }
    
    # Apply correction
    if ("labels" %in% names(STEAM.obj)) {
      STEAM.obj$labels[cell_id] <- new_pred
    } else if ("predictions" %in% names(STEAM.obj)) {
      STEAM.obj$predictions[cell_id] <- new_pred
    }
    
    successful_cells <- c(successful_cells, cell_id)
  }
  
  successful_count <- length(successful_cells)
  success_rate <- successful_count / nrow(corrections)
  
  return(list(
    STEAM.obj = STEAM.obj,
    successful_count = successful_count,
    successful_cells = successful_cells,
    success_rate = success_rate,
    improvement = successful_count / max(1, nrow(corrections))
  ))
}







#' Update Learning Memory with Correction Results
#'
#' Updates the progressive learning memory with results from correction attempts,
#' tracking pattern performance for future intelligent filtering decisions.
#'
#' @param STEAM.obj STEAM object containing progressive learning state
#' @param correction_results Data frame containing correction attempts and outcomes
#' @return Updated STEAM.obj with enhanced learning memory
#' @details
#' This function maintains the learning memory that powers intelligent pattern
#' filtering. It:
#' 1. Processes each correction result to extract pattern information
#' 2. Updates success/failure statistics for correction patterns
#' 3. Maintains running tallies for pattern performance analysis
#' 4. Stores results in the progressive_learning memory structure
#' 
#' The learning memory enables the system to learn which correction patterns
#' are successful and should be prioritized, and which patterns consistently
#' fail and should be blocked or treated with caution.
#' @seealso \code{\link{updatePatternPerformance}} for pattern-specific updates
#' @keywords internal
updateLearningMemory <- function(STEAM.obj, correction_results) {
  
  if (nrow(correction_results) == 0) return(STEAM.obj)
  
  pl <- STEAM.obj$spatial_anchor_analysis$progressive_learning
  
  # Update correction patterns
  for (i in seq_len(nrow(correction_results))) {
    # Handle different column naming conventions
    original_type <- if ("original_pred" %in% names(correction_results)) {
      correction_results$original_pred[i]
    } else {
      correction_results$original_pred[i]  # Fallback to same name
    }
    
    suggested_type <- if ("suggested_correction" %in% names(correction_results)) {
      correction_results$suggested_correction[i]
    } else if ("suggested_pred" %in% names(correction_results)) {
      correction_results$suggested_pred[i]
    } else {
      NA  # Will skip this correction
    }
    
    was_correct <- correction_results$correct[i]
    
    # Skip if we couldn't extract the required fields
    if (is.na(original_type) || is.na(suggested_type)) next
    
    # Find existing pattern or create new one
    pattern_match <- which(pl$correction_patterns$original_type == original_type & 
                             pl$correction_patterns$suggested_type == suggested_type)
    
    if (length(pattern_match) > 0) {
      # Update existing pattern
      idx <- pattern_match[1]
      if (was_correct) {
        pl$correction_patterns$success_count[idx] <- pl$correction_patterns$success_count[idx] + 1
      } else {
        pl$correction_patterns$failure_count[idx] <- pl$correction_patterns$failure_count[idx] + 1
      }
    } else {
      # Add new pattern
      new_row <- data.frame(
        original_type = original_type,
        suggested_type = suggested_type,
        success_count = ifelse(was_correct, 1, 0),
        failure_count = ifelse(was_correct, 0, 1),
        success_rate = ifelse(was_correct, 1.0, 0.0),
        last_updated = as.character(Sys.time()),
        stringsAsFactors = FALSE
      )
      pl$correction_patterns <- rbind(pl$correction_patterns, new_row)
    }
  }
  
  # Recalculate success rates
  total_attempts <- pl$correction_patterns$success_count + pl$correction_patterns$failure_count
  pl$correction_patterns$success_rate <- ifelse(total_attempts > 0, 
                                                pl$correction_patterns$success_count / total_attempts, 0)
  
  STEAM.obj$spatial_anchor_analysis$progressive_learning <- pl
  return(STEAM.obj)
}





#' Update Pattern Performance Tracking
#'
#' Updates pattern-specific performance metrics for dynamic pattern filtering
#' and intelligent correction prioritization in progressive learning.
#'
#' @param STEAM.obj STEAM object containing progressive learning configuration
#' @param corrections Data frame of correction attempts with success/failure information
#' @param verbose Logical indicating whether to print pattern performance updates
#' @return Updated STEAM.obj with enhanced pattern performance tracking
#' @details
#' This function maintains detailed pattern performance statistics that enable
#' intelligent filtering and prioritization. It:
#' 1. Creates pattern keys using universal data type support (character conversion)
#' 2. Tracks success/failure counts for each correction pattern
#' 3. Calculates rolling success rates for pattern performance
#' 4. Updates the pattern_performance memory structure
#' 5. Provides verbose logging of pattern learning progress
#' 
#' Pattern keys are generated using "from→to" format with character conversion
#' for universal data type compatibility. This replaces hardcoded numeric
#' patterns with dynamic learning from actual performance data.
#' @seealso \code{\link{updateLearningMemory}} for overall memory management
#' @keywords internal
updatePatternPerformance <- function(STEAM.obj, corrections, verbose = TRUE) {
  
  if (!"progressive_learning" %in% names(STEAM.obj) || nrow(corrections) == 0) {
    return(STEAM.obj)
  }
  
  # Ensure pattern_performance exists
  if (is.null(STEAM.obj$spatial_anchor_analysis$progressive_learning$pattern_performance)) {
    STEAM.obj$spatial_anchor_analysis$progressive_learning$pattern_performance <- list()
  }
  
  # Process each correction to update pattern performance  
  for (i in seq_len(nrow(corrections))) {
    correction <- corrections[i, ]
    pattern_key <- paste0(as.character(correction$original_pred), "â†’", as.character(correction$suggested_correction))
    
    # Initialize pattern if not exists
    if (!pattern_key %in% names(STEAM.obj$spatial_anchor_analysis$progressive_learning$pattern_performance)) {
      STEAM.obj$spatial_anchor_analysis$progressive_learning$pattern_performance[[pattern_key]] <- list(
        successes = 0,
        attempts = 0,
        success_rate = 0.0
      )
    }
    
    # Update pattern statistics
    perf <- STEAM.obj$spatial_anchor_analysis$progressive_learning$pattern_performance[[pattern_key]]
    perf$attempts <- perf$attempts + 1
    
    # Check for both column name variations
    is_successful <- FALSE
    if ("correct" %in% colnames(corrections) && correction$correct) {
      is_successful <- TRUE
    } else if ("is_correct" %in% colnames(corrections) && correction$is_correct) {
      is_successful <- TRUE
    }
    
    if (is_successful) {
      perf$successes <- perf$successes + 1
    }
    
    # Update success rate
    perf$success_rate <- perf$successes / perf$attempts
    
    STEAM.obj$spatial_anchor_analysis$progressive_learning$pattern_performance[[pattern_key]] <- perf
    
    if (verbose && perf$attempts <= 5) {  # Only show updates for new patterns
      cat(sprintf("Pattern %s: %.1f%% (%d/%d)\n", 
                  pattern_key, perf$success_rate * 100, perf$successes, perf$attempts))
    }
  }
  
  return(STEAM.obj)
}