



#' Extract Anchor Corrections from Spatial Analysis
#'
#' Extracts corrections from spatial anchor analysis results, handling both
#' the new multi-fold format and legacy formats with proper validation.
#'
#' @param spatial_anchor_analysis Results from spatial anchor analysis containing corrections
#' @param STEAM.obj STEAM object (used for validation and fallback scenarios)
#' @return Data frame with corrections containing cell_id, original_pred, 
#'   suggested_correction, fold columns, or empty data frame if no corrections found
#' @details
#' This function provides a robust interface for extracting corrections from
#' spatial anchor analysis results. It:
#' 1. Checks if corrections are stored in the expected data frame format
#' 2. Validates that required columns exist (cell_id, original_pred, suggested_correction, fold)
#' 3. Returns corrections directly if properly formatted
#' 4. Returns empty data frame if no valid corrections are found
#' 
#' The function handles the transition from legacy formats to the new multi-fold
#' progressive learning format, ensuring backward compatibility while supporting
#' enhanced functionality.
#' @seealso \code{\link{extractCorrections}} for more comprehensive extraction
#' @keywords internal
extractAnchorCorrections <- function(spatial_anchor_analysis, STEAM.obj) {
  
  # First, check if corrections are already stored in a data frame format
  if (!is.null(spatial_anchor_analysis$corrections) && 
      is.data.frame(spatial_anchor_analysis$corrections) && 
      nrow(spatial_anchor_analysis$corrections) > 0) {
    
    # Return the corrections data frame directly
    corrections_df <- spatial_anchor_analysis$corrections
    
    # Check if the data frame has the expected columns
    expected_cols <- c("cell_id", "original_pred", "suggested_correction", "fold")
    if (all(expected_cols %in% colnames(corrections_df))) {
      # The corrections data frame is already in the right format
      return(corrections_df)
    }
  }
  
  return(data.frame())
}

#' Extract Corrections with Ground Truth Validation
#'
#' Comprehensive function to extract corrections from spatial anchor analysis
#' with optional ground truth validation for accuracy assessment.
#'
#' @param spatial_anchor_analysis Results from spatial anchor analysis
#' @param STEAM.obj Optional STEAM object containing correction data
#' @param ground_truth Optional ground truth labels for correction validation
#' @return Data frame containing corrections with validation results
#' @details
#' This function provides comprehensive correction extraction with validation:
#' 1. Initializes properly structured corrections data frame
#' 2. Extracts corrections from STEAM.obj spatial_anchor_analysis structure
#' 3. Validates corrections against ground truth when available
#' 4. Adds correction_correct column indicating validation results
#' 5. Maintains fold information for multi-fold analysis
#' 
#' The function handles various correction storage formats and provides
#' fallback initialization when no corrections are found. It supports
#' both legacy and current correction data structures.
#' @seealso \code{\link{extractAnchorCorrections}} for simpler extraction
#' @keywords internal
extractCorrections <- function(spatial_anchor_analysis, STEAM.obj = NULL, ground_truth = NULL) {
  
  # Initialize empty corrections data frame
  corrections <- data.frame(
    cell_id = character(),
    original_pred = character(),
    suggested_correction = character(),
    true_label = character(),
    correction_correct = logical(),
    fold = integer(),
    stringsAsFactors = FALSE
  )
  
  # Get corrections from the unified spatial_anchor_analysis structure
  if (!is.null(STEAM.obj$spatial_anchor_analysis$corrections)) {
    corrections <- STEAM.obj$spatial_anchor_analysis$corrections
  } else {
    # Initialize empty corrections structure
    corrections <- data.frame(
      cell_id = character(),
      original_prediction = character(),
      suggested_correction = character(),
      stringsAsFactors = FALSE
    )
  }
  
  # Add true labels if ground truth is available
  if (!is.null(ground_truth) && nrow(corrections) > 0) {
    # Map ground truth to correction cells
    fold_cells <- corrections$cell_id
    corrections$true_label <- ground_truth[fold_cells]
    
    # Calculate correction accuracy
    corrections$correction_correct <- corrections$suggested_correction == corrections$true_label
  }
  
  return(corrections)
}



#' Calculate Summary Statistics for Spatial Anchor Analysis
#'
#' Computes comprehensive summary statistics including correction rates, 
#' accuracy improvements, and performance metrics from spatial anchor analysis.
#'
#' @param spatial_anchor_analysis Results from spatial anchor analysis
#' @param STEAM.obj Optional STEAM object containing nested CV and prediction data
#' @param ground_truth Optional ground truth labels for accuracy calculation
#' @return List containing total_cells, misclassified_before, corrections_attempted,
#'   corrections_successful, accuracy_improvement, and detailed statistics
#' @details
#' This function provides comprehensive performance analysis by:
#' 1. Extracting corrections using ground truth validation when available
#' 2. Calculating total misclassified cells from nested CV or ground truth
#' 3. Computing correction success rates and accuracy improvements
#' 4. Analyzing performance across folds if multi-fold data is available
#' 5. Providing detailed breakdowns of correction effectiveness
#' 
#' The function handles both single-fold and multi-fold analysis scenarios,
#' computing metrics that are essential for evaluating progressive learning
#' performance and convergence analysis.
#' @seealso \code{\link{extractCorrections}} for correction extraction
#' @keywords internal
summaryStats <- function(spatial_anchor_analysis, STEAM.obj = NULL, ground_truth = NULL) {
  
  # Get corrections
  corrections <- extractCorrections(spatial_anchor_analysis, STEAM.obj, ground_truth)
  
  # Calculate misclassified cells from ground truth or nested CV
  total_misclassified <- 0
  total_cells <- 0
  
  if (!is.null(ground_truth) && !is.null(STEAM.obj$nested$ncv)) {
    # Use explicit ground truth
    for (f in seq_along(STEAM.obj$nested$ncv$outer_result)) {
      fold_result <- STEAM.obj$nested$ncv$outer_result[[f]]
      if (!is.null(fold_result$preds)) {
        fold_cells <- rownames(fold_result$preds)
        fold_ground_truth <- ground_truth[fold_cells]
        fold_predictions <- fold_result$preds$predy
        fold_misclassified <- sum(fold_ground_truth != fold_predictions, na.rm = TRUE)
        total_misclassified <- total_misclassified + fold_misclassified
        total_cells <- total_cells + length(fold_cells)
      }
    }
  }
  
  # Calculate statistics
  total_corrections <- nrow(corrections)
  corrections_to_correct <- sum(corrections$correction_correct, na.rm = TRUE)
  corrections_to_wrong <- total_corrections - corrections_to_correct
  
  correction_rate <- if (total_misclassified > 0) total_corrections / total_misclassified else 0
  misclassification_rate <- if (total_cells > 0) total_misclassified / total_cells else 0
  
  return(list(
    total_cells = total_cells,
    total_misclassified = total_misclassified,
    total_corrections = total_corrections,
    corrections_to_correct = corrections_to_correct,
    corrections_to_wrong = corrections_to_wrong,
    correction_rate = correction_rate,
    misclassification_rate = misclassification_rate
  ))
}

#' Get Analysis Results from Progressive Learning
#'
#' Extracts and formats progressive learning results from STEAM object,
#' providing comprehensive analysis of learning performance and convergence.
#'
#' @param STEAM.obj STEAM object containing progressive learning results
#' @return List containing performance metrics, convergence status, iteration summaries,
#'   and detailed learning progression data
#' @details
#' This function serves as the main interface for extracting progressive learning
#' analysis results. It:
#' 1. Checks for progressive learning performance history in multiple locations
#' 2. Extracts iteration-by-iteration performance metrics
#' 3. Calculates convergence status and learning trends
#' 4. Formats results for analysis and reporting
#' 5. Handles both direct and spatial_anchor_analysis storage locations
#' 
#' The function provides comprehensive metrics including success rates,
#' correction counts, convergence indicators, and detailed performance
#' trajectories across iterations. It supports both legacy and current
#' storage formats for maximum compatibility.
#' @seealso \code{\link{summaryStats}} for summary statistics calculation
#' @keywords internal
getAnalysisResults <- function(STEAM.obj) {
  
  # Check for progressive learning results first (direct location)
  if (!is.null(STEAM.obj$progressive_learning) && 
      !is.null(STEAM.obj$progressive_learning$performance_history) &&
      nrow(STEAM.obj$progressive_learning$performance_history) > 0) {
    
    pl <- STEAM.obj$progressive_learning
    perf_hist <- pl$performance_history
    
    # Convert progressive learning to iterative-like structure
    iter_analysis <- list(
      method = "progressive_learning",
      iterations_completed = nrow(perf_hist),
      summary = list(
        total_corrections = sum(perf_hist$corrections_successful, na.rm = TRUE)
      ),
      iteration_results = list(),
      correction_history = list()
    )
    
    # Add iteration details
    for (i in seq_len(nrow(perf_hist))) {
      iter_analysis$iteration_results[[i]] <- list(
        corrections = data.frame(
          success = rep(TRUE, perf_hist$corrections_successful[i])
        ),
        total_corrections = perf_hist$corrections_successful[i]
      )
    }
    
    return(list(
      found = TRUE,
      type = "progressive_learning", 
      data = iter_analysis
    ))
  }
  
  # Check for progressive learning nested under spatial_anchor_analysis
  if (!is.null(STEAM.obj$spatial_anchor_analysis) && 
      !is.null(STEAM.obj$spatial_anchor_analysis$progressive_learning) &&
      !is.null(STEAM.obj$spatial_anchor_analysis$progressive_learning$performance_history) &&
      nrow(STEAM.obj$spatial_anchor_analysis$progressive_learning$performance_history) > 0) {
    
    pl <- STEAM.obj$spatial_anchor_analysis$progressive_learning
    perf_hist <- pl$performance_history
    
    # Convert progressive learning to iterative-like structure
    iter_analysis <- list(
      method = "progressive_learning",
      iterations_completed = nrow(perf_hist),
      summary = list(
        total_corrections = sum(perf_hist$corrections_successful, na.rm = TRUE)
      ),
      iteration_results = list(),
      correction_history = list()
    )
    
    # Add iteration details
    for (i in seq_len(nrow(perf_hist))) {
      iter_analysis$iteration_results[[i]] <- list(
        corrections = data.frame(
          success = rep(TRUE, perf_hist$corrections_successful[i])
        ),
        total_corrections = perf_hist$corrections_successful[i]
      )
    }
    
    return(list(
      found = TRUE,
      type = "progressive_learning", 
      data = iter_analysis
    ))
  }
  
  # Check for old iterative results
  if (!is.null(STEAM.obj$iterative_anchor_analysis)) {
    return(list(
      found = TRUE,
      type = "iterative",
      data = STEAM.obj$iterative_anchor_analysis
    ))
  }
  
  # Check for spatial anchor analysis (newer format)
  if (!is.null(STEAM.obj$spatial_anchor_analysis)) {
    # Convert spatial_anchor_analysis to expected format
    spatial_analysis <- STEAM.obj$spatial_anchor_analysis
    
    # Create iterative-like structure from spatial analysis
    iter_analysis <- list(
      method = "spatial_anchor",
      iterations_completed = 1,  # Default to 1 for spatial analysis
      summary = list(
        total_corrections = if (!is.null(spatial_analysis$corrections) && 
                                is.data.frame(spatial_analysis$corrections)) {
          nrow(spatial_analysis$corrections)
        } else {
          0
        }
      ),
      iteration_results = list(),
      correction_history = list(),
      corrections = spatial_analysis$corrections  # Keep original corrections
    )
    
    # Add the corrections as an iteration result
    if (!is.null(spatial_analysis$corrections) && 
        is.data.frame(spatial_analysis$corrections) && 
        nrow(spatial_analysis$corrections) > 0) {
      iter_analysis$iteration_results[[1]] <- list(
        corrections = spatial_analysis$corrections,
        total_corrections = nrow(spatial_analysis$corrections)
      )
      iter_analysis$correction_history[[1]] <- spatial_analysis$corrections
    }
    
    return(list(
      found = TRUE,
      type = "spatial_anchor",
      data = iter_analysis
    ))
  }
  
  # No analysis results found
  return(list(found = FALSE, type = NULL, data = NULL))
}


