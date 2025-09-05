
#' Identify Uncertain Cells for Progressive Learning
#'
#' Identifies cells with low confidence scores that are candidates for correction
#' in progressive learning. Uses multimodal confidence calculation based on 
#' spatial consistency, expression similarity, and model probabilities.
#'
#' @param STEAM.obj STEAM object containing spatial coordinates, expression data, 
#'   and prediction information
#' @param candidate_cells Character vector of cell IDs to evaluate for uncertainty
#' @param confidence_threshold Numeric threshold below which cells are considered uncertain
#' @param expression_matrix Optional expression matrix (currently unused, for future expansion)
#' @param verbose Logical indicating whether to print progress information
#' @return Character vector of cell IDs that have confidence scores below the threshold
#' @details
#' This function:
#' 1. Extracts predictions for candidate cells from the STEAM object
#' 2. Calculates multimodal confidence scores using calculateConfidenceScores()
#' 3. Identifies cells below the confidence threshold as "uncertain"
#' 4. Returns these uncertain cells for potential correction
#' 
#' The confidence calculation considers spatial neighborhood consistency, expression
#' similarity with neighbors, and model prediction probabilities when available.
#' @seealso \code{\link{calculateConfidenceScores}} for confidence calculation details
#' @keywords internal
identifyUncertainCells <- function(STEAM.obj, candidate_cells, confidence_threshold, 
                                   expression_matrix = NULL, verbose = TRUE) {
  
  # Get predictions for candidate cells
  all_predictions <- extractPreds(STEAM.obj)
  candidate_predictions <- all_predictions[candidate_cells]
  
  # Calculate confidence scores
  confidence_scores <- calculateConfidenceScores(STEAM.obj, candidate_predictions)
  
  # Identify uncertain cells (below threshold)
  uncertain_mask <- confidence_scores < confidence_threshold
  uncertain_cells <- names(candidate_predictions)[uncertain_mask]
  
  if (verbose) {
    cat(sprintf("Identified %d uncertain cells from %d candidates (confidence < %.3f)\n", 
                length(uncertain_cells), length(candidate_cells), confidence_threshold))
  }
  
  return(uncertain_cells)
}



