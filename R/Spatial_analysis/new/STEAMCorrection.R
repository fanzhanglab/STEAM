#' STEAM Correction for Misclassified Cells (Multi-fold)
#'
#' This implements automatic processing of all available folds with auto-calculated fractions.
#'
#' @param STEAM.obj STEAM object with spatial coordinates and predictions  
#' @param k Number of nearest neighbors to consider (default: 8)
#' @param min_consensus Minimum proportion of neighbors that must agree for correction (default: 0.6)
#' @param max_iterations Maximum number of iterations (default: 10)
#' @param verbose Whether to print progress (default: TRUE)
#' @return Updated STEAM object with naive corrections applied

library(RANN)

# Helper function to extract analysis results from STEAM object
getAnalysisResults <- function(STEAM.obj) {
  
  if (!is.null(STEAM.obj$spatial_anchor_analysis$correction_results)) {
    correction_results <- STEAM.obj$spatial_anchor_analysis$correction_results
    
    total_corrections <- 0
    for (fold_name in names(correction_results)) {
      fold_data <- correction_results[[fold_name]]
      if (!is.null(fold_data$total_corrections)) {
        total_corrections <- total_corrections + fold_data$total_corrections
      }
    }
    
    return(list(
      found = TRUE,
      type = "correction_results",
      data = list(
        iterations_completed = length(correction_results),  # Number of folds
        summary = list(
          total_corrections = total_corrections
        ),
        correction_results = correction_results
      )
    ))
  }

  return(list(found = FALSE, type = "none", data = NULL))
}

STEAMCorrection <- function(STEAM.obj, 
                           k = 8, 
                           min_consensus = 0.6, 
                           max_iterations = 10,
                           verbose = TRUE) {

  if (is.null(STEAM.obj$nested$ncv$outer_result)) {
    stop("No cross-validation results found in STEAM object. Cannot automatically process folds.")
  }
  
  num_folds <- length(STEAM.obj$nested$ncv$outer_result)
  
  if (verbose) {
    cat("=== STEAM Correction System ===\n")
    cat(sprintf("Automatically processing all %d folds\n", num_folds))
    cat(sprintf("Using k=%d neighbors, min_consensus=%.2f\n", k, min_consensus))
  }

  all_fold_results <- list()
  
  for (fold in 1:num_folds) {
    if (verbose) {
      cat(sprintf("\n=== PROCESSING FOLD %d/%d ===\n", fold, num_folds))
    }

  coordinates <- STEAM.obj$spatial
  if (is.null(coordinates)) {
    stop("No spatial coordinates found in STEAM object")
  }
  coordinates <- as.matrix(coordinates)
  
  fold_result <- STEAM.obj$nested$ncv$outer_result[[fold]]
  if (is.null(fold_result$preds) || !"predy" %in% names(fold_result$preds)) {
    stop(sprintf("No predictions found in fold %d", fold))
  }
  
  current_preds <- as.character(fold_result$preds$predy)
  if (is.null(names(fold_result$preds$predy))) {
    names(current_preds) <- rownames(fold_result$preds)
  } else {
    names(current_preds) <- names(fold_result$preds$predy)
  }
  
  ground_truth <- NULL
  if (!is.null(fold_result$preds$testy)) {
    ground_truth <- as.character(fold_result$preds$testy)
    if (is.null(names(fold_result$preds$testy))) {
      names(ground_truth) <- rownames(fold_result$preds)
    } else {
      names(ground_truth) <- names(fold_result$preds$testy)
    }
  } else {
    stop("No ground truth found - cannot identify misclassified cells")
  }
  
  common_cells <- intersect(names(current_preds), names(ground_truth))
  if (length(common_cells) == 0) {
    stop("No matching cells between predictions and ground truth")
  }
  
  initial_accuracy <- sum(current_preds[common_cells] == ground_truth[common_cells]) / length(common_cells)
  if (verbose) {
    cat(sprintf("Fold %d initial accuracy: %.1f%%\n", fold, initial_accuracy * 100))
  }
  
  misclassified_cells <- common_cells[current_preds[common_cells] != ground_truth[common_cells]]
  if (length(misclassified_cells) == 0) {
    if (verbose) cat("No misclassified cells found\n")
    return(STEAM.obj)
  }
  
  if (verbose) {
    cat(sprintf("Found %d misclassified cells to potentially correct\n", length(misclassified_cells)))
  }
  
  num_misclassified <- length(misclassified_cells)
  if (num_misclassified <= 50) {
    fractions <- c(0.2, 0.3, 0.4)
  } else if (num_misclassified <= 200) {
    fractions <- c(0.1, 0.15, 0.2, 0.25)
  } else {
    fractions <- c(0.05, 0.1, 0.15, 0.2, 0.25)
  }
  
  if (verbose) {
    cat(sprintf("Auto-calculated fractions: %s\n", paste(fractions, collapse=", ")))
  }
  
  fold_common_cells <- intersect(rownames(coordinates), common_cells)
  fold_coordinates <- coordinates[fold_common_cells, , drop = FALSE]
  fold_preds <- current_preds[fold_common_cells]
  fold_ground_truth <- ground_truth[fold_common_cells]
  
  knn_result <- nn2(fold_coordinates, k = min(k, nrow(fold_coordinates) - 1))
  
  correction_history <- list()
  iteration <- 0
  current_misclassified <- names(fold_preds)[fold_preds != fold_ground_truth]
  total_corrections <- 0
  
  correction_history[[1]] <- list(
    iteration = 0,
    corrections_applied = 0,
    total_corrections = 0,
    accuracy = initial_accuracy,
    remaining_misclassified = length(current_misclassified),
    predictions_snapshot = fold_preds
  )
  
  while (iteration < max_iterations && length(current_misclassified) > 0) {
    iteration <- iteration + 1
    
    if (verbose) {
      cat(sprintf("\n--- Iteration %d ---\n", iteration))
      cat(sprintf("Currently misclassified cells: %d\n", length(current_misclassified)))
    }
    
    correction_suggestions <- generateCorrections(
      fold_preds, fold_coordinates, knn_result, min_consensus, current_misclassified
    )
    
    if (length(correction_suggestions) == 0) {
      if (verbose) cat("No valid corrections found - stopping\n")
      break
    }
    
    if (verbose) {
      cat(sprintf("Generated %d correction suggestions\n", length(correction_suggestions)))
    }
    
    corrections_this_iteration <- 0
    
    for (fraction in fractions) {
      num_to_apply <- max(1, round(fraction * length(current_misclassified)))
      num_to_apply <- min(num_to_apply, length(correction_suggestions))
      
      if (num_to_apply > 0) {
        corrections_to_apply <- correction_suggestions[1:num_to_apply]
        consensus_values <- sapply(corrections_to_apply, function(x) x$consensus)
        
        if (verbose) {
          cat(sprintf("  Fraction %.2f: Applying %d corrections (consensus range: %.2f-%.2f)\n", 
                     fraction, length(corrections_to_apply), min(consensus_values), max(consensus_values)))
        }
        
        applied_count <- 0
        for (correction in corrections_to_apply) {
          if (fold_preds[correction$cell_id] == correction$current_prediction) {
            fold_preds[correction$cell_id] <- correction$new_prediction
            applied_count <- applied_count + 1
            
            if (verbose) {
              is_correct <- fold_ground_truth[correction$cell_id] == correction$new_prediction
              status <- if (is_correct) "✓" else "✗"
              cat(sprintf("    %s: %s → %s (consensus: %.2f) %s\n", 
                         correction$cell_id, correction$current_prediction, 
                         correction$new_prediction, correction$consensus, status))
            }
          }
        }
        corrections_this_iteration <- corrections_this_iteration + applied_count
        
        if (applied_count > 0) {
          correction_suggestions <- correction_suggestions[-(1:num_to_apply)]
        }
      }
    }
    
    current_misclassified <- names(fold_preds)[fold_preds != fold_ground_truth]
    total_corrections <- total_corrections + corrections_this_iteration
    
    current_accuracy <- sum(fold_preds == fold_ground_truth) / length(fold_ground_truth)
    
    if (verbose) {
      cat(sprintf("Iteration %d: Applied %d corrections\n", iteration, corrections_this_iteration))
    }
    
    correction_history[[iteration + 1]] <- list(
      iteration = iteration,
      corrections_applied = corrections_this_iteration,
      total_corrections = total_corrections,
      accuracy = current_accuracy,
      remaining_misclassified = length(current_misclassified),
      predictions_snapshot = fold_preds
    )
    
    if (corrections_this_iteration == 0) {
      if (verbose) cat("No valid corrections found - stopping\n")
      break
    }
  }
  
  final_accuracy <- sum(fold_preds == fold_ground_truth) / length(fold_ground_truth)
  
  if (verbose) {
    cat(sprintf("\n=== Final Summary ===\n"))
    cat(sprintf("Total corrections made: %d\n", total_corrections))
    cat(sprintf("Original accuracy: %.3f\n", initial_accuracy))
    cat(sprintf("Final accuracy: %.3f\n", final_accuracy))
    cat(sprintf("Improvement: %.3f\n", final_accuracy - initial_accuracy))
    cat(sprintf("Misclassified cells reduced: %d → %d\n", 
               sum(current_preds[common_cells] != ground_truth[common_cells]),
               sum(fold_preds != fold_ground_truth)))
  }
  
    all_fold_results[[paste0("fold_", fold)]] <- list(
      corrected_labels = fold_preds,
      correction_history = correction_history,
      initial_accuracy = initial_accuracy,
      final_accuracy = final_accuracy,
      total_corrections = total_corrections,
      fold_number = fold
    )
  }
  
  if (verbose) {
    cat(sprintf("\n=== MULTI-FOLD SUMMARY ===\n"))
    for (i in seq_along(all_fold_results)) {
      fold_result <- all_fold_results[[i]]
      cat(sprintf("Fold %d: %.1f%% → %.1f%% (%.1f%% improvement, %d corrections)\n",
                 fold_result$fold_number,
                 fold_result$initial_accuracy * 100,
                 fold_result$final_accuracy * 100,
                 (fold_result$final_accuracy - fold_result$initial_accuracy) * 100,
                 fold_result$total_corrections))
    }
  }
  
  if (!is.null(STEAM.obj$spatial_anchor_analysis)) {
    STEAM.obj$spatial_anchor_analysis$correction_results <- all_fold_results
  } else {
    STEAM.obj$spatial_anchor_analysis <- list(
      correction_results = all_fold_results
    )
  }
  
  STEAM.obj$spatial_anchor_analysis$parameters <- list(
    k = k,
    min_consensus = min_consensus,
    max_iterations = max_iterations
  )

  return(STEAM.obj)
}

generateCorrections <- function(predictions, coordinates, knn_result, min_consensus, target_cells = NULL) {
  
  if (is.null(target_cells)) {
    target_cells <- names(predictions)
  }
  
  corrections <- list()
  
  for (cell_name in target_cells) {
    if (!cell_name %in% names(predictions)) next
    
    cell_idx <- which(rownames(coordinates) == cell_name)
    if (length(cell_idx) == 0) next
    
    neighbor_indices <- knn_result$nn.idx[cell_idx, ]
    neighbor_names <- rownames(coordinates)[neighbor_indices]

    valid_neighbors <- intersect(neighbor_names, names(predictions))
    if (length(valid_neighbors) < 3) next  # Need at least 3 neighbors
    
    neighbor_preds <- predictions[valid_neighbors]
    neighbor_preds <- neighbor_preds[!is.na(neighbor_preds)]
    
    if (length(neighbor_preds) == 0) next
    
    pred_counts <- table(neighbor_preds)
    majority_pred <- names(pred_counts)[which.max(pred_counts)]
    consensus <- max(pred_counts) / length(neighbor_preds)

    current_pred <- predictions[cell_name]
    if (consensus >= min_consensus && majority_pred != current_pred) {
      corrections <- append(corrections, list(list(
        cell_id = cell_name,
        current_prediction = current_pred,
        new_prediction = majority_pred,
        consensus = consensus
      )))
    }
  }
  
  if (length(corrections) > 0) {
    consensus_values <- sapply(corrections, function(x) x$consensus)
    corrections <- corrections[order(consensus_values, decreasing = TRUE)]
  }
  
  return(corrections)
}
