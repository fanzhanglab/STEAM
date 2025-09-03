#' Helper Functions for Iterative Analysis Without Ground Truth
#' 
#' This file contains all helper functions needed for the no-ground-truth
#' iterative analysis approach, implementing alternative validation methods
#' based on spatial coherence, expression consistency, and model confidence.
#' 
#' Functions are organized into the following sections:
#' 1. Core Metrics & Validation
#' 2. Uncertainty Identification
#' 3. Correction Generation & Evaluation
#' 4. Validation & Application
#' 5. Utility Functions

# =============================================================================
# 1. CORE METRICS & VALIDATION
# =============================================================================

#' Compute Baseline Validation Metrics
#' 
#' @param STEAM.obj STEAM object
#' @param expression_matrix Expression matrix
#' @return List of baseline validation metrics
#' @keywords internal
baselineMetrics <- function(STEAM.obj, expression_matrix) {
  
  # Extract predictions and spatial data
  all_predictions <- extractPreds(STEAM.obj)
  coords <- STEAM.obj$spatial
  
  # Align data
  common_cells <- intersect(names(all_predictions), rownames(coords))
  if (!is.null(expression_matrix)) {
    expr_cells <- if (ncol(expression_matrix) > nrow(expression_matrix)) {
      colnames(expression_matrix)  # genes in rows
    } else {
      rownames(expression_matrix)  # genes in columns
    }
    common_cells <- intersect(common_cells, expr_cells)
  }
  
  if (length(common_cells) < 50) {
    stop("Too few cells for baseline validation metrics")
  }
  
  # Spatial coherence: measure neighborhood homogeneity
  spatial_coherence <- spatialCoherence(
    predictions = all_predictions[common_cells],
    coordinates = coords[common_cells, , drop = FALSE],
    k = 15
  )
  
  # Expression coherence: if available
  expression_coherence <- if (!is.null(expression_matrix)) {
    expressionCoherence(
      predictions = all_predictions[common_cells],
      expression_matrix = expression_matrix,
      common_cells = common_cells,
      genes_transpose = TRUE
    )
  } else {
    0.5  # neutral baseline
  }
  
  # Model confidence metrics
  confidence_metrics <- confidenceMetrics(STEAM.obj, common_cells)
  
  # Global consistency
  global_entropy <- predictionEntropy(all_predictions[common_cells])
  
  return(list(
    spatial_coherence = spatial_coherence,
    expression_coherence = expression_coherence,
    mean_confidence = confidence_metrics$mean_confidence,
    confidence_entropy = confidence_metrics$confidence_entropy,
    global_entropy = global_entropy,
    n_cells = length(common_cells),
    baseline_timestamp = Sys.time()
  ))
}

#' Calculate Spatial Coherence
#' 
#' @param predictions Cell predictions
#' @param coordinates Spatial coordinates
#' @param k Number of neighbors
#' @return Spatial coherence score (0-1)
#' @keywords internal
spatialCoherence <- function(predictions, coordinates, k = 15) {
  
  n_cells <- length(predictions)
  if (n_cells < k) k <- max(1, n_cells - 1)
  
  # Build KNN
  knn <- buildKnn(coordinates, k)
  
  # Calculate coherence for each cell
  coherence_scores <- sapply(seq_len(n_cells), function(i) {
    neighbors <- knn$idx_list[[i]]
    if (length(neighbors) == 0) return(0.5)
    
    neighbor_preds <- predictions[neighbors]
    agreement <- mean(neighbor_preds == predictions[i], na.rm = TRUE)
    return(agreement)
  })
  
  return(mean(coherence_scores, na.rm = TRUE))
}

#' Calculate Expression Coherence
#' 
#' @param predictions Cell predictions
#' @param expression_matrix Expression matrix
#' @param common_cells Common cell IDs
#' @param genes_transpose Whether genes are in columns
#' @return Expression coherence score (0-1)
#' @keywords internal
expressionCoherence <- function(predictions, expression_matrix, 
                                         common_cells, genes_transpose = TRUE) {
  
  # Prepare expression data
  if (genes_transpose) {
    if (ncol(expression_matrix) > nrow(expression_matrix)) {
      expr_mat <- t(expression_matrix)  # transpose if genes in rows
    } else {
      expr_mat <- expression_matrix
    }
  } else {
    expr_mat <- expression_matrix
  }
  
  # Get common cells
  expr_cells <- intersect(common_cells, rownames(expr_mat))
  if (length(expr_cells) < 20) return(0.5)
  
  expr_data <- expr_mat[expr_cells, , drop = FALSE]
  cell_preds <- predictions[expr_cells]
  
  # Calculate within-class vs between-class expression similarity
  unique_classes <- unique(cell_preds)
  if (length(unique_classes) < 2) return(0.5)
  
  within_class_sim <- 0
  between_class_sim <- 0
  n_within <- 0
  n_between <- 0
  
  for (class1 in unique_classes) {
    cells_class1 <- names(cell_preds)[cell_preds == class1]
    if (length(cells_class1) < 2) next
    
    # Within-class similarity
    class1_expr <- expr_data[cells_class1, , drop = FALSE]
    if (nrow(class1_expr) >= 2) {
      cor_matrix <- cor(t(class1_expr), use = "pairwise.complete.obs")
      diag(cor_matrix) <- NA
      within_vals <- cor_matrix[!is.na(cor_matrix)]
      if (length(within_vals) > 0) {
        within_class_sim <- within_class_sim + sum(within_vals, na.rm = TRUE)
        n_within <- n_within + length(within_vals)
      }
    }
    
    # Between-class similarity
    for (class2 in unique_classes) {
      if (class1 >= class2) next  # avoid double counting
      
      cells_class2 <- names(cell_preds)[cell_preds == class2]
      if (length(cells_class2) < 1) next
      
      class2_expr <- expr_data[cells_class2, , drop = FALSE]
      
      # Calculate cross-class correlations
      for (i in seq_len(min(10, nrow(class1_expr)))) {  # sample to avoid computation explosion
        for (j in seq_len(min(10, nrow(class2_expr)))) {
          cor_val <- cor(class1_expr[i, ], class2_expr[j, ], use = "pairwise.complete.obs")
          if (!is.na(cor_val)) {
            between_class_sim <- between_class_sim + cor_val
            n_between <- n_between + 1
          }
        }
      }
    }
  }
  
  # Calculate coherence as relative within-class similarity
  if (n_within == 0 || n_between == 0) return(0.5)
  
  mean_within <- within_class_sim / n_within
  mean_between <- between_class_sim / n_between
  
  # Normalize to 0-1 scale
  coherence <- (mean_within - mean_between + 2) / 4  # +2 to handle negative correlations
  return(max(0, min(1, coherence)))
}

#' Calculate Confidence Metrics
#' 
#' @param STEAM.obj STEAM object
#' @param common_cells Common cell IDs
#' @return List with confidence metrics
#' @keywords internal
confidenceMetrics <- function(STEAM.obj, common_cells) {
  
  # Extract probabilities from nested CV
  all_probs <- list()
  all_cells <- character()
  
  for (f in seq_along(STEAM.obj$nested$ncv$outer_result)) {
    fold_result <- STEAM.obj$nested$ncv$outer_result[[f]]
    if (!is.null(fold_result$probs)) {
      fold_cells <- intersect(rownames(fold_result$probs), common_cells)
      if (length(fold_cells) > 0) {
        all_probs[[f]] <- fold_result$probs[fold_cells, , drop = FALSE]
        all_cells <- c(all_cells, fold_cells)
      }
    }
  }
  
  if (length(all_probs) == 0) {
    return(list(mean_confidence = 0.5, confidence_entropy = 1.0))
  }
  
  # Combine probabilities
  combined_probs <- do.call(rbind, all_probs)
  
  # Calculate max confidence (highest probability per cell)
  max_confidences <- apply(combined_probs, 1, max, na.rm = TRUE)
  mean_confidence <- mean(max_confidences, na.rm = TRUE)
  
  # Calculate entropy of confidence distribution
  confidence_entropy <- -sum(max_confidences * log(max_confidences + 1e-10), na.rm = TRUE) / length(max_confidences)
  
  return(list(
    mean_confidence = mean_confidence,
    confidence_entropy = confidence_entropy,
    n_cells_with_probs = nrow(combined_probs)
  ))
}

#' Calculate Prediction Entropy
#' 
#' @param predictions Cell predictions
#' @return Entropy of prediction distribution
#' @keywords internal
predictionEntropy <- function(predictions) {
  
  pred_table <- table(predictions)
  pred_probs <- pred_table / sum(pred_table)
  
  entropy <- -sum(pred_probs * log(pred_probs + 1e-10))
  return(entropy)
}

# =============================================================================
# 2. UNCERTAINTY IDENTIFICATION
# =============================================================================

#' Identify Uncertain Cells Using Alternative Criteria
#' 
#' @param STEAM.obj STEAM object
#' @param expression_matrix Expression matrix
#' @param iteration Current iteration
#' @param verbose Print messages
#' @return Vector of uncertain cell IDs
#' @keywords internal
identifyUncertain <- function(STEAM.obj, expression_matrix, 
                                                iteration, verbose = FALSE) {
  
  # Extract all predictions
  all_predictions <- extractPreds(STEAM.obj)
  coords <- STEAM.obj$spatial
  common_cells <- intersect(names(all_predictions), rownames(coords))
  
  if (length(common_cells) < 20) {
    if (verbose) message("Too few cells for uncertainty analysis")
    return(character())
  }
  
  uncertain_cells <- character()
  
  # 1. Low model confidence cells
  confidence_uncertain <- lowConfidence(STEAM.obj, common_cells)
  
  # 2. Spatially inconsistent cells
  spatial_uncertain <- spatialInconsistent(
    predictions = all_predictions[common_cells],
    coordinates = coords[common_cells, , drop = FALSE],
    k = 15
  )
  
  # 3. Expression-prediction mismatch (if expression available)
  expression_uncertain <- if (!is.null(expression_matrix)) {
    expressionMismatch(
      predictions = all_predictions[common_cells],
      expression_matrix = expression_matrix,
      common_cells = common_cells
    )
  } else {
    character()
  }
  
  # Combine criteria with voting
  all_uncertain_candidates <- unique(c(confidence_uncertain, spatial_uncertain, expression_uncertain))
  
  # Require at least 2 criteria to be met (more conservative)
  uncertain_votes <- sapply(all_uncertain_candidates, function(cell) {
    votes <- 0
    if (cell %in% confidence_uncertain) votes <- votes + 1
    if (cell %in% spatial_uncertain) votes <- votes + 1
    if (cell %in% expression_uncertain) votes <- votes + 1
    return(votes)
  })
  
  # In later iterations, be more lenient (accept 1 vote), in early iterations require 2+
  min_votes <- if (iteration <= 2) 2 else 1
  uncertain_cells <- all_uncertain_candidates[uncertain_votes >= min_votes]
  
  # Cap at 2% of total cells to prevent cascade errors
  max_uncertain <- max(1, round(0.02 * length(common_cells)))
  if (length(uncertain_cells) > max_uncertain) {
    # Prioritize by number of votes
    vote_order <- order(uncertain_votes[uncertain_cells], decreasing = TRUE)
    uncertain_cells <- uncertain_cells[vote_order[1:max_uncertain]]
  }
  
  if (verbose) {
    message(sprintf("Identified %d uncertain cells: %d confidence, %d spatial, %d expression",
                   length(uncertain_cells), length(confidence_uncertain), 
                   length(spatial_uncertain), length(expression_uncertain)))
  }
  
  return(uncertain_cells)
}

#' Identify Low Confidence Cells
#' 
#' @param STEAM.obj STEAM object
#' @param common_cells Common cell IDs
#' @return Vector of low confidence cell IDs
#' @keywords internal
lowConfidence <- function(STEAM.obj, common_cells) {
  
  # Extract probabilities
  all_probs <- list()
  for (f in seq_along(STEAM.obj$nested$ncv$outer_result)) {
    fold_result <- STEAM.obj$nested$ncv$outer_result[[f]]
    if (!is.null(fold_result$probs)) {
      fold_cells <- intersect(rownames(fold_result$probs), common_cells)
      if (length(fold_cells) > 0) {
        all_probs[[f]] <- fold_result$probs[fold_cells, , drop = FALSE]
      }
    }
  }
  
  if (length(all_probs) == 0) return(character())
  
  combined_probs <- do.call(rbind, all_probs)
  max_confidences <- apply(combined_probs, 1, max, na.rm = TRUE)
  
  # Identify cells with confidence below threshold
  confidence_threshold <- quantile(max_confidences, 0.25, na.rm = TRUE)  # bottom quartile
  low_confidence_cells <- names(max_confidences)[max_confidences < confidence_threshold]
  
  return(low_confidence_cells)
}

#' Identify Spatially Inconsistent Cells
#' 
#' @param predictions Cell predictions
#' @param coordinates Spatial coordinates
#' @param k Number of neighbors
#' @return Vector of spatially inconsistent cell IDs
#' @keywords internal
spatialInconsistent <- function(predictions, coordinates, k = 15) {
  
  n_cells <- length(predictions)
  if (n_cells < k) k <- max(1, n_cells - 1)
  
  # Build KNN
  knn <- buildKnn(coordinates, k)
  
  # Find cells with low spatial agreement
  spatial_agreements <- sapply(seq_len(n_cells), function(i) {
    neighbors <- knn$idx_list[[i]]
    if (length(neighbors) == 0) return(1.0)
    
    neighbor_preds <- predictions[neighbors]
    agreement <- mean(neighbor_preds == predictions[i], na.rm = TRUE)
    return(agreement)
  })
  
  # Identify cells with agreement below threshold
  agreement_threshold <- quantile(spatial_agreements, 0.25, na.rm = TRUE)  # bottom quartile
  inconsistent_cells <- names(predictions)[spatial_agreements < agreement_threshold]
  
  return(inconsistent_cells)
}

#' Identify Expression-Prediction Mismatched Cells
#' 
#' @param predictions Cell predictions
#' @param expression_matrix Expression matrix
#' @param common_cells Common cell IDs
#' @return Vector of mismatched cell IDs
#' @keywords internal
expressionMismatch <- function(predictions, expression_matrix, common_cells) {
  
  # Prepare expression data
  if (ncol(expression_matrix) > nrow(expression_matrix)) {
    expr_mat <- t(expression_matrix)  # transpose if genes in rows
  } else {
    expr_mat <- expression_matrix
  }
  
  expr_cells <- intersect(common_cells, rownames(expr_mat))
  if (length(expr_cells) < 20) return(character())
  
  expr_data <- expr_mat[expr_cells, , drop = FALSE]
  cell_preds <- predictions[expr_cells]
  
  # Calculate expression similarity within vs across predicted classes
  mismatched_cells <- character()
  
  for (cell in expr_cells) {
    cell_class <- cell_preds[cell]
    same_class_cells <- names(cell_preds)[cell_preds == cell_class & names(cell_preds) != cell]
    
    if (length(same_class_cells) < 5) next  # need enough cells for comparison
    
    # Sample to avoid computation explosion
    same_class_sample <- sample(same_class_cells, min(10, length(same_class_cells)))
    
    # Calculate similarity to same class
    same_class_cors <- sapply(same_class_sample, function(other_cell) {
      cor(expr_data[cell, ], expr_data[other_cell, ], use = "pairwise.complete.obs")
    })
    mean_same_class_cor <- mean(same_class_cors, na.rm = TRUE)
    
    # Calculate similarity to different classes (sample)
    other_class_cells <- names(cell_preds)[cell_preds != cell_class]
    if (length(other_class_cells) >= 5) {
      other_class_sample <- sample(other_class_cells, min(10, length(other_class_cells)))
      other_class_cors <- sapply(other_class_sample, function(other_cell) {
        cor(expr_data[cell, ], expr_data[other_cell, ], use = "pairwise.complete.obs")
      })
      mean_other_class_cor <- mean(other_class_cors, na.rm = TRUE)
      
      # Flag if more similar to other classes than own class
      if (!is.na(mean_same_class_cor) && !is.na(mean_other_class_cor) &&
          mean_other_class_cor > mean_same_class_cor + 0.1) {  # threshold
        mismatched_cells <- c(mismatched_cells, cell)
      }
    }
  }
  
  return(mismatched_cells)
}

# =============================================================================
# 3. CORRECTION GENERATION & EVALUATION
# =============================================================================

#' Generate Corrections with Alternative Validation
#' 
#' @param STEAM.obj STEAM object
#' @param uncertain_cells Uncertain cell IDs
#' @param expression_matrix Expression matrix
#' @param params Current parameters
#' @param baseline_metrics Baseline validation metrics
#' @param iteration Current iteration
#' @param verbose Print messages
#' @return Data frame of corrections
#' @keywords internal
generateCorrections <- function(STEAM.obj, uncertain_cells, 
                                                       expression_matrix, params,
                                                       baseline_metrics, iteration, 
                                                       verbose = FALSE) {
  
  if (length(uncertain_cells) == 0) return(NULL)
  
  # Extract data
  all_predictions <- extractPreds(STEAM.obj)
  coords <- STEAM.obj$spatial
  
  # Prepare expression data if available
  expr_data <- if (!is.null(expression_matrix)) {
    prepareExpression(expression_matrix, coords)
  } else {
    NULL
  }
  
  corrections <- data.frame()
  
  # Process uncertain cells for potential corrections
  for (cell_id in uncertain_cells) {
    if (!(cell_id %in% rownames(coords))) next
    
    correction <- evaluateCorrection(
      cell_id = cell_id,
      all_predictions = all_predictions,
      coordinates = coords,
      expression_data = expr_data,
      params = params,
      baseline_metrics = baseline_metrics,
      iteration = iteration
    )
    
    if (!is.null(correction)) {
      corrections <- rbind(corrections, correction)
    }
  }
  
  # Apply holdout validation on subset of corrections
  if (nrow(corrections) > 5) {
    validated_corrections <- holdoutValidation(
      corrections = corrections,
      STEAM.obj = STEAM.obj,
      expression_matrix = expression_matrix,
      baseline_metrics = baseline_metrics,
      verbose = verbose
    )
    
    return(validated_corrections)
  }
  
  return(corrections)
}

#' Evaluate Cell Correction with Alternative Methods
#' 
#' @param cell_id Cell ID to evaluate
#' @param all_predictions All current predictions
#' @param coordinates Spatial coordinates
#' @param expression_data Expression data (prepared)
#' @param params Current parameters
#' @param baseline_metrics Baseline metrics
#' @param iteration Current iteration
#' @return Correction data frame or NULL
#' @keywords internal
evaluateCorrection <- function(cell_id, all_predictions, coordinates,
                                                expression_data, params, baseline_metrics,
                                                iteration) {
  
  current_pred <- all_predictions[cell_id]
  if (is.na(current_pred)) return(NULL)
  
  # Get spatial neighbors
  cell_coords <- coordinates[cell_id, , drop = FALSE]
  k <- min(15, nrow(coordinates) - 1)
  
  # Simple distance-based neighbors - fix matrix broadcasting
  coord_matrix <- as.matrix(coordinates)
  cell_coord_vector <- as.vector(cell_coords)
  distances <- sqrt(rowSums((coord_matrix - matrix(cell_coord_vector, nrow = nrow(coord_matrix), ncol = ncol(coord_matrix), byrow = TRUE))^2))
  neighbor_indices <- order(distances)[2:(k+1)]  # exclude self
  neighbor_ids <- rownames(coordinates)[neighbor_indices]
  
  # Get neighbor predictions
  neighbor_preds <- all_predictions[neighbor_ids]
  neighbor_preds <- neighbor_preds[!is.na(neighbor_preds)]
  
  if (length(neighbor_preds) < 3) return(NULL)
  
  # Find candidate correction (most common neighbor class)
  neighbor_table <- table(neighbor_preds)
  candidate_class <- names(which.max(neighbor_table))
  spatial_consensus <- max(neighbor_table) / length(neighbor_preds)
  
  # Don't correct if already matches dominant class
  if (candidate_class == current_pred) return(NULL)
  
  # Apply decayed consensus threshold
  consensus_threshold <- params$consensus_multiplier * 0.7  # base threshold
  if (spatial_consensus < consensus_threshold) return(NULL)
  
  # Calculate validation scores
  spatial_improvement <- spatialImprovement(
    cell_id = cell_id,
    current_pred = current_pred,
    candidate_pred = candidate_class,
    coordinates = coordinates,
    all_predictions = all_predictions
  )
  
  expression_improvement <- if (!is.null(expression_data)) {
    expressionImprovement(
      cell_id = cell_id,
      current_pred = current_pred,
      candidate_pred = candidate_class,
      expression_data = expression_data,
      all_predictions = all_predictions
    )
  } else {
    0  # neutral
  }
  
  # Combined improvement score
  total_improvement <- 0.6 * spatial_improvement + 0.4 * expression_improvement
  
  # Apply decayed minimum gain threshold
  min_gain_threshold <- params$min_gain_multiplier * 0.1  # base threshold
  if (total_improvement < min_gain_threshold) return(NULL)
  
  # Create correction record
  return(data.frame(
    cell_id = cell_id,
    original_prediction = current_pred,
    suggested_correction = candidate_class,
    spatial_consensus = spatial_consensus,
    spatial_improvement = spatial_improvement,
    expression_improvement = expression_improvement,
    total_improvement = total_improvement,
    iteration = iteration,
    n_neighbors = length(neighbor_preds),
    validation_method = "alternative",
    stringsAsFactors = FALSE
  ))
}

#' Calculate Spatial Improvement
#' 
#' @param cell_id Cell ID
#' @param current_pred Current prediction
#' @param candidate_pred Candidate prediction
#' @param coordinates Spatial coordinates
#' @param all_predictions All predictions
#' @return Spatial improvement score
#' @keywords internal
spatialImprovement <- function(cell_id, current_pred, candidate_pred,
                                        coordinates, all_predictions) {
  
  # Get neighbors
  cell_coords <- coordinates[cell_id, , drop = FALSE]
  k <- min(15, nrow(coordinates) - 1)
  distances <- sqrt(rowSums((coordinates - rep(cell_coords, each = nrow(coordinates)))^2))
  neighbor_indices <- order(distances)[2:(k+1)]
  neighbor_ids <- rownames(coordinates)[neighbor_indices]
  neighbor_preds <- all_predictions[neighbor_ids]
  neighbor_preds <- neighbor_preds[!is.na(neighbor_preds)]
  
  if (length(neighbor_preds) == 0) return(0)
  
  # Current spatial agreement
  current_agreement <- mean(neighbor_preds == current_pred, na.rm = TRUE)
  
  # Candidate spatial agreement
  candidate_agreement <- mean(neighbor_preds == candidate_pred, na.rm = TRUE)
  
  # Return improvement
  return(candidate_agreement - current_agreement)
}

#' Calculate Expression Improvement
#' 
#' @param cell_id Cell ID
#' @param current_pred Current prediction
#' @param candidate_pred Candidate prediction
#' @param expression_data Expression data
#' @param all_predictions All predictions
#' @return Expression improvement score
#' @keywords internal
expressionImprovement <- function(cell_id, current_pred, candidate_pred,
                                           expression_data, all_predictions) {
  
  if (is.null(expression_data) || !(cell_id %in% expression_data$common_cells)) {
    return(0)
  }
  
  expr_mat <- expression_data$expression
  cell_expr <- expr_mat[cell_id, ]
  
  # Find cells of current and candidate classes
  current_class_cells <- names(all_predictions)[all_predictions == current_pred & 
                                                names(all_predictions) %in% rownames(expr_mat) &
                                                names(all_predictions) != cell_id]
  candidate_class_cells <- names(all_predictions)[all_predictions == candidate_pred & 
                                                  names(all_predictions) %in% rownames(expr_mat)]
  
  if (length(current_class_cells) < 3 || length(candidate_class_cells) < 3) {
    return(0)
  }
  
  # Sample to avoid computation explosion
  current_sample <- sample(current_class_cells, min(10, length(current_class_cells)))
  candidate_sample <- sample(candidate_class_cells, min(10, length(candidate_class_cells)))
  
  # Calculate similarities
  current_cors <- sapply(current_sample, function(other_cell) {
    cor(cell_expr, expr_mat[other_cell, ], use = "pairwise.complete.obs")
  })
  
  candidate_cors <- sapply(candidate_sample, function(other_cell) {
    cor(cell_expr, expr_mat[other_cell, ], use = "pairwise.complete.obs")
  })
  
  current_mean_cor <- mean(current_cors, na.rm = TRUE)
  candidate_mean_cor <- mean(candidate_cors, na.rm = TRUE)
  
  if (is.na(current_mean_cor) || is.na(candidate_mean_cor)) return(0)
  
  # Return improvement (candidate should be more similar)
  return(candidate_mean_cor - current_mean_cor)
}

# =============================================================================
# 4. VALIDATION & APPLICATION
# =============================================================================

#' Apply Holdout Validation
#' 
#' @param corrections Proposed corrections
#' @param STEAM.obj STEAM object
#' @param expression_matrix Expression matrix
#' @param baseline_metrics Baseline metrics
#' @param verbose Print messages
#' @return Validated corrections
#' @keywords internal
holdoutValidation <- function(corrections, STEAM.obj, expression_matrix,
                                   baseline_metrics, verbose = FALSE) {
  
  if (nrow(corrections) < 5) return(corrections)
  
  # Use 20% of corrections for validation
  n_validation <- max(1, round(0.2 * nrow(corrections)))
  validation_indices <- sample(nrow(corrections), n_validation)
  validation_corrections <- corrections[validation_indices, ]
  remaining_corrections <- corrections[-validation_indices, ]
  
  # Apply validation corrections temporarily
  temp_predictions <- extractPreds(STEAM.obj)
  temp_predictions[validation_corrections$cell_id] <- validation_corrections$suggested_correction
  
  # Calculate validation metrics
  validation_metrics <- validationMetrics(
    predictions = temp_predictions,
    coordinates = STEAM.obj$spatial,
    expression_matrix = expression_matrix,
    baseline_metrics = baseline_metrics
  )
  
  # Check if validation improves metrics
  improvement_score <- calcImprovement(baseline_metrics, validation_metrics)
  
  if (verbose) {
    message(sprintf("Holdout validation: tested %d corrections, improvement: %+.4f",
                   n_validation, improvement_score))
  }
  
  # If validation successful, return all corrections; otherwise be more conservative
  if (improvement_score > 0) {
    return(corrections)
  } else {
    # Return only the most confident corrections
    confidence_threshold <- quantile(remaining_corrections$total_improvement, 0.75)
    return(remaining_corrections[remaining_corrections$total_improvement >= confidence_threshold, ])
  }
}

#' Compute Validation Metrics Temporarily
#' 
#' @param predictions Temporary predictions
#' @param coordinates Spatial coordinates
#' @param expression_matrix Expression matrix
#' @param baseline_metrics Baseline metrics
#' @return Validation metrics
#' @keywords internal
validationMetrics <- function(predictions, coordinates, expression_matrix,
                                          baseline_metrics) {
  
  common_cells <- intersect(names(predictions), rownames(coordinates))
  
  spatial_coherence <- spatialCoherence(
    predictions = predictions[common_cells],
    coordinates = coordinates[common_cells, , drop = FALSE],
    k = 15
  )
  
  expression_coherence <- if (!is.null(expression_matrix)) {
    expressionCoherence(
      predictions = predictions[common_cells],
      expression_matrix = expression_matrix,
      common_cells = common_cells,
      genes_transpose = TRUE
    )
  } else {
    baseline_metrics$expression_coherence
  }
  
  global_entropy <- predictionEntropy(predictions[common_cells])
  
  return(list(
    spatial_coherence = spatial_coherence,
    expression_coherence = expression_coherence,
    global_entropy = global_entropy,
    n_cells = length(common_cells)
  ))
}

#' Apply Alternative Corrections
#' 
#' @param STEAM.obj STEAM object
#' @param corrections Corrections to apply
#' @param iteration Current iteration
#' @return Updated STEAM object
#' @keywords internal
applyCorrections <- function(STEAM.obj, corrections, iteration) {
  
  if (is.null(corrections) || nrow(corrections) == 0) {
    return(STEAM.obj)
  }
  
  # Store corrections in the STEAM object
  if (is.null(STEAM.obj$spatial_anchor_analysis$corrections)) {
    STEAM.obj$spatial_anchor_analysis$corrections <- corrections
  } else {
    STEAM.obj$spatial_anchor_analysis$corrections <- rbind(
      STEAM.obj$spatial_anchor_analysis$corrections,
      corrections
    )
  }
  
  # Update prediction tracking
  corrected_cells <- corrections$cell_id
  new_predictions <- corrections$suggested_correction
  names(new_predictions) <- corrected_cells
  
  # Store in temporary prediction update structure
  if (is.null(STEAM.obj$temp_prediction_updates)) {
    STEAM.obj$temp_prediction_updates <- list()
  }
  STEAM.obj$temp_prediction_updates[[iteration]] <- new_predictions
  
  return(STEAM.obj)
}

#' Compute Updated Validation Metrics
#' 
#' @param STEAM.obj STEAM object
#' @param expression_matrix Expression matrix
#' @param corrections Recent corrections
#' @param baseline_metrics Baseline metrics
#' @return Updated validation metrics
#' @keywords internal
updateMetrics <- function(STEAM.obj, expression_matrix, 
                                             corrections, baseline_metrics) {
  
  # Get current predictions including all corrections
  current_predictions <- extractPreds(STEAM.obj)
  
  # Apply all corrections made so far
  if (!is.null(STEAM.obj$spatial_anchor_analysis$corrections)) {
    all_corrections <- STEAM.obj$spatial_anchor_analysis$corrections
    current_predictions[all_corrections$cell_id] <- all_corrections$suggested_correction
  }
  
  # Calculate updated metrics
  return(validationMetrics(
    predictions = current_predictions,
    coordinates = STEAM.obj$spatial,
    expression_matrix = expression_matrix,
    baseline_metrics = baseline_metrics
  ))
}

#' Calculate Alternative Improvement
#' 
#' @param prev_metrics Previous metrics
#' @param current_metrics Current metrics
#' @return Overall improvement score
#' @keywords internal
calcImprovement <- function(prev_metrics, current_metrics) {
  
  # Weight different improvements
  spatial_weight <- 0.5
  expression_weight <- 0.3
  entropy_weight <- 0.2
  
  spatial_improvement <- current_metrics$spatial_coherence - prev_metrics$spatial_coherence
  expression_improvement <- current_metrics$expression_coherence - prev_metrics$expression_coherence
  entropy_improvement <- prev_metrics$global_entropy - current_metrics$global_entropy  # lower entropy is better
  
  total_improvement <- (
    spatial_weight * spatial_improvement +
    expression_weight * expression_improvement +
    entropy_weight * entropy_improvement
  )
  
  return(total_improvement)
}

# =============================================================================
# 5. UTILITY FUNCTIONS
# =============================================================================

#' Extract All Predictions from STEAM Object
#' 
#' @param STEAM.obj STEAM object
#' @return Named vector of all predictions
#' @keywords internal
extractPreds <- function(STEAM.obj) {
  
  all_predictions <- character()
  
  # Extract from nested CV results
  for (f in seq_along(STEAM.obj$nested$ncv$outer_result)) {
    fold_result <- STEAM.obj$nested$ncv$outer_result[[f]]
    if (!is.null(fold_result$preds)) {
      preds_df <- fold_result$preds
      fold_preds <- setNames(preds_df$predy, rownames(preds_df))
      all_predictions <- c(all_predictions, fold_preds)
    }
  }
  
  # Remove duplicates (keep first occurrence)
  all_predictions <- all_predictions[!duplicated(names(all_predictions))]
  
  return(all_predictions)
}

#' Build K-Nearest Neighbors Index
#' 
#' @param coordinates Spatial coordinates matrix
#' @param k Number of neighbors
#' @return List with neighbor indices
#' @keywords internal
buildKnn <- function(coordinates, k) {
  
  n_cells <- nrow(coordinates)
  coord_matrix <- as.matrix(coordinates)
  
  # Initialize result list
  idx_list <- vector("list", n_cells)
  
  # Calculate neighbors for each cell
  for (i in seq_len(n_cells)) {
    cell_coord <- coord_matrix[i, , drop = FALSE]
    
    # Calculate distances to all other cells
    distances <- sqrt(rowSums((coord_matrix - matrix(rep(cell_coord, n_cells), 
                                                   nrow = n_cells, byrow = TRUE))^2))
    
    # Get k nearest neighbors (excluding self)
    neighbor_indices <- order(distances)[2:(min(k + 1, n_cells))]
    idx_list[[i]] <- neighbor_indices
  }
  
  return(list(idx_list = idx_list))
}

#' Prepare Expression Data
#' 
#' @param expression_matrix Expression matrix
#' @param coordinates Spatial coordinates
#' @return List with prepared expression data
#' @keywords internal
prepareExpression <- function(expression_matrix, coordinates) {
  
  # Determine matrix orientation
  if (ncol(expression_matrix) > nrow(expression_matrix)) {
    expr_mat <- t(expression_matrix)  # transpose if genes in rows
  } else {
    expr_mat <- expression_matrix
  }
  
  # Find common cells
  common_cells <- intersect(rownames(expr_mat), rownames(coordinates))
  
  if (length(common_cells) < 10) {
    warning("Very few common cells between expression and spatial data")
    return(NULL)
  }
  
  return(list(
    expression = expr_mat[common_cells, , drop = FALSE],
    common_cells = common_cells,
    n_genes = ncol(expr_mat),
    n_cells = length(common_cells)
  ))
}
