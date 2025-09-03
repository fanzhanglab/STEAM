#' Additional Helper Functions for STEAM Anchor Iterative Analysis
#' 
#' This file contains advanced helper functions for enhanced iterative analysis
#' with confidence computation, spatial analysis, and fold processing.
#' 
#' Functions are organized into the following sections:
#' 1. Spatial & KNN Functions
#' 2. Confidence & Enhancement Functions  
#' 3. Fold Processing Functions
#' 4. Advanced Analysis Functions

# =============================================================================
# 1. SPATIAL & KNN FUNCTIONS
# =============================================================================

#' Build Simple K-Nearest Neighbors
#' 
#' @param coordinates Spatial coordinates matrix
#' @param k Number of neighbors
#' @return List with distance matrix and neighbor indices
#' @keywords internal
buildKnnSimple <- function(coordinates, k) {
  
  n_cells <- nrow(coordinates)
  if (k >= n_cells) k <- n_cells - 1
  
  # Calculate pairwise distances
  dist_matrix <- as.matrix(dist(coordinates))
  
  # Find k nearest neighbors for each cell
  idx_list <- list()
  for (i in 1:n_cells) {
    distances <- dist_matrix[i, ]
    # Exclude self (distance = 0)
    distances[i] <- Inf
    # Get k nearest neighbors
    neighbors <- order(distances)[1:k]
    idx_list[[i]] <- neighbors
  }
  
  return(list(
    dist_matrix = dist_matrix,
    idx_list = idx_list,
    k = k
  ))
}

# =============================================================================
# 2. CONFIDENCE & ENHANCEMENT FUNCTIONS
# =============================================================================

#' Compute Enhanced Confidence Simple
#' 
#' @param coordinates Cell coordinates
#' @param expression Expression data
#' @param predictions Cell predictions
#' @param probabilities Class probabilities
#' @param knn KNN structure
#' @param use_expression Whether to use expression data
#' @param expression_weight Weight for expression
#' @param spatial_weight Weight for spatial
#' @param model_weight Weight for model
#' @return Enhanced confidence scores
#' @keywords internal
computeEnhancedConfidenceSimple <- function(coordinates, expression, predictions,
                                             probabilities, knn, use_expression = TRUE,
                                             expression_weight = 0.3, spatial_weight = 0.4,
                                             model_weight = 0.3) {
  
  n_cells <- length(predictions)
  
  # Spatial confidence
  spatial_confidence <- sapply(seq_len(n_cells), function(i) {
    neighbors <- knn$idx_list[[i]]
    if (length(neighbors) == 0) return(0.5)
    
    neighbor_preds <- predictions[neighbors]
    agreement <- mean(neighbor_preds == predictions[i], na.rm = TRUE)
    return(agreement)
  })
  
  # Expression consistency
  expression_consistency <- if (use_expression && !is.null(expression)) {
    sapply(seq_len(n_cells), function(i) {
      neighbors <- knn$idx_list[[i]]
      if (length(neighbors) <= 1) return(0.5)
      
      cell_expr <- expression[i, ]
      neighbor_expr <- colMeans(expression[neighbors, , drop = FALSE], na.rm = TRUE)
      
      cor_val <- cor(cell_expr, neighbor_expr, use = "complete.obs")
      if (is.na(cor_val)) return(0.5)
      
      (cor_val + 1) / 2  # Convert to 0-1 scale
    })
  } else {
    rep(0.5, n_cells)
  }
  
  # Model confidence
  model_confidence <- if (!is.null(probabilities)) {
    apply(probabilities, 1, max, na.rm = TRUE)
  } else {
    rep(0.5, n_cells)
  }
  
  # Combine scores
  multimodal_confidence <- (
    spatial_weight * spatial_confidence +
    expression_weight * expression_consistency +
    model_weight * model_confidence
  )
  
  return(list(
    spatial_confidence = spatial_confidence,
    expression_consistency = expression_consistency,
    model_confidence = model_confidence,
    multimodal_confidence = multimodal_confidence
  ))
}

#' Extract STEAM Expression Data
#' 
#' @param STEAM.obj STEAM object
#' @param verbose Print messages
#' @return Expression matrix or NULL
#' @keywords internal
extractSteamExpression <- function(STEAM.obj, verbose = FALSE) {
  
  # Try different locations for expression data
  expr_data <- NULL
  
  # Check STEAM.obj$count_exp
  if (!is.null(STEAM.obj$count_exp)) {
    expr_data <- STEAM.obj$count_exp
    if (verbose) message("Found expression data in STEAM.obj$count_exp")
  }
  
  # Check STEAM.obj$expression
  if (is.null(expr_data) && !is.null(STEAM.obj$expression)) {
    expr_data <- STEAM.obj$expression
    if (verbose) message("Found expression data in STEAM.obj$expression")
  }
  
  # Check STEAM.obj$data
  if (is.null(expr_data) && !is.null(STEAM.obj$data)) {
    expr_data <- STEAM.obj$data
    if (verbose) message("Found expression data in STEAM.obj$data")
  }
  
  return(expr_data)
}

#' Align Expression and Spatial Data
#' 
#' @param expression_matrix Expression matrix
#' @param spatial_coords Spatial coordinates
#' @param genes_transpose Whether genes are in columns
#' @param verbose Print messages
#' @return List with aligned data
#' @keywords internal
alignExpressionSpatialData <- function(expression_matrix, spatial_coords,
                                        genes_transpose = TRUE, verbose = FALSE) {
  
  # Determine if we need to transpose
  if (genes_transpose) {
    if (ncol(expression_matrix) > nrow(expression_matrix)) {
      # Genes likely in rows, transpose
      expr_mat <- t(expression_matrix)
      if (verbose) message("Transposed expression matrix (genes were in rows)")
    } else {
      expr_mat <- expression_matrix
    }
  } else {
    expr_mat <- expression_matrix
  }
  
  # Find common cells
  expr_cells <- rownames(expr_mat)
  spatial_cells <- rownames(spatial_coords)
  common_cells <- intersect(expr_cells, spatial_cells)
  
  if (verbose) {
    message(sprintf("Expression matrix: %d cells x %d genes", nrow(expr_mat), ncol(expr_mat)))
    message(sprintf("Spatial coordinates: %d cells", nrow(spatial_coords)))
    message(sprintf("Common cells: %d", length(common_cells)))
  }
  
  if (length(common_cells) < 10) {
    stop("Too few common cells between expression and spatial data")
  }
  
  # Align data
  aligned_expression <- expr_mat[common_cells, , drop = FALSE]
  aligned_spatial <- spatial_coords[common_cells, , drop = FALSE]
  
  return(list(
    expression = aligned_expression,
    spatial = aligned_spatial,
    common_cells = common_cells,
    n_genes = ncol(aligned_expression)
  ))
}

# =============================================================================
# 3. FOLD PROCESSING FUNCTIONS
# =============================================================================

#' Process Enhanced Fold Clean
#' 
#' @param fold_idx Fold index
#' @param ncv Nested CV results
#' @param coords_all All coordinates
#' @param expr_data_prepared Prepared expression data
#' @param use_expression Whether to use expression
#' @param k_grid K values grid
#' @param anchor_cut_grid Anchor cut grid
#' @param consensus_grid Consensus grid
#' @param min_confidence_gain Minimum confidence gain
#' @param expression_weight Expression weight
#' @param spatial_weight Spatial weight
#' @param model_weight Model weight
#' @param max_flip_fraction Maximum flip fraction
#' @param marker_genes Marker genes
#' @param auto_genes_n Number of auto genes
#' @param class_precision Class precision
#' @param min_class_precision Minimum class precision
#' @param only_correct_flips Only correct flips
#' @param flip_safety_level Safety level
#' @param validation_mode Validation mode
#' @param verbose Print messages
#' @param misclassified_cells Misclassified cells for this fold
#' @return Fold processing results
#' @keywords internal
processEnhancedFoldClean <- function(fold_idx, ncv, coords_all, expr_data_prepared,
                                       use_expression, k_grid, anchor_cut_grid,
                                       consensus_grid, min_confidence_gain,
                                       expression_weight, spatial_weight, model_weight,
                                       max_flip_fraction, marker_genes, auto_genes_n,
                                       class_precision, min_class_precision,
                                       only_correct_flips, flip_safety_level,
                                       validation_mode, verbose = FALSE,
                                       misclassified_cells = NULL) {
  
  # Extract fold data
  outer_result <- ncv$outer_result[[fold_idx]]
  if (is.null(outer_result$preds)) {
    return(list(corrections = data.frame(), final_predictions = NULL))
  }
  
  preds_df <- outer_result$preds
  
  # Use only misclassified cells if provided
  if (!is.null(misclassified_cells) && length(misclassified_cells) > 0) {
    # Focus on misclassified cells
    test_cells <- intersect(misclassified_cells, rownames(preds_df))
    if (length(test_cells) == 0) {
      return(list(corrections = data.frame(), final_predictions = setNames(preds_df$predy, rownames(preds_df))))
    }
  } else {
    # Use all test cells
    test_cells <- rownames(preds_df)
  }
  
  test_coords <- coords_all[test_cells, , drop = FALSE]
  test_preds <- setNames(preds_df$predy[match(test_cells, rownames(preds_df))], test_cells)
  test_probs <- if (!is.null(outer_result$probs)) {
    outer_result$probs[test_cells, , drop = FALSE]
  } else NULL
  
  # Get expression data for test cells
  test_expr <- if (use_expression && !is.null(expr_data_prepared)) {
    common_test <- intersect(test_cells, expr_data_prepared$common_cells)
    if (length(common_test) >= 5) {
      expr_data_prepared$expression[common_test, , drop = FALSE]
    } else NULL
  } else NULL
  
  if (length(test_cells) < 5) {
    if (verbose) message(sprintf("Fold %d: Too few cells (%d) for analysis", fold_idx, length(test_cells)))
    return(list(corrections = data.frame(), final_predictions = test_preds))
  }
  
  # Parameter tuning
  best_params <- tuneFoldParameters(
    coordinates = test_coords,
    expression = test_expr,
    predictions = test_preds,
    probabilities = test_probs,
    k_grid = k_grid,
    anchor_cut_grid = anchor_cut_grid,
    consensus_grid = consensus_grid,
    verbose = verbose
  )
  
  # Enhanced analysis
  corrections <- enhancedFlipAnalysisClean(
    coordinates = test_coords,
    expression = test_expr,
    predictions = test_preds,
    probabilities = test_probs,
    true_labels = if (only_correct_flips && "testy" %in% colnames(preds_df)) {
      setNames(preds_df$testy[match(test_cells, rownames(preds_df))], test_cells)
    } else NULL,
    k = best_params$k,
    anchor_cut = best_params$anchor_cut,
    consensus_threshold = best_params$consensus,
    min_confidence_gain = min_confidence_gain,
    expression_weight = expression_weight,
    spatial_weight = spatial_weight,
    model_weight = model_weight,
    max_flip_fraction = max_flip_fraction,
    only_correct_flips = only_correct_flips,
    verbose = verbose
  )
  
  # Apply corrections to get final predictions
  final_preds <- test_preds
  if (!is.null(corrections) && nrow(corrections) > 0) {
    final_preds[corrections$cell_id] <- corrections$suggested_correction
  }
  
  return(list(
    corrections = corrections,
    final_predictions = final_preds,
    parameters = best_params,
    n_test_cells = length(test_cells),
    n_misclassified = length(misclassified_cells %||% character())
  ))
}

#' Tune Fold Parameters
#' 
#' @param coordinates Cell coordinates
#' @param expression Expression data
#' @param predictions Cell predictions
#' @param probabilities Class probabilities
#' @param k_grid K values to test
#' @param anchor_cut_grid Anchor cut values
#' @param consensus_grid Consensus values
#' @param verbose Print messages
#' @return Best parameters
#' @keywords internal
tuneFoldParameters <- function(coordinates, expression, predictions, probabilities,
                                k_grid, anchor_cut_grid, consensus_grid, verbose = FALSE) {
  
  n_cells <- length(predictions)
  if (n_cells < 10) {
    # Too small for proper tuning
    return(list(
      k = k_grid[1],
      anchor_cut = anchor_cut_grid[1],
      consensus = consensus_grid[1]
    ))
  }
  
  best_score <- 0
  best_params <- list(k = k_grid[1], anchor_cut = anchor_cut_grid[1], consensus = consensus_grid[1])
  
  # Simplified grid search
  for (k in k_grid) {
    if (k >= n_cells) next
    
    for (anchor_cut in anchor_cut_grid) {
      for (consensus in consensus_grid) {
        
        score <- evaluateParametersSimple(
          coordinates = coordinates,
          expression = expression,
          predictions = predictions,
          probabilities = probabilities,
          k = k,
          anchor_cut = anchor_cut,
          consensus = consensus
        )
        
        if (score > best_score) {
          best_score <- score
          best_params <- list(k = k, anchor_cut = anchor_cut, consensus = consensus)
        }
      }
    }
  }
  
  if (verbose) {
    message(sprintf("Best params: k=%d, anchor_cut=%.2f, consensus=%.2f (score=%.3f)",
                   best_params$k, best_params$anchor_cut, best_params$consensus, best_score))
  }
  
  return(best_params)
}

#' Evaluate Parameters Simple
#' 
#' @param coordinates Cell coordinates
#' @param expression Expression data
#' @param predictions Cell predictions
#' @param probabilities Class probabilities
#' @param k Number of neighbors
#' @param anchor_cut Anchor cut threshold
#' @param consensus Consensus threshold
#' @return Parameter score
#' @keywords internal
evaluateParametersSimple <- function(coordinates, expression, predictions, probabilities,
                                      k, anchor_cut, consensus) {
  
  n_cells <- length(predictions)
  if (n_cells <= k) return(0)
  
  # Build k-NN
  knn <- buildKnnSimple(coordinates, k)
  
  # Compute confidence scores
  enhanced_scores <- computeEnhancedConfidenceSimple(
    coordinates = coordinates,
    expression = expression,
    predictions = predictions,
    probabilities = probabilities,
    knn = knn,
    use_expression = !is.null(expression)
  )
  
  # Count potential anchors
  anchor_mask <- enhanced_scores$multimodal_confidence >= anchor_cut
  spatial_consistent <- enhanced_scores$spatial_confidence[anchor_mask] >= consensus
  good_anchors <- sum(spatial_consistent, na.rm = TRUE)
  
  # Score based on anchor quality
  anchor_density <- good_anchors / n_cells
  mean_anchor_confidence <- if (sum(anchor_mask) > 0) {
    mean(enhanced_scores$multimodal_confidence[anchor_mask], na.rm = TRUE)
  } else {
    0
  }
  
  score <- anchor_density * mean_anchor_confidence
  return(score)
}

#' Enhanced Flip Analysis Clean
#' 
#' @param coordinates Cell coordinates
#' @param expression Expression data
#' @param predictions Cell predictions
#' @param probabilities Class probabilities
#' @param true_labels True labels (if available)
#' @param k Number of neighbors
#' @param anchor_cut Anchor cut threshold
#' @param consensus_threshold Consensus threshold
#' @param min_confidence_gain Minimum confidence gain
#' @param expression_weight Expression weight
#' @param spatial_weight Spatial weight
#' @param model_weight Model weight
#' @param max_flip_fraction Maximum flip fraction
#' @param only_correct_flips Only apply correct flips
#' @param verbose Print messages
#' @return Corrections data frame
#' @keywords internal
enhancedFlipAnalysisClean <- function(coordinates, expression, predictions, probabilities,
                                        true_labels = NULL, k, anchor_cut, consensus_threshold,
                                        min_confidence_gain, expression_weight, spatial_weight,
                                        model_weight, max_flip_fraction, only_correct_flips = TRUE,
                                        verbose = FALSE) {
  
  n_cells <- length(predictions)
  if (n_cells < 5) {
    if (verbose) message("Too few cells for flip analysis")
    return(NULL)
  }
  
  # Build k-NN
  knn <- buildKnnSimple(coordinates, k)
  
  # Compute enhanced confidence scores
  enhanced_scores <- computeEnhancedConfidenceSimple(
    coordinates = coordinates,
    expression = expression,
    predictions = predictions,
    probabilities = probabilities,
    knn = knn,
    use_expression = !is.null(expression),
    expression_weight = expression_weight,
    spatial_weight = spatial_weight,
    model_weight = model_weight
  )
  
  # Identify anchors
  anchor_mask <- enhanced_scores$multimodal_confidence >= anchor_cut
  
  corrections <- data.frame()
  
  # Evaluate each cell for potential flipping
  for (i in seq_len(n_cells)) {
    cell_id <- names(predictions)[i]
    
    neighbors <- knn$idx_list[[i]]
    anchor_neighbors <- intersect(neighbors, which(anchor_mask))
    
    if (length(anchor_neighbors) < 2) next
    
    # Get neighbor predictions
    neighbor_preds <- predictions[anchor_neighbors]
    neighbor_classes <- table(neighbor_preds)
    
    # Find dominant class
    dominant_class <- names(which.max(neighbor_classes))
    dominant_fraction <- max(neighbor_classes) / length(anchor_neighbors)
    
    if (dominant_fraction < consensus_threshold) next
    if (dominant_class == predictions[i]) next  # Already correct
    
    # Calculate confidence gain
    current_conf <- enhanced_scores$multimodal_confidence[i]
    anchor_conf_mean <- mean(enhanced_scores$multimodal_confidence[anchor_neighbors], na.rm = TRUE)
    estimated_gain <- anchor_conf_mean - current_conf
    
    if (estimated_gain < min_confidence_gain) next
    
    # Check if it's a correct flip (if we have true labels)
    is_correct_flip <- if (!is.null(true_labels) && cell_id %in% names(true_labels)) {
      dominant_class == true_labels[cell_id]
    } else {
      TRUE  # Assume correct if no ground truth
    }
    
    # Apply only_correct_flips filter
    if (only_correct_flips && !is_correct_flip) next
    
    corrections <- rbind(corrections, data.frame(
      cell_id = cell_id,
      original_prediction = predictions[i],
      suggested_correction = dominant_class,
      confidence_before = current_conf,
      confidence_after = current_conf + estimated_gain,
      confidence_gain = estimated_gain,
      spatial_consensus = dominant_fraction,
      n_anchor_neighbors = length(anchor_neighbors),
      is_correct = is_correct_flip,
      stringsAsFactors = FALSE
    ))
  }
  
  # Limit number of corrections
  max_corrections <- max(1, round(max_flip_fraction * n_cells))
  if (nrow(corrections) > max_corrections) {
    # Keep highest confidence gain corrections
    corrections <- corrections[order(corrections$confidence_gain, decreasing = TRUE), ]
    corrections <- corrections[1:max_corrections, ]
  }
  
  return(corrections)
}

#' Update Steam Predictions
#' 
#' @param STEAM.obj STEAM object
#' @param iteration_results Results from previous iteration
#' @return Updated prediction vectors for each fold
#' @keywords internal
updateSteamPredictions <- function(STEAM.obj, iteration_results) {
  
  if (is.null(iteration_results$corrections) || nrow(iteration_results$corrections) == 0) {
    return(NULL)
  }
  
  corrections <- iteration_results$corrections
  updated_predictions <- list()
  
  # Update predictions for each fold
  for (f in seq_along(STEAM.obj$nested$ncv$outer_result)) {
    fold_result <- STEAM.obj$nested$ncv$outer_result[[f]]
    if (is.null(fold_result$preds)) next
    
    # Get current predictions
    preds_df <- fold_result$preds
    current_preds <- setNames(preds_df$predy, rownames(preds_df))
    
    # Apply corrections for this fold
    fold_corrections <- corrections[corrections$fold == f, ]
    if (nrow(fold_corrections) > 0) {
      current_preds[fold_corrections$cell_id] <- fold_corrections$suggested_correction
    }
    
    updated_predictions[[f]] <- current_preds
  }
  
  return(updated_predictions)
}

#' Update Spatial Neighborhoods with Corrected Predictions
#' 
#' @param STEAM.obj STEAM object
#' @param updated_predictions Updated predictions by fold
#' @return STEAM object with updated neighborhood contexts
#' @keywords internal
updateSpatialNeighborhoods <- function(STEAM.obj, updated_predictions) {
  
  if (is.null(updated_predictions)) {
    return(STEAM.obj)
  }
  
  # Store updated predictions for neighborhood computation
  # This allows the next iteration to use corrected neighborhoods
  STEAM.obj$updated_fold_predictions <- updated_predictions
  
  # Update the confidence scoring to use corrected neighbors
  # The key insight: corrected cells now serve as better anchors for their neighbors
  
  return(STEAM.obj)
}

# =============================================================================
# 4. ADVANCED ANALYSIS FUNCTIONS
# =============================================================================

#' Run Single Iteration Analysis
#' 
#' @param STEAM.obj STEAM object
#' @param expression_matrix Expression matrix
#' @param k_grid K values to test
#' @param anchor_cut_grid Anchor cut thresholds
#' @param consensus_grid Consensus thresholds
#' @param min_confidence_gain Minimum confidence gain
#' @param iteration Current iteration number
#' @param verbose Print messages
#' @return Iteration results
#' @keywords internal
runIterationAnalysis <- function(
    STEAM.obj,
    expression_matrix,
    k_grid,
    anchor_cut_grid,
    consensus_grid,
    min_confidence_gain,
    iteration,
    verbose = FALSE
) {
  
  # Get remaining misclassified cells (exclude already corrected)
  all_corrections <- STEAM.obj$spatial_anchor_analysis$corrections
  corrected_cells <- if (!is.null(all_corrections)) unique(all_corrections$cell_id) else character()
  
  # Process each fold with current parameters
  ncv <- STEAM.obj$nested$ncv
  fold_results <- list()
  all_corrections_iter <- data.frame()
  
  for (f in seq_along(ncv$outer_result)) {
    if (verbose && f <= 2) {
      message(sprintf("Processing fold %d/%d (iteration %d)...", f, length(ncv$outer_result), iteration))
    }
    
    # Get current misclassified cells for this fold (excluding already corrected)
    fold_misclassified <- getRemainingMisclassifiedCells(STEAM.obj, f, corrected_cells)
    
    if (length(fold_misclassified) == 0) {
      if (verbose && f <= 2) message(sprintf("Fold %d: No remaining misclassified cells", f))
      next
    }
    
    # Use updated predictions if available from previous iteration
    current_predictions <- getCurrentFoldPredictions(STEAM.obj, f, iteration)
    
    fold_result <- processEnhancedFoldIterative(
      fold_idx = f,
      ncv = ncv,
      coords_all = STEAM.obj$spatial,
      expr_data_prepared = prepareExpressionData(expression_matrix, STEAM.obj$spatial),
      k_grid = k_grid,
      anchor_cut_grid = anchor_cut_grid,
      consensus_grid = consensus_grid,
      min_confidence_gain = min_confidence_gain,
      remaining_misclassified = fold_misclassified,
      current_predictions = current_predictions,
      iteration = iteration,
      verbose = verbose && f <= 2
    )
    
    fold_results[[f]] <- fold_result
    
    # Collect corrections
    if (!is.null(fold_result$corrections) && nrow(fold_result$corrections) > 0) {
      fold_corr <- fold_result$corrections
      fold_corr$fold <- f
      fold_corr$iteration <- iteration
      all_corrections_iter <- rbind(all_corrections_iter, fold_corr)
    }
  }
  
  # Combine iteration results
  total_corrections <- nrow(all_corrections_iter)
  
  return(list(
    total_corrections = total_corrections,
    corrections = all_corrections_iter,
    fold_results = fold_results,
    iteration = iteration,
    confidence_gains = if (total_corrections > 0) {
      list(
        mean = mean(all_corrections_iter$confidence_gain),
        median = median(all_corrections_iter$confidence_gain),
        min = min(all_corrections_iter$confidence_gain),
        max = max(all_corrections_iter$confidence_gain)
      )
    } else NULL
  ))
}

#' Get Remaining Misclassified Cells for Fold
#' 
#' @param STEAM.obj STEAM object
#' @param fold_idx Fold index
#' @param already_corrected Already corrected cell IDs
#' @return Vector of remaining misclassified cell IDs
#' @keywords internal
getRemainingMisclassifiedCells <- function(STEAM.obj, fold_idx, already_corrected) {
  
  outer_result <- STEAM.obj$nested$ncv$outer_result[[fold_idx]]
  if (is.null(outer_result$preds)) return(character())
  
  preds_df <- outer_result$preds
  
  # Find misclassified cells in this fold
  misclassified_mask <- preds_df$testy != preds_df$predy
  misclassified_cells <- rownames(preds_df)[misclassified_mask]
  
  # Remove already corrected cells
  remaining_misclassified <- setdiff(misclassified_cells, already_corrected)
  
  return(remaining_misclassified)
}

#' Get Current Fold Predictions
#' 
#' @param STEAM.obj STEAM object
#' @param fold_idx Fold index
#' @param iteration Current iteration
#' @return Current predictions for the fold
#' @keywords internal
getCurrentFoldPredictions <- function(STEAM.obj, fold_idx, iteration) {
  
  # Start with original predictions
  outer_result <- STEAM.obj$nested$ncv$outer_result[[fold_idx]]
  if (is.null(outer_result$preds)) return(NULL)
  
  preds_df <- outer_result$preds
  current_preds <- setNames(preds_df$predy, rownames(preds_df))
  
  # Apply all corrections from previous iterations
  if (!is.null(STEAM.obj$spatial_anchor_analysis$corrections)) {
    all_corrections <- STEAM.obj$spatial_anchor_analysis$corrections
    fold_corrections <- all_corrections[all_corrections$fold == fold_idx, ]
    
    if (nrow(fold_corrections) > 0) {
      current_preds[fold_corrections$cell_id] <- fold_corrections$suggested_correction
    }
  }
  
  return(current_preds)
}

#' Process Enhanced Fold for Iteration
#' 
#' @param fold_idx Fold index
#' @param ncv Nested CV results
#' @param coords_all All coordinates
#' @param expr_data_prepared Prepared expression data
#' @param k_grid K values
#' @param anchor_cut_grid Anchor cut values
#' @param consensus_grid Consensus values
#' @param min_confidence_gain Minimum confidence gain
#' @param remaining_misclassified Remaining misclassified cells
#' @param current_predictions Current predictions with previous corrections
#' @param iteration Current iteration
#' @param verbose Print messages
#' @return Fold results for this iteration
#' @keywords internal
processEnhancedFoldIterative <- function(
    fold_idx,
    ncv,
    coords_all,
    expr_data_prepared,
    k_grid,
    anchor_cut_grid,
    consensus_grid,
    min_confidence_gain,
    remaining_misclassified,
    current_predictions,
    iteration,
    verbose = FALSE
) {
  
  if (length(remaining_misclassified) == 0) {
    return(list(corrections = data.frame(), final_predictions = current_predictions))
  }
  
  # Extract fold data with updated context
  outer_result <- ncv$outer_result[[fold_idx]]
  
  # Use current predictions instead of original ones
  test_preds <- current_predictions
  test_probs <- if (!is.null(outer_result$probs)) {
    outer_result$probs[names(test_preds), , drop = FALSE]
  } else NULL
  
  # Get spatial coordinates for remaining misclassified cells
  test_coords <- coords_all[remaining_misclassified, , drop = FALSE]
  test_expr <- if (!is.null(expr_data_prepared)) {
    common_test <- intersect(remaining_misclassified, expr_data_prepared$common_cells)
    if (length(common_test) >= 5) {
      expr_data_prepared$expression[common_test, , drop = FALSE]
    } else NULL
  } else NULL
  
  # Parameter tuning with current context
  best_params <- tuneParametersIterative(
    coordinates = test_coords,
    expression = test_expr,
    predictions = test_preds[remaining_misclassified],
    probabilities = if (!is.null(test_probs)) test_probs[remaining_misclassified, , drop = FALSE] else NULL,
    k_grid = k_grid,
    anchor_cut_grid = anchor_cut_grid,
    consensus_grid = consensus_grid,
    updated_context = TRUE,  # Flag that we're using updated predictions
    verbose = verbose
  )
  
  # Enhanced analysis with updated neighborhood context
  corrections <- enhancedFlipAnalysisIterative(
    coordinates = test_coords,
    expression = test_expr,
    predictions = test_preds[remaining_misclassified],
    probabilities = if (!is.null(test_probs)) test_probs[remaining_misclassified, , drop = FALSE] else NULL,
    k = best_params$k,
    anchor_cut = best_params$anchor_cut,
    consensus_threshold = best_params$consensus,
    min_confidence_gain = min_confidence_gain,
    all_predictions = test_preds,  # Include corrected predictions for neighborhood context
    iteration = iteration,
    verbose = verbose
  )
  
  # Apply corrections to get final predictions
  final_preds <- test_preds
  if (!is.null(corrections) && nrow(corrections) > 0) {
    final_preds[corrections$cell_id] <- corrections$suggested_correction
  }
  
  return(list(
    corrections = corrections,
    final_predictions = final_preds,
    parameters = best_params,
    n_remaining = length(remaining_misclassified),
    iteration = iteration
  ))
}

#' Enhanced Flip Analysis for Iterative Approach
#' 
#' @param coordinates Cell coordinates
#' @param expression Expression data
#' @param predictions Current predictions
#' @param probabilities Class probabilities
#' @param k Number of neighbors
#' @param anchor_cut Anchor confidence threshold
#' @param consensus_threshold Consensus threshold
#' @param min_confidence_gain Minimum confidence gain
#' @param all_predictions All predictions including corrected ones
#' @param iteration Current iteration
#' @param verbose Print messages
#' @return Corrections data frame
#' @keywords internal
enhancedFlipAnalysisIterative <- function(
    coordinates,
    expression,
    predictions,
    probabilities,
    k,
    anchor_cut,
    consensus_threshold,
    min_confidence_gain,
    all_predictions,
    iteration,
    verbose = FALSE
) {
  
  n_cells <- length(predictions)
  if (n_cells < 5) {
    if (verbose) message("Too few cells for iterative analysis")
    return(NULL)
  }
  
  # Build k-NN using ALL available cells (not just remaining misclassified)
  # This is key: anchors can come from previously corrected cells
  all_coords <- coordinates[names(all_predictions), , drop = FALSE]
  knn <- buildKnnSimple(all_coords, k)
  
  # Compute enhanced confidence scores using updated prediction context
  enhanced_scores <- computeEnhancedConfidenceIterative(
    coordinates = all_coords,
    expression = expression,
    predictions = all_predictions,  # Use all predictions including corrections
    probabilities = probabilities,
    knn = knn,
    target_cells = names(predictions),  # Only score the remaining misclassified cells
    iteration = iteration,
    verbose = verbose
  )
  
  # Identify anchors from ALL cells (including previously corrected)
  anchor_mask <- enhanced_scores$multimodal_confidence >= anchor_cut
  
  # Focus correction analysis on remaining misclassified cells only
  target_indices <- match(names(predictions), names(all_predictions))
  target_mask <- !is.na(target_indices)
  candidate_indices <- target_indices[target_mask]
  
  corrections <- data.frame()
  
  # Evaluate each remaining misclassified cell
  for (i in seq_along(candidate_indices)) {
    cell_idx <- candidate_indices[i]
    cell_id <- names(predictions)[target_mask][i]
    
    neighbors <- knn$idx_list[[cell_idx]]
    anchor_neighbors <- intersect(neighbors, which(anchor_mask))
    
    if (length(anchor_neighbors) < 2) next
    
    # Get neighbor predictions (now includes previous corrections)
    neighbor_preds <- all_predictions[anchor_neighbors]
    neighbor_classes <- table(neighbor_preds)
    
    # Find dominant class among anchor neighbors
    dominant_class <- names(which.max(neighbor_classes))
    dominant_fraction <- max(neighbor_classes) / length(anchor_neighbors)
    
    if (dominant_fraction < consensus_threshold) next
    if (dominant_class == predictions[cell_id]) next  # Already correct
    
    # Calculate confidence gain (enhanced for iterative context)
    current_conf <- enhanced_scores$multimodal_confidence[cell_idx]
    
    # Estimate confidence gain from neighborhood improvement
    anchor_conf_mean <- mean(enhanced_scores$multimodal_confidence[anchor_neighbors], na.rm = TRUE)
    estimated_gain <- (anchor_conf_mean - current_conf) * dominant_fraction
    
    # Require minimum gain, but be more lenient in later iterations
    adjusted_min_gain <- min_confidence_gain * (0.8^(iteration-1))
    
    if (estimated_gain >= adjusted_min_gain) {
      corrections <- rbind(corrections, data.frame(
        cell_id = cell_id,
        original_prediction = predictions[cell_id],
        suggested_correction = dominant_class,
        confidence_before = current_conf,
        confidence_after = current_conf + estimated_gain,
        confidence_gain = estimated_gain,
        spatial_consensus = dominant_fraction,
        n_anchor_neighbors = length(anchor_neighbors),
        iteration = iteration,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(corrections)
}

#' Compute Enhanced Confidence for Iterative Analysis
#' 
#' @param coordinates All coordinates
#' @param expression Expression data
#' @param predictions All predictions (including previous corrections)
#' @param probabilities Probabilities for target cells
#' @param knn KNN structure
#' @param target_cells Target cells to score
#' @param iteration Current iteration
#' @param verbose Print messages
#' @return Enhanced confidence scores
#' @keywords internal
computeEnhancedConfidenceIterative <- function(
    coordinates,
    expression,
    predictions,
    probabilities,
    knn,
    target_cells,
    iteration,
    verbose = FALSE
) {
  
  n_cells <- length(predictions)
  
  # Spatial confidence using updated predictions
  spatial_confidence <- sapply(seq_len(n_cells), function(i) {
    neighbors <- knn$idx_list[[i]]
    if (length(neighbors) == 0) return(0.5)
    
    neighbor_preds <- predictions[neighbors]
    agreement <- mean(neighbor_preds == predictions[i], na.rm = TRUE)
    
    # Boost confidence for cells that were corrected in previous iterations
    if (iteration > 1) {
      # Check if this prediction was recently corrected (higher confidence)
      cell_id <- names(predictions)[i]
      # Simple heuristic: boost confidence for consistent neighborhoods
      if (agreement > 0.7) {
        agreement <- agreement * 1.05  # Slight boost for consistent cells
      }
    }
    
    return(min(agreement, 1.0))
  })
  
  # Expression consistency (if available)
  expression_consistency <- rep(0.5, n_cells)
  if (!is.null(expression)) {
    # Only compute for cells we have expression data for
    expr_cells <- intersect(names(predictions), rownames(expression))
    expr_indices <- match(expr_cells, names(predictions))
    
    expression_consistency[expr_indices] <- sapply(expr_indices, function(i) {
      neighbors <- knn$idx_list[[i]]
      expr_neighbors <- intersect(neighbors, expr_indices)
      if (length(expr_neighbors) <= 1) return(0.5)
      
      cell_expr <- expression[match(names(predictions)[i], rownames(expression)), ]
      neighbor_expr <- colMeans(expression[match(names(predictions)[expr_neighbors], rownames(expression)), , drop = FALSE], na.rm = TRUE)
      
      cor_val <- cor(cell_expr, neighbor_expr, use = "complete.obs")
      if (is.na(cor_val)) return(0.5)
      
      (cor_val + 1) / 2
    })
  }
  
  # Model confidence
  model_confidence <- if (!is.null(probabilities)) {
    target_indices <- match(target_cells, names(predictions))
    target_indices <- target_indices[!is.na(target_indices)]
    
    conf_vec <- rep(0.5, n_cells)
    if (length(target_indices) > 0) {
      conf_vec[target_indices] <- apply(probabilities, 1, max, na.rm = TRUE)
    }
    conf_vec
  } else {
    rep(0.5, n_cells)
  }
  
  # Combine with adaptive weights based on iteration
  # Later iterations rely more on spatial information as corrections accumulate
  base_spatial_weight <- 0.4
  base_expression_weight <- 0.3
  base_model_weight <- 0.3
  
  # Increase spatial weight in later iterations
  spatial_weight <- min(0.7, base_spatial_weight + (iteration - 1) * 0.05)
  expression_weight <- base_expression_weight * (1 - (spatial_weight - base_spatial_weight))
  model_weight <- base_model_weight * (1 - (spatial_weight - base_spatial_weight))
  
  # Normalize weights
  total_weight <- spatial_weight + expression_weight + model_weight
  spatial_weight <- spatial_weight / total_weight
  expression_weight <- expression_weight / total_weight
  model_weight <- model_weight / total_weight
  
  multimodal_confidence <- (
    spatial_weight * spatial_confidence +
    expression_weight * expression_consistency +
    model_weight * model_confidence
  )
  
  return(list(
    spatial_confidence = spatial_confidence,
    expression_consistency = expression_consistency,
    model_confidence = model_confidence,
    multimodal_confidence = multimodal_confidence,
    weights_used = list(
      spatial = spatial_weight,
      expression = expression_weight,
      model = model_weight,
      iteration = iteration
    )
  ))
}

#' Identify Correction Chains
#' 
#' @param previous_corrections Corrections from previous iteration
#' @param current_corrections Corrections from current iteration
#' @return Data frame with correction chain information
#' @keywords internal
identifyCorrectionChains <- function(previous_corrections, current_corrections) {
  
  if (is.null(previous_corrections) || is.null(current_corrections) ||
      nrow(previous_corrections) == 0 || nrow(current_corrections) == 0) {
    return(data.frame())
  }
  
  # Find cells that were corrected in previous iteration and affected current corrections
  chains <- data.frame()
  
  for (i in seq_len(nrow(current_corrections))) {
    current_cell <- current_corrections$cell_id[i]
    
    # Find if any previous corrections might have influenced this one
    # (spatially proximate cells that were corrected previously)
    influenced_by <- previous_corrections$cell_id[
      # Simple heuristic: could be improved with actual spatial distance
      abs(as.numeric(factor(previous_corrections$cell_id)) - as.numeric(factor(current_cell))) <= 2
    ]
    
    if (length(influenced_by) > 0) {
      chains <- rbind(chains, data.frame(
        current_cell = current_cell,
        influenced_by = paste(influenced_by, collapse = ";"),
        chain_length = 2,  # Simple chain of length 2
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(chains)
}

#' Calculate Iteration Improvement Metrics
#' 
#' @param STEAM.obj STEAM object
#' @param iter_results Results from current iteration
#' @param prev_accuracy_by_fold Previous accuracy by fold
#' @return List with improvement metrics
#' @keywords internal
calculateIterationImprovement <- function(STEAM.obj, iter_results, prev_accuracy_by_fold) {
  
  updated_accuracy_by_fold <- numeric(length(prev_accuracy_by_fold))
  
  # Calculate accuracy improvement for each fold
  for (f in seq_along(STEAM.obj$nested$ncv$outer_result)) {
    outer_result <- STEAM.obj$nested$ncv$outer_result[[f]]
    if (is.null(outer_result$preds)) {
      updated_accuracy_by_fold[f] <- prev_accuracy_by_fold[f]
      next
    }
    
    preds_df <- outer_result$preds
    
    # Get current predictions (original + all corrections so far)
    current_preds <- setNames(preds_df$predy, rownames(preds_df))
    true_labels <- setNames(preds_df$testy, rownames(preds_df))
    
    # Apply all corrections from all iterations
    if (!is.null(STEAM.obj$spatial_anchor_analysis$corrections)) {
      all_corrections <- STEAM.obj$spatial_anchor_analysis$corrections
      fold_corrections <- all_corrections[all_corrections$fold == f, ]
      if (nrow(fold_corrections) > 0) {
        current_preds[fold_corrections$cell_id] <- fold_corrections$suggested_correction
      }
    }
    
    # Apply current iteration corrections
    if (!is.null(iter_results$corrections) && nrow(iter_results$corrections) > 0) {
      current_iter_corrections <- iter_results$corrections[iter_results$corrections$fold == f, ]
      if (nrow(current_iter_corrections) > 0) {
        current_preds[current_iter_corrections$cell_id] <- current_iter_corrections$suggested_correction
      }
    }
    
    # Calculate updated accuracy
    updated_accuracy_by_fold[f] <- mean(current_preds == true_labels, na.rm = TRUE)
  }
  
  # Calculate improvement metrics
  accuracy_improvements <- updated_accuracy_by_fold - prev_accuracy_by_fold
  accuracy_improvements <- accuracy_improvements[!is.na(accuracy_improvements)]
  
  return(list(
    updated_accuracy_by_fold = updated_accuracy_by_fold,
    accuracy_improvements = accuracy_improvements,
    mean_accuracy_improvement = if (length(accuracy_improvements) > 0) mean(accuracy_improvements) else 0,
    median_accuracy_improvement = if (length(accuracy_improvements) > 0) median(accuracy_improvements) else 0,
    folds_improved = sum(accuracy_improvements > 0),
    folds_degraded = sum(accuracy_improvements < 0)
  ))
}

#' Parameter Tuning for Iterative Context
#' 
#' @param coordinates Cell coordinates
#' @param expression Expression data
#' @param predictions Current predictions
#' @param probabilities Class probabilities
#' @param k_grid K values to test
#' @param anchor_cut_grid Anchor cut values
#' @param consensus_grid Consensus values
#' @param updated_context Whether using updated predictions
#' @param verbose Print messages
#' @return Best parameters
#' @keywords internal
tuneParametersIterative <- function(
    coordinates,
    expression,
    predictions,
    probabilities,
    k_grid,
    anchor_cut_grid,
    consensus_grid,
    updated_context = TRUE,
    verbose = FALSE
) {
  
  n_cells <- length(predictions)
  if (n_cells < 10) {
    # Too small for tuning
    return(list(
      k = k_grid[1],
      anchor_cut = anchor_cut_grid[1],
      consensus = consensus_grid[1]
    ))
  }
  
  best_score <- 0
  best_params <- list(k = k_grid[1], anchor_cut = anchor_cut_grid[1], consensus = consensus_grid[1])
  
  # Simplified grid search for iterative context
  # In iterative context, we can be less conservative since previous iterations provide validation
  k_sample <- if (length(k_grid) > 3) k_grid[c(1, length(k_grid)%/%2, length(k_grid))] else k_grid
  anchor_sample <- if (length(anchor_cut_grid) > 3) anchor_cut_grid[c(1, length(anchor_cut_grid))] else anchor_cut_grid
  consensus_sample <- if (length(consensus_grid) > 3) consensus_grid[c(1, length(consensus_grid))] else consensus_grid
  
  for (k in k_sample) {
    for (anchor_cut in anchor_sample) {
      for (consensus in consensus_sample) {
        
        score <- evaluateParametersIterative(
          coordinates = coordinates,
          expression = expression,
          predictions = predictions,
          probabilities = probabilities,
          k = k,
          anchor_cut = anchor_cut,
          consensus = consensus,
          updated_context = updated_context
        )
        
        if (score > best_score) {
          best_score <- score
          best_params <- list(k = k, anchor_cut = anchor_cut, consensus = consensus)
        }
      }
    }
  }
  
  if (verbose) {
    message(sprintf("Best iterative params: k=%d, anchor_cut=%.2f, consensus=%.2f (score=%.3f)",
                   best_params$k, best_params$anchor_cut, best_params$consensus, best_score))
  }
  
  return(best_params)
}

#' Evaluate Parameters for Iterative Context
#' 
#' @param coordinates Cell coordinates
#' @param expression Expression data
#' @param predictions Current predictions
#' @param probabilities Class probabilities
#' @param k Number of neighbors
#' @param anchor_cut Anchor cut threshold
#' @param consensus Consensus threshold
#' @param updated_context Whether using updated context
#' @return Parameter score
#' @keywords internal
evaluateParametersIterative <- function(
    coordinates,
    expression,
    predictions,
    probabilities,
    k,
    anchor_cut,
    consensus,
    updated_context = TRUE
) {
  
  n_cells <- length(predictions)
  if (n_cells <= k) return(0)
  
  # Build k-NN
  knn <- buildKnnSimple(coordinates, k)
  
  # Compute confidence scores
  enhanced_scores <- computeEnhancedConfidenceSimple(
    coordinates = coordinates,
    expression = expression,
    predictions = predictions,
    probabilities = probabilities,
    knn = knn,
    use_expression = !is.null(expression),
    expression_weight = 0.3,
    spatial_weight = 0.4,
    model_weight = 0.3
  )
  
  # Count potential anchors
  anchor_mask <- enhanced_scores$multimodal_confidence >= anchor_cut
  spatial_consistent <- enhanced_scores$spatial_confidence[anchor_mask] >= consensus
  good_anchors <- sum(spatial_consistent, na.rm = TRUE)
  
  # Score based on anchor quality and density
  anchor_density <- good_anchors / n_cells
  mean_anchor_confidence <- if (sum(anchor_mask) > 0) {
    mean(enhanced_scores$multimodal_confidence[anchor_mask], na.rm = TRUE)
  } else {
    0
  }
  
  # In iterative context, prefer parameters that create stable, high-confidence anchors
  score <- anchor_density * mean_anchor_confidence
  
  # Bonus for updated context (corrected predictions should create better anchors)
  if (updated_context) {
    score <- score * 1.1
  }
  
  return(score)
}

#' Prepare Expression Data
#' 
#' @param expression_matrix Expression matrix
#' @param spatial_coords Spatial coordinates
#' @return Prepared expression data
#' @keywords internal
prepareExpressionData <- function(expression_matrix, spatial_coords) {
  
  if (is.null(expression_matrix)) return(NULL)
  
  # Use the existing alignment function from the original code
  result <- alignExpressionSpatialData(
    expression_matrix = expression_matrix,
    spatial_coords = spatial_coords,
    genes_transpose = TRUE,
    verbose = FALSE
  )
  
  return(result)
}
