

#' Calculate Adaptive Convergence Threshold
#'
#' Dynamically calculates convergence threshold based on performance stability,
#' success rates, and iteration progress. Used to determine when progressive learning
#' has converged and can stop iterating.
#'
#' @param STEAM.obj STEAM object containing progressive learning state
#' @return Numeric value representing the adaptive convergence threshold
#' @details
#' The function considers:
#' - Performance stability (variance in recent success rates)
#' - Average success rate (higher rates allow tighter thresholds)
#' - Iteration progress (later iterations use tighter thresholds)
#' - Bounds checking to ensure reasonable threshold values
#' @keywords internal
calculateAdaptiveConvergenceThreshold <- function(STEAM.obj) {
  
  pl <- STEAM.obj$spatial_anchor_analysis$progressive_learning
  conv <- pl$convergence
  
  if (!conv$enabled) {
    return(conv$base_threshold)
  }
  
  # Get performance history
  perf_history <- pl$performance_history
  
  if (nrow(perf_history) < 2) {
    return(conv$base_threshold)
  }
  
  # Analyze recent performance characteristics
  recent_window <- min(3, nrow(perf_history))
  recent_perf <- tail(perf_history, recent_window)
  
  # Performance stability factor
  stability_factor <- 1.0
  if (nrow(recent_perf) >= 2) {
    perf_variance <- var(recent_perf$success_rate)
    if (perf_variance < 0.005) {
      # Very stable performance - can use tighter convergence threshold
      stability_factor <- 0.5
    } else if (perf_variance > 0.05) {
      # Volatile performance - use looser convergence threshold  
      stability_factor <- 2.0
    }
  }
  
  # Success rate factor - higher success rates can use tighter thresholds
  avg_success_rate <- mean(recent_perf$success_rate, na.rm = TRUE)
  success_factor <- if (avg_success_rate > 0.7) {
    0.8  # Tighter threshold for good performance
  } else if (avg_success_rate < 0.5) {
    1.5  # Looser threshold for poor performance
  } else {
    1.0  # Standard threshold
  }
  
  # Iteration progress factor - later iterations can use tighter thresholds
  current_iter <- pl$current_iteration
  max_iter <- pl$max_iterations
  iteration_factor <- 1.0 - (current_iter / max_iter) * 0.3  # Up to 30% reduction
  
  # Calculate adaptive threshold
  adaptive_threshold <- conv$base_threshold * 
    stability_factor * 
    success_factor * 
    iteration_factor * 
    conv$adaptive_multiplier
  
  # Bounds checking
  adaptive_threshold <- max(0.001, min(0.05, adaptive_threshold))
  
  return(adaptive_threshold)
}


#' Calculate Adaptive Confidence Threshold
#'
#' Dynamically adjusts the confidence threshold for progressive learning based on
#' current performance metrics, convergence state, and iteration progress.
#'
#' @param STEAM.obj STEAM object containing progressive learning configuration
#' @param iteration_performance List containing success_rate and other performance metrics
#' @return Numeric value representing the new adaptive confidence threshold
#' @details
#' Adjusts threshold based on:
#' - Success rate (high success → lower threshold for more corrections)
#' - Performance trend (converging, improving, or declining)
#' - Convergence analysis using adaptive convergence threshold
#' - Iteration progress for dynamic aggressiveness
#' @keywords internal
calculateAdaptiveThreshold <- function(STEAM.obj, iteration_performance) {
  
  pl <- STEAM.obj$spatial_anchor_analysis$progressive_learning
  conv <- pl$convergence
  current_threshold <- pl$confidence_threshold
  success_rate <- iteration_performance$success_rate
  
  # Get performance history for convergence analysis
  perf_history <- pl$performance_history
  current_iteration <- pl$current_iteration
  
  # Calculate adaptive convergence threshold
  adaptive_conv_threshold <- calculateAdaptiveConvergenceThreshold(STEAM.obj)
  
  # Base threshold adjustment (existing logic)
  base_adjustment <- 0
  if (success_rate >= 0.8) {
    # High success - be more aggressive (lower threshold)
    base_adjustment <- -0.05
  } else if (success_rate <= 0.4) {
    # Low success - be more conservative (higher threshold)  
    base_adjustment <- 0.05
  } else {
    # Moderate success - proportional adjustment
    base_adjustment <- (0.6 - success_rate) * 0.02
  }
  
  # Apply convergence-aware modifications
  if (conv$enabled && nrow(perf_history) >= 2) {
    
    # Performance trend analysis
    recent_window <- min(conv$performance_window, nrow(perf_history))
    recent_perf <- tail(perf_history$success_rate, recent_window)
    
    # Calculate performance trend
    if (length(recent_perf) >= 2) {
      trend_slope <- (tail(recent_perf, 1) - head(recent_perf, 1)) / (length(recent_perf) - 1)
      
      # Convergence-based adjustment
      if (abs(trend_slope) < adaptive_conv_threshold) {
        # Performance is converging - can be more aggressive
        convergence_bonus <- -0.02 * conv$aggressiveness_factor
        base_adjustment <- base_adjustment + convergence_bonus
        
        # Stability bonus for consistent performance
        if (var(recent_perf) < 0.01) {  # Low variance = stable
          base_adjustment <- base_adjustment - 0.01
        }
        
      } else if (trend_slope < -adaptive_conv_threshold) {
        # Performance declining - be more conservative
        divergence_penalty <- 0.03 * conv$aggressiveness_factor  
        base_adjustment <- base_adjustment + divergence_penalty
        
      } else if (trend_slope > adaptive_conv_threshold) {
        # Performance improving - can be slightly more aggressive
        improvement_bonus <- -0.015 * conv$aggressiveness_factor
        base_adjustment <- base_adjustment + improvement_bonus
      }
    }
    
    # Adaptive aggressiveness based on iteration progress
    iteration_factor <- min(1.0, current_iteration / pl$max_iterations)
    base_adjustment <- base_adjustment * (1 + iteration_factor * 0.2)
  }
  
  # Apply the adjustment
  new_threshold <- current_threshold + base_adjustment
  new_threshold <- max(pl$min_threshold, min(0.95, new_threshold))
  
  return(new_threshold)
}

#' Calculate Multimodal Confidence Scores
#'
#' Computes confidence scores for cell type predictions using spatial consistency,
#' expression similarity, and model probabilities. Supports adaptive weighting
#' across iterations.
#'
#' @param coordinates Matrix of spatial coordinates for cells
#' @param expression Expression matrix (cells x genes) 
#' @param predictions Vector of current cell type predictions
#' @param probabilities Matrix of prediction probabilities from model
#' @param knn K-nearest neighbor structure with neighbor indices
#' @param target_cells Optional vector of specific cells to compute confidence for
#' @param iteration Current iteration number (affects weight adaptation)
#' @param use_expression Whether to include expression-based confidence
#' @param expression_weight Weight for expression consistency component
#' @param spatial_weight Weight for spatial consistency component  
#' @param model_weight Weight for model probability component
#' @param verbose Whether to print debug information
#' @return List containing spatial_confidence, expression_consistency, model_confidence,
#'         multimodal_confidence vectors and optionally weights_used
#' @details
#' The function adapts weights across iterations, increasing spatial importance
#' as corrections accumulate. Expression consistency uses correlation between
#' cell and neighbor expression profiles.
#' @keywords internal
calculateConfidence <- function(
    coordinates,
    expression,
    predictions,
    probabilities,
    knn,
    target_cells = NULL,
    iteration = 1,
    use_expression = TRUE,
    expression_weight = 0.3,
    spatial_weight = 0.4,
    model_weight = 0.3,
    verbose = FALSE
) {
  
  n_cells <- length(predictions)
  
  # Spatial confidence using predictions (potentially updated from previous iterations)
  spatial_confidence <- sapply(seq_len(n_cells), function(i) {
    neighbors <- knn$idx_list[[i]]
    if (length(neighbors) == 0) return(0.5)
    
    neighbor_preds <- predictions[neighbors]
    agreement <- mean(neighbor_preds == predictions[i], na.rm = TRUE)
    
    # Boost confidence for cells in consistent neighborhoods in later iterations
    if (iteration > 1 && agreement > 0.7) {
      agreement <- agreement * 1.05  # Slight boost for consistent cells
    }
    
    return(min(agreement, 1.0))
  })
  
  # Expression consistency
  expression_consistency <- rep(0.5, n_cells)
  if (use_expression && !is.null(expression)) {
    # Handle both matrix formats and target cell subsetting
    expr_cells <- intersect(names(predictions), rownames(expression))
    expr_indices <- match(expr_cells, names(predictions))
    
    if (length(expr_indices) > 0) {
      expression_consistency[expr_indices] <- sapply(expr_indices, function(i) {
        neighbors <- knn$idx_list[[i]]
        expr_neighbors <- intersect(neighbors, expr_indices)
        if (length(expr_neighbors) <= 1) return(0.5)
        
        # Handle both row indexing methods
        cell_expr <- if (names(predictions)[i] %in% rownames(expression)) {
          expression[match(names(predictions)[i], rownames(expression)), ]
        } else {
          expression[i, ]
        }
        
        neighbor_expr <- if (names(predictions)[i] %in% rownames(expression)) {
          colMeans(expression[match(names(predictions)[expr_neighbors], rownames(expression)), , drop = FALSE], na.rm = TRUE)
        } else {
          colMeans(expression[expr_neighbors, , drop = FALSE], na.rm = TRUE)
        }
        
        cor_val <- cor(cell_expr, neighbor_expr, use = "complete.obs")
        if (is.na(cor_val)) return(0.5)
        
        (cor_val + 1) / 2  # Convert to 0-1 scale
      })
    }
  }
  
  # Model confidence
  model_confidence <- if (!is.null(probabilities)) {
    if (!is.null(target_cells)) {
      # Only compute for target cells
      target_indices <- match(target_cells, names(predictions))
      target_indices <- target_indices[!is.na(target_indices)]
      
      conf_vec <- rep(0.5, n_cells)
      if (length(target_indices) > 0) {
        conf_vec[target_indices] <- apply(probabilities, 1, max, na.rm = TRUE)
      }
      conf_vec
    } else {
      apply(probabilities, 1, max, na.rm = TRUE)
    }
  } else {
    rep(0.5, n_cells)
  }
  
  # Adaptive weight adjustment for iterations
  if (iteration > 1) {
    # Later iterations rely more on spatial information as corrections accumulate
    base_spatial_weight <- spatial_weight
    base_expression_weight <- expression_weight
    base_model_weight <- model_weight
    
    # Increase spatial weight in later iterations
    spatial_weight <- min(0.7, base_spatial_weight + (iteration - 1) * 0.05)
    expression_weight <- base_expression_weight * (1 - (spatial_weight - base_spatial_weight))
    model_weight <- base_model_weight * (1 - (spatial_weight - base_spatial_weight))
  }
  
  # Normalize weights
  total_weight <- spatial_weight + expression_weight + model_weight
  spatial_weight <- spatial_weight / total_weight
  expression_weight <- expression_weight / total_weight
  model_weight <- model_weight / total_weight
  
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
    multimodal_confidence = multimodal_confidence,
    weights_used = if (iteration > 1) {
      list(
        spatial = spatial_weight,
        expression = expression_weight,
        model = model_weight,
        iteration = iteration
      )
    } else NULL
  ))
}



#' Calculate Confidence Scores for Predictions
#'
#' High-level function to calculate confidence scores for cell predictions using
#' available data modalities (spatial, expression, model probabilities).
#'
#' @param STEAM.obj STEAM object containing spatial coordinates, expression data, etc.
#' @param predictions Named vector of cell type predictions
#' @return Named vector of confidence scores (0-1) for each cell
#' @details
#' Attempts to use multimodal confidence calculation if spatial data is available,
#' falls back to prediction probabilities, or uses random baseline if no suitable
#' data is found. Automatically builds k-NN structure for spatial analysis.
#' @keywords internal
calculateConfidenceScores <- function(STEAM.obj, predictions) {
  
  # Try to build k-NN structure for enhanced confidence calculation
  coords <- STEAM.obj$spatial
  expression <- if ("count_exp" %in% names(STEAM.obj)) STEAM.obj$count_exp else NULL
  probabilities <- if ("prediction_probs" %in% names(STEAM.obj)) STEAM.obj$prediction_probs else NULL
  
  # Build k-NN if we have spatial data
  if (!is.null(coords)) {
    # Get current iteration from progressive learning state
    current_iteration <- if ("progressive_learning" %in% names(STEAM.obj) && 
                             "current_iteration" %in% names(STEAM.obj$spatial_anchor_analysis$progressive_learning)) {
      STEAM.obj$spatial_anchor_analysis$progressive_learning$current_iteration
    } else 1
    
    # Build k-NN structure
    coords_matrix <- as.matrix(coords)
    common_cells <- intersect(names(predictions), rownames(coords_matrix))
    
    if (length(common_cells) > 5) {
      # Create simple k-NN structure
      knn_k <- min(8, length(common_cells) - 1)
      knn <- list(idx_list = list())
      
      for (i in seq_along(common_cells)) {
        cell <- common_cells[i]
        cell_coords <- coords_matrix[cell, ]
        
        # Calculate distances to all other cells
        coord_diff <- sweep(coords_matrix[common_cells, , drop = FALSE], 2, cell_coords, "-")
        distances <- sqrt(rowSums(coord_diff^2))
        names(distances) <- common_cells
        
        # Get nearest neighbors (excluding self)
        neighbors <- names(sort(distances))[2:(knn_k + 1)]
        neighbor_indices <- match(neighbors, common_cells)
        neighbor_indices <- neighbor_indices[!is.na(neighbor_indices)]
        
        knn$idx_list[[i]] <- neighbor_indices
      }
      
      # Use enhanced confidence calculation
      confidence_result <- calculateConfidence(
        coordinates = coords_matrix[common_cells, , drop = FALSE],
        expression = expression,
        predictions = predictions[common_cells],
        probabilities = probabilities,
        knn = knn,
        iteration = current_iteration,
        verbose = FALSE
      )
      
      # Extract multimodal confidence scores
      confidence <- rep(0.5, length(predictions))
      names(confidence) <- names(predictions)
      confidence[common_cells] <- confidence_result$multimodal_confidence
      
      return(confidence)
    }
  }
  
  # Fallback: Use prediction probabilities if available
  if (!is.null(probabilities) && (is.matrix(probabilities) || is.data.frame(probabilities))) {
    confidence <- apply(probabilities, 1, max, na.rm = TRUE)
    if (length(confidence) == length(predictions)) {
      names(confidence) <- names(predictions)
      return(confidence)
    }
  }
  
  # Final fallback: random confidence around 0.5
  confidence <- runif(length(predictions), 0.4, 0.6)
  names(confidence) <- names(predictions)
  return(confidence)
}




#' Calculate Correction Momentum
#'
#' Calculates a momentum factor for correction confidence based on proximity
#' to previously successful corrections. Cells near successful corrections
#' receive a confidence boost.
#'
#' @param cell_id ID of the cell to calculate momentum for
#' @param STEAM.obj STEAM object containing correction history
#' @param coordinates Matrix of spatial coordinates
#' @param iteration Current iteration number
#' @return Numeric momentum multiplier (1.0 = no boost, >1.0 = confidence boost)
#' @details
#' Uses exponential decay based on distance to nearest successful correction.
#' Only applies in iteration 2 and beyond. Maximum boost is 25%.
#' @keywords internal
calculateCorrectionMomentum <- function(cell_id, STEAM.obj, coordinates, iteration) {
  
  if (iteration <= 1) return(1.0)
  
  # Get successful corrections from previous iterations
  all_corrections <- STEAM.obj$spatial_anchor_analysis$corrections
  if (is.null(all_corrections)) return(1.0)
  
  successful_corrections <- all_corrections[all_corrections$correct == TRUE & 
                                              all_corrections$iteration < iteration, ]
  
  if (nrow(successful_corrections) == 0) return(1.0)
  
  # Calculate distance to nearest successful correction
  if (!(cell_id %in% rownames(coordinates))) return(1.0)
  
  cell_coord <- coordinates[cell_id, ]
  success_coords <- coordinates[successful_corrections$cell_id, , drop = FALSE]
  
  if (nrow(success_coords) == 0) return(1.0)
  
  # Find minimum distance to successful corrections
  distances <- apply(success_coords, 1, function(row) {
    sqrt(sum((as.numeric(cell_coord) - as.numeric(row))^2))
  })
  
  min_dist <- min(distances)
  
  # Boost confidence for cells near successful corrections
  momentum <- exp(-min_dist / 50)  # Adjust decay parameter as needed
  return(1 + momentum * 0.25)  # Up to 25% boost
}


#' Calculate Progressive Learning Iteration Metrics
#'
#' Computes performance metrics for a progressive learning iteration including
#' success trends, correction efficiency, and coverage rates.
#'
#' @param iteration_result Result object from progressive learning iteration
#' @param iteration Current iteration number
#' @param verbose Whether to print detailed metrics
#' @return List containing success_trend correlation and coverage_rate
#' @details
#' Analyzes performance trends over iterations using correlation analysis,
#' calculates correction efficiency, and measures what proportion of cells
#' have been processed. Used for convergence analysis and performance monitoring.
#' @keywords internal
calculateIterationMetrics <- function(iteration_result, iteration, verbose = TRUE) {
  
  if (verbose && iteration > 1) {
    cat("\n=== Progressive Learning Metrics ===\n")
    
    # Calculate trend if we have multiple data points
    if (length(iteration_result$progressive_state$success_trends) >= 2) {
      success_trend <- cor(
        seq_along(iteration_result$progressive_state$success_trends),
        iteration_result$progressive_state$success_trends,
        use = "complete.obs"
      )
      
      trend_label <- if (is.na(success_trend)) "Stable" else
        if (success_trend > 0.1) "Improving" else
          if (success_trend < -0.1) "Declining" else "Stable"
      
      cat(sprintf("Success Trend: %s (r=%.3f)\n", trend_label, ifelse(is.na(success_trend), 0, success_trend)))
    }
    
    # Efficiency metrics
    avg_corrections <- mean(sapply(iteration_result$progressive_state$performance_history$corrections_attempted, function(x) x))
    cat(sprintf("Correction Efficiency: %.1f corrections/iteration\n", avg_corrections))
    
    # Coverage rate
    total_cells <- length(extractPreds(iteration_result$STEAM.obj))
    coverage_rate <- length(iteration_result$progressive_state$corrected_cells) / total_cells
    cat(sprintf("Coverage Rate: %.1f%% of cells processed\n", coverage_rate * 100))
  }
  
  return(list(
    success_trend = if(length(iteration_result$progressive_state$success_trends) >= 2) {
      cor(seq_along(iteration_result$progressive_state$success_trends),
          iteration_result$progressive_state$success_trends, use = "complete.obs")
    } else 0,
    coverage_rate = length(iteration_result$progressive_state$corrected_cells) / 
      length(extractPreds(iteration_result$STEAM.obj))
  ))
}


#' Calculate Progressive Learning Performance Metrics
#'
#' Analyzes convergence trends, learning efficiency, and overall performance
#' across progressive learning iterations.
#'
#' @param STEAM.obj STEAM object containing progressive learning history
#' @param iteration Current iteration number
#' @param verbose Whether to print metric summaries
#' @return List containing convergence_slope, convergence_stability, corrections_per_second
#' @details
#' Uses linear regression to analyze convergence trends from recent performance
#' history. Calculates learning efficiency as corrections per second. Provides
#' convergence stability based on variance in success rates.
#' @keywords internal
calculateProgressiveMetrics <- function(STEAM.obj, iteration, verbose = TRUE) {
  
  pl <- STEAM.obj$spatial_anchor_analysis$progressive_learning
  metrics <- list()
  
  if (nrow(pl$performance_history) > 1) {
    # Calculate convergence trend
    recent_data <- tail(pl$performance_history, min(5, nrow(pl$performance_history)))
    if (nrow(recent_data) > 1) {
      trend_model <- lm(success_rate ~ iteration, data = recent_data)
      metrics$convergence_slope <- coef(trend_model)[2]
      metrics$convergence_stability <- 1 - var(recent_data$success_rate)
    }
    
    # Calculate learning efficiency
    total_corrections <- sum(pl$performance_history$corrections_successful)
    total_time <- sum(pl$performance_history$execution_time, na.rm = TRUE)
    metrics$corrections_per_second <- if (total_time > 0) total_corrections / total_time else 0
    
    if (verbose) {
      cat("Progressive Metrics:\n")
      if (!is.null(metrics$convergence_slope) && !is.na(metrics$convergence_slope)) {
        cat(sprintf("  Convergence trend: %.4f\n", metrics$convergence_slope))
      }
      if (!is.null(metrics$corrections_per_second)) {
        cat(sprintf("  Learning efficiency: %.2f corrections/sec\n", metrics$corrections_per_second))
      }
    }
  }
  
  return(metrics)
}





#' Estimate Improved Success Rate for Filtered Corrections
#'
#' Predicts the expected success rate for a set of filtered corrections based
#' on learned pattern performance from previous iterations.
#'
#' @param filtered_corrections Data frame of corrections after pattern filtering
#' @param learning_memory Optional learning memory containing pattern performance
#' @return Numeric success rate estimate as percentage (0-100)
#' @details
#' Uses weighted average of learned pattern success rates when available.
#' Falls back to conservative baseline estimates for unknown patterns.
#' Incorporates priority weights in calculation.
#' @keywords internal
estimate_improved_success_rate <- function(filtered_corrections, learning_memory = NULL) {
  
  if (nrow(filtered_corrections) == 0) return(0)
  
  # If no learning memory, use conservative estimate
  if (is.null(learning_memory) || is.null(learning_memory$pattern_performance)) {
    return(75)  # Conservative baseline estimate
  }
  
  pattern_perf <- learning_memory$pattern_performance
  
  # Calculate weighted average success rate based on actual learned performance
  total_weight <- 0
  weighted_success <- 0
  
  for (i in seq_len(nrow(filtered_corrections))) {
    correction <- filtered_corrections[i, ]
    pattern_key <- paste0(as.character(correction$original_pred), "→", as.character(correction$suggested_correction))
    
    if (pattern_key %in% names(pattern_perf)) {
      # Use learned success rate
      pattern_success <- pattern_perf[[pattern_key]]$success_rate
      weight <- correction$priority_weight
      
      weighted_success <- weighted_success + (pattern_success * weight)
      total_weight <- total_weight + weight
    } else {
      # Unknown pattern - use moderate estimate
      weight <- correction$priority_weight
      weighted_success <- weighted_success + (0.5 * weight)  # 50% for unknown patterns
      total_weight <- total_weight + weight
    }
  }
  
  if (total_weight > 0) {
    return((weighted_success / total_weight) * 100)
  } else {
    return(75)  # Fallback estimate
  }
}


#' Estimate Success Rate for Filtered Corrections (Universal Data Types)
#'
#' Estimates the expected success rate for correction patterns using learned
#' performance data or similarity-based heuristics. Supports any R data type.
#'
#' @param filtered_corrections Data frame containing corrections to estimate
#' @return Numeric success rate estimate as percentage (0-100)
#' @details
#' This function replaces hardcoded numeric patterns with universal data type support.
#' Uses dynamic pattern learning when available, or falls back to similarity-based
#' estimation. For numeric types, uses distance-based heuristics. All pattern
#' keys are created using character conversion for universal compatibility.
#' @keywords internal
estimateFilteredSuccessRate <- function(filtered_corrections) {
  
  if (nrow(filtered_corrections) == 0) return(0)
  
  # Use dynamic pattern success rates based on learning memory if available
  # This replaces hardcoded numeric patterns with universal data type support
  success_estimate <- 75  # Default conservative estimate
  
  # If we have learning memory, use it to estimate success rates
  if (exists("STEAM.obj") && "progressive_learning" %in% names(STEAM.obj) && 
      !is.null(STEAM.obj$spatial_anchor_analysis$progressive_learning$pattern_performance)) {
    
    pattern_perf <- STEAM.obj$spatial_anchor_analysis$progressive_learning$pattern_performance
    total_weight <- 0
    weighted_success <- 0
    
    for (i in seq_len(nrow(filtered_corrections))) {
      correction <- filtered_corrections[i, ]
      # Create pattern key that works with any data type
      from_val <- as.character(correction$original_pred)
      to_val <- if ("suggested_correction" %in% colnames(correction)) {
        as.character(correction$suggested_correction) 
      } else {
        as.character(correction$suggested_pred)
      }
      pattern_key <- paste0(from_val, "→", to_val)
      
      if (pattern_key %in% names(pattern_perf)) {
        perf <- pattern_perf[[pattern_key]]
        weight <- perf$attempts
        success_rate <- perf$success_rate
        
        total_weight <- total_weight + weight
        weighted_success <- weighted_success + (success_rate * weight)
      } else {
        # Use default for unknown patterns
        total_weight <- total_weight + 1
        weighted_success <- weighted_success + 0.75
      }
    }
    
    if (total_weight > 0) {
      success_estimate <- (weighted_success / total_weight) * 100
    }
  } else {
    # Fallback: estimate based on correction characteristics
    total_weight <- 0
    weighted_success <- 0
    
    for (i in seq_len(nrow(filtered_corrections))) {
      correction <- filtered_corrections[i, ]
      from_val <- correction$original_pred
      to_val <- if ("suggested_correction" %in% colnames(correction)) {
        correction$suggested_correction 
      } else {
        correction$suggested_pred
      }
      
      # Universal similarity-based success estimation
      estimated_success <- 0.75  # Base rate
      
      # For numeric types, use distance-based estimation
      if (is.numeric(from_val) && is.numeric(to_val)) {
        distance <- abs(as.numeric(from_val) - as.numeric(to_val))
        if (distance == 1) {
          estimated_success <- 0.85
        } else if (distance >= 3) {
          estimated_success <- 0.65
        }
      }
      
      total_weight <- total_weight + 1
      weighted_success <- weighted_success + estimated_success
    }
    
    if (total_weight > 0) {
      success_estimate <- (weighted_success / total_weight) * 100
    }
  }
  
  return(success_estimate)
}




