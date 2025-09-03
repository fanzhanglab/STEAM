#' Run Iterative Analysis Without Ground Truth
#' 
#' @param STEAM.obj STEAM object
#' @param expression_matrix Expression matrix
#' @param k_grid K values grid
#' @param anchor_cut_grid Anchor cut grid
#' @param consensus_grid Consensus grid
#' @param min_confidence_gain Minimum confidence gain
#' @param expression_weight Expression weight
#' @param spatial_weight Spatial weight
#' @param model_weight Model weight
#' @param max_flip_fraction Maximum flip fraction
#' @param max_iterations Maximum iterations
#' @param convergence_threshold Convergence threshold
#' @param confidence_decay Confidence decay factor
#' @param update_neighborhoods Whether to update neighborhoods
#' @param track_correction_chains Whether to track chains
#' @param verbose Print messages
#' @return STEAM object with iterative results
#' @keywords internal
runIterative <- function(STEAM.obj, expression_matrix, k_grid,
                                        anchor_cut_grid, consensus_grid, min_confidence_gain,
                                        expression_weight, spatial_weight, model_weight,
                                        max_flip_fraction, max_iterations, convergence_threshold,
                                        confidence_decay, update_neighborhoods, track_correction_chains, verbose) {
  
  # Start timing for iterative analysis  
  t0 <- Sys.time()
  
  if (verbose) message("Using alternative validation (no ground truth available)")
  
  # Initialize validation metrics
  baseline_metrics <- baselineMetrics(STEAM.obj, expression_matrix)
  iteration_metrics <- list()
  iteration_metrics[[1]] <- baseline_metrics
  
  correction_history <- list()
  converged <- FALSE
  
  for (iter in 2:max_iterations) {
    if (verbose) message(sprintf("\n--- Iteration %d (No Ground Truth) ---", iter))
    
    # More conservative parameter decay without ground truth
    current_params <- list(
      anchor_cut_multiplier = confidence_decay^(iter-1),
      consensus_multiplier = confidence_decay^(iter-1),
      min_gain_multiplier = 0.95^(iter-1)  # Slower decay for minimum gain
    )
    
    # Identify uncertain cells using alternative criteria
    uncertain_cells <- identifyUncertain(
      STEAM.obj = STEAM.obj,
      expression_matrix = expression_matrix,
      iteration = iter,
      verbose = verbose
    )
    
    if (length(uncertain_cells) == 0) {
      if (verbose) message("No uncertain cells identified. Converged.")
      converged <- TRUE
      break
    }
    
    # Generate corrections with alternative validation
    iter_corrections <- generateCorrections(
      STEAM.obj = STEAM.obj,
      uncertain_cells = uncertain_cells,
      expression_matrix = expression_matrix,
      params = current_params,
      baseline_metrics = baseline_metrics,
      iteration = iter,
      verbose = verbose
    )
    
    if (is.null(iter_corrections) || nrow(iter_corrections) == 0) {
      if (verbose) message("No valid corrections generated. Converged.")
      converged <- TRUE
      break
    }
    
    # Apply corrections and update STEAM object
    STEAM.obj <- applyCorrections(STEAM.obj, iter_corrections, iter)
    correction_history[[iter]] <- iter_corrections
    
    # Compute updated validation metrics
    updated_metrics <- updateMetrics(
      STEAM.obj = STEAM.obj,
      expression_matrix = expression_matrix,
      corrections = iter_corrections,
      baseline_metrics = baseline_metrics
    )
    
    iteration_metrics[[iter]] <- updated_metrics
    
    # Check improvement using alternative metrics
    improvement <- calcImprovement(
      prev_metrics = iteration_metrics[[iter-1]],
      current_metrics = updated_metrics
    )
    
    if (verbose) {
      message(sprintf("Iteration %d: %d corrections, improvement: %+.4f",
                     iter, nrow(iter_corrections), improvement))
    }
    
    # Convergence check
    if (improvement < convergence_threshold) {
      if (verbose) message("Alternative validation convergence reached.")
      converged <- TRUE
      break
    }
    
    # Update baseline
    baseline_metrics <- updated_metrics
  }
  
  # Store results in same structure as ground truth version
  STEAM.obj$iterative_anchor_analysis <- list(
    method = "iterative_no_ground_truth",
    iterations_completed = length(iteration_metrics),
    converged = converged,
    validation_approach = "alternative_metrics",
    
    iteration_metrics = iteration_metrics,
    correction_history = correction_history,
    
    summary = list(
      total_corrections = {
        valid_corrections <- correction_history[!sapply(correction_history, is.null)]
        if (length(valid_corrections) > 0) {
          sum(sapply(valid_corrections, function(x) if (is.data.frame(x)) nrow(x) else 0))
        } else {
          0
        }
      },
      final_spatial_coherence = iteration_metrics[[length(iteration_metrics)]]$spatial_coherence,
      final_expression_coherence = iteration_metrics[[length(iteration_metrics)]]$expression_coherence,
      overall_improvement = calcImprovement(
        prev_metrics = iteration_metrics[[1]],
        current_metrics = iteration_metrics[[length(iteration_metrics)]]
      )
    )
  )
  
  return(STEAM.obj)
}

#' Run Iterative Analysis With Ground Truth
#' 
#' @param STEAM.obj STEAM object
#' @param expression_matrix Expression matrix
#' @param ground_truth Ground truth labels for validation
#' @param k_grid K values grid
#' @param anchor_cut_grid Anchor cut grid
#' @param consensus_grid Consensus grid
#' @param min_confidence_gain Minimum confidence gain
#' @param expression_weight Expression weight
#' @param spatial_weight Spatial weight
#' @param model_weight Model weight
#' @param max_flip_fraction Maximum flip fraction
#' @param max_iterations Maximum iterations
#' @param convergence_threshold Convergence threshold
#' @param confidence_decay Confidence decay factor
#' @param min_confidence_gain_decay Minimum gain decay
#' @param update_neighborhoods Whether to update neighborhoods
#' @param track_correction_chains Whether to track chains
#' @param verbose Print messages
#' @return STEAM object with iterative results
#' @keywords internal
runGroundTruth <- function(STEAM.obj, expression_matrix, ground_truth, 
                                          k_grid, anchor_cut_grid, consensus_grid,
                                          min_confidence_gain, expression_weight, spatial_weight,
                                          model_weight, max_flip_fraction, max_iterations,
                                          convergence_threshold, confidence_decay, 
                                          min_confidence_gain_decay, update_neighborhoods,
                                          track_correction_chains, verbose) {
  
  # Start timing for iterative analysis
  t0 <- Sys.time()
  
  if (verbose) message("Using ground truth validation")
  
  # Validate ground truth parameter
  if (is.null(ground_truth)) {
    stop("Ground truth labels must be provided when using ground truth validation mode")
  }
  
  # Extract original parameters from single-pass analysis
  original_params <- STEAM.obj$spatial_anchor_analysis$parameters
  
  # Initialize tracking
  iteration_results <- list()
  correction_chains <- list()
  cumulative_corrections <- data.frame()
  prev_accuracy_by_fold <- numeric()
  
  # Store initial results
  iteration_results[[1]] <- STEAM.obj$spatial_anchor_analysis$results
  
  # Calculate initial accuracy by fold (BASELINE - never changes)
  baseline_accuracy_by_fold <- numeric()
  for (f in seq_along(STEAM.obj$nested$ncv$outer_result)) {
    fold_result <- STEAM.obj$nested$ncv$outer_result[[f]]
    if (!is.null(fold_result$preds)) {
      preds_df <- fold_result$preds
      fold_cells <- rownames(preds_df)
      fold_ground_truth <- ground_truth[fold_cells]
      fold_predictions <- preds_df$predy
      initial_acc <- mean(fold_ground_truth == fold_predictions, na.rm = TRUE)
      baseline_accuracy_by_fold[f] <- initial_acc
    }
  }
  
  # prev_accuracy_by_fold should always be the baseline for comparison
  prev_accuracy_by_fold <- baseline_accuracy_by_fold
  
  # Add baseline improvement metrics for iteration 1
  if (length(baseline_accuracy_by_fold) > 0) {
    iteration_results[[1]]$improvement_metrics <- list(
      updated_accuracy_by_fold = baseline_accuracy_by_fold,
      accuracy_improvements = rep(0, length(baseline_accuracy_by_fold)),
      mean_accuracy_improvement = 0,
      median_accuracy_improvement = 0,
      folds_improved = 0,
      folds_degraded = 0
    )
  }
  
  converged <- FALSE
  
  # Iterative loop
  for (iter in 2:max_iterations) {
    if (verbose) message(sprintf("\n--- Iteration %d (With Ground Truth) ---", iter))
    
    # Apply parameter decay
    current_anchor_cut <- original_params$anchor_cut_grid * (confidence_decay^(iter-1))
    current_consensus <- original_params$consensus_grid * (confidence_decay^(iter-1))
    current_min_gain <- original_params$min_confidence_gain * (min_confidence_gain_decay^(iter-1))
    
    # Ensure reasonable bounds
    current_anchor_cut <- pmax(current_anchor_cut, 0.2)
    current_consensus <- pmax(current_consensus, 0.4)
    current_min_gain <- max(current_min_gain, 0.02)
    
    if (verbose) {
      message(sprintf("Decayed thresholds: anchor_cut=%.3f-%.3f, consensus=%.3f-%.3f, min_gain=%.3f",
                     min(current_anchor_cut), max(current_anchor_cut),
                     min(current_consensus), max(current_consensus),
                     current_min_gain))
    }
    
    # Update predictions in STEAM object with previous iteration's corrections
    updated_predictions <- updatePreds(STEAM.obj, iteration_results[[iter-1]])
    
    # Recompute spatial neighborhoods if requested
    if (update_neighborhoods) {
      # Rebuild KNN with updated predictions - this allows corrections to propagate
      if (verbose) message("Updating spatial neighborhoods with corrected predictions...")
      
      # This step is key - neighborhoods now reflect corrected labels
      STEAM.obj <- updateNeighbors(STEAM.obj, updated_predictions)
    }
    
    # Run analysis on remaining misclassified cells
    iter_results <- runIteration(
      STEAM.obj = STEAM.obj,
      expression_matrix = expression_matrix,
      ground_truth = ground_truth,
      k_grid = original_params$k_grid,
      anchor_cut_grid = current_anchor_cut,
      consensus_grid = current_consensus,
      min_confidence_gain = current_min_gain,
      iteration = iter,
      verbose = verbose
    )
    
    if (is.null(iter_results) || iter_results$total_corrections == 0) {
      if (verbose) message("No new corrections found. Convergence reached.")
      converged <- TRUE
      break
    }
    
    # Track correction chains
    if (track_correction_chains && iter > 2) {
      chains <- correctionChains(
        previous_corrections = iteration_results[[iter-1]]$corrections,
        current_corrections = iter_results$corrections
      )
      correction_chains[[iter]] <- chains
    }
    
    # Store iteration results
    iteration_results[[iter]] <- iter_results
    
    # Calculate improvement metrics using BASELINE accuracy (not previous iteration)
    improvement_metrics <- iterationImprovement(
      STEAM.obj = STEAM.obj,
      iter_results = iter_results,
      prev_accuracy_by_fold = baseline_accuracy_by_fold,  # Always compare vs baseline
      ground_truth = ground_truth
    )
    
    # Add improvement metrics to iteration results
    iteration_results[[iter]]$improvement_metrics <- improvement_metrics
    
    if (verbose) {
      message(sprintf("Iteration %d: %d new corrections, %.4f mean accuracy improvement vs baseline",
                     iter, iter_results$total_corrections,
                     improvement_metrics$mean_accuracy_improvement))
    }
    
    # Check convergence
    if (improvement_metrics$mean_accuracy_improvement < convergence_threshold) {
      if (verbose) message(sprintf("Convergence reached: improvement %.6f < threshold %.6f",
                                  improvement_metrics$mean_accuracy_improvement, 
                                  convergence_threshold))
      converged <- TRUE
      break
    }
    
    # DON'T update prev_accuracy_by_fold - it should always be the baseline!
    # prev_accuracy_by_fold should remain baseline_accuracy_by_fold
    cumulative_corrections <- rbind(cumulative_corrections, iter_results$corrections)
  }
  
  # Store comprehensive iterative results
  STEAM.obj$iterative_anchor_analysis <- list(
    method = "iterative_with_ground_truth",
    iterations_completed = length(iteration_results),
    converged = converged,
    total_runtime = difftime(Sys.time(), t0, units = "secs"),
    
    # Results by iteration
    iteration_results = iteration_results,
    correction_chains = if (track_correction_chains) correction_chains else NULL,
    
    # Summary statistics
    summary = list(
      total_unique_corrections = length(unique(cumulative_corrections$cell_id)),
      final_total_corrections = sum(sapply(iteration_results, function(x) x$total_corrections)),
      iterations_with_corrections = sum(sapply(iteration_results, function(x) x$total_corrections > 0)),
      mean_corrections_per_iteration = mean(sapply(iteration_results, function(x) x$total_corrections))
    ),
    
    # Parameters used
    parameters = list(
      max_iterations = max_iterations,
      convergence_threshold = convergence_threshold,
      confidence_decay = confidence_decay,
      min_confidence_gain_decay = min_confidence_gain_decay,
      update_neighborhoods = update_neighborhoods,
      track_correction_chains = track_correction_chains
    )
  )
  
  return(STEAM.obj)
}

#' Enhanced Spatial Anchor Analysis with Optional Iterative Refinement
#'
#' @description
#' Self-contained enhanced spatial anchor analysis that focuses exclusively on 
#' correcting misclassified cells identified during nested cross-validation.
#' This function integrates both spatial and expression data using multimodal 
#' confidence scoring and sophisticated flipping logic to improve classification 
#' accuracy. Uses iterative refinement by default for optimal results.
#' 
#' The function automatically extracts expression data from STEAM.obj$count_exp
#' if no expression_matrix is provided (expression data is required for analysis). 
#' By default (only_correct_flips = TRUE), only corrections that result in the 
#' correct label are applied and counted. Wrong corrections are ignored.
#'
#' @param STEAM.obj STEAM object with nested CV results, spatial coordinates, and 
#'   expression data (in STEAM.obj$count_exp)
#' @param expression_matrix Gene expression matrix (cells x genes or genes x cells). 
#'   If NULL, will automatically extract from STEAM.obj$count_exp (required)
#' @param genes_transpose Whether genes are in columns (TRUE) or rows (FALSE)
#' @param ground_truth Ground truth labels for validation. If NULL, uses alternative validation 
#'   without ground truth. Should be a named vector where names match cell identifiers
#' @param iterative Whether to use iterative refinement (default: TRUE for multiple iterations)
#' @param max_iterations Maximum number of iterations (default: 3, optimal from parameter sweep)
#' @param convergence_threshold Minimum improvement required to continue (default: 0.001)
#' @param confidence_decay Decay factor for confidence thresholds over iterations (default: 0.95)
#' @param min_confidence_gain_decay Decay for minimum confidence gain (default: 0.9)
#' @param update_neighborhoods Whether to recompute neighborhoods after corrections (default: TRUE)
#' @param track_correction_chains Whether to track correction cascades (default: TRUE)
#' @param k_grid Grid of k values for spatial neighbors. If NULL, defaults to 
#'   optimized c(5,8,12) for best precision/recall balance
#' @param anchor_cut_grid Grid of confidence thresholds for anchors. If NULL, defaults to
#'   optimized c(0.75,0.82,0.90) for best precision/recall balance
#' @param consensus_grid Grid of consensus thresholds. If NULL, defaults to
#'   optimized c(0.82,0.87,0.92) for best precision/recall balance
#' @param min_confidence_gain Minimum confidence gain required to flip a cell (default: 0.20, optimized)
#' @param expression_weight Weight for expression-based scoring (0-1, default: 0.4, optimized)
#' @param spatial_weight Weight for spatial consensus (0-1, default: 0.4)
#' @param model_weight Weight for model confidence (0-1, default: 0.2, optimized)
#' @param max_flip_fraction Maximum fraction of cells that can be flipped (default: 0.05, optimized)
#' @param marker_genes Specific genes to use (NULL for auto-selection)
#' @param auto_genes_n Number of genes for auto-selection (default: 50)
#' @param class_precision Per-class precision scores (optional)
#' @param min_class_precision Minimum class precision threshold (default: 0.70)
#' @param flip_safety_level Safety level to prevent wrong flips: "low", "medium", "high" (default: "high")
#' @param validation_mode Either "cross_validation" or "production"
#' @param parallel Parallelization mode: "none", "multicore", or "multisession"
#' @param workers Number of parallel workers
#' @param verbose Whether to print detailed progress messages
#'
#' @return Enhanced STEAM object with integrated analysis results
#'
#' @export
STEAM_anchor <- function(
    STEAM.obj,
    expression_matrix = NULL,
    genes_transpose = TRUE,
    ground_truth = NULL,  # Explicit ground truth parameter
    iterative = TRUE,  # Always use iterative mode for better results
    max_iterations = 3,  # Optimal value from parameter sweep (was NULL)
    convergence_threshold = 0.001,
    confidence_decay = 0.95,
    min_confidence_gain_decay = 0.9,
    update_neighborhoods = TRUE,
    track_correction_chains = TRUE,
    k_grid = NULL,  # Will be set to optimized values if NULL
    anchor_cut_grid = NULL,  # Will be set to optimized values if NULL
    consensus_grid = NULL,  # Will be set to optimized values if NULL
    min_confidence_gain = 0.20,  # Optimized value (was 0.10)
    expression_weight = 0.4,  # Optimized value (was 0.3)
    spatial_weight = 0.4,  # Same as before
    model_weight = 0.2,  # Optimized value (was 0.3)
    max_flip_fraction = 0.05,  # Optimized value (was 0.15)
    marker_genes = NULL,
    auto_genes_n = 50,
    class_precision = NULL,
    min_class_precision = 0.70,
    flip_safety_level = c("high", "medium", "low"),
    validation_mode = c("cross_validation", "production"),
    parallel = c("none", "multicore", "multisession"),
    workers = max(1, parallel::detectCores() - 1),
    verbose = TRUE
) {
  
  t0 <- Sys.time()
  validation_mode <- match.arg(validation_mode)
  parallel <- match.arg(parallel)
  flip_safety_level <- match.arg(flip_safety_level)
  
  if (verbose) {
    if (iterative) {
      message("=== STEAM Anchor Analysis (Iterative Mode - Default) ===")
    } else {
      message("=== STEAM Anchor Analysis (Single-Pass Mode) ===")
    }
    message("Analyzing misclassified cells for spatial correction...")
  }
  
  # Determine ground truth availability from explicit parameter
  has_ground_truth <- !is.null(ground_truth)
  
  if (verbose) {
    message(sprintf("Ground truth provided: %s", has_ground_truth))
    if (iterative) {
      message(sprintf("Iterative mode: up to %d iterations", max_iterations))
    }
  }
  
  # Extract expression matrix from STEAM.obj (always required)
  if (is.null(expression_matrix)) {
    if (verbose) message("Extracting expression matrix from STEAM.obj...")
    expression_matrix <- extractExpression(STEAM.obj, verbose = verbose)
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
  
  # Skip single-pass analysis and go directly to iterative refinement
  if (verbose) message("\nStarting iterative analysis...")
  
  if (has_ground_truth) {
    # Use ground truth validation approach
    STEAM.obj <- runGroundTruth(
      STEAM.obj = STEAM.obj,
      expression_matrix = expression_matrix,
      ground_truth = ground_truth,
      k_grid = k_grid,
      anchor_cut_grid = anchor_cut_grid,
      consensus_grid = consensus_grid,
      min_confidence_gain = min_confidence_gain,
      expression_weight = expression_weight,
      spatial_weight = spatial_weight,
      model_weight = model_weight,
      max_flip_fraction = max_flip_fraction,
      max_iterations = max_iterations,
      convergence_threshold = convergence_threshold,
      confidence_decay = confidence_decay,
      min_confidence_gain_decay = min_confidence_gain_decay,
      update_neighborhoods = update_neighborhoods,
      track_correction_chains = track_correction_chains,
      verbose = verbose
    )
  } else {
    # Use alternative validation approach
    STEAM.obj <- runIterative(
      STEAM.obj = STEAM.obj,
      expression_matrix = expression_matrix,
      k_grid = k_grid,
      anchor_cut_grid = anchor_cut_grid,
      consensus_grid = consensus_grid,
      min_confidence_gain = min_confidence_gain,
      expression_weight = expression_weight,
      spatial_weight = spatial_weight,
      model_weight = model_weight,
      max_flip_fraction = max_flip_fraction,
      max_iterations = max_iterations,
      convergence_threshold = convergence_threshold * 5,  # More lenient
      confidence_decay = confidence_decay * 0.95,  # Slower decay
      update_neighborhoods = update_neighborhoods,
      track_correction_chains = track_correction_chains,
      verbose = verbose
    )
  }
  
  # Consolidate iterative results into main spatial_anchor_analysis structure
  if (!is.null(STEAM.obj$iterative_anchor_analysis)) {
    if (verbose) message("Consolidating iterative results...")
    
    # Get all corrections from iterative analysis
    iterative_results <- STEAM.obj$iterative_anchor_analysis$iteration_results
    all_iterative_corrections <- do.call(rbind, lapply(iterative_results, function(x) x$corrections))
    
    if (!is.null(all_iterative_corrections) && nrow(all_iterative_corrections) > 0) {
      # Remove duplicates (keep last correction for each cell in each fold)
      all_iterative_corrections$unique_id <- paste(all_iterative_corrections$cell_id, all_iterative_corrections$fold, sep = "_")
      unique_corrections <- all_iterative_corrections[!duplicated(all_iterative_corrections$unique_id, fromLast = TRUE), ]
      unique_corrections$unique_id <- NULL  # Remove temporary column
      
      # Update the main spatial_anchor_analysis with final consolidated results
      STEAM.obj$spatial_anchor_analysis$corrections <- unique_corrections
      STEAM.obj$spatial_anchor_analysis$method <- "enhanced_spatial_anchor_iterative"
      STEAM.obj$spatial_anchor_analysis$iterations_completed <- STEAM.obj$iterative_anchor_analysis$iterations_completed
      STEAM.obj$spatial_anchor_analysis$iterative_summary <- STEAM.obj$iterative_anchor_analysis$summary
      
      if (verbose) {
        message(sprintf("Consolidated %d unique corrections from %d iterations (%d total before deduplication)", 
                       nrow(unique_corrections),
                       STEAM.obj$iterative_anchor_analysis$iterations_completed,
                       nrow(all_iterative_corrections)))
      }
    }
  }
  
  if (verbose) {
    total_runtime <- difftime(Sys.time(), t0, units = "secs")
    message(sprintf("Total STEAM_anchor runtime: %.1f seconds", as.numeric(total_runtime)))
  }
  
  return(STEAM.obj)
}
  
#' Detect Ground Truth Availability (DEPRECATED)
#'
#' @description 
#' This function is deprecated. Ground truth should now be provided explicitly
#' as a parameter to STEAM_anchor() instead of being automatically detected.
#' 
#' @param STEAM.obj STEAM object
#' @return Logical indicating if ground truth is available
#' @keywords internal
#' @deprecated Use explicit ground_truth parameter in STEAM_anchor() instead
detectTruth <- function(STEAM.obj) {  # Check if we have true labels available in nested CV results
  if (!is.null(STEAM.obj$nested) && !is.null(STEAM.obj$nested$ncv)) {
    for (f in seq_along(STEAM.obj$nested$ncv$outer_result)) {
      outer_result <- STEAM.obj$nested$ncv$outer_result[[f]]
      if (!is.null(outer_result$preds) && "testy" %in% colnames(outer_result$preds)) {
        return(TRUE)
      }
    }
  }
  
  # Check train CV results
  if (!is.null(STEAM.obj$train) && !is.null(STEAM.obj$train$ncv)) {
    for (f in seq_along(STEAM.obj$train$ncv$outer_result)) {
      outer_result <- STEAM.obj$train$ncv$outer_result[[f]]
      if (!is.null(outer_result$preds) && "testy" %in% colnames(outer_result$preds)) {
        return(TRUE)
      }
    }
  }
  
  return(FALSE)
}

#' Run Single-Pass Analysis (Original STEAM_anchor Logic)
#' 
#' @param STEAM.obj STEAM object
#' @param ... All parameters passed through
#' @return STEAM object with single-pass results
#' @keywords internal
runSingle <- function(STEAM.obj, ...) {
  
  # Start timing for single-pass analysis
  t0 <- Sys.time()
  
  # This is essentially your original STEAM_anchor function
  # Extract all the parameters and run the original analysis
  args <- list(...)
  
  # Input validation
  if (is.null(STEAM.obj$nested) || is.null(STEAM.obj$nested$ncv)) {
    stop("STEAM.obj lacks nested CV results. Run RunSTEAM(mode='nested') first.")
  }
  if (is.null(STEAM.obj$spatial)) {
    stop("STEAM.obj has no spatial coordinates.")
  }
  
  # Extract parameters
  expression_matrix <- args$expression_matrix
  ground_truth <- args$ground_truth
  genes_transpose <- args$genes_transpose %||% TRUE
  k_grid <- args$k_grid %||% c(5, 8, 12)
  anchor_cut_grid <- args$anchor_cut_grid %||% c(0.75, 0.82, 0.90)
  consensus_grid <- args$consensus_grid %||% c(0.82, 0.87, 0.92)
  min_confidence_gain <- args$min_confidence_gain %||% 0.20
  expression_weight <- args$expression_weight %||% 0.4
  spatial_weight <- args$spatial_weight %||% 0.4
  model_weight <- args$model_weight %||% 0.2
  max_flip_fraction <- args$max_flip_fraction %||% 0.05
  marker_genes <- args$marker_genes
  auto_genes_n <- args$auto_genes_n %||% 50
  class_precision <- args$class_precision
  min_class_precision <- args$min_class_precision %||% 0.70
  flip_safety_level <- args$flip_safety_level %||% "high"
  validation_mode <- args$validation_mode %||% "cross_validation"
  parallel <- args$parallel %||% "none"
  workers <- args$workers %||% max(1, parallel::detectCores() - 1)
  verbose <- args$verbose %||% TRUE
  
  # Hardcode safety feature - always only apply correct flips
  only_correct_flips <- TRUE
  
  # Extract data
  # Check for prediction results - support both simple CV and nested CV (same as plot_misclassified_cells)
  coords_all <- STEAM.obj$spatial
  
  # Always identify misclassified cells - this is now the default behavior
  if (verbose) {
    if (is.null(ground_truth)) {
      message("Identifying low-confidence cells for alternative validation (no ground truth provided)...")
    } else {
      message("Identifying misclassified cells for spatial correction...")
    }
  }
  
  # Extract misclassified cells using EXACT same logic as plot_misclassified_cells
  misclassified_by_fold <- list()
  total_misclassified <- 0
  
  # ---- SIMPLE CV path (same as plot_misclassified_cells) ----
  if (!is.null(STEAM.obj$test) && !is.null(STEAM.obj$test$predictions)) {
    if (verbose) message("Using simple CV predictions for misclassification detection...")
    
    test_lbls <- as.character(STEAM.obj$test$test.data.labels)
    preds <- as.character(STEAM.obj$test$predictions)
    test_coords <- STEAM.obj$test$test.data.coords
    test_ids <- rownames(as.data.frame(test_coords))
    
    if (is.null(ground_truth)) {
      # For alternative validation in simple CV: use prediction confidence if available
      if (!is.null(STEAM.obj$test$probs)) {
        model_confidence <- apply(STEAM.obj$test$probs, 1, max, na.rm = TRUE)
        low_confidence_mask <- model_confidence < 0.8
        mis_cells <- test_ids[low_confidence_mask]
        
        if (length(mis_cells) > 0) {
          misclassified_by_fold[[1]] <- mis_cells
          total_misclassified <- length(mis_cells)
          if (verbose) message(sprintf("Simple CV: Found %d low-confidence cells for alternative validation", length(mis_cells)))
        }
      } else {
        if (verbose) message("Simple CV: No prediction probabilities available for alternative validation")
      }
    } else {
      # Use same logic as plot_misclassified_cells: test_lbls != preds
      mis_mask <- test_lbls != preds
      mis_cells <- test_ids[mis_mask & !is.na(mis_mask)]
      
      if (length(mis_cells) > 0) {
        misclassified_by_fold[[1]] <- mis_cells
        total_misclassified <- length(mis_cells)
        if (verbose) message(sprintf("Simple CV: Found %d misclassified cells", length(mis_cells)))
      }
    }
    
  # ---- NESTED CV path (same as plot_misclassified_cells) ----
  } else if (!is.null(STEAM.obj$nested) && !is.null(STEAM.obj$nested$ncv$outer_result)) {
    if (verbose) message("Using nested CV predictions for misclassification detection...")
    
    ncv <- STEAM.obj$nested$ncv
    n_folds <- length(ncv$outer_result)
    
    for (f in seq_along(ncv$outer_result)) {
      preds_df <- ncv$outer_result[[f]]$preds
      if (is.null(preds_df) || nrow(preds_df) == 0) next
      
      # Get cell names for this fold
      fold_cells <- rownames(preds_df)
      
      if (is.null(ground_truth)) {
        # For alternative validation: identify low-confidence cells instead of misclassified cells
        fold_probs <- ncv$outer_result[[f]]$probs
        if (!is.null(fold_probs)) {
          # Use cells with low model confidence as candidates for correction
          model_confidence <- apply(fold_probs, 1, max, na.rm = TRUE)
          low_confidence_mask <- model_confidence < 0.8  # Threshold for low confidence
          mis_cells <- fold_cells[low_confidence_mask]
          
          if (length(mis_cells) > 0) {
            misclassified_by_fold[[f]] <- mis_cells
            total_misclassified <- total_misclassified + length(mis_cells)
            if (verbose) message(sprintf("Fold %d: Found %d low-confidence cells for alternative validation", f, length(mis_cells)))
          }
        } else {
          if (verbose) message(sprintf("Fold %d: No prediction probabilities available, skipping...", f))
        }
        next
      }
      
      # EXACT same logic as plot_misclassified_cells: use testy vs predy from nested CV
      # This ensures the misclassified cells identified here match those shown in plot_misclassified_cells
      if ("testy" %in% colnames(preds_df) && "predy" %in% colnames(preds_df)) {
        # Use the nested CV's ground truth (testy) vs predictions (predy)
        mis_mask <- preds_df$testy != preds_df$predy
        mis_cells <- fold_cells[mis_mask & !is.na(mis_mask)]
        
        if (length(mis_cells) > 0) {
          misclassified_by_fold[[f]] <- mis_cells
          total_misclassified <- total_misclassified + length(mis_cells)
          if (verbose) message(sprintf("Fold %d: Found %d misclassified cells (using nested CV labels)", f, length(mis_cells)))
        }
      } else {
        # Fallback: use provided ground truth if testy column is not available
        fold_ground_truth <- ground_truth[fold_cells]
        fold_predictions <- preds_df$predy
        
        # Identify misclassified cells
        mis_mask <- fold_ground_truth != fold_predictions
        mis_cells <- fold_cells[mis_mask & !is.na(mis_mask)]
        
        if (length(mis_cells) > 0) {
          misclassified_by_fold[[f]] <- mis_cells
          total_misclassified <- total_misclassified + length(mis_cells)
          if (verbose) message(sprintf("Fold %d: Found %d misclassified cells (using provided ground truth)", f, length(mis_cells)))
        }
      }
    }
    
  } else {
    stop("No predictions found. Run RunSTEAM with mode='simple' or mode='nested' first.")
  }
  
  if (total_misclassified == 0) {
    if (is.null(ground_truth)) {
      warning("No low-confidence cells found across any fold for alternative validation!")
    } else {
      warning("No misclassified cells found across any fold. All predictions were correct!")
    }
    # Return the object unchanged
    STEAM.obj$spatial_anchor_analysis <- list(
      method = "enhanced_spatial_anchor",
      results = list(
        total_misclassified = 0,
        total_corrections = 0,
        message = "No misclassified cells to analyze"
      ),
      runtime = difftime(Sys.time(), t0, units = "secs")
    )
    return(STEAM.obj)
  }
  
  # Expression data handling - always extract from STEAM.obj if not provided
  expr_data_prepared <- NULL
  if (is.null(expression_matrix)) {
    expression_matrix <- extractExpression(STEAM.obj, verbose = verbose)
    if (is.null(expression_matrix)) {
      stop("No expression data found in STEAM.obj. Expression data is required for STEAM_anchor analysis. Ensure STEAM.obj$count_exp is available.")
    }
    if (verbose) message("Using expression data from STEAM.obj")
  }
  
  if (verbose) {
    message("=== Pre-alignment debugging ===")
    message("Expression matrix dimensions: ", nrow(expression_matrix), " x ", ncol(expression_matrix))
    message("Expression rownames (first 3): ", paste(head(rownames(expression_matrix), 3), collapse = ", "))
    message("Expression colnames (first 3): ", paste(head(colnames(expression_matrix), 3), collapse = ", "))
    message("Spatial coords dimensions: ", nrow(coords_all), " x ", ncol(coords_all))
    message("Spatial rownames (first 3): ", paste(head(rownames(coords_all), 3), collapse = ", "))
    message("===============================")
  }
  
  # Align expression and spatial data
  alignment_result <- alignData(
    expression_matrix = expression_matrix,
    spatial_coords = coords_all,
    genes_transpose = genes_transpose,
    verbose = verbose
  )
  
  if (length(alignment_result$common_cells) < 50) {
    stop("Too few common cells between expression and spatial data.")
  }
  
  expr_data_prepared <- alignment_result
  use_expression <- TRUE
  
  # Normalize weights
  total_weight <- expression_weight + spatial_weight + model_weight
  expression_weight <- expression_weight / total_weight
  spatial_weight <- spatial_weight / total_weight
  model_weight <- model_weight / total_weight
  
  if (verbose) {
    message(sprintf("Processing %d CV folds with multimodal analysis...", n_folds))
    message(sprintf("Weights: expression=%.2f, spatial=%.2f, model=%.2f",
                   expression_weight, spatial_weight, model_weight))
  }
  
  # Process each fold (using original logic)
  fold_results <- list()
  all_corrections <- data.frame()
  all_final_predictions <- list()
  
  for (f in seq_len(n_folds)) {
    if (verbose && f <= 3) {
      message(sprintf("Processing fold %d/%d...", f, n_folds))
    }
    
    fold_result <- processFold(
      fold_idx = f,
      ncv = ncv,
      coords_all = coords_all,
      expr_data_prepared = expr_data_prepared,
      ground_truth = ground_truth,
      use_expression = use_expression,
      k_grid = k_grid,
      anchor_cut_grid = anchor_cut_grid,
      consensus_grid = consensus_grid,
      min_confidence_gain = min_confidence_gain,
      expression_weight = expression_weight,
      spatial_weight = spatial_weight,
      model_weight = model_weight,
      max_flip_fraction = max_flip_fraction,
      marker_genes = marker_genes,
      auto_genes_n = auto_genes_n,
      class_precision = class_precision,
      min_class_precision = min_class_precision,
      only_correct_flips = only_correct_flips,
      flip_safety_level = flip_safety_level,
      validation_mode = validation_mode,
      verbose = verbose && f <= 2,
      misclassified_cells = misclassified_by_fold[[f]]
    )
    
    fold_results[[f]] <- fold_result
    
    # Collect corrections and predictions
    if (!is.null(fold_result$corrections) && nrow(fold_result$corrections) > 0) {
      fold_corr <- fold_result$corrections
      fold_corr$fold <- f
      all_corrections <- rbind(all_corrections, fold_corr)
    }
    
    if (!is.null(fold_result$final_predictions)) {
      all_final_predictions[[f]] <- fold_result$final_predictions
    }
  }
  
  # Store single-pass results
  total_corrections <- nrow(all_corrections)
  
  STEAM.obj$spatial_anchor_analysis <- list(
    method = "enhanced_spatial_anchor",
    parameters = list(
      k_grid = k_grid,
      anchor_cut_grid = anchor_cut_grid,
      consensus_grid = consensus_grid,
      min_confidence_gain = min_confidence_gain,
      expression_weight = expression_weight,
      spatial_weight = spatial_weight,
      model_weight = model_weight,
      max_flip_fraction = max_flip_fraction,
      use_expression = use_expression,
      genes_transpose = genes_transpose
    ),
    results = list(
      total_misclassified = total_misclassified,
      total_corrections = total_corrections,
      correction_type = "correct_only",
      correction_rate = if (total_misclassified > 0) total_corrections / total_misclassified else 0,
      n_folds_processed = n_folds
    ),
    corrections = all_corrections,
    final_by_fold = all_final_predictions,
    fold_results = fold_results,
    runtime = difftime(Sys.time(), t0, units = "secs")
  )
  
  return(STEAM.obj)
}

# Load helper functions for iterative analysis
source("iterativeHelpers.R")


# Load helper functions for iterative analysis
source("iterativeHelpers.R")

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
enhancedConfidence <- function(coordinates, expression, predictions,
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
extractExpression <- function(STEAM.obj, verbose = FALSE) {
  
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
alignData <- function(expression_matrix, spatial_coords,
                                        genes_transpose = TRUE, verbose = FALSE) {
  
  # Determine if we need to transpose
  if (genes_transpose) {
    # Check if genes are in rows by comparing rownames vs colnames with spatial coords
    expr_rows_are_cells <- length(intersect(rownames(expression_matrix), rownames(spatial_coords))) > 0
    expr_cols_are_cells <- length(intersect(colnames(expression_matrix), rownames(spatial_coords))) > 0
    
    if (verbose) {
      message("Expression rownames match spatial: ", expr_rows_are_cells)
      message("Expression colnames match spatial: ", expr_cols_are_cells)
    }
    
    if (expr_cols_are_cells && !expr_rows_are_cells) {
      # Cells are in columns (genes in rows), need to transpose
      expr_mat <- t(expression_matrix)
      if (verbose) message("Transposed expression matrix (genes were in rows, cells in columns)")
    } else if (expr_rows_are_cells && !expr_cols_are_cells) {
      # Cells are already in rows
      expr_mat <- expression_matrix
      if (verbose) message("Expression matrix orientation is correct (cells in rows)")
    } else {
      # Fallback to original dimension-based logic
      if (nrow(expression_matrix) > ncol(expression_matrix)) {
        # More rows than columns, likely genes in rows
        expr_mat <- t(expression_matrix)
        if (verbose) message("Transposed expression matrix (more rows than columns, assuming genes in rows)")
      } else {
        expr_mat <- expression_matrix
        if (verbose) message("Kept expression matrix as-is (more columns than rows)")
      }
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
    
    # Diagnostic information for debugging cell ID mismatches
    if (length(common_cells) == 0) {
      message("\n=== DEBUGGING CELL ID MISMATCH ===")
      message("Sample expression cell IDs (first 5):")
      message(paste(head(expr_cells, 5), collapse = ", "))
      message("Sample spatial cell IDs (first 5):")
      message(paste(head(spatial_cells, 5), collapse = ", "))
      
      # Check for potential formatting differences
      expr_sample <- head(expr_cells, 10)
      spatial_sample <- head(spatial_cells, 10)
      
      message("\nPotential fixes to try:")
      message("1. Check if cell IDs have different prefixes/suffixes")
      message("2. Check if one uses dots vs underscores vs dashes")
      message("3. Check if case sensitivity is an issue")
      message("4. Check if spatial coords need different rowname extraction")
      
      # Try some common transformations
      expr_no_prefix <- gsub("^[^0-9]*", "", expr_sample)
      spatial_no_prefix <- gsub("^[^0-9]*", "", spatial_sample)
      
      if (length(intersect(expr_no_prefix, spatial_no_prefix)) > 0) {
        message("5. Try removing prefixes - some matches found!")
      }
      
      message("=====================================\n")
    }
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

#' Process Enhanced Fold Clean
#' 
#' @param fold_idx Fold index
#' @param ncv Nested CV results
#' @param coords_all All coordinates
#' @param expr_data_prepared Prepared expression data
#' @param ground_truth Ground truth labels
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
processFold <- function(fold_idx, ncv, coords_all, expr_data_prepared,
                                       ground_truth, use_expression, k_grid, anchor_cut_grid,
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
  best_params <- tuneParams(
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
  corrections <- flipAnalysis(
    coordinates = test_coords,
    expression = test_expr,
    predictions = test_preds,
    probabilities = test_probs,
    true_labels = if (only_correct_flips && !is.null(ground_truth)) {
      ground_truth[test_cells]
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

#' Prepare Expression Data
#' 
#' @param expression_matrix Expression matrix
#' @param spatial_coords Spatial coordinates
#' @return Prepared expression data
#' @keywords internal
prepareExpression <- function(expression_matrix, spatial_coords) {
  
  if (is.null(expression_matrix)) return(NULL)
  
  # Use the existing alignment function from the original code
  result <- alignData(
    expression_matrix = expression_matrix,
    spatial_coords = spatial_coords,
    genes_transpose = TRUE,
    verbose = FALSE
  )
  
  return(result)
}

# Utility operator
`%||%` <- function(x, y) if (is.null(x)) y else x

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
tuneParams <- function(coordinates, expression, predictions, probabilities,
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
        
        score <- evalParams(
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
evalParams <- function(coordinates, expression, predictions, probabilities,
                                      k, anchor_cut, consensus) {
  
  n_cells <- length(predictions)
  if (n_cells <= k) return(0)
  
  # Build k-NN
  knn <- buildKnn(coordinates, k)
  
  # Compute confidence scores
  enhanced_scores <- enhancedConfidence(
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
flipAnalysis <- function(coordinates, expression, predictions, probabilities,
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
  knn <- buildKnn(coordinates, k)
  
  # Compute enhanced confidence scores
  enhanced_scores <- enhancedConfidence(
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
updatePreds <- function(STEAM.obj, iteration_results) {
  
  if (is.null(iteration_results$corrections) || nrow(iteration_results$corrections) == 0) {
    return(NULL)
  }
  
  corrections <- iteration_results$corrections
  updated_predictions <- list()
  
  # Update predictions for each fold
  for (f in seq_along(STEAM.obj$nested$ncv$outer_result)) {
    fold_result <- STEAM.obj$nested$ncv$outer_result[[f]]
    if (is.null(fold_result$preds)) next
    
    # Get current predictions    # New explicit approach
    STEAM.obj <- STEAM_anchor(
      STEAM.obj = my_steam_obj,
      ground_truth = my_ground_truth_labels,  # Explicit parameter
      iterative = TRUE,
      verbose = TRUE
    )
    
    # No ground truth (alternative validation)
    STEAM.obj <- STEAM_anchor(
      STEAM.obj = my_steam_obj,
      ground_truth = NULL,  # Explicit - will use alternative validation
      iterative = TRUE,
      verbose = TRUE
    )
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
updateNeighbors <- function(STEAM.obj, updated_predictions) {
  
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

#' Run Single Iteration Analysis
#' 
#' @param STEAM.obj STEAM object
#' @param expression_matrix Expression matrix
#' @param ground_truth Ground truth labels
#' @param k_grid K values to test
#' @param anchor_cut_grid Anchor cut thresholds
#' @param consensus_grid Consensus thresholds
#' @param min_confidence_gain Minimum confidence gain
#' @param iteration Current iteration number
#' @param verbose Print messages
#' @return Iteration results
#' @keywords internal
runIteration <- function(
    STEAM.obj,
    expression_matrix,
    ground_truth,
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
    fold_misclassified <- remainingCells(STEAM.obj, f, corrected_cells, ground_truth)
    
    if (length(fold_misclassified) == 0) {
      if (verbose && f <= 2) message(sprintf("Fold %d: No remaining misclassified cells", f))
      next
    }
    
    # Use updated predictions if available from previous iteration
    current_predictions <- foldPreds(STEAM.obj, f, iteration)
    
    fold_result <- processFoldIterative(
      fold_idx = f,
      ncv = ncv,
      coords_all = STEAM.obj$spatial,
      expr_data_prepared = prepareExpression(expression_matrix, STEAM.obj$spatial),
      ground_truth = ground_truth,
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
#' @param ground_truth Ground truth labels
#' @return Vector of remaining misclassified cell IDs
#' @keywords internal
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

#' Get Current Fold Predictions
#' 
#' @param STEAM.obj STEAM object
#' @param fold_idx Fold index
#' @param iteration Current iteration
#' @return Current predictions for the fold
#' @keywords internal
foldPreds <- function(STEAM.obj, fold_idx, iteration) {
  
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
#' @param ground_truth Ground truth labels
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
processFoldIterative <- function(
    fold_idx,
    ncv,
    coords_all,
    expr_data_prepared,
    ground_truth,
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
  best_params <- tuneIterative(
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
  corrections <- flipIterative(
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
flipIterative <- function(
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
  knn <- buildKnn(all_coords, k)
  
  # Compute enhanced confidence scores using updated prediction context
  enhanced_scores <- confidenceIterative(
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
confidenceIterative <- function(
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
correctionChains <- function(previous_corrections, current_corrections) {
  
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
#' @param ground_truth Ground truth labels
#' @return List with improvement metrics
#' @keywords internal
iterationImprovement <- function(STEAM.obj, iter_results, prev_accuracy_by_fold, ground_truth) {
  
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
    fold_cells <- rownames(preds_df)
    true_labels <- ground_truth[fold_cells]
    
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
tuneIterative <- function(
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
        
        score <- evalIterative(
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
evalIterative <- function(
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
  knn <- buildKnn(coordinates, k)
  
  # Compute confidence scores
  enhanced_scores <- enhancedConfidence(
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

#' Extract Anchor Corrections from STEAM Analysis Results
#' 
#' @param spatial_anchor_analysis Results from STEAM anchor analysis
#' @param STEAM.obj STEAM object for context
#' @param ground_truth Ground truth labels (optional)
#' @return Data frame of corrections
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
  
  # Check for iterative analysis results first
  if (!is.null(STEAM.obj$iterative_anchor_analysis)) {
    iter_analysis <- STEAM.obj$iterative_anchor_analysis
    
    if (iter_analysis$method == "iterative_with_ground_truth") {
      # Extract from iteration results
      for (i in seq_along(iter_analysis$iteration_results)) {
        if (!is.null(iter_analysis$iteration_results[[i]]$corrections)) {
          iter_corrections <- iter_analysis$iteration_results[[i]]$corrections
          iter_corrections$iteration <- i
          corrections <- rbind(corrections, iter_corrections[, intersect(colnames(corrections), colnames(iter_corrections))])
        }
      }
    } else {
      # Extract from correction history (no ground truth method)
      for (i in seq_along(iter_analysis$correction_history)) {
        if (!is.null(iter_analysis$correction_history[[i]])) {
          iter_corrections <- iter_analysis$correction_history[[i]]
          iter_corrections$iteration <- i
          # Map column names
          if ("suggested_correction" %in% colnames(iter_corrections)) {
            iter_corrections$suggested_correction <- iter_corrections$suggested_correction
          }
          if ("original_prediction" %in% colnames(iter_corrections)) {
            iter_corrections$original_pred <- iter_corrections$original_prediction
          }
          corrections <- rbind(corrections, iter_corrections[, intersect(colnames(corrections), colnames(iter_corrections))])
        }
      }
    }
  }
  
  # If no iterative results, try single-pass results
  if (nrow(corrections) == 0 && !is.null(spatial_anchor_analysis$corrections)) {
    single_corrections <- spatial_anchor_analysis$corrections
    
    # Map to standard format
    if (nrow(single_corrections) > 0) {
      std_corrections <- data.frame(
        cell_id = single_corrections$cell_id %||% rownames(single_corrections),
        original_pred = single_corrections$original_pred %||% single_corrections$from_label,
        suggested_correction = single_corrections$suggested_correction %||% single_corrections$to_label,
        true_label = single_corrections$true_label %||% NA,
        correction_correct = single_corrections$correction_correct %||% 
                           single_corrections$is_correct %||% TRUE,
        fold = single_corrections$fold %||% 1,
        stringsAsFactors = FALSE
      )
      corrections <- rbind(corrections, std_corrections)
    }
  }
  
  # Add true labels if missing and ground truth is available
  if (!is.null(ground_truth) && any(is.na(corrections$true_label))) {
    # Use explicit ground truth instead of extracting from nested CV
    fold_cells <- corrections$cell_id
    corrections$true_label[is.na(corrections$true_label)] <- ground_truth[fold_cells[is.na(corrections$true_label)]]
    
    # Update correction_correct if we have true labels
    if (!all(is.na(corrections$true_label))) {
      corrections$correction_correct <- corrections$suggested_correction == corrections$true_label
    }
  } else if (!is.null(STEAM.obj) && any(is.na(corrections$true_label))) {
    # Fallback: Extract true labels from nested CV results (deprecated approach)
    for (f in seq_along(STEAM.obj$nested$ncv$outer_result)) {
      fold_result <- STEAM.obj$nested$ncv$outer_result[[f]]
      if (!is.null(fold_result$preds) && "testy" %in% colnames(fold_result$preds)) {
        true_labels <- setNames(fold_result$preds$testy, rownames(fold_result$preds))
        
        # Update missing true labels
        fold_mask <- corrections$fold == f & is.na(corrections$true_label)
        corrections$true_label[fold_mask] <- true_labels[corrections$cell_id[fold_mask]]
      }
    }
    
    # Update correction_correct if we have true labels
    if (!all(is.na(corrections$true_label))) {
      corrections$correction_correct <- corrections$suggested_correction == corrections$true_label
    }
  }
  
  return(corrections)
}

#' Extract Anchor Summary Statistics
#' 
#' @param spatial_anchor_analysis Results from STEAM anchor analysis
#' @param STEAM.obj STEAM object for context
#' @param ground_truth Ground truth labels (optional)
#' @return List of summary statistics
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
  } else if (!is.null(STEAM.obj$nested$ncv)) {
    # Fallback: use testy column (deprecated)
    for (f in seq_along(STEAM.obj$nested$ncv$outer_result)) {
      fold_result <- STEAM.obj$nested$ncv$outer_result[[f]]
      if (!is.null(fold_result$preds) && all(c("testy", "predy") %in% colnames(fold_result$preds))) {
        fold_misclassified <- sum(fold_result$preds$testy != fold_result$preds$predy, na.rm = TRUE)
        total_misclassified <- total_misclassified + fold_misclassified
        total_cells <- total_cells + nrow(fold_result$preds)
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
