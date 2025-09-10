#' STEAM Correction for Misclassified Cells
#'
#' Applies hierarchical neighbor-based correction to misclassified cells
#' across all folds in a STEAM object. Uses multi-level voting (primary,
#' secondary, tertiary) with configurable consensus thresholds.
#'
#' @param STEAM.obj A STEAM object with nested CV results and predictions.
#' @param k Number of nearest neighbors to consider (default: 8).
#' @param min_consensus Minimum proportion of neighbors that must agree for
#'   1st-choice corrections (default: 0.6).
#' @param max_iterations Maximum number of correction iterations (default: 10).
#' @param max_voting_levels Maximum number of voting levels to try (default: 3).
#' @param min_secondary_consensus Minimum consensus for 2nd/3rd choice
#'   corrections (default: 0.3).
#' @param verbose Logical; print progress to console (default: TRUE).
#'
#' @return An updated STEAM object with hierarchical corrections stored in
#'   \code{$spatial_anchor_analysis$correction_results}.
#'
#' @export
STEAMCorrection <- function(STEAM.obj,
                            k = 8,
                            min_consensus = 0.6,
                            max_iterations = 10,
                            max_voting_levels = 3,
                            min_secondary_consensus = 0.3,
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
    stuck_cells <- c()  # Track cells that haven't been corrected for multiple iterations

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

      # Track stuck cells (misclassified for multiple iterations)
      if (iteration >= 2) {
        # Cells that were misclassified in previous iteration and still misclassified
        previous_misclassified <- correction_history[[iteration]]$remaining_misclassified_cells
        if (!is.null(previous_misclassified)) {
          newly_stuck <- intersect(current_misclassified, previous_misclassified)
          stuck_cells <- union(stuck_cells, newly_stuck)
        }
      }

      correction_suggestions <- generateCorrections(
        fold_preds, fold_coordinates, knn_result, min_consensus, current_misclassified,
        max_voting_levels, min_secondary_consensus, stuck_cells, iteration
      )

      if (length(correction_suggestions) == 0) {
        if (verbose) cat("No valid corrections found - stopping\n")
        break
      }

      if (verbose) {
        cat(sprintf("Generated %d correction suggestions\n", length(correction_suggestions)))
        if (length(stuck_cells) > 0) {
          cat(sprintf("Stuck cells (multiple failed attempts): %d\n", length(stuck_cells)))
        }

        # Show voting level breakdown
        if (length(correction_suggestions) > 0) {
          voting_levels <- sapply(correction_suggestions, function(x) x$voting_level)
          level_counts <- table(voting_levels)
          for (level in sort(unique(voting_levels))) {
            cat(sprintf("  %s suggestions: %d\n",
                        ifelse(level == 1, "1st choice", ifelse(level == 2, "2nd choice", paste0(level, "th choice"))),
                        level_counts[as.character(level)]))
          }
        }
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
                status <- if (is_correct) "OK" else "FAIL"
                stuck_indicator <- if (correction$is_stuck_cell) " [STUCK]" else ""
                cat(sprintf("    %s: %s -> %s (%s, %.2f consensus) %s%s\n",
                            correction$cell_id, correction$current_prediction,
                            correction$new_prediction, correction$voting_rank,
                            correction$consensus, status, stuck_indicator))


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
        remaining_misclassified_cells = current_misclassified,
        predictions_snapshot = fold_preds,
        stuck_cells = stuck_cells
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
      cat(sprintf("Misclassified cells reduced: %d -> %d\n",
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
      cat(sprintf("Fold %d: %.1f%% -> %.1f%% (%.1f%% improvement, %d corrections)\n",
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




#' Extract Correction Analysis Results (Internal)
#'
#' Helper function to extract correction results from a STEAM object
#' after running \code{\link{STEAMCorrection}}.
#'
#' @param STEAM.obj A STEAM object that may contain
#'   \code{$spatial_anchor_analysis$correction_results}.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{found} (logical): Whether results were found.
#'     \item \code{type} (character): Type of results.
#'     \item \code{data} (list): Details including number of folds,
#'       total corrections, and per-fold results.
#'   }
#'
#' @keywords internal
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






#' Generate Correction Suggestions (Internal)
#'
#' Internal helper that suggests label corrections for misclassified
#' cells using neighbor voting with hierarchical consensus thresholds.
#'
#' @param predictions Vector of current predictions.
#' @param coordinates Matrix or data.frame of spatial coordinates.
#' @param knn_result Result from \code{RANN::nn2()}.
#' @param min_consensus Minimum consensus threshold for first-choice corrections.
#' @param target_cells Character vector of cells to process (default: all).
#' @param max_voting_levels Maximum number of voting levels to consider.
#' @param min_secondary_consensus Consensus threshold for secondary/tertiary votes.
#' @param stuck_cells Cells repeatedly misclassified in earlier iterations.
#' @param iteration Current iteration number.
#'
#' @return A list of suggested corrections with fields:
#'   \code{cell_id}, \code{current_prediction}, \code{new_prediction},
#'   \code{consensus}, \code{voting_level}, \code{voting_rank},
#'   \code{is_stuck_cell}.
#'
#' @keywords internal
generateCorrections <- function(predictions, coordinates, knn_result, min_consensus, target_cells = NULL,
                                max_voting_levels = 3, min_secondary_consensus = 0.3,
                                stuck_cells = NULL, iteration = 1) {

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

    # Get vote counts sorted by frequency
    pred_counts <- table(neighbor_preds)
    pred_counts_sorted <- sort(pred_counts, decreasing = TRUE)

    current_pred <- predictions[cell_name]
    is_stuck_cell <- !is.null(stuck_cells) && cell_name %in% stuck_cells

    # Try multiple voting levels
    max_levels <- min(max_voting_levels, length(pred_counts_sorted))
    if (max_levels > 0) {
      for (voting_level in seq_len(max_levels)) {
        candidate_pred <- names(pred_counts_sorted)[voting_level]
        candidate_consensus <- pred_counts_sorted[voting_level] / length(neighbor_preds)

        # Skip if same as current prediction
        if (candidate_pred == current_pred) next

        # Determine consensus threshold based on voting level and iteration
        if (voting_level == 1) {
          # Primary choice: use standard threshold
          threshold <- min_consensus
        } else {
          # Secondary+ choices: use lower threshold, especially for stuck cells or later iterations
          if (is_stuck_cell || iteration > 3) {
            threshold <- min_secondary_consensus * 0.7  # Even more lenient for difficult cases
          } else {
            threshold <- min_secondary_consensus
          }
        }

        # Accept this correction if consensus meets threshold
        if (candidate_consensus >= threshold) {
          corrections <- append(corrections, list(list(
            cell_id = cell_name,
            current_prediction = current_pred,
            new_prediction = candidate_pred,
            consensus = candidate_consensus,
            voting_level = voting_level,
            voting_rank = paste0(voting_level,
                                 if (voting_level == 1) "st" else if (voting_level == 2) "nd" else if (voting_level == 3) "rd" else "th",
                                 " choice"),
            is_stuck_cell = is_stuck_cell
          )))
          break  # Take first valid correction for this cell
        }
      }
    }
  }

  if (length(corrections) > 0) {
    # Sort by consensus within voting level (primary choices first, then by consensus)
    corrections <- corrections[order(
      sapply(corrections, function(x) x$voting_level),     # Primary sort: voting level
      -sapply(corrections, function(x) x$consensus),       # Secondary sort: consensus (descending)
      decreasing = FALSE
    )]
  }

  return(corrections)
}

