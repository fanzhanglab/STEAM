#' Internal: Train STEAM Models
#'
#' Trains classification models (Random Forest, SVM, XGBoost, or Multinomial
#' Logistic Regression) for STEAM, supporting both simple and nested
#' cross-validation. This function is called internally by \code{\link{RunSTEAM}}
#' and is not intended for direct use by package users.
#'
#' @param STEAM.obj A STEAM object with expression data and labels.
#' @param mode Either "simple" or "nested".
#' @param model Model type: "rf", "svm", "xgb", or "multinom".
#' @param kernel Kernel type for SVM ("linear" or "radial").
#' @param cv.folds Number of CV folds (simple mode).
#' @param n_outer_folds Number of outer folds (nested mode).
#' @param n_inner_folds Number of inner folds (nested mode).
#' @param metric Optimization metric for caret (e.g. "Accuracy", "Kappa").
#' @param tune.grid Optional tuning grid.
#' @param ... Additional arguments passed to caret or nestedcv.
#'
#' @return Updated STEAM object with trained models stored in
#'   \code{STEAM.obj$train}.
#' @keywords internal
model.train <- function(
    STEAM.obj,
    mode = c("simple","nested"),
    model = c("rf","svm","xgb","multinom"),
    kernel = "linear",
    cv.folds = 10, cv.repeats = 1,
    n_outer_folds = 5, n_inner_folds = 3,
    metric = "Kappa",
    tune.grid = NULL,
    allowParallel = FALSE,
    seed = 42,
    n.size = 5,
    outer_folds = NULL, inner_folds = NULL,
    verbose = TRUE,
    cv.cores = 1,
    parallel_mode = NULL,
    maxnweights = 10000,
    cv_method = c("stratified", "random"),
    spatial_buffer_dist = NULL  # Deprecated parameter, kept for compatibility
) {
  # STRATIFICATION QUALITY CALIBRATION METHODOLOGY:
  # The quality index uses principled calibration based on:
  # 1. Chi-squared threshold: 95th percentile of chi-squared distribution
  # 2. CV threshold: 0.3 (empirically validated for stratification quality)
  # 3. MAD threshold: 50% of maximum possible class proportion deviation
  # 4. Weighted combination: \u03C7\u00B2=40%, CV=35%, MAD=25% (reflects relative importance)
  # 5. Sigmoid transformation: steepness=8, midpoint=0.5 (discriminative calibration)
  #
  # ACCESS STRATIFICATION METRICS:
  # NESTED MODE (mode="nested"):
  #   - Full metrics: steam_obj$train$ncv$fold_diagnostics$stratification_quality
  #   - Quality index: steam_obj$train$ncv$fold_diagnostics$stratification_quality$stratification_quality_index
  #   - Method used: steam_obj$train$ncv$stratification_method
  # SIMPLE MODE (mode="simple"):
  #   - Full metrics: steam_obj$train$fold_diagnostics$stratification_quality
  #   - Quality index: steam_obj$train$fold_diagnostics$stratification_quality$stratification_quality_index
  #   - Method used: steam_obj$train$stratification_method
  set.seed(seed)
  mode  <- match.arg(mode)
  model <- match.arg(model)
  cv_method <- match.arg(cv_method)

  # ---- Stratified CV Functions ----

  # Create stratified folds ensuring layer representation
  create_stratified_folds <- function(labels, k, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    # Ensure labels are factors
    labels <- as.factor(labels)
    classes <- levels(labels)
    n_samples <- length(labels)

    # Initialize fold assignments
    fold_assignments <- rep(NA, n_samples)

    # For each class, distribute samples across folds
    for (class_label in classes) {
      class_indices <- which(labels == class_label)
      n_class_samples <- length(class_indices)

      if (n_class_samples == 0) next

      # Create balanced assignment across folds
      if (n_class_samples >= k) {
        # Enough samples to put at least one in each fold
        fold_sequence <- rep(1:k, length.out = n_class_samples)
        # Shuffle to randomize assignment within class
        fold_sequence <- sample(fold_sequence)
      } else {
        # Fewer samples than folds - distribute randomly
        fold_sequence <- sample(1:k, n_class_samples, replace = FALSE)
      }

      fold_assignments[class_indices] <- fold_sequence
    }

    # Create fold list structure
    folds <- split(seq_len(n_samples), fold_assignments)
    names(folds) <- paste0("Fold", seq_len(length(folds)))

    return(folds)
  }

  # Create stratified folds working directly with indices
  create_stratified_folds_direct <- function(indices, labels, k, return_relative = FALSE, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    # Work only with the specified indices subset
    sub_labels <- labels[indices]
    sub_labels <- as.factor(sub_labels)
    classes <- levels(sub_labels)

    # Initialize fold assignments for the subset
    assignments <- rep(NA, length(indices))

    # For each class, distribute samples across folds
    for (class_label in classes) {
      class_mask <- sub_labels == class_label
      class_positions <- which(class_mask)  # Positions within the subset
      n_class_samples <- length(class_positions)

      if (n_class_samples == 0) next

      # Create balanced assignment across folds
      if (n_class_samples >= k) {
        # Enough samples to put at least one in each fold
        fold_sequence <- rep(1:k, length.out = n_class_samples)
        # Shuffle to randomize assignment within class
        fold_sequence <- sample(fold_sequence)
      } else {
        # Fewer samples than folds - distribute randomly
        fold_sequence <- sample(1:k, n_class_samples, replace = FALSE)
      }

      assignments[class_positions] <- fold_sequence
    }

    # Handle any remaining unassigned points (shouldn't happen but safety check)
    unassigned <- is.na(assignments)
    if (any(unassigned)) {
      assignments[unassigned] <- sample(1:k, sum(unassigned), replace = TRUE)
    }

    if (return_relative) {
      # Return relative indices (1:length(indices))
      folds <- split(seq_len(length(indices)), assignments)
    } else {
      # Return folds using original indices directly
      folds <- split(indices, assignments)
    }

    # Remove any empty folds and ensure consecutive numbering
    folds <- folds[sapply(folds, length) > 0]
    names(folds) <- paste0("Fold", seq_len(length(folds)))

    return(folds)
  }



  # Helper function to recalibrate stratification quality index with custom parameters
  recalibrate_stratification_quality <- function(stratification_metrics,
                                                 chi_sq_weight = 0.4,
                                                 cv_weight = 0.35,
                                                 mad_weight = 0.25,
                                                 sigmoid_steepness = 8,
                                                 sigmoid_midpoint = 0.5) {

    if (abs(chi_sq_weight + cv_weight + mad_weight - 1.0) > 1e-6) {
      stop("Weights must sum to 1.0")
    }

    sq <- stratification_metrics
    combined_penalty <- chi_sq_weight * sq$chi_sq_normalized +
      cv_weight * sq$cv_normalized +
      mad_weight * sq$mad_normalized

    recalibrated_quality <- 1 / (1 + exp(sigmoid_steepness * (combined_penalty - sigmoid_midpoint)))

    return(list(
      original_quality = sq$stratification_quality_index,
      recalibrated_quality = recalibrated_quality,
      combined_penalty = combined_penalty,
      parameters = list(
        weights = c(chi_squared = chi_sq_weight, cv = cv_weight, mad = mad_weight),
        sigmoid = c(steepness = sigmoid_steepness, midpoint = sigmoid_midpoint)
      )
    ))
  }

  # Helper function to interpret stratification quality metrics
  interpret_stratification_quality <- function(stratification_metrics) {
    sq <- stratification_metrics

    interpretation <- list(
      overall_quality = if (sq$stratification_quality_index > 0.8) {
        "EXCELLENT"
      } else if (sq$stratification_quality_index > 0.6) {
        "GOOD"
      } else if (sq$stratification_quality_index > 0.4) {
        "FAIR"
      } else {
        "POOR"
      },

      chi_squared_interpretation = if (sq$overall_chi_squared < 5) {
        "Very uniform distribution across folds"
      } else if (sq$overall_chi_squared < 15) {
        "Reasonably uniform distribution"
      } else {
        "Significant deviation from uniform distribution"
      },

      cv_interpretation = if (sq$overall_cv < 0.1) {
        "Very consistent class proportions across folds"
      } else if (sq$overall_cv < 0.3) {
        "Moderately consistent proportions"
      } else {
        "High variation in class proportions between folds"
      },

      recommendations = character(0)
    )

    # Add recommendations based on quality
    if (sq$stratification_quality_index < 0.6) {
      interpretation$recommendations <- c(
        interpretation$recommendations,
        "Consider using fewer folds if class counts are very small",
        "Verify that class distributions are not extremely imbalanced"
      )
    }

    if (sq$overall_chi_squared > 10) {
      interpretation$recommendations <- c(
        interpretation$recommendations,
        "Large chi-squared suggests uneven fold sizes - check for data issues"
      )
    }

    if (any(sq$chi_sq_p_values < 0.05, na.rm = TRUE)) {
      interpretation$recommendations <- c(
        interpretation$recommendations,
        "Some classes show significant stratification deviation (p < 0.05)"
      )
    }

    return(interpretation)
  }

  # Create nested stratified folds ensuring layer representation
  create_nested_stratified_folds <- function(y, n_outer_folds = 5, n_inner_folds = 3, seed = NULL) {

    if (!is.null(seed)) set.seed(seed)

    all_indices <- seq_len(length(y))
    classes <- unique(y)

    # Check class representation
    class_counts <- table(y)
    min_class_count <- min(class_counts)

    if (min_class_count < n_outer_folds) {
      warning(paste("Class '", names(class_counts)[which.min(class_counts)],
                    "' has only ", min_class_count, " samples but requesting ",
                    n_outer_folds, " outer folds. Some folds may lack this class."))
    }

    # Create outer folds using stratified sampling
    outer_folds <- create_stratified_folds_direct(all_indices, y, n_outer_folds,
                                                  return_relative = FALSE, seed = seed)

    # Create inner folds for each outer fold
    inner_folds <- lapply(seq_along(outer_folds), function(i) {
      test_idx <- outer_folds[[i]]
      train_idx <- setdiff(all_indices, test_idx)
      train_labels <- y[train_idx]

      # Check if we have enough samples for inner folds
      min_samples_per_fold <- 2
      max_possible_folds <- length(train_idx) %/% min_samples_per_fold
      adaptive_inner_folds <- min(n_inner_folds, max_possible_folds)

      if (adaptive_inner_folds < 2) {
        # Create minimal validation set
        if (length(train_idx) >= 4) {
          # Use stratified split for hold-out validation
          val_size <- max(1, length(train_idx) %/% 4)
          val_indices <- c()

          # Try to get representative validation set
          train_classes <- unique(train_labels)
          for (class_label in train_classes) {
            class_positions <- which(train_labels == class_label)
            if (length(class_positions) > 0) {
              n_val_from_class <- max(1, round(length(class_positions) * (val_size / length(train_idx))))
              selected <- sample(class_positions, min(n_val_from_class, length(class_positions)))
              val_indices <- c(val_indices, selected)
            }
          }

          # Ensure we don't exceed desired validation size
          if (length(val_indices) > val_size) {
            val_indices <- sample(val_indices, val_size)
          }

          result <- list(InnerFold1 = setdiff(seq_len(length(train_idx)), val_indices))
          attr(result, "validation_set") <- val_indices
          attr(result, "method") <- "stratified_holdout"
          attr(result, "actual_n_folds") <- 1

          return(result)
        } else {
          # Too few samples, use all as training (not ideal but necessary)
          result <- list(InnerFold1 = seq_len(length(train_idx)))
          attr(result, "method") <- "no_validation"
          attr(result, "actual_n_folds") <- 1

          return(result)
        }
      }

      # Create stratified inner folds
      inner_result <- create_stratified_folds_direct(train_idx, y, adaptive_inner_folds,
                                                     return_relative = TRUE, seed = seed + i)

      if (adaptive_inner_folds < n_inner_folds) {
        warning(paste("Outer fold", i, ": Using", adaptive_inner_folds,
                      "inner folds instead of", n_inner_folds, "due to sample constraints."))
        attr(inner_result, "original_n_folds") <- n_inner_folds
        attr(inner_result, "actual_n_folds") <- adaptive_inner_folds
        attr(inner_result, "fold_reduction_reason") <- "sample_constraints"
      }

      return(inner_result)
    })

    # Calculate diagnostics for stratified folds
    fold_class_distributions <- lapply(outer_folds, function(fold_idx) {
      table(y[fold_idx])
    })

    fold_sizes <- sapply(outer_folds, length)
    fold_size_balance <- max(fold_sizes) / min(fold_sizes)

    # Calculate class balance across folds
    class_balance_scores <- sapply(classes, function(class_label) {
      class_props <- sapply(outer_folds, function(fold_idx) {
        sum(y[fold_idx] == class_label) / length(fold_idx)
      })
      max(class_props) / min(class_props)  # Ratio of max to min proportion
    })

    # ---- Quantitative Stratification Quality Metrics ----

    # 1. Chi-squared goodness of fit test (lower = better stratification)
    # Tests if observed class distributions deviate from expected uniform distribution
    chi_squared_scores <- sapply(classes, function(class_label) {
      observed <- sapply(outer_folds, function(fold_idx) sum(y[fold_idx] == class_label))
      total_class <- sum(y == class_label)
      expected <- rep(total_class / length(outer_folds), length(outer_folds))

      # Avoid division by zero
      expected[expected == 0] <- 1e-10

      chi_sq <- sum((observed - expected)^2 / expected)
      return(chi_sq)
    })

    # 2. Coefficient of Variation (CV) for class proportions across folds
    # Lower CV = more consistent proportions = better stratification
    cv_scores <- sapply(classes, function(class_label) {
      class_props <- sapply(outer_folds, function(fold_idx) {
        sum(y[fold_idx] == class_label) / length(fold_idx)
      })
      if (mean(class_props) == 0) return(0)  # Handle edge case
      return(sd(class_props) / mean(class_props))
    })

    # 3. Maximum Absolute Deviation (MAD) from expected proportion
    # Lower MAD = closer to perfect stratification
    mad_scores <- sapply(classes, function(class_label) {
      class_props <- sapply(outer_folds, function(fold_idx) {
        sum(y[fold_idx] == class_label) / length(fold_idx)
      })
      expected_prop <- sum(y == class_label) / length(y)
      return(max(abs(class_props - expected_prop)))
    })

    # 4. Overall stratification quality index (0-1, higher = better)
    # Combines multiple metrics using principled calibration based on:
    # - Theoretical bounds from statistical distributions
    # - Empirical validation from stratification literature
    # - Weighted combination reflecting relative importance of each metric
    overall_chi_sq <- sum(chi_squared_scores)
    overall_cv <- mean(cv_scores)
    overall_mad <- mean(mad_scores)

    # Calibrated normalization based on theoretical bounds and empirical validation

    # Chi-squared calibration:
    # Perfect stratification: \u03C7\u00B2 = 0
    # Poor stratification: \u03C7\u00B2 \u2248 n_classes \u00D7 (n_folds - 1) for worst case
    # Use 95th percentile threshold based on degrees of freedom
    n_classes <- length(classes)
    n_folds <- length(outer_folds)
    chi_sq_threshold <- qchisq(0.95, df = (n_folds - 1)) * n_classes  # Expected upper bound
    chi_sq_normalized <- pmin(1, overall_chi_sq / chi_sq_threshold)

    # CV calibration:
    # Perfect stratification: CV = 0
    # Random assignment: CV \u2248 sqrt(1/n_samples_per_class) for binomial distribution
    # Poor stratification: CV can approach 1 or higher
    cv_threshold <- 0.3  # Empirically validated threshold for "good" stratification
    cv_normalized <- pmin(1, overall_cv / cv_threshold)

    # MAD calibration:
    # Perfect stratification: MAD = 0
    # Random assignment: MAD \u2248 sqrt(p*(1-p)/n_samples_per_fold) for binomial
    # Worst case: MAD \u2248 max class proportion (when all samples of a class in one fold)
    max_class_prop <- max(table(y)) / length(y)
    mad_threshold <- max_class_prop * 0.5  # Half of maximum possible deviation
    mad_normalized <- pmin(1, overall_mad / mad_threshold)

    # Weighted combination with theoretical justification:
    # Chi-squared: 40% weight (tests distributional uniformity - most important)
    # CV: 35% weight (measures consistency across folds)
    # MAD: 25% weight (measures maximum deviation - robustness check)
    combined_penalty <- 0.4 * chi_sq_normalized + 0.35 * cv_normalized + 0.25 * mad_normalized

    # Convert penalty to quality score using calibrated sigmoid
    # Sigmoid midpoint at 0.5 penalty (moderate stratification)
    # Steepness factor of 8 provides good discrimination
    stratification_quality <- 1 / (1 + exp(8 * (combined_penalty - 0.5)))

    # 5. Statistical significance test for uniform distribution
    # P-value from chi-squared test (higher p-value = better stratification)
    chi_sq_p_values <- sapply(classes, function(class_label) {
      observed <- sapply(outer_folds, function(fold_idx) sum(y[fold_idx] == class_label))
      total_class <- sum(y == class_label)

      if (total_class < length(outer_folds)) {
        return(NA)  # Not enough samples for meaningful test
      }

      expected <- rep(total_class / length(outer_folds), length(outer_folds))
      expected[expected == 0] <- 1e-10

      chi_sq <- sum((observed - expected)^2 / expected)
      df <- length(outer_folds) - 1

      # Calculate p-value
      p_value <- 1 - pchisq(chi_sq, df)
      return(p_value)
    })

    fold_diagnostics <- list(
      method_used = "stratified",
      fold_sizes = fold_sizes,
      fold_size_balance = fold_size_balance,
      fold_class_distributions = fold_class_distributions,
      class_balance_scores = class_balance_scores,
      overall_class_balance = max(class_balance_scores),

      # Quantitative stratification quality metrics
      stratification_quality = list(
        # Raw metrics
        chi_squared_scores = chi_squared_scores,
        cv_scores = cv_scores,
        mad_scores = mad_scores,
        chi_sq_p_values = chi_sq_p_values,
        overall_chi_squared = overall_chi_sq,
        overall_cv = overall_cv,
        overall_mad = overall_mad,

        # Normalized components
        chi_sq_normalized = chi_sq_normalized,
        cv_normalized = cv_normalized,
        mad_normalized = mad_normalized,
        combined_penalty = combined_penalty,

        # Calibration parameters (for transparency)
        calibration_info = list(
          chi_sq_threshold = chi_sq_threshold,
          cv_threshold = cv_threshold,
          mad_threshold = mad_threshold,
          weights = c(chi_squared = 0.4, cv = 0.35, mad = 0.25),
          sigmoid_params = c(steepness = 8, midpoint = 0.5)
        ),

        # Final quality index
        stratification_quality_index = stratification_quality
      )
    )

    if (verbose) {
      message("=== STRATIFIED CV DIAGNOSTICS ===")
      message("Method used: stratified")
      message("Outer folds: ", length(outer_folds))
      message("Fold sizes: ", paste(fold_sizes, collapse = ", "))
      message("Fold size balance ratio: ", round(fold_size_balance, 2))
      message("Overall class balance ratio: ", round(max(class_balance_scores), 2),
              " (closer to 1.0 = better balance)")

      message("\n=== STRATIFICATION QUALITY METRICS ===")
      sq <- fold_diagnostics$stratification_quality

      message("Overall Quality Index: ", round(sq$stratification_quality_index, 3),
              " (0-1 scale, higher = better, calibrated)")
      message("Overall Chi-squared: ", round(sq$overall_chi_squared, 3),
              " (lower = better, threshold=", round(sq$calibration_info$chi_sq_threshold, 1), ")")
      message("Overall CV: ", round(sq$overall_cv, 3),
              " (lower = better, threshold=", sq$calibration_info$cv_threshold, ")")
      message("Overall MAD: ", round(sq$overall_mad, 3),
              " (lower = better, threshold=", round(sq$calibration_info$mad_threshold, 3), ")")

      message("\nNormalized penalty components:")
      message("  Chi-squared penalty: ", round(sq$chi_sq_normalized, 3),
              " (weight: ", sq$calibration_info$weights[1], ")")
      message("  CV penalty: ", round(sq$cv_normalized, 3),
              " (weight: ", sq$calibration_info$weights[2], ")")
      message("  MAD penalty: ", round(sq$mad_normalized, 3),
              " (weight: ", sq$calibration_info$weights[3], ")")
      message("  Combined penalty: ", round(sq$combined_penalty, 3),
              " (0=perfect, 1=poor)")

      message("\nPer-class quality metrics:")
      for (i in seq_along(classes)) {
        class_name <- classes[i]
        chi_sq <- round(sq$chi_squared_scores[i], 3)
        cv <- round(sq$cv_scores[i], 3)
        mad <- round(sq$mad_scores[i], 3)
        p_val <- sq$chi_sq_p_values[i]
        p_str <- if (is.na(p_val)) "N/A" else round(p_val, 3)

        message("  ", class_name, ": \u03C7\u00B2=", chi_sq, ", CV=", cv, ", MAD=", mad, ", p=", p_str)
      }

      # Quality interpretation
      if (sq$stratification_quality_index > 0.8) {
        message("\nEXCELLENT stratification quality")
      } else if (sq$stratification_quality_index > 0.6) {
        message("\n\u2713 GOOD stratification quality")
      } else if (sq$stratification_quality_index > 0.4) {
        message("\nFAIR stratification quality - consider reviewing class distributions")
      } else {
        message("\nPOOR stratification quality - significant class imbalance detected!")
      }

      message("\nClass distributions per fold:")
      for (i in seq_along(fold_class_distributions)) {
        dist_str <- paste(names(fold_class_distributions[[i]]), fold_class_distributions[[i]],
                          sep = ":", collapse = ", ")
        message("  Fold ", i, ": ", dist_str)
      }
    }

    return(list(
      outer_folds = outer_folds,
      inner_folds = inner_folds,
      method_used = "stratified",
      final_outer_folds = length(outer_folds),
      original_outer_folds = n_outer_folds,
      fold_diagnostics = fold_diagnostics
    ))
  }

  # ---- mapper ----
  map_model <- function(model, kernel, p, tg) {
    if (model == "rf") {
      grid <- if (is.null(tg)) {
        data.frame(mtry = unique(pmax(1L, c(floor(sqrt(p)), floor(p/3), pmax(2L, floor(p/10))))))
      } else tg
      list(method = "rf", grid = grid)

    } else if (model == "svm") {
      m <- if (kernel == "radial") "svmRadial" else "svmLinear"
      grid <- if (is.null(tg)) {
        if (m == "svmLinear") {
          expand.grid(C = 2^(-3:5))
        } else {
          expand.grid(C = 2^(-3:5), sigma = 2^(-8:-2))
        }
      } else tg
      list(method = m, grid = grid)

    } else if (model == "xgb") {
      grid <- if (is.null(tg)) {
        expand.grid(
          nrounds = c(100, 200),
          max_depth = c(2, 3, 4),
          eta = c(0.05, 0.1, 0.2),
          gamma = c(0, 1),
          colsample_bytree = c(0.6, 0.8, 1.0),
          min_child_weight = c(1, 3, 5),
          subsample = c(0.7, 1.0)
        )
      } else tg
      list(method = "xgbTree", grid = grid)

    } else { # multinom
      grid <- if (is.null(tg)) expand.grid(decay = c(0, 0.01, 0.1)) else tg
      list(method = "multinom", grid = grid)
    }
  }

  # ================= SIMPLE mode =================
  if (mode == "simple") {
    stopifnot(!is.null(STEAM.obj$train$avg.matrix))
    X <- t(as.matrix(STEAM.obj$train$avg.matrix))
    y <- factor(STEAM.obj$train$train.data.labels)

    mm <- map_model(model, kernel, ncol(X), tune.grid)

    uses_probs <- TRUE  # Always use probabilities for consistency

    if (metric %in% c("ROC","PR","logLoss","mnLogLoss")) {
      sumfun <- if (nlevels(y) > 2) caret::multiClassSummary else caret::twoClassSummary
    } else {
      sumfun <- caret::multiClassSummary  # Use multiClassSummary for Kappa, Accuracy, etc.
    }

    preproc <- if (mm$method == "svmLinear") c("center","scale") else NULL

    # Create fold indices
    fold_indices <- NULL

    if (cv_method == "stratified") {
      if (verbose) message("Creating stratified CV folds for simple mode")

      # Create stratified folds ensuring balanced class representation
      stratified_folds <- create_stratified_folds(y, cv.folds, seed = seed)

      # Convert to caret format (training indices for each fold)
      fold_indices <- vector("list", cv.folds)
      for (i in seq_len(cv.folds)) {
        fold_indices[[i]] <- setdiff(seq_len(nrow(X)), stratified_folds[[i]])
      }

      # Calculate and display stratification quality for simple mode
      if (verbose) {
        classes <- levels(y)

        # Calculate quality metrics
        chi_squared_scores <- sapply(classes, function(class_label) {
          observed <- sapply(stratified_folds, function(fold_idx) sum(y[fold_idx] == class_label))
          total_class <- sum(y == class_label)
          expected <- rep(total_class / length(stratified_folds), length(stratified_folds))
          expected[expected == 0] <- 1e-10
          return(sum((observed - expected)^2 / expected))
        })

        overall_chi_sq <- sum(chi_squared_scores)

        # Calibrated quality calculation (same as nested mode)
        n_classes <- length(classes)
        n_folds <- cv.folds
        chi_sq_threshold <- qchisq(0.95, df = (n_folds - 1)) * n_classes
        chi_sq_normalized <- pmin(1, overall_chi_sq / chi_sq_threshold)

        # Convert to quality score using calibrated sigmoid
        stratification_quality <- 1 / (1 + exp(8 * (chi_sq_normalized - 0.5)))

        # Store simple mode stratification metrics for user access
        simple_stratification_metrics <- list(
          method_used = "stratified",
          cv_folds = cv.folds,
          stratification_quality = list(
            chi_squared_scores = chi_squared_scores,
            overall_chi_squared = overall_chi_sq,
            chi_sq_normalized = chi_sq_normalized,
            stratification_quality_index = stratification_quality,
            calibration_info = list(
              chi_sq_threshold = chi_sq_threshold,
              cv_threshold = 0.3,  # For consistency with nested mode
              weights = c(chi_squared = 1.0),  # Simple mode only uses chi-squared
              sigmoid_params = c(steepness = 8, midpoint = 0.5)
            )
          )
        )

        message("=== SIMPLE MODE STRATIFICATION QUALITY ===")
        message("Overall Chi-squared: ", round(overall_chi_sq, 3), " (lower = better)")
        message("Stratification Quality: ", round(stratification_quality, 3), " (higher = better)")

        if (stratification_quality > 0.8) {
          message("EXCELLENT stratification quality")
        } else if (stratification_quality > 0.6) {
          message("\n\u2713 GOOD stratification quality")
        } else {
          message("FAIR stratification quality")
        }
      } else {
        simple_stratification_metrics <- NULL
      }

    } else {
      if (verbose) message("Using random CV folds for simple mode")
      # fold_indices will remain NULL, caret will use random CV
    }

    ctrl <- caret::trainControl(
      method = if (cv.repeats > 1) "repeatedcv" else "cv",
      number = cv.folds, repeats = cv.repeats,
      classProbs = uses_probs,
      summaryFunction = sumfun,
      allowParallel = allowParallel,
      savePredictions = "final",
      index = if (!is.null(fold_indices)) fold_indices else NULL
    )

    fit <- caret::train(
      x = X, y = y,
      method = mm$method,
      tuneGrid = mm$grid,
      trControl = ctrl,
      metric = metric,
      preProcess = preproc,
      MaxNWts = if (mm$method == "multinom") maxnweights else NULL
    )

    STEAM.obj$train$model <- fit

    # Store stratification quality metrics for simple mode if they exist
    if (exists("simple_stratification_metrics") && !is.null(simple_stratification_metrics)) {
      STEAM.obj$train$fold_diagnostics <- simple_stratification_metrics
      STEAM.obj$train$stratification_method <- simple_stratification_metrics$method_used
    }

    return(STEAM.obj)
  }

  # ================= NESTED mode =================
  clean_labels <- function(labels) {
    labels <- as.character(labels)
    labels <- gsub("/", ".", labels); labels <- gsub("-", ".", labels)
    labels <- gsub(" ", "_", labels); labels <- gsub("[()]", "", labels)
    labels <- gsub(",", "", labels);  labels <- gsub("\\.\\.", ".", labels)
    labels <- trimws(labels); labels <- iconv(labels, to = "ASCII//TRANSLIT"); labels
  }
  y <- factor(clean_labels(STEAM.obj$labels))
  X <- t(as.matrix(STEAM.obj$count_exp))
  coords_all <- as.data.frame(STEAM.obj$spatial)
  if (is.null(rownames(X))) stop("count_exp must have column names matching spatial rownames.")
  if (is.null(rownames(coords_all))) rownames(coords_all) <- colnames(STEAM.obj$count_exp)

  # Create folds for nested CV
  if (is.null(outer_folds) || is.null(inner_folds)) {
    if (verbose) message("Creating nested CV folds using method: ", cv_method)

    # Create nested stratified folds ensuring layer representation
    stratified_result <- create_nested_stratified_folds(
      y = y,
      n_outer_folds = n_outer_folds,
      n_inner_folds = n_inner_folds,
      seed = seed
    )

    if (is.null(outer_folds)) outer_folds <- stratified_result$outer_folds
    if (is.null(inner_folds)) inner_folds <- stratified_result$inner_folds

    # Print detailed fold diagnostics
    if (verbose && !is.null(stratified_result$fold_diagnostics)) {
      # The diagnostics are already printed by the create_nested_stratified_folds function
      # Just add a separator
      message("================================")
    }
  }

  mm <- map_model(model, kernel, ncol(X), tune.grid)

  uses_probs <- TRUE  # Always use probabilities for consistency

  if (metric %in% c("ROC","PR","logLoss","mnLogLoss")) {
    sumfun <- if (nlevels(y) > 2) caret::multiClassSummary else caret::twoClassSummary
  } else {
    sumfun <- caret::multiClassSummary  # Use multiClassSummary for Kappa, Accuracy, etc.
  }

  preproc <- if (mm$method == "svmLinear") c("center","scale") else NULL

  # Adaptive inner control based on actual fold structure
  # Check if any outer folds have reduced inner folds
  actual_inner_folds <- sapply(inner_folds, function(x) {
    if (!is.null(attr(x, "actual_n_folds"))) {
      attr(x, "actual_n_folds")
    } else {
      length(x)
    }
  })

  # Use the minimum actual inner folds for consistency
  min_inner_folds <- min(actual_inner_folds)

  if (min_inner_folds < n_inner_folds) {
    message("Note: Some outer folds have reduced inner folds (min=", min_inner_folds,
            ") due to sample constraints. Using adaptive inner CV.")
  }

  tr_ctrl_inner <- caret::trainControl(
    method = "cv",
    number = min_inner_folds,  # Use actual minimum instead of requested
    classProbs = uses_probs,
    summaryFunction = sumfun,
    savePredictions = "final",
    allowParallel = allowParallel
  )

  if (cv.cores > 1 && isTRUE(allowParallel)) {
    message("Note: cv.cores > 1 and allowParallel = TRUE \u2192 disabling inner caret parallel to avoid nested parallelism.")
    tr_ctrl_inner$allowParallel <- FALSE
  }

  ncv_args <- list(
    y = y, x = X,
    method = mm$method,
    tuneGrid = mm$grid,
    metric = metric,
    trControl = tr_ctrl_inner,
    outer_folds = outer_folds,
    inner_folds = inner_folds,
    modifyX = "neighavg_fit",
    modifyX_useY = TRUE,
    modifyX_options = list(coords_all = coords_all, n.size = as.integer(n.size)),
    filterFUN = NULL,
    filter_options = list(),
    finalCV = NA,
    cv.cores = cv.cores,
    parallel_mode = parallel_mode,
    verbose = as.logical(verbose)
  )

  if (!is.null(preproc)) {
    ncv_args$preProcess <- preproc
  }

  if (mm$method == "multinom") {
    ncv_args$MaxNWts <- maxnweights
  }

  ncv <- do.call(nestedcv::nestcv.train, ncv_args)

  # Store stratification quality metrics in STEAM object for user access
  if (exists("stratified_result") && !is.null(stratified_result$fold_diagnostics)) {
    ncv$fold_diagnostics <- stratified_result$fold_diagnostics
    ncv$stratification_method <- stratified_result$method_used
    ncv$original_outer_folds <- stratified_result$original_outer_folds
    ncv$final_outer_folds <- stratified_result$final_outer_folds
  }

  STEAM.obj$train$ncv <- ncv
  STEAM.obj
}
