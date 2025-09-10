# Statistical Comparison Framework for STEAM
# Comprehensive statistical testing, effect size calculation, and multiple testing correction

# Input Validation Helper Function
validate_steam_results <- function(results_list, metric = "Kappa", verbose = TRUE) {

  if (verbose) {
    cat("=== VALIDATING STEAM RESULTS INPUT ===\n")
  }

  validation_results <- list()
  all_valid <- TRUE

  # Check if results_list is a list
  if (!is.list(results_list)) {
    stop("results_list must be a list of STEAM training results")
  }

  # Check each result
  for (i in seq_along(results_list)) {
    method_name <- names(results_list)[i]
    if (is.null(method_name) || method_name == "") {
      method_name <- paste0("Method_", i)
    }

    result <- results_list[[i]]
    validation <- validate_single_steam_result(result, method_name, metric, verbose)
    validation_results[[method_name]] <- validation

    if (!validation$valid) {
      all_valid <- FALSE
    }
  }

  if (verbose) {
    if (all_valid) {
      cat("All results are properly formatted!\n")
    } else {
      cat("Some results have formatting issues. Please fix before proceeding.\n")
    }
    cat("=== VALIDATION COMPLETE ===\n\n")
  }

  return(list(
    all_valid = all_valid,
    individual_results = validation_results,
    summary = list(
      n_methods = length(results_list),
      valid_methods = sum(sapply(validation_results, function(x) x$valid)),
      metric_available = all(sapply(validation_results, function(x) x$metric_available))
    )
  ))
}

# Validate single STEAM result
validate_single_steam_result <- function(steam_result, method_name, metric, verbose = TRUE) {

  if (verbose) {
    cat("\nChecking method:", method_name, "\n")
  }

  validation <- list(
    method_name = method_name,
    valid = FALSE,
    structure_type = "unknown",
    metric_available = FALSE,
    n_folds = 0,
    available_metrics = character(0),
    issues = character(0)
  )

  # Check basic structure
  if (!"train" %in% names(steam_result)) {
    validation$issues <- c(validation$issues, "Missing 'train' component")
    if (verbose) cat("Missing 'train' component\n")
    return(validation)
  }

  train_obj <- steam_result$train

  # Check for nested CV structure
  if ("ncv" %in% names(train_obj)) {
    validation$structure_type <- "nested"
    ncv <- train_obj$ncv

    if (verbose) cat("Nested CV structure detected\n")

    # Check for outer_result
    if ("outer_result" %in% names(ncv)) {
      outer_result <- ncv$outer_result

      if (is.list(outer_result) && length(outer_result) > 0) {
        validation$n_folds <- length(outer_result)

        # Check first fold structure
        first_fold <- outer_result[[1]]
        if (is.list(first_fold)) {
          validation$available_metrics <- names(first_fold)
          validation$metric_available <- metric %in% validation$available_metrics

          if (verbose) {
            cat("Found", validation$n_folds, "outer folds\n")
            cat("Available metrics:", paste(validation$available_metrics, collapse = ", "), "\n")

            if (validation$metric_available) {
              cat("Requested metric '", metric, "' is available\n")
            } else {
              cat("Requested metric '", metric, "' not found\n")
            }
          }
        } else {
          validation$issues <- c(validation$issues, "outer_result elements are not lists")
          if (verbose) cat("outer_result elements should be lists with named metrics\n")
        }
      } else {
        validation$issues <- c(validation$issues, "outer_result is empty or not a list")
        if (verbose) cat("outer_result is empty or not a list\n")
      }

      # Check for summary as fallback
    } else if ("summary" %in% names(ncv)) {
      summary_obj <- ncv$summary

      if (is.list(summary_obj)) {
        validation$available_metrics <- names(summary_obj)
        validation$metric_available <- metric %in% validation$available_metrics
        validation$n_folds <- 5  # Default assumption

        if (verbose) {
          cat("No outer_result found, using summary\n")
          cat("Available metrics:", paste(validation$available_metrics, collapse = ", "), "\n")

          if (validation$metric_available) {
            cat("Requested metric '", metric, "' is available\n")
          } else {
            cat("Requested metric '", metric, "' not found in summary\n")
          }
        }
      } else {
        validation$issues <- c(validation$issues, "summary is not a list")
        if (verbose) cat("summary should be a list with named metrics\n")
      }
    } else {
      validation$issues <- c(validation$issues, "Neither outer_result nor summary found in ncv")
      if (verbose) cat("Neither outer_result nor summary found in ncv\n")
    }

    # Check for simple CV structure
  } else if ("model" %in% names(train_obj)) {
    validation$structure_type <- "simple"
    model_obj <- train_obj$model

    if (verbose) cat("Simple CV structure detected\n")

    # Check for resample data
    if ("resample" %in% names(model_obj)) {
      resample_data <- model_obj$resample

      if (is.data.frame(resample_data) && nrow(resample_data) > 0) {
        validation$n_folds <- nrow(resample_data)
        validation$available_metrics <- names(resample_data)
        validation$metric_available <- metric %in% validation$available_metrics

        if (verbose) {
          cat("Found resample data with", validation$n_folds, "folds\n")
          cat("Available metrics:", paste(validation$available_metrics, collapse = ", "), "\n")

          if (validation$metric_available) {
            cat("Requested metric '", metric, "' is available\n")
          } else {
            cat("Requested metric '", metric, "' not found\n")
          }
        }
      } else {
        validation$issues <- c(validation$issues, "resample is empty or not a data.frame")
        if (verbose) cat("resample should be a data.frame with fold results\n")
      }

      # Check for results as fallback
    } else if ("results" %in% names(model_obj)) {
      results_obj <- model_obj$results

      if (is.data.frame(results_obj) && nrow(results_obj) > 0) {
        validation$available_metrics <- names(results_obj)
        validation$metric_available <- metric %in% validation$available_metrics
        validation$n_folds <- 10  # Default assumption

        if (verbose) {
          cat("No resample found, using results summary\n")
          cat("Available metrics:", paste(validation$available_metrics, collapse = ", "), "\n")

          if (validation$metric_available) {
            cat("Requested metric '", metric, "' is available\n")
          } else {
            cat("Requested metric '", metric, "' not found in results\n")
          }
        }
      } else {
        validation$issues <- c(validation$issues, "results is empty or not a data.frame")
        if (verbose) cat("results should be a data.frame\n")
      }
    } else {
      validation$issues <- c(validation$issues, "Neither resample nor results found in model")
      if (verbose) cat("Neither resample nor results found in model\n")
    }

  } else {
    validation$issues <- c(validation$issues, "Neither ncv nor model found in train")
    if (verbose) cat("Neither ncv nor model found in train\n")
  }

  # Final validation
  validation$valid <- (length(validation$issues) == 0) &&
    validation$metric_available &&
    (validation$n_folds >= 3)

  if (verbose) {
    if (validation$valid) {
      cat("Method '", method_name, "' passed validation\n")
    } else {
      cat("Method '", method_name, "' failed validation\n")
      if (length(validation$issues) > 0) {
        cat("Issues found:\n")
        for (issue in validation$issues) {
          cat("  - ", issue, "\n")
        }
      }
    }
  }

  return(validation)
}

# Main Statistical Comparison Function
statistical_comparison <- function(
    results_list,
    method_names = NULL,
    metric = "Kappa",
    comparison_type = c("pairwise", "one_vs_all", "nested_models"),
    statistical_tests = c("paired_t", "wilcoxon", "mcnemar", "friedman"),
    effect_size_methods = c("cohens_d", "hedges_g", "cliff_delta", "eta_squared"),
    multiple_correction = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr", "none"),
    alpha = 0.05,
    bootstrap_iterations = 1000,
    confidence_level = 0.95,
    verbose = TRUE,
    return_raw = FALSE
) {

  # Input validation and setup
  comparison_type <- match.arg(comparison_type)
  multiple_correction <- match.arg(multiple_correction)

  if (is.null(method_names)) {
    method_names <- names(results_list)
    if (is.null(method_names)) {
      method_names <- paste0("Method_", seq_along(results_list))
    }
  }

  if (verbose) {
    message("=== STATISTICAL COMPARISON FRAMEWORK ===")
    message("Methods: ", paste(method_names, collapse = ", "))
    message("Metric: ", metric)
    message("Comparison type: ", comparison_type)
    message("Multiple correction: ", multiple_correction)
  }

  # Validate inputs
  validation <- validate_steam_results(results_list, metric, verbose)
  if (!validation$all_valid) {
    stop("Input validation failed. Please check your STEAM results format.")
  }

  # Extract performance metrics from results
  performance_data <- extract_performance_metrics(results_list, metric, method_names)

  # Perform statistical comparisons
  statistical_results <- perform_statistical_tests(
    performance_data,
    statistical_tests,
    comparison_type,
    verbose
  )

  # Calculate effect sizes
  effect_size_results <- calculate_effect_sizes(
    performance_data,
    effect_size_methods,
    comparison_type,
    confidence_level,
    bootstrap_iterations,
    verbose
  )

  # Apply multiple testing correction
  corrected_results <- apply_multiple_testing_correction(
    statistical_results,
    multiple_correction,
    alpha,
    verbose
  )

  # Generate comprehensive summary
  summary_results <- generate_comparison_summary(
    performance_data,
    statistical_results,
    effect_size_results,
    corrected_results,
    method_names,
    metric,
    alpha
  )

  # Create visualization data
  viz_data <- prepare_visualization_data(
    performance_data,
    statistical_results,
    effect_size_results,
    method_names
  )

  # Compile final results
  final_results <- list(
    summary = summary_results,
    performance_data = performance_data,
    statistical_tests = statistical_results,
    effect_sizes = effect_size_results,
    corrected_results = corrected_results,
    visualization_data = viz_data,
    metadata = list(
      method_names = method_names,
      metric = metric,
      comparison_type = comparison_type,
      multiple_correction = multiple_correction,
      alpha = alpha,
      n_methods = length(method_names),
      n_comparisons = nrow(corrected_results$pairwise_comparisons)
    )
  )

  # Add raw data if requested
  if (return_raw) {
    final_results$raw_results <- results_list
  }

  class(final_results) <- c("statistical_comparison", "list")

  if (verbose) {
    print_comparison_summary(final_results)
  }

  return(final_results)
}

# Performance Metric Extraction
extract_performance_metrics <- function(results_list, metric, method_names) {

  performance_matrix <- NULL
  fold_data <- list()

  for (i in seq_along(results_list)) {
    result <- results_list[[i]]
    method_name <- method_names[i]

    # Extract performance based on result structure
    if ("ncv" %in% names(result$train)) {
      # Nested CV results
      perf_values <- extract_nested_cv_performance(result, metric)
    } else if ("model" %in% names(result$train)) {
      # Simple CV results
      perf_values <- extract_simple_cv_performance(result, metric)
    } else {
      stop(paste("Unknown result structure for method:", method_name))
    }

    fold_data[[method_name]] <- perf_values

    # Build performance matrix
    if (is.null(performance_matrix)) {
      performance_matrix <- matrix(NA, nrow = length(perf_values), ncol = length(results_list))
      colnames(performance_matrix) <- method_names
    }

    performance_matrix[, i] <- perf_values
  }

  # Calculate summary statistics
  summary_stats <- data.frame(
    method = method_names,
    mean = apply(performance_matrix, 2, mean, na.rm = TRUE),
    median = apply(performance_matrix, 2, median, na.rm = TRUE),
    sd = apply(performance_matrix, 2, sd, na.rm = TRUE),
    min = apply(performance_matrix, 2, min, na.rm = TRUE),
    max = apply(performance_matrix, 2, max, na.rm = TRUE),
    q25 = apply(performance_matrix, 2, quantile, 0.25, na.rm = TRUE),
    q75 = apply(performance_matrix, 2, quantile, 0.75, na.rm = TRUE),
    n = apply(performance_matrix, 2, function(x) sum(!is.na(x))),
    stringsAsFactors = FALSE
  )

  return(list(
    matrix = performance_matrix,
    fold_data = fold_data,
    summary = summary_stats
  ))
}

# Extract performance from nested CV results
extract_nested_cv_performance <- function(result, metric) {
  ncv_result <- result$train$ncv

  if ("outer_result" %in% names(ncv_result)) {
    # Extract from outer fold results
    outer_results <- ncv_result$outer_result
    if (is.list(outer_results)) {
      perf_values <- sapply(outer_results, function(x) {
        if (metric %in% names(x)) {
          return(x[[metric]])
        } else if ("performance" %in% names(x) && metric %in% names(x$performance)) {
          return(x$performance[[metric]])
        } else {
          return(NA)
        }
      })
    } else {
      perf_values <- rep(NA, 5)  # Default to 5 folds
    }
  } else if ("summary" %in% names(ncv_result)) {
    # Extract from summary statistics
    if (metric %in% names(ncv_result$summary)) {
      perf_values <- rep(ncv_result$summary[[metric]], 5)
    } else {
      perf_values <- rep(NA, 5)
    }
  } else {
    perf_values <- rep(NA, 5)
  }

  return(perf_values)
}

# Extract performance from simple CV results
extract_simple_cv_performance <- function(result, metric) {
  model_result <- result$train$model

  if ("resample" %in% names(model_result)) {
    # Extract from resampling results
    resample_data <- model_result$resample
    if (metric %in% names(resample_data)) {
      perf_values <- resample_data[[metric]]
    } else {
      perf_values <- rep(NA, nrow(resample_data))
    }
  } else if ("results" %in% names(model_result)) {
    # Extract from results summary
    results_data <- model_result$results
    if (metric %in% names(results_data)) {
      perf_values <- rep(results_data[[metric]], 10)  # Default to 10 folds
    } else {
      perf_values <- rep(NA, 10)
    }
  } else {
    perf_values <- rep(NA, 10)
  }

  return(perf_values)
}

# Statistical Testing Framework
perform_statistical_tests <- function(performance_data, statistical_tests, comparison_type, verbose) {

  perf_matrix <- performance_data$matrix
  method_names <- colnames(perf_matrix)
  n_methods <- ncol(perf_matrix)

  test_results <- list()

  if (verbose) {
    message("\n=== STATISTICAL TESTS ===")
  }

  # Pairwise comparisons
  if (comparison_type %in% c("pairwise", "one_vs_all")) {

    pairwise_results <- data.frame()

    for (i in 1:(n_methods-1)) {
      for (j in (i+1):n_methods) {

        method1 <- method_names[i]
        method2 <- method_names[j]

        data1 <- perf_matrix[, i]
        data2 <- perf_matrix[, j]

        # Remove NA values
        complete_pairs <- complete.cases(data1, data2)
        data1_clean <- data1[complete_pairs]
        data2_clean <- data2[complete_pairs]

        if (length(data1_clean) < 3) {
          warning(paste("Insufficient data for comparison:", method1, "vs", method2))
          next
        }

        # Perform requested tests
        test_row <- data.frame(
          method1 = method1,
          method2 = method2,
          n_pairs = length(data1_clean),
          stringsAsFactors = FALSE
        )

        # Paired t-test
        if ("paired_t" %in% statistical_tests) {
          t_test <- tryCatch({
            t.test(data1_clean, data2_clean, paired = TRUE)
          }, error = function(e) NULL)

          if (!is.null(t_test)) {
            test_row$t_statistic <- t_test$statistic
            test_row$t_df <- t_test$parameter
            test_row$t_pvalue <- t_test$p.value
            test_row$t_ci_lower <- t_test$conf.int[1]
            test_row$t_ci_upper <- t_test$conf.int[2]
          }
        }

        # Wilcoxon signed-rank test
        if ("wilcoxon" %in% statistical_tests) {
          wilcox_test <- tryCatch({
            wilcox.test(data1_clean, data2_clean, paired = TRUE)
          }, error = function(e) NULL)

          if (!is.null(wilcox_test)) {
            test_row$wilcox_statistic <- wilcox_test$statistic
            test_row$wilcox_pvalue <- wilcox_test$p.value
          }
        }

        # McNemar test (for binary outcomes)
        if ("mcnemar" %in% statistical_tests) {
          # Convert to binary outcomes (above/below median)
          median_perf <- median(c(data1_clean, data2_clean), na.rm = TRUE)
          binary1 <- data1_clean > median_perf
          binary2 <- data2_clean > median_perf

          mcnemar_table <- table(binary1, binary2)
          if (all(dim(mcnemar_table) == c(2, 2))) {
            mcnemar_test <- tryCatch({
              mcnemar.test(mcnemar_table)
            }, error = function(e) NULL)

            if (!is.null(mcnemar_test)) {
              test_row$mcnemar_statistic <- mcnemar_test$statistic
              test_row$mcnemar_pvalue <- mcnemar_test$p.value
            }
          }
        }

        pairwise_results <- rbind(pairwise_results, test_row)
      }
    }

    test_results$pairwise <- pairwise_results
  }

  # Friedman test (if multiple methods and sufficient data)
  if ("friedman" %in% statistical_tests && n_methods >= 3) {

    # Prepare data for Friedman test
    complete_rows <- complete.cases(perf_matrix)
    if (sum(complete_rows) >= 3) {

      friedman_data <- perf_matrix[complete_rows, ]

      friedman_test <- tryCatch({
        friedman.test(as.matrix(friedman_data))
      }, error = function(e) NULL)

      if (!is.null(friedman_test)) {
        test_results$friedman <- list(
          statistic = friedman_test$statistic,
          df = friedman_test$parameter,
          p_value = friedman_test$p.value,
          method = friedman_test$method
        )

        if (verbose) {
          message("Friedman test: chi^2 = ", round(friedman_test$statistic, 3),
                  ", df = ", friedman_test$parameter,
                  ", p = ", round(friedman_test$p.value, 4))
        }
      }
    }
  }

  return(test_results)
}

# Effect Size Calculation Framework
calculate_effect_sizes <- function(performance_data, effect_size_methods, comparison_type,
                                   confidence_level, bootstrap_iterations, verbose) {

  perf_matrix <- performance_data$matrix
  method_names <- colnames(perf_matrix)
  n_methods <- ncol(perf_matrix)

  effect_size_results <- list()

  if (verbose) {
    message("\n=== EFFECT SIZE CALCULATIONS ===")
  }

  # Pairwise effect sizes
  pairwise_effects <- data.frame()

  for (i in 1:(n_methods-1)) {
    for (j in (i+1):n_methods) {

      method1 <- method_names[i]
      method2 <- method_names[j]

      data1 <- perf_matrix[, i]
      data2 <- perf_matrix[, j]

      # Remove NA values
      complete_pairs <- complete.cases(data1, data2)
      data1_clean <- data1[complete_pairs]
      data2_clean <- data2[complete_pairs]

      if (length(data1_clean) < 3) next

      effect_row <- data.frame(
        method1 = method1,
        method2 = method2,
        n_pairs = length(data1_clean),
        stringsAsFactors = FALSE
      )

      # Cohen's d
      if ("cohens_d" %in% effect_size_methods) {
        cohens_d <- calculate_cohens_d(data1_clean, data2_clean)
        effect_row$cohens_d <- cohens_d$estimate

        # Bootstrap confidence interval
        if (bootstrap_iterations > 0) {
          cohens_d_ci <- bootstrap_cohens_d(data1_clean, data2_clean,
                                            bootstrap_iterations, confidence_level)
          effect_row$cohens_d_ci_lower <- cohens_d_ci[1]
          effect_row$cohens_d_ci_upper <- cohens_d_ci[2]
        }
      }

      # Hedges' g (bias-corrected Cohen's d)
      if ("hedges_g" %in% effect_size_methods) {
        hedges_g <- calculate_hedges_g(data1_clean, data2_clean)
        effect_row$hedges_g <- hedges_g$estimate

        # Bootstrap confidence interval
        if (bootstrap_iterations > 0) {
          hedges_g_ci <- bootstrap_hedges_g(data1_clean, data2_clean,
                                            bootstrap_iterations, confidence_level)
          effect_row$hedges_g_ci_lower <- hedges_g_ci[1]
          effect_row$hedges_g_ci_upper <- hedges_g_ci[2]
        }
      }

      # Cliff's Delta (non-parametric effect size)
      if ("cliff_delta" %in% effect_size_methods) {
        cliff_delta <- calculate_cliff_delta(data1_clean, data2_clean)
        effect_row$cliff_delta <- cliff_delta$estimate
        effect_row$cliff_delta_magnitude <- cliff_delta$magnitude

        # Bootstrap confidence interval
        if (bootstrap_iterations > 0) {
          cliff_delta_ci <- bootstrap_cliff_delta(data1_clean, data2_clean,
                                                  bootstrap_iterations, confidence_level)
          effect_row$cliff_delta_ci_lower <- cliff_delta_ci[1]
          effect_row$cliff_delta_ci_upper <- cliff_delta_ci[2]
        }
      }

      pairwise_effects <- rbind(pairwise_effects, effect_row)
    }
  }

  effect_size_results$pairwise <- pairwise_effects

  # Eta-squared (for overall effect when multiple methods)
  if ("eta_squared" %in% effect_size_methods && n_methods >= 3) {
    eta_squared <- calculate_eta_squared(perf_matrix)
    effect_size_results$eta_squared <- eta_squared

    if (verbose) {
      message("Eta-squared (overall effect): ", round(eta_squared$estimate, 3))
    }
  }

  return(effect_size_results)
}

# Cohen's d calculation
calculate_cohens_d <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)

  mean_diff <- mean(x) - mean(y)
  pooled_sd <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))

  d <- mean_diff / pooled_sd

  return(list(estimate = d, interpretation = interpret_cohens_d(d)))
}

# Hedges' g calculation (bias-corrected Cohen's d)
calculate_hedges_g <- function(x, y) {
  cohens_d_result <- calculate_cohens_d(x, y)
  d <- cohens_d_result$estimate

  n1 <- length(x)
  n2 <- length(y)
  df <- n1 + n2 - 2

  # Bias correction factor
  correction_factor <- 1 - (3 / (4 * df - 1))
  g <- d * correction_factor

  return(list(estimate = g, interpretation = interpret_cohens_d(g)))
}

# Cliff's Delta calculation
calculate_cliff_delta <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)

  # Calculate all pairwise comparisons
  comparisons <- outer(x, y, function(a, b) sign(a - b))

  # Cliff's delta is the proportion of positive minus negative comparisons
  delta <- sum(comparisons) / (n1 * n2)

  magnitude <- interpret_cliff_delta(delta)

  return(list(estimate = delta, magnitude = magnitude))
}

# Eta-squared calculation
calculate_eta_squared <- function(data_matrix) {
  complete_data <- data_matrix[complete.cases(data_matrix), ]

  if (nrow(complete_data) < 3 || ncol(complete_data) < 2) {
    return(list(estimate = NA, interpretation = "Insufficient data"))
  }

  # Convert to long format for ANOVA
  long_data <- data.frame(
    value = as.vector(complete_data),
    method = rep(colnames(complete_data), each = nrow(complete_data)),
    subject = rep(seq_len(nrow(complete_data)), ncol(complete_data))
  )

  # Perform repeated measures ANOVA
  aov_result <- tryCatch({
    aov(value ~ method + Error(factor(subject)), data = long_data)
  }, error = function(e) NULL)

  if (is.null(aov_result)) {
    return(list(estimate = NA, interpretation = "ANOVA failed"))
  }

  # Extract sum of squares
  aov_summary <- summary(aov_result)
  if (length(aov_summary) >= 2) {
    ss_method <- aov_summary$`Error: Within`[[1]]["method", "Sum Sq"]
    ss_error <- aov_summary$`Error: Within`[[1]]["Residuals", "Sum Sq"]

    eta_squared <- ss_method / (ss_method + ss_error)

    return(list(
      estimate = eta_squared,
      interpretation = interpret_eta_squared(eta_squared)
    ))
  }

  return(list(estimate = NA, interpretation = "Could not calculate"))
}

# Bootstrap confidence intervals for Cohen's d
bootstrap_cohens_d <- function(x, y, n_bootstrap, confidence_level) {
  bootstrap_estimates <- replicate(n_bootstrap, {
    n1 <- length(x)
    n2 <- length(y)

    x_boot <- sample(x, n1, replace = TRUE)
    y_boot <- sample(y, n2, replace = TRUE)

    calculate_cohens_d(x_boot, y_boot)$estimate
  })

  alpha <- 1 - confidence_level
  ci <- quantile(bootstrap_estimates, c(alpha/2, 1 - alpha/2), na.rm = TRUE)

  return(ci)
}

# Bootstrap confidence intervals for Hedges' g
bootstrap_hedges_g <- function(x, y, n_bootstrap, confidence_level) {
  bootstrap_estimates <- replicate(n_bootstrap, {
    n1 <- length(x)
    n2 <- length(y)

    x_boot <- sample(x, n1, replace = TRUE)
    y_boot <- sample(y, n2, replace = TRUE)

    calculate_hedges_g(x_boot, y_boot)$estimate
  })

  alpha <- 1 - confidence_level
  ci <- quantile(bootstrap_estimates, c(alpha/2, 1 - alpha/2), na.rm = TRUE)

  return(ci)
}

# Bootstrap confidence intervals for Cliff's Delta
bootstrap_cliff_delta <- function(x, y, n_bootstrap, confidence_level) {
  bootstrap_estimates <- replicate(n_bootstrap, {
    n1 <- length(x)
    n2 <- length(y)

    x_boot <- sample(x, n1, replace = TRUE)
    y_boot <- sample(y, n2, replace = TRUE)

    calculate_cliff_delta(x_boot, y_boot)$estimate
  })

  alpha <- 1 - confidence_level
  ci <- quantile(bootstrap_estimates, c(alpha/2, 1 - alpha/2), na.rm = TRUE)

  return(ci)
}

# Multiple Testing Correction
apply_multiple_testing_correction <- function(statistical_results, correction_method, alpha, verbose) {

  corrected_results <- list()

  if (verbose) {
    message("\n=== MULTIPLE TESTING CORRECTION ===")
    message("Method: ", correction_method)
    message("Alpha level: ", alpha)
  }

  # Extract p-values for correction
  if ("pairwise" %in% names(statistical_results)) {
    pairwise <- statistical_results$pairwise

    # Collect all p-values
    p_values <- c()
    test_names <- c()

    if ("t_pvalue" %in% names(pairwise)) {
      p_values <- c(p_values, pairwise$t_pvalue)
      test_names <- c(test_names, paste(pairwise$method1, "vs", pairwise$method2, "(t-test)"))
    }

    if ("wilcox_pvalue" %in% names(pairwise)) {
      p_values <- c(p_values, pairwise$wilcox_pvalue)
      test_names <- c(test_names, paste(pairwise$method1, "vs", pairwise$method2, "(Wilcoxon)"))
    }

    if ("mcnemar_pvalue" %in% names(pairwise)) {
      p_values <- c(p_values, pairwise$mcnemar_pvalue)
      test_names <- c(test_names, paste(pairwise$method1, "vs", pairwise$method2, "(McNemar)"))
    }

    # Remove NA values
    valid_indices <- !is.na(p_values)
    p_values_clean <- p_values[valid_indices]
    test_names_clean <- test_names[valid_indices]

    if (length(p_values_clean) > 0) {
      # Apply correction
      if (correction_method != "none") {
        corrected_p <- p.adjust(p_values_clean, method = correction_method)
      } else {
        corrected_p <- p_values_clean
      }

      # Determine significance
      significant <- corrected_p < alpha

      correction_results <- data.frame(
        comparison = test_names_clean,
        raw_p = p_values_clean,
        corrected_p = corrected_p,
        significant = significant,
        stringsAsFactors = FALSE
      )

      corrected_results$pairwise_comparisons <- correction_results

      if (verbose) {
        n_significant <- sum(significant)
        message("Significant comparisons (after correction): ", n_significant, "/", length(corrected_p))
      }
    }
  }

  # Add Friedman test result if available
  if ("friedman" %in% names(statistical_results)) {
    friedman_p <- statistical_results$friedman$p_value
    friedman_significant <- friedman_p < alpha

    corrected_results$friedman <- list(
      p_value = friedman_p,
      significant = friedman_significant,
      alpha = alpha
    )
  }

  corrected_results$correction_method <- correction_method
  corrected_results$alpha <- alpha

  return(corrected_results)
}

# Generate Comprehensive Summary
generate_comparison_summary <- function(performance_data, statistical_results, effect_size_results,
                                        corrected_results, method_names, metric, alpha) {

  summary_list <- list()

  # Performance summary
  perf_summary <- performance_data$summary
  perf_summary$rank <- rank(-perf_summary$mean, ties.method = "min")

  summary_list$performance_ranking <- perf_summary[order(perf_summary$rank), ]

  # Statistical significance summary
  if ("pairwise_comparisons" %in% names(corrected_results)) {
    sig_comparisons <- corrected_results$pairwise_comparisons[
      corrected_results$pairwise_comparisons$significant, ]

    summary_list$significant_differences <- sig_comparisons
    summary_list$n_significant_comparisons <- nrow(sig_comparisons)
  }

  # Effect size summary
  if ("pairwise" %in% names(effect_size_results)) {
    effect_summary <- effect_size_results$pairwise

    # Categorize effect sizes
    if ("cohens_d" %in% names(effect_summary)) {
      effect_summary$cohens_d_category <- sapply(effect_summary$cohens_d,
                                                 interpret_cohens_d)
    }

    if ("cliff_delta" %in% names(effect_summary)) {
      effect_summary$cliff_delta_category <- effect_summary$cliff_delta_magnitude
    }

    summary_list$effect_size_summary <- effect_summary
  }

  # Overall assessment
  best_method <- method_names[which.min(perf_summary$rank)]

  summary_list$overall_assessment <- list(
    best_performing_method = best_method,
    best_mean_performance = perf_summary$mean[perf_summary$method == best_method],
    metric_used = metric,
    total_methods_compared = length(method_names),
    alpha_level = alpha
  )

  return(summary_list)
}

# Prepare Visualization Data
prepare_visualization_data <- function(performance_data, statistical_results, effect_size_results, method_names) {

  viz_data <- list()

  # Boxplot data
  viz_data$boxplot_data <- data.frame(
    value = as.vector(performance_data$matrix),
    method = rep(method_names, each = nrow(performance_data$matrix))
  )

  # Pairwise comparison matrix
  if ("pairwise" %in% names(statistical_results)) {
    pairwise <- statistical_results$pairwise
    n_methods <- length(method_names)

    # P-value matrix
    p_matrix <- matrix(1, n_methods, n_methods)
    rownames(p_matrix) <- method_names
    colnames(p_matrix) <- method_names

    # Effect size matrix
    effect_matrix <- matrix(0, n_methods, n_methods)
    rownames(effect_matrix) <- method_names
    colnames(effect_matrix) <- method_names

    for (i in seq_len(nrow(pairwise))) {
      m1 <- pairwise$method1[i]
      m2 <- pairwise$method2[i]

      # Use t-test p-value if available, otherwise Wilcoxon
      p_val <- if ("t_pvalue" %in% names(pairwise)) {
        pairwise$t_pvalue[i]
      } else if ("wilcox_pvalue" %in% names(pairwise)) {
        pairwise$wilcox_pvalue[i]
      } else {
        1
      }

      p_matrix[m1, m2] <- p_val
      p_matrix[m2, m1] <- p_val

      # Effect size (Cohen's d if available)
      if ("pairwise" %in% names(effect_size_results) &&
          "cohens_d" %in% names(effect_size_results$pairwise)) {
        effect_val <- effect_size_results$pairwise$cohens_d[i]
        effect_matrix[m1, m2] <- effect_val
        effect_matrix[m2, m1] <- -effect_val
      }
    }

    viz_data$p_value_matrix <- p_matrix
    viz_data$effect_size_matrix <- effect_matrix
  }

  return(viz_data)
}

# Interpretation Functions
interpret_cohens_d <- function(d) {
  abs_d <- abs(d)
  if (abs_d < 0.2) return("negligible")
  if (abs_d < 0.5) return("small")
  if (abs_d < 0.8) return("medium")
  return("large")
}

interpret_cliff_delta <- function(delta) {
  abs_delta <- abs(delta)
  if (abs_delta < 0.147) return("negligible")
  if (abs_delta < 0.33) return("small")
  if (abs_delta < 0.474) return("medium")
  return("large")
}

interpret_eta_squared <- function(eta_sq) {
  if (eta_sq < 0.01) return("small")
  if (eta_sq < 0.06) return("medium")
  return("large")
}

# Print Summary Function
print_comparison_summary <- function(results) {

  cat("\n")
  cat(repeat_char("=", 60), "\n")
  cat("STATISTICAL COMPARISON SUMMARY\n")
  cat(repeat_char("=", 60), "\n")

  # Performance ranking
  cat("\nPERFORMANCE RANKING:\n")
  ranking <- results$summary$performance_ranking
  for (i in seq_len(nrow(ranking))) {
    cat(sprintf("%d. %s: %.4f +/- %.4f\n",
                ranking$rank[i],
                ranking$method[i],
                ranking$mean[i],
                ranking$sd[i]))
  }

  # Significant differences
  if ("significant_differences" %in% names(results$summary)) {
    cat("\nSIGNIFICANT DIFFERENCES:\n")
    sig_diffs <- results$summary$significant_differences
    if (nrow(sig_diffs) > 0) {
      for (i in seq_len(nrow(sig_diffs))) {
        cat(sprintf("%s (p = %.4f)\n",
                    sig_diffs$comparison[i],
                    sig_diffs$corrected_p[i]))
      }
    } else {
      cat("No significant differences found.\n")
    }
  }

  # Overall assessment
  assessment <- results$summary$overall_assessment
  cat("\nOVERALL ASSESSMENT:\n")
  cat(sprintf("Best performing method: %s\n", assessment$best_performing_method))
  cat(sprintf("Best mean performance: %.4f\n", assessment$best_mean_performance))
  cat(sprintf("Total methods compared: %d\n", assessment$total_methods_compared))
  cat(sprintf("Alpha level: %.3f\n", assessment$alpha_level))

  cat("\n", repeat_char("=", 60), "\n")
}

# Visualization Helper Functions
plot_comparison_results <- function(results, plot_type = c("boxplot", "heatmap", "ranking")) {
  plot_type <- match.arg(plot_type)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Please install it.")
  }

  if (plot_type == "boxplot") {
    viz_data <- results$visualization_data$boxplot_data

    p <- ggplot(viz_data, aes(x = .data$method, y = .data$value, fill = .data$method)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Method Comparison",
           x = "Method",
           y = paste("Performance (", results$metadata$metric, ")"),
           fill = "Method")

    return(p)
  }

  # Additional plot types can be implemented here
}

# Helper operator for string repetition
`%R%` <- function(x, n) {
  paste(rep(x, n), collapse = "")
}

# Move %R% definition to top of file to ensure visibility
repeat_char <- function(char, n) {
  paste(rep(char, n), collapse = "")
}
