## DLPFC RESULTS

# Load required packages for statistical analysis (suppress startup messages)
if (!require("effectsize", quietly = TRUE)) install.packages("effectsize")
if (!require("multcomp", quietly = TRUE)) install.packages("multcomp") 
if (!require("emmeans", quietly = TRUE)) install.packages("emmeans")
if (!require("broom", quietly = TRUE)) install.packages("broom")
if (!require("PMCMRplus", quietly = TRUE)) install.packages("PMCMRplus")

suppressPackageStartupMessages({
  library(effectsize)
  library(multcomp)
  library(emmeans)
  library(broom)
  library(PMCMRplus)
})

DLPFCs <- readRDS("~/Desktop/STEAM/DLPFC/DLPFC12.RDS")

set.seed(0)
## ---------- 1) helpers ----------

# Master function for comprehensive metric extraction
extract_comprehensive_metrics <- function(steam_result, sample_id) {
  message("Extracting comprehensive metrics for sample ", sample_id)
  
  # Initialize results list
  results <- list(
    per_fold_metrics = NULL,
    pooled_metrics = NULL,
    raw_metrics = NULL
  )
  
  # Try different paths to find the nested CV results
  outer_result <- NULL
  if (!is.null(steam_result$nested$ncv$outer_result)) {
    outer_result <- steam_result$nested$ncv$outer_result
  } else if (!is.null(steam_result$train$ncv$outer_result)) {
    outer_result <- steam_result$train$ncv$outer_result
  } else if (!is.null(steam_result$nested)) {
    outer_result <- steam_result$nested
  }
  
  if (!is.null(outer_result)) {
    # Extract per-fold metrics
    results$per_fold_metrics <- compute_outer_fold_metrics(outer_result)
    
    # Extract pooled metrics
    results$pooled_metrics <- pooled_metrics_row(steam_result)
    
    # Store raw metrics for debugging
    results$raw_metrics <- outer_result
    
    message("Successfully extracted metrics - Folds: ", nrow(results$per_fold_metrics), 
            ", Metrics: ", ncol(results$per_fold_metrics) - 1)
  } else {
    message("Warning: Could not find nested CV results for sample ", sample_id)
    # Create empty results with expected structure
    results$pooled_metrics <- data.frame(
      Sample = sample_id,
      Accuracy = NA,
      Balanced_Accuracy = NA,
      Kappa = NA,
      stringsAsFactors = FALSE
    )
  }
  
  return(results)
}

compute_outer_fold_metrics <- function(outer_res) {
  stopifnot(!is.null(outer_res), length(outer_res) > 0)
  
  do.call(rbind, lapply(seq_along(outer_res), function(i) {
    of <- outer_res[[i]]
    preds <- of$preds
    
    # coerce to character, drop NAs, then align to a common level set
    pr <- as.character(preds$predy)
    tr <- as.character(preds$testy)
    keep <- !(is.na(pr) | is.na(tr))
    pr <- pr[keep]; tr <- tr[keep]
    
    # union of labels → guarantees a square confusion matrix
    lev <- sort(unique(c(pr, tr)))
    prf <- factor(pr, levels = lev)
    trf <- factor(tr, levels = lev)
    
    # accuracy/kappa via caret on factors (also safe now)
    cm_full <- caret::confusionMatrix(prf, trf)
    acc   <- unname(cm_full$overall["Accuracy"])
    kappa <- unname(cm_full$overall["Kappa"])
    
    # Calculate F1 Score (macro average)
    # Extract per-class precision, recall, and F1 from confusion matrix
    cm_byclass <- cm_full$byClass
    if (is.matrix(cm_byclass)) {
      # Multi-class case
      f1_scores <- cm_byclass[, "F1"]
      macro_f1 <- mean(f1_scores, na.rm = TRUE)
    } else if (is.vector(cm_byclass) && "F1" %in% names(cm_byclass)) {
      # Binary case
      macro_f1 <- cm_byclass["F1"]
    } else {
      # Fallback calculation if F1 not available
      precision <- cm_byclass[, "Precision"]
      recall <- cm_byclass[, "Recall"]
      f1_scores <- 2 * (precision * recall) / (precision + recall)
      macro_f1 <- mean(f1_scores, na.rm = TRUE)
    }
    
    # external indices
    ari <- mclust::adjustedRandIndex(tr, pr)
    nmi <- aricode::NMI(tr, pr)
    ami <- aricode::AMI(tr, pr)
    
    # Simplified core metrics only - focus on statistical analysis
    result <- data.frame(
      Outer_Fold    = i,
      Accuracy      = as.numeric(acc),
      Kappa         = as.numeric(kappa),
      Macro_F1      = as.numeric(macro_f1),
      ARI           = as.numeric(ari),
      NMI           = as.numeric(nmi),
      AMI           = as.numeric(ami),
      stringsAsFactors = FALSE
    )
    
    return(result)
  }))
}


# Enhanced statistical summary using established packages
summarize_outer_folds <- function(outer_df, confidence_level = 0.95) {
  stopifnot(is.data.frame(outer_df), "Outer_Fold" %in% names(outer_df))
  num_cols <- setdiff(names(outer_df), "Outer_Fold")
  
  # Calculate statistics for each metric using broom for tidy output
  stats_list <- lapply(num_cols, function(col) {
    x <- outer_df[[col]]
    x <- x[!is.na(x)]  # Remove NAs
    n <- length(x)
    
    if (n == 0) {
      return(data.frame(
        Metric = col, Mean = NA, SD = NA, SE = NA, 
        CI_Lower = NA, CI_Upper = NA, Cohen_d = NA, Cohen_d_CI_Lower = NA, Cohen_d_CI_Upper = NA,
        stringsAsFactors = FALSE
      ))
    }
    
    # Validate data is numeric
    if (!is.numeric(x)) {
      warning("Non-numeric data encountered in column: ", col)
      return(data.frame(
        Metric = col, Mean = NA, SD = NA, SE = NA, 
        CI_Lower = NA, CI_Upper = NA, Cohen_d = NA, Cohen_d_CI_Lower = NA, Cohen_d_CI_Upper = NA,
        n_folds = n,
        stringsAsFactors = FALSE
      ))
    }
    
    # Check for constant values (all identical)
    if (sd(x) == 0 || var(x) < .Machine$double.eps) {
      warning("Constant values detected in column: ", col, ". Cannot compute statistics.")
      return(data.frame(
        Metric = col, 
        Mean = mean(x), 
        SD = 0, 
        SE = 0, 
        CI_Lower = mean(x), 
        CI_Upper = mean(x), 
        Cohen_d = NA,  # Cannot compute effect size for constant data
        Cohen_d_CI_Lower = NA, 
        Cohen_d_CI_Upper = NA,
        n_folds = n,
        stringsAsFactors = FALSE
      ))
    }
    
    # Basic statistics using broom::tidy for t.test
    t_test_result <- tryCatch({
      t.test(x, conf.level = confidence_level)
    }, error = function(e) {
      warning("t.test failed for column: ", col, ". Error: ", e$message)
      return(NULL)
    })
    
    if (is.null(t_test_result)) {
      return(data.frame(
        Metric = col, Mean = mean(x), SD = sd(x), SE = sd(x)/sqrt(n), 
        CI_Lower = NA, CI_Upper = NA, Cohen_d = NA, Cohen_d_CI_Lower = NA, Cohen_d_CI_Upper = NA,
        n_folds = n,
        stringsAsFactors = FALSE
      ))
    }
    
    tidy_result <- broom::tidy(t_test_result)
    
    # Effect size using effectsize package - Cohen's d against null (0)
    cohens_d_result <- tryCatch({
      # Only compute if we have sufficient variation and reasonable values
      if (sd(x) > 0 && is.finite(mean(x)) && is.finite(sd(x))) {
        effectsize::cohens_d(x, mu = 0)
      } else {
        list(Cohens_d = NA, CI_low = NA, CI_high = NA)
      }
    }, error = function(e) {
      # Fallback calculation if effectsize fails
      if (sd(x) > 0 && is.finite(mean(x)) && is.finite(sd(x))) {
        list(Cohens_d = mean(x) / sd(x), CI_low = NA, CI_high = NA)
      } else {
        list(Cohens_d = NA, CI_low = NA, CI_high = NA)
      }
    })
    
    # Extract Cohen's d and its confidence interval
    if (is.list(cohens_d_result) && "Cohens_d" %in% names(cohens_d_result)) {
      cohen_d_val <- cohens_d_result$Cohens_d
      cohen_d_ci_low <- cohens_d_result$CI_low
      cohen_d_ci_high <- cohens_d_result$CI_high
    } else if (is.numeric(cohens_d_result)) {
      cohen_d_val <- cohens_d_result[1]
      cohen_d_ci_low <- NA
      cohen_d_ci_high <- NA
    } else {
      # Safe fallback calculation
      if (sd(x) > 0 && is.finite(mean(x)) && is.finite(sd(x))) {
        cohen_d_val <- mean(x) / sd(x)
      } else {
        cohen_d_val <- NA
      }
      cohen_d_ci_low <- NA
      cohen_d_ci_high <- NA
    }
    
    # Calculate SE properly - avoid division by zero
    se_val <- if (is.finite(tidy_result$statistic) && tidy_result$statistic != 0) {
      abs(tidy_result$estimate / tidy_result$statistic)
    } else {
      sd(x) / sqrt(n)  # Standard formula
    }
    
    data.frame(
      Metric = col,
      Mean = tidy_result$estimate,
      SD = sd(x),
      SE = se_val,
      CI_Lower = tidy_result$conf.low,
      CI_Upper = tidy_result$conf.high,
      Cohen_d = cohen_d_val,
      Cohen_d_CI_Lower = cohen_d_ci_low,
      Cohen_d_CI_Upper = cohen_d_ci_high,
      n_folds = n,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, stats_list)
}


pooled_metrics_row <- function(out_nested) {
  # Try multiple paths to find metrics
  m <- NULL
  if (!is.null(out_nested[["nested"]][["metrics"]])) {
    m <- out_nested[["nested"]][["metrics"]]
  } else if (!is.null(out_nested[["train"]][["ncv"]][["summary"]])) {
    m <- out_nested[["train"]][["ncv"]][["summary"]]
  } else if (!is.null(out_nested[["train"]][["ncv"]][["final_result"]])) {
    m <- out_nested[["train"]][["ncv"]][["final_result"]]
  }
  
  if (is.null(m)) {
    warning("pooled_metrics_row: No metrics found in nested results. Available names: ", 
            paste(names(out_nested), collapse = ", "))
    # Return empty row with all expected columns
    return(data.frame(
      Accuracy = NA_real_, Kappa = NA_real_, ARI = NA_real_, NMI = NA_real_, AMI = NA_real_,
      Macro_Precision = NA_real_, Macro_Recall = NA_real_, Macro_F1 = NA_real_,
      Test_Accuracy = NA_real_, Train_Accuracy = NA_real_,
      Best_Tune_mtry = NA_real_,
      Processing_Time = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  # Initialize result data frame with proper structure
  result <- data.frame(
    Accuracy = NA_real_,
    Kappa = NA_real_,
    ARI = NA_real_,
    NMI = NA_real_,
    AMI = NA_real_,
    Macro_Precision = NA_real_,
    Macro_Recall = NA_real_,
    Macro_F1 = NA_real_,
    Test_Accuracy = NA_real_,
    Train_Accuracy = NA_real_,
    Best_Tune_mtry = NA_real_,
    Processing_Time = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # Handle different summary structure types
  if (is.list(m) && "metrics" %in% names(m)) {
    # Extract from nested metrics structure
    metrics <- m$metrics
    if (is.numeric(metrics) && length(metrics) >= 3) {
      # Assuming order: Multiclass AUC, Accuracy, Balanced accuracy
      result$Accuracy <- as.numeric(metrics[2])  # Accuracy is typically 2nd
      # Extract other metrics if available
    }
  }
  
  # Core metrics with flexible names - handle both nested and direct access
  if (is.list(m)) {
    # Try direct access first
    result$Accuracy <- as.numeric(m$Accuracy %||% m$accuracy %||% 
                                 (if("metrics" %in% names(m) && is.numeric(m$metrics) && length(m$metrics) >= 2) m$metrics[2] else NA_real_))
    result$Kappa <- as.numeric(m$Kappa %||% m$kappa %||% NA_real_)
    result$ARI <- as.numeric(m$ARI %||% m$ari %||% NA_real_)
    result$NMI <- as.numeric(m$NMI %||% m$nmi %||% NA_real_)
    result$AMI <- as.numeric(m$AMI %||% m$ami %||% NA_real_)
    
    # Handle multiclass AUC if it's the first metric
    if ("metrics" %in% names(m) && is.numeric(m$metrics) && length(m$metrics) >= 3) {
      if (is.na(result$Accuracy)) result$Accuracy <- as.numeric(m$metrics[2])
      # Could add AUC if needed: result$Multiclass_AUC <- as.numeric(m$metrics[1])
    }
    
    # PRF metrics
    p <- m$precision %||% m$Precision
    r <- m$recall %||% m$Recall  
    f <- m$f1_score %||% m$F1 %||% m$f1
    
    macro_p <- if (!is.null(p)) mean(p, na.rm = TRUE) else as.numeric(m$Macro_Precision %||% NA_real_)
    macro_r <- if (!is.null(r)) mean(r, na.rm = TRUE) else as.numeric(m$Macro_Recall %||% NA_real_)
    macro_f <- if (!is.null(f)) {
      mean(f, na.rm = TRUE)
    } else if (is.finite(macro_p) && is.finite(macro_r) && (macro_p + macro_r) > 0) {
      2 * (macro_p * macro_r) / (macro_p + macro_r)
    } else as.numeric(m$Macro_F1 %||% NA_real_)
    
    result$Macro_Precision <- macro_p
    result$Macro_Recall <- macro_r  
    result$Macro_F1 <- macro_f
    
    # Additional metrics
    result$Test_Accuracy <- as.numeric(m$Test_Accuracy %||% m$test_accuracy %||% result$Accuracy)
    result$Train_Accuracy <- as.numeric(m$Train_Accuracy %||% m$train_accuracy %||% NA_real_)
    
    # Model tuning information
    if (!is.null(m$best_tune) || !is.null(m$bestTune)) {
      bt <- m$best_tune %||% m$bestTune
      if (is.data.frame(bt)) {
        if ("mtry" %in% names(bt)) result$Best_Tune_mtry <- as.numeric(bt$mtry[1])
      } else if (is.list(bt)) {
        result$Best_Tune_mtry <- as.numeric(bt$mtry %||% NA_real_)
      }
    }
    
    # Processing time if available
    result$Processing_Time <- as.numeric(m$processing_time %||% m$elapsed_time %||% NA_real_)
  }
  
  return(result)
}

# Helper function for null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Pairwise comparisons using established packages
pairwise_sample_comparisons <- function(all_per_fold, metric = "Accuracy", 
                                       correction_method = "bonferroni") {
  
  # Prepare data for analysis
  long_data <- data.frame()
  for (sample_id in names(all_per_fold)) {
    if (metric %in% names(all_per_fold[[sample_id]])) {
      sample_df <- data.frame(
        Sample = sample_id,
        Metric_Value = all_per_fold[[sample_id]][[metric]],
        Fold = all_per_fold[[sample_id]]$Outer_Fold,
        stringsAsFactors = FALSE
      )
      long_data <- rbind(long_data, sample_df)
    }
  }
  
  if (nrow(long_data) == 0 || length(unique(long_data$Sample)) < 2) {
    warning("Insufficient data for pairwise comparisons")
    return(NULL)
  }
  
  # Method 1: Using emmeans for pairwise comparisons
  tryCatch({
    # Fit a simple model
    model <- lm(Metric_Value ~ Sample, data = long_data)
    
    # Get estimated marginal means
    emm <- emmeans::emmeans(model, ~ Sample)
    
    # Pairwise comparisons with correction
    pairwise_result <- pairs(emm, adjust = correction_method)
    pairwise_tidy <- broom::tidy(pairwise_result)
    
    # Add effect sizes using effectsize
    effect_sizes <- data.frame()
    samples <- unique(long_data$Sample)
    
    for (i in 1:(length(samples)-1)) {
      for (j in (i+1):length(samples)) {
        s1 <- samples[i]
        s2 <- samples[j]
        
        x1 <- long_data$Metric_Value[long_data$Sample == s1]
        x2 <- long_data$Metric_Value[long_data$Sample == s2]
        
        # Cohen's d for independent samples using effectsize package
        cohens_d_result <- tryCatch({
          effectsize::cohens_d(x1, x2)
        }, error = function(e) {
          # Fallback calculation
          pooled_sd <- sqrt(((length(x1)-1)*var(x1) + (length(x2)-1)*var(x2))/(length(x1)+length(x2)-2))
          list(Cohens_d = (mean(x1) - mean(x2)) / pooled_sd, CI_low = NA, CI_high = NA)
        })
        
        effect_size_row <- data.frame(
          contrast = paste(s1, "-", s2),
          Effect_Size = if(is.list(cohens_d_result)) cohens_d_result$Cohens_d else cohens_d_result,
          Effect_Size_CI_Low = if(is.list(cohens_d_result)) cohens_d_result$CI_low else NA,
          Effect_Size_CI_High = if(is.list(cohens_d_result)) cohens_d_result$CI_high else NA,
          stringsAsFactors = FALSE
        )
        effect_sizes <- rbind(effect_sizes, effect_size_row)
      }
    }
    
    # Merge results
    final_results <- merge(pairwise_tidy, effect_sizes, by = "contrast", all.x = TRUE)
    
    # Add interpretation using effectsize package
    final_results$Interpretation <- effectsize::interpret_cohens_d(abs(final_results$Effect_Size))
    final_results$Metric <- metric
    final_results$Correction_Method <- correction_method
    final_results$Significant <- final_results$adj.p.value < 0.05
    
    return(final_results)
    
  }, error = function(e) {
    warning("emmeans approach failed, using PMCMRplus fallback: ", e$message)
    
    # Method 2: Fallback using PMCMRplus
    tryCatch({
      # Kruskal-Wallis with pairwise comparisons using PMCMRplus
      kw_result <- PMCMRplus::kwAllPairsConoverTest(
        x = long_data$Metric_Value, 
        g = factor(long_data$Sample),
        p.adjust.method = correction_method
      )
      
      # Convert to data frame
      pairwise_p <- kw_result$p.value
      comparison_results <- data.frame()
      
      sample_names <- rownames(pairwise_p)
      for (i in 1:nrow(pairwise_p)) {
        for (j in 1:ncol(pairwise_p)) {
          if (!is.na(pairwise_p[i, j])) {
            comparison_results <- rbind(comparison_results, data.frame(
              contrast = paste(sample_names[i], "-", colnames(pairwise_p)[j]),
              adj.p.value = pairwise_p[i, j],
              Significant = pairwise_p[i, j] < 0.05,
              Metric = metric,
              Correction_Method = correction_method,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
      
      return(comparison_results)
      
    }, error = function(e2) {
      warning("Both emmeans and PMCMRplus failed: ", e2$message)
      return(NULL)
    })
  })
}

# Cross-validation fold comparison (paired t-tests across folds)
cv_fold_comparisons <- function(all_per_fold, metrics = c("Accuracy", "Kappa"), 
                               correction_method = "bonferroni") {
  results_list <- list()
  
  for (metric in metrics) {
    comparison_result <- pairwise_sample_comparisons(
      all_per_fold, metric, correction_method
    )
    
    if (!is.null(comparison_result)) {
      results_list[[metric]] <- comparison_result
    }
  }
  
  if (length(results_list) > 0) {
    combined_results <- do.call(rbind, results_list)
    
    # Overall multiple testing correction across all metrics
    # Check if we have p.value (from emmeans/broom) or adj.p.value (from PMCMRplus)
    p_col <- if("p.value" %in% names(combined_results)) "p.value" else "adj.p.value"
    combined_results$p_global_adjusted <- p.adjust(
      combined_results[[p_col]], 
      method = correction_method
    )
    combined_results$Significant_Global <- combined_results$p_global_adjusted < 0.05
    
    return(combined_results)
  } else {
    return(NULL)
  }
}

# Performance testing using established statistical packages
statistical_performance_summary <- function(all_per_fold, all_summary) {
  # Combine all samples for overall analysis
  all_metrics <- do.call(rbind, all_per_fold)
  
  # Test each metric against null hypothesis
  metric_cols <- c("Accuracy", "Kappa", "ARI", "NMI", "AMI")
  
  # Define null hypotheses for each metric
  null_values <- list(
    Accuracy = 1/length(unique(all_metrics$Sample)),  # Random accuracy for multiclass
    Kappa = 0,      # No agreement
    ARI = 0,        # Random clustering
    NMI = 0,        # No mutual information
    AMI = 0         # No adjusted mutual information
  )
  
  significance_tests <- data.frame()
  
  for (metric in metric_cols) {
    if (metric %in% names(all_metrics)) {
      values <- all_metrics[[metric]]
      values <- values[!is.na(values)]
      
      if (length(values) > 1) {
        null_val <- null_values[[metric]] %||% 0
        
        # One-sample t-test using broom for tidy output
        test_result <- t.test(values, mu = null_val)
        tidy_test <- broom::tidy(test_result)
        
        # Effect size using effectsize package
        cohens_d_result <- tryCatch({
          effectsize::cohens_d(values, mu = null_val)
        }, error = function(e) {
          # Fallback calculation
          (mean(values) - null_val) / sd(values)
        })
        
        # Extract effect size and interpretation
        if (is.list(cohens_d_result)) {
          effect_size <- cohens_d_result$Cohens_d
          effect_interpretation <- effectsize::interpret_cohens_d(abs(effect_size))
          effect_ci_low <- cohens_d_result$CI_low
          effect_ci_high <- cohens_d_result$CI_high
        } else {
          effect_size <- cohens_d_result
          effect_interpretation <- effectsize::interpret_cohens_d(abs(effect_size))
          effect_ci_low <- NA
          effect_ci_high <- NA
        }
        
        # Compile results
        test_row <- data.frame(
          Metric = metric,
          Mean_Performance = tidy_test$estimate,
          Null_Hypothesis = null_val,
          t_statistic = tidy_test$statistic,
          p_value = tidy_test$p.value,
          CI_Lower = tidy_test$conf.low,
          CI_Upper = tidy_test$conf.high,
          Effect_Size = effect_size,
          Effect_Size_CI_Low = effect_ci_low,
          Effect_Size_CI_High = effect_ci_high,
          Interpretation = as.character(effect_interpretation),
          stringsAsFactors = FALSE
        )
        
        significance_tests <- rbind(significance_tests, test_row)
      }
    }
  }
  
  # Multiple testing correction using p.adjust
  if (nrow(significance_tests) > 0) {
    significance_tests$p_adjusted <- p.adjust(
      significance_tests$p_value, 
      method = "bonferroni"
    )
    significance_tests$Significant <- significance_tests$p_adjusted < 0.05
  }
  
  return(significance_tests)
}

# Advanced statistical analysis using additional packages
advanced_cv_analysis <- function(all_per_fold, metrics = c("Accuracy", "Kappa")) {
  
  # Prepare data for advanced analysis
  long_data <- data.frame()
  for (sample_id in names(all_per_fold)) {
    for (metric in metrics) {
      if (metric %in% names(all_per_fold[[sample_id]])) {
        sample_df <- data.frame(
          Sample = sample_id,
          Metric = metric,
          Value = all_per_fold[[sample_id]][[metric]],
          Fold = all_per_fold[[sample_id]]$Outer_Fold,
          stringsAsFactors = FALSE
        )
        long_data <- rbind(long_data, sample_df)
      }
    }
  }
  
  results <- list()
  
  # 1. Mixed-effects analysis (accounting for fold dependencies)
  if (requireNamespace("lme4", quietly = TRUE)) {
    tryCatch({
      library(lme4)
      mixed_model <- lme4::lmer(Value ~ Sample + Metric + (1|Fold), data = long_data)
      results$mixed_effects <- broom.mixed::tidy(mixed_model)
      results$mixed_effects_anova <- car::Anova(mixed_model)
    }, error = function(e) {
      message("Mixed-effects analysis failed: ", e$message)
    })
  }
  
  # 2. Bayesian analysis using rstanarm or brms
  if (requireNamespace("rstanarm", quietly = TRUE)) {
    tryCatch({
      library(rstanarm)
      bayesian_model <- rstanarm::stan_lmer(
        Value ~ Sample + Metric + (1|Fold), 
        data = long_data,
        cores = 2, 
        iter = 1000  # Reduced for faster computation
      )
      results$bayesian_summary <- summary(bayesian_model)
      results$bayesian_intervals <- posterior_interval(bayesian_model)
    }, error = function(e) {
      message("Bayesian analysis failed: ", e$message)
    })
  }
  
  # 3. Non-parametric alternatives using coin package
  if (requireNamespace("coin", quietly = TRUE)) {
    tryCatch({
      library(coin)
      # Friedman test for repeated measures (CV folds)
      results$friedman_tests <- list()
      for (metric in metrics) {
        metric_data <- long_data[long_data$Metric == metric, ]
        if (nrow(metric_data) > 0) {
          friedman_result <- coin::friedman_test(
            Value ~ Sample | Fold, 
            data = metric_data
          )
          results$friedman_tests[[metric]] <- friedman_result
        }
      }
    }, error = function(e) {
      message("Non-parametric analysis failed: ", e$message)
    })
  }
  
  # 4. Effect size confidence intervals using MBESS
  if (requireNamespace("MBESS", quietly = TRUE)) {
    tryCatch({
      library(MBESS)
      results$effect_size_cis <- list()
      samples <- unique(long_data$Sample)
      
      for (i in 1:(length(samples)-1)) {
        for (j in (i+1):length(samples)) {
          s1 <- samples[i]
          s2 <- samples[j]
          
          for (metric in metrics) {
            x1 <- long_data$Value[long_data$Sample == s1 & long_data$Metric == metric]
            x2 <- long_data$Value[long_data$Sample == s2 & long_data$Metric == metric]
            
            if (length(x1) > 1 && length(x2) > 1) {
              ci_result <- MBESS::ci.smd(
                ncp = (mean(x1) - mean(x2)) / sqrt(((length(x1)-1)*var(x1) + (length(x2)-1)*var(x2))/(length(x1)+length(x2)-2)),
                n.1 = length(x1),
                n.2 = length(x2)
              )
              
              key <- paste(s1, "vs", s2, metric, sep = "_")
              results$effect_size_cis[[key]] <- ci_result
            }
          }
        }
      }
    }, error = function(e) {
      message("MBESS effect size CIs failed: ", e$message)
    })
  }
  
  return(results)
}

# Comprehensive statistical report generator
generate_statistical_report <- function(per_fold_table, summary_table, pooled_table,
                                       pairwise_accuracy, pairwise_kappa, 
                                       cv_comparisons, performance_tests,
                                       output_file = NULL) {
  
  report <- list()
  
  # 1. Summary Statistics with Confidence Intervals
  report$summary_stats <- summary_table
  
  # 2. Effect Size Interpretation Summary
  effect_size_summary <- aggregate(Cohen_d ~ Metric, data = summary_table, 
                                   FUN = function(x) mean(x, na.rm = TRUE))
  names(effect_size_summary)[2] <- "Mean_Effect_Size"
  effect_size_summary$Interpretation <- ifelse(abs(effect_size_summary$Mean_Effect_Size) < 0.2, "Negligible",
                                              ifelse(abs(effect_size_summary$Mean_Effect_Size) < 0.5, "Small", 
                                              ifelse(abs(effect_size_summary$Mean_Effect_Size) < 0.8, "Medium", "Large")))
  report$effect_sizes <- effect_size_summary
  
  # 3. Significant Pairwise Comparisons
  if (!is.null(pairwise_accuracy)) {
    report$significant_accuracy_comparisons <- pairwise_accuracy[pairwise_accuracy$Significant, ]
  }
  
  if (!is.null(pairwise_kappa)) {
    report$significant_kappa_comparisons <- pairwise_kappa[pairwise_kappa$Significant, ]
  }
  
  # 4. Cross-validation significance
  if (!is.null(cv_comparisons)) {
    report$significant_cv_comparisons <- cv_comparisons[cv_comparisons$Significant_Global, ]
  }
  
  # 5. Performance vs null hypothesis
  report$performance_significance <- performance_tests
  
  # 6. Overall model performance summary
  overall_summary <- data.frame(
    Metric = c("Samples_Analyzed", "Total_CV_Folds", "Significant_Comparisons", "Large_Effect_Sizes"),
    Value = c(
      length(unique(per_fold_table$Sample)),
      nrow(per_fold_table),
      sum(c(nrow(report$significant_accuracy_comparisons %||% data.frame()), 
            nrow(report$significant_kappa_comparisons %||% data.frame()))),
      sum(abs(summary_table$Cohen_d) > 0.8, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  report$overall_summary <- overall_summary
  
  # Save to file if requested
  if (!is.null(output_file)) {
    saveRDS(report, paste0(output_file, "_statistical_report.RDS"))
    cat("Statistical report saved to:", paste0(output_file, "_statistical_report.RDS"), "\n")
  }
  
  return(report)
}

# Quick interpretation function
interpret_results <- function(statistical_report) {
  cat("=== STEAM STATISTICAL ANALYSIS INTERPRETATION ===\n\n")
  
  cat("1. OVERALL PERFORMANCE:\n")
  sig_performance <- statistical_report$performance_significance[statistical_report$performance_significance$Significant, ]
  if (nrow(sig_performance) > 0) {
    for (i in 1:nrow(sig_performance)) {
      cat(sprintf("   - %s: Mean = %.3f, significantly better than null (p = %.3e, Effect = %s)\n",
                  sig_performance$Metric[i], sig_performance$Mean_Performance[i], 
                  sig_performance$p_adjusted[i], sig_performance$Interpretation[i]))
    }
  } else {
    cat("   - No metrics significantly better than null hypothesis\n")
  }
  
  cat("\n2. EFFECT SIZES:\n")
  large_effects <- statistical_report$effect_sizes[abs(statistical_report$effect_sizes$Mean_Effect_Size) > 0.8, ]
  if (nrow(large_effects) > 0) {
    for (i in 1:nrow(large_effects)) {
      cat(sprintf("   - %s: Large effect size (d = %.3f)\n", 
                  large_effects$Metric[i], large_effects$Mean_Effect_Size[i]))
    }
  } else {
    cat("   - No large effect sizes detected\n")
  }
  
  cat("\n3. SAMPLE COMPARISONS:\n")
  total_sig <- nrow(statistical_report$significant_accuracy_comparisons %||% data.frame()) + 
               nrow(statistical_report$significant_kappa_comparisons %||% data.frame())
  
  cat(sprintf("   - %d significant pairwise differences found (Bonferroni corrected)\n", total_sig))
  
  if (!is.null(statistical_report$significant_cv_comparisons) && 
      nrow(statistical_report$significant_cv_comparisons) > 0) {
    cat(sprintf("   - %d cross-metric comparisons remain significant after global correction\n",
                nrow(statistical_report$significant_cv_comparisons)))
  }
  
  cat("\n4. RECOMMENDATIONS:\n")
  mean_acc <- mean(statistical_report$summary_stats[statistical_report$summary_stats$Metric == "Accuracy", "Mean"], na.rm = TRUE)
  
  if (mean_acc > 0.8) {
    cat("   - Model shows strong performance across samples\n")
  } else if (mean_acc > 0.6) {
    cat("   - Model shows moderate performance - consider optimization\n")
  } else {
    cat("   - Model performance may need significant improvement\n")
  }
  
  if (total_sig > 0) {
    cat("   - Significant sample differences detected - investigate sample-specific factors\n")
  } else {
    cat("   - Consistent performance across samples\n")
  }
}

# Function to explore STEAM result structure and extract all available metrics
explore_steam_results <- function(out_nested, sample_id = "unknown") {
  cat("=== STEAM Results Structure for Sample", sample_id, "===\n")
  
  # Top level structure
  cat("Top level names:", paste(names(out_nested), collapse = ", "), "\n")
  
  # Nested structure
  if ("nested" %in% names(out_nested)) {
    cat("Nested names:", paste(names(out_nested$nested), collapse = ", "), "\n")
  }
  
  # Train structure  
  if ("train" %in% names(out_nested)) {
    cat("Train names:", paste(names(out_nested$train), collapse = ", "), "\n")
    
    if ("ncv" %in% names(out_nested$train)) {
      cat("NCV names:", paste(names(out_nested$train$ncv), collapse = ", "), "\n")
      
      # Outer results structure
      if ("outer_result" %in% names(out_nested$train$ncv)) {
        cat("Number of outer folds:", length(out_nested$train$ncv$outer_result), "\n")
        if (length(out_nested$train$ncv$outer_result) > 0) {
          cat("Outer fold 1 names:", paste(names(out_nested$train$ncv$outer_result[[1]]), collapse = ", "), "\n")
        }
      }
      
      # Summary structure
      if ("summary" %in% names(out_nested$train$ncv)) {
        cat("Summary names:", paste(names(out_nested$train$ncv$summary), collapse = ", "), "\n")
        cat("Summary values:\n")
        print(out_nested$train$ncv$summary)
      }
    }
  }
  
  cat("=== End Structure Exploration ===\n\n")
}

# Enhanced function to extract comprehensive metrics from STEAM results
extract_comprehensive_metrics <- function(out_nested, sample_id) {
  # First explore the structure
  explore_steam_results(out_nested, sample_id)
  
  # Extract metrics using multiple approaches
  results <- list()
  
  # 1. Outer fold metrics (detailed per fold)
  if (!is.null(out_nested$train$ncv$outer_result)) {
    results$per_fold_metrics <- compute_outer_fold_metrics(out_nested$train$ncv$outer_result)
  }
  
  # 2. Summary metrics (aggregated)
  if (!is.null(out_nested$train$ncv$summary)) {
    results$summary_metrics <- out_nested$train$ncv$summary
  }
  
  # 3. Pooled metrics (our computed version)
  results$pooled_metrics <- pooled_metrics_row(out_nested)
  
  # 4. Best model information
  if (!is.null(out_nested$train$ncv$final_result)) {
    results$final_model <- out_nested$train$ncv$final_result
  }
  
  # 5. Fold diagnostics if available
  if (!is.null(out_nested$train$ncv$fold_diagnostics)) {
    results$fold_diagnostics <- out_nested$train$ncv$fold_diagnostics
  }
  
  return(results)
}



## ---------- 2) driver over samples ----------

sample_ids <- c("151507","151508","151509","151510",
                "151669","151670","151671","151672",
                "151673","151674","151675","151676")

# storage
all_per_fold   <- list()
all_summary    <- list()
all_pooled     <- list()
all_out_nested <- vector("list", length(sample_ids))
names(all_out_nested) <- sample_ids

for (sid in sample_ids) {
  message("=== Running sample ", sid, " ===")
  DLPFC <- DLPFCs[[sid]]
  
  # 1) Load STEAM
  STEAM.obj <- LoadSTEAM(
    Seurat.obj   = DLPFC,
    label.column = "Labels",
    assay        = "SCT"
  )
  
  # 2) Run nested CV (tune args as you like)
  out_nested <- RunSTEAM(
    STEAM.obj, mode = "nested",
    n_outer_folds = 5, n_inner_folds = 3,
    cv.cores = 10,               # outer folds in parallel
    allowParallel = FALSE,       # keep caret inner loop sequential
    n.size = 5, model = "rf",
    metric = "Kappa",
    maxnweights = 15000
  )
  all_out_nested[[sid]] <- out_nested
  
  # 3) Extract comprehensive metrics
  comprehensive_results <- extract_comprehensive_metrics(out_nested, sid)
  
  # 4) Per-outer-fold metrics (enhanced)
  if (!is.null(comprehensive_results$per_fold_metrics)) {
    outer_df <- comprehensive_results$per_fold_metrics
    outer_df$Sample <- sid
    all_per_fold[[sid]] <- outer_df
    
    # Mean ± SD across outer folds  
    outer_summary <- summarize_outer_folds(outer_df)
    outer_summary$Sample <- sid
    all_summary[[sid]] <- outer_summary
  } else {
    # Fallback to original method
    if (!is.null(out_nested$train$ncv$outer_result)) {
      outer_df <- compute_outer_fold_metrics(out_nested$train$ncv$outer_result)
      outer_df$Sample <- sid
      all_per_fold[[sid]] <- outer_df
      
      outer_summary <- summarize_outer_folds(outer_df)
      outer_summary$Sample <- sid
      all_summary[[sid]] <- outer_summary
    }
  }
  
  # 5) Pooled metrics (enhanced extraction)
  pooled_df <- comprehensive_results$pooled_metrics
  pooled_df$Sample <- sid
  all_pooled[[sid]] <- pooled_df
}

## ---------- 3) bind results and perform statistical analysis ----------

per_fold_table   <- do.call(rbind, all_per_fold)      # one row per outer fold per sample
summary_table   <- do.call(rbind, all_summary)       # mean ± CI ± effect size per metric per sample
pooled_table   <- do.call(rbind, all_pooled)        # pooled metrics per sample (single row)

# Optional: add model label (useful if you compare methods)
pooled_table$Model  <- "STEAM_RF"
per_fold_table$Model <- "STEAM_RF" 
summary_table$Model  <- "STEAM_RF"

## ---------- 4) Statistical Analysis Using Established Packages ----------

# Pairwise comparisons using emmeans/PMCMRplus
pairwise_accuracy <- pairwise_sample_comparisons(all_per_fold, "Accuracy", "bonferroni")
pairwise_kappa <- pairwise_sample_comparisons(all_per_fold, "Kappa", "bonferroni")

# Cross-validation fold comparisons across multiple metrics
cv_comparisons <- cv_fold_comparisons(all_per_fold, 
                                     metrics = c("Accuracy", "Kappa", "ARI", "NMI", "AMI"),
                                     correction_method = "bonferroni")

# Performance testing using broom + effectsize packages
performance_tests <- statistical_performance_summary(all_per_fold, all_summary)

# Advanced statistical analysis (mixed-effects, Bayesian, non-parametric)
advanced_results <- advanced_cv_analysis(all_per_fold, metrics = c("Accuracy", "Kappa"))

# Validation and testing function
validate_comprehensive_extraction <- function(steam_result, sample_id = "test") {
  cat("=== STEAM Result Structure Validation ===\n")
  cat("Sample ID:", sample_id, "\n\n")
  
  # Test comprehensive extraction
  comprehensive_results <- extract_comprehensive_metrics(steam_result, sample_id)
  
  # Report results
  cat("Per-fold metrics:\n")
  if (!is.null(comprehensive_results$per_fold_metrics)) {
    print(head(comprehensive_results$per_fold_metrics))
    cat("Dimensions:", dim(comprehensive_results$per_fold_metrics), "\n")
    cat("Available metrics:", colnames(comprehensive_results$per_fold_metrics), "\n\n")
  } else {
    cat("No per-fold metrics extracted\n\n")
  }
  
  cat("Pooled metrics:\n")
  if (!is.null(comprehensive_results$pooled_metrics)) {
    print(comprehensive_results$pooled_metrics)
    cat("Available pooled metrics:", colnames(comprehensive_results$pooled_metrics), "\n\n")
  } else {
    cat("No pooled metrics extracted\n\n")
  }
  
  # Return for further inspection
  invisible(comprehensive_results)
}

## ---------- 5) Generate comprehensive statistical report ----------

statistical_report <- generate_statistical_report(
  per_fold_table, summary_table, pooled_table,
  pairwise_accuracy, pairwise_kappa, 
  cv_comparisons, performance_tests,
  output_file = "STEAM_analysis"  # Will save as STEAM_analysis_statistical_report.RDS
)

## ---------- 6) View results and interpretation ----------

# Core results with statistical enhancements
per_fold_table    # Individual fold performance
summary_table     # Statistical summaries with confidence intervals and effect sizes
pooled_table      # Pooled metrics

# Comprehensive interpretation
interpret_results(statistical_report)

# Detailed statistical tables
cat("\n=== DETAILED STATISTICAL RESULTS ===\n")
print("Summary Statistics with Confidence Intervals:")
print(statistical_report$summary_stats)

print("Effect Size Summary:")
print(statistical_report$effect_sizes)

if (!is.null(statistical_report$significant_accuracy_comparisons) && 
    nrow(statistical_report$significant_accuracy_comparisons) > 0) {
  print("Significant Accuracy Comparisons:")
  print(statistical_report$significant_accuracy_comparisons)
}

if (!is.null(statistical_report$performance_significance)) {
  print("Performance vs Null Hypothesis:")
  print(statistical_report$performance_significance)
}





## PLOT - Enhanced with Confidence Intervals

# Load required packages
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("scales", quietly = TRUE)) install.packages("scales")
library(ggplot2)
library(scales)

per_fold_table$Model  <- factor(per_fold_table$Model)
per_fold_table$Sample <- factor(per_fold_table$Sample, levels = unique(per_fold_table$Sample))

# Set up custom color palette for samples
# Generate enough colors for all samples (12 total)
n_samples <- length(levels(per_fold_table$Sample))
custom_colors <- scales::hue_pal()(n_samples)
names(custom_colors) <- levels(per_fold_table$Sample)

# Prepare summary data for confidence interval plots
summary_plot_data <- summary_table[!is.na(summary_table$Sample), ]
summary_plot_data$Sample <- factor(summary_plot_data$Sample, levels = unique(per_fold_table$Sample))

# Create a mapping of metric names to clean labels
metric_labels <- c(
  "Accuracy" = "Accuracy",
  "Kappa" = "Kappa", 
  "ARI" = "Adjusted Rand Index",
  "Macro_F1" = "F1 Score"
)

# Enhanced Kappa plot with confidence intervals
samples_kappa <- ggplot() +
  # Add boxplots from raw fold data
  geom_boxplot(data = per_fold_table, 
               aes(x = Sample, y = Kappa, fill = Sample), 
               outlier.shape = NA, width = 0.6, alpha = 0.7, linewidth = 0.3) +
  # Add individual points
  geom_jitter(data = per_fold_table,
              aes(x = Sample, y = Kappa), 
              color = "black", width = 0.15, alpha = 0.4, size = 1.6) +
  # Add confidence interval error bars from summary data
  geom_errorbar(data = summary_plot_data[summary_plot_data$Metric == "Kappa", ],
                aes(x = Sample, ymin = CI_Lower, ymax = CI_Upper),
                width = 0.3, size = 0.6, color = "black", alpha = 0.8) +
  # Add mean points
  geom_point(data = summary_plot_data[summary_plot_data$Metric == "Kappa", ],
             aes(x = Sample, y = Mean),
             color = "black", size = 3, shape = 18) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 25, color = "black", angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 25, color = "black"),
    axis.text.y  = element_text(size = 20, color = "black"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Kappa", 
       title = "Kappa Performance with 95% Confidence Intervals",
       subtitle = "Black diamonds = means, Black bars = 95% CI")

samples_kappa




# Enhanced Accuracy plot with confidence intervals
samples_accuracy <- ggplot() +
  # Add boxplots from raw fold data
  geom_boxplot(data = per_fold_table, 
               aes(x = Sample, y = Accuracy, fill = Sample), 
               outlier.shape = NA, width = 0.6, alpha = 0.7, linewidth = 0.3) +
  # Add individual points
  geom_jitter(data = per_fold_table,
              aes(x = Sample, y = Accuracy), 
              color = "black", width = 0.15, alpha = 0.4, size = 1.6) +
  # Add confidence interval error bars
  geom_errorbar(data = summary_plot_data[summary_plot_data$Metric == "Accuracy", ],
                aes(x = Sample, ymin = CI_Lower, ymax = CI_Upper),
                width = 0.3, size = 0.6, color = "black", alpha = 0.8) +
  # Add mean points
  geom_point(data = summary_plot_data[summary_plot_data$Metric == "Accuracy", ],
             aes(x = Sample, y = Mean),
             color = "black", size = 3, shape = 18) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 25, color = "black", angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 25, color = "black"),
    axis.text.y  = element_text(size = 20, color = "black"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Accuracy",
       title = "Accuracy Performance with 95% Confidence Intervals",
       subtitle = "Black diamonds = means, Black bars = 95% CI")

samples_accuracy


# Enhanced ARI plot with confidence intervals
samples_ari <- ggplot() +
  # Add boxplots from raw fold data
  geom_boxplot(data = per_fold_table, 
               aes(x = Sample, y = ARI, fill = Sample), 
               outlier.shape = NA, width = 0.6, alpha = 0.7, linewidth = 0.3) +
  # Add individual points
  geom_jitter(data = per_fold_table,
              aes(x = Sample, y = ARI), 
              color = "black", width = 0.15, alpha = 0.4, size = 1.6) +
  # Add confidence interval error bars
  geom_errorbar(data = summary_plot_data[summary_plot_data$Metric == "ARI", ],
                aes(x = Sample, ymin = CI_Lower, ymax = CI_Upper),
                width = 0.3, size = 0.6, color = "black", alpha = 0.8) +
  # Add mean points
  geom_point(data = summary_plot_data[summary_plot_data$Metric == "ARI", ],
             aes(x = Sample, y = Mean),
             color = "black", size = 3, shape = 18) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 25, color = "black", angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 25, color = "black"),
    axis.text.y  = element_text(size = 20, color = "black"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Adjusted Rand Index",
       title = "ARI Performance with 95% Confidence Intervals", 
       subtitle = "Black diamonds = means, Black bars = 95% CI")

samples_ari


# Enhanced F1 Score plot with confidence intervals (using pooled_table data as boxplot representation)
# Note: F1 Score comes from pooled metrics, converting to boxplot-style representation
samples_f1 <- ggplot() +
  # Create boxplot-style bars using geom_col to mimic boxplot appearance
  geom_col(data = pooled_table, 
           aes(x = Sample, y = Macro_F1, fill = Sample),
           alpha = 0.7, width = 0.6) +
  # Add central line to mimic boxplot median
  geom_segment(data = pooled_table,
               aes(x = as.numeric(Sample) - 0.25, xend = as.numeric(Sample) + 0.25, 
                   y = Macro_F1, yend = Macro_F1),
               color = "black", size = 0.6) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 25, color = "black", angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 25, color = "black"),
    axis.text.y  = element_text(size = 20, color = "black"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "F1 Score",
       title = "F1 Score Performance Across Samples",
       subtitle = "Macro F1 from pooled results")

samples_f1

# Function to create all plots with confidence intervals
create_metric_plots <- function(per_fold_data, summary_data, pooled_data = NULL, metrics = c("Accuracy", "Kappa", "ARI", "Macro_F1")) {
  
  # Ensure data is properly formatted
  per_fold_data$Sample <- factor(per_fold_data$Sample, levels = unique(per_fold_data$Sample))
  summary_data$Sample <- factor(summary_data$Sample, levels = unique(per_fold_data$Sample))
  
  # Set up custom color palette for samples
  # Generate enough colors for all samples
  n_samples <- length(levels(per_fold_data$Sample))
  custom_colors <- scales::hue_pal()(n_samples)
  names(custom_colors) <- levels(per_fold_data$Sample)
  
  plots <- list()
  
  # Create plot for each metric
  for (metric in metrics) {
    if (metric == "Macro_F1" && !is.null(pooled_data)) {
      # Special handling for F1 Score from pooled data - boxplot style
      plots[[metric]] <- ggplot() +
        # Create boxplot-style bars
        geom_col(data = pooled_data, 
                 aes(x = Sample, y = Macro_F1, fill = Sample),
                 alpha = 0.7, width = 0.6) +
        # Add central line to mimic boxplot median
        geom_segment(data = pooled_data,
                     aes(x = as.numeric(Sample) - 0.25, xend = as.numeric(Sample) + 0.25, 
                         y = Macro_F1, yend = Macro_F1),
                     color = "black", size = 0.6) +
        scale_fill_manual(values = custom_colors) +
        theme_minimal() +
        theme(
          axis.title.x = element_blank(),
          axis.text.x  = element_text(size = 12, color = "black", angle = 45, hjust = 1),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 14, color = "black"),
          axis.text.y  = element_text(size = 10, color = "black"),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.8),
          legend.position = "none",
          plot.title = element_text(size = 16, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray50")
        ) +
        coord_cartesian(ylim = c(0, 1)) +
        labs(y = "F1 Score",
             title = "F1 Score Performance Across Samples",
             subtitle = "Boxplot-style representation")
    } else {
      # Regular plots for metrics with confidence intervals
      plots[[metric]] <- ggplot() +
        # Boxplots from raw fold data
        geom_boxplot(data = per_fold_data, 
                     aes(x = Sample, y = .data[[metric]], fill = Sample), 
                     outlier.shape = NA, width = 0.6, alpha = 0.7, linewidth = 0.3) +
        # Individual points
        geom_jitter(data = per_fold_data,
                    aes(x = Sample, y = .data[[metric]]), 
                    color = "black", width = 0.15, alpha = 0.4, size = 1.6) +
        # Confidence interval error bars
        geom_errorbar(data = summary_data[summary_data$Metric == metric, ],
                      aes(x = Sample, ymin = CI_Lower, ymax = CI_Upper),
                      width = 0.3, size = 0.6, color = "black", alpha = 0.8) +
        # Mean points
        geom_point(data = summary_data[summary_data$Metric == metric, ],
                   aes(x = Sample, y = Mean),
                   color = "black", size = 3, shape = 18) +
        scale_fill_manual(values = custom_colors) +
        theme_minimal() +
        theme(
          axis.title.x = element_blank(),
          axis.text.x  = element_text(size = 12, color = "black", angle = 45, hjust = 1),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 14, color = "black"),
          axis.text.y  = element_text(size = 10, color = "black"),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.8),
          legend.position = "none",
          plot.title = element_text(size = 16, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray50")
        ) +
        coord_cartesian(ylim = c(0, 1)) +
        labs(y = metric_labels[metric] %||% metric,
             title = paste(metric, "Performance Across Samples"),
             subtitle = "Black diamonds = means, Black bars = 95% CI")
    }
  }
  
  return(plots)
}

# Create all plots for the four selected metrics
metric_plots <- create_metric_plots(per_fold_table, summary_plot_data, pooled_table, 
                                   metrics = c("Accuracy", "Kappa", "ARI", "Macro_F1"))

# Display individual plots
metric_plots$Accuracy
metric_plots$Kappa  
metric_plots$ARI
metric_plots$Macro_F1

# Optional: Create a combined plot using patchwork (if available)
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  
  combined_plot <- metric_plots$Accuracy + metric_plots$Kappa + 
                   metric_plots$ARI + metric_plots$Macro_F1 +
                   plot_layout(ncol = 2, nrow = 2)
  
  print("Combined plot created - use 'combined_plot' to display")
  combined_plot
} else {
  message("Install 'patchwork' package to create combined plots: install.packages('patchwork')")
}

# Function to save plots
save_metric_plots <- function(plots, prefix = "STEAM_metrics", width = 12, height = 8, dpi = 300) {
  for (metric_name in names(plots)) {
    filename <- paste0(prefix, "_", metric_name, ".png")
    ggsave(filename, plots[[metric_name]], width = width, height = height, dpi = dpi)
    message("Saved: ", filename)
  }
  
  # Save combined plot if available
  if (exists("combined_plot")) {
    combined_filename <- paste0(prefix, "_combined.png")
    ggsave(combined_filename, combined_plot, width = 16, height = 12, dpi = dpi)
    message("Saved: ", combined_filename)
  }
}

# Uncomment to save plots
# save_metric_plots(metric_plots)
