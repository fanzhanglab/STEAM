library(caret)
library(pROC)
library(ggplot2)
library(reshape2)

model.train <- function(matrix, labels, model, seed, k = 5, rf_grid = NULL, svm_grid = NULL) {
  set.seed(seed)

  # Prepare data frame
  df <- data.frame(labels = factor(labels), t(matrix))
  df[-1] <- lapply(df[-1], as.numeric)

  # Cross-validation setup with accuracy metric (adaptable for multiclass)
  train_control <- trainControl(
    method = "cv",
    number = k,
    search = "grid",
    classProbs = TRUE,
    summaryFunction = multiClassSummary,  # Summary for multiclass metrics
    savePredictions = "final"
  )

  # Create "train" directory if it doesn't exist
  if (!dir.exists("train")) {
    dir.create("train")
  }

  # Function to save ROC plots for each class in multiclass setting
  save_multiclass_roc_plot <- function(predictions, labels, model_name) {
    levels_list <- levels(labels)
    for (level in levels_list) {
      roc_obj <- roc(as.numeric(predictions == level), as.numeric(labels == level))
      roc_auc <- auc(roc_obj)

      # Plot ROC curve for the current class
      roc_plot <- ggroc(roc_obj) +
        ggtitle(paste(model_name, "ROC Curve - Class", level, "(AUC:", round(roc_auc, 2), ")")) +
        xlab("False Positive Rate") +
        ylab("True Positive Rate")

      # Save plot
      ggsave(filename = paste0("train/", model_name, "_ROC_Curve_Class_", level, ".png"), plot = roc_plot)
    }
  }

  if (model == "rf") {
    # Default grid if none provided
    if (is.null(rf_grid)) {
      rf_grid <- expand.grid(mtry = c(1, floor(sqrt(ncol(df) - 1)), floor((ncol(df) - 1) / 3)))
    }

    rf_model <- train(
      labels ~ .,
      data = df,
      method = "rf",
      trControl = train_control,
      tuneGrid = rf_grid,
      metric = "Accuracy"
    )

    message("Finished RandomForest model building with cross-validation and grid search")

    # Save ROC plots for each class
    save_multiclass_roc_plot(rf_model$pred$pred, rf_model$pred$obs, "RandomForest")

    # Save model and metrics
    saveRDS(rf_model, file = "train/rf_model.rds")
    write.csv(rf_model$results, file = "train/rf_metrics.csv", row.names = FALSE)

    return(list(model = rf_model, metrics = rf_model$results))

  } else if (model == "svm") {
    # Default grid if none provided
    if (is.null(svm_grid)) {
      svm_grid <- expand.grid(C = c(0.1, 1, 10), sigma = c(0.01, 0.05, 0.1))
    }

    svm_model <- train(
      labels ~ .,
      data = df,
      method = "svmRadial",
      trControl = train_control,
      tuneGrid = svm_grid,
      preProcess = c("center", "scale"),
      metric = "Accuracy"
    )

    message("Finished SVM model building with cross-validation and grid search")

    # Save ROC plots for each class
    save_multiclass_roc_plot(svm_model$pred$pred, svm_model$pred$obs, "SVM")

    # Save model and metrics
    saveRDS(svm_model, file = "train/svm_model.rds")
    write.csv(svm_model$results, file = "train/svm_metrics.csv", row.names = FALSE)

    return(list(model = svm_model, metrics = svm_model$results))

  } else {
    stop("Model type not supported. Please choose 'rf' (for RandomForest) or 'svm' (for SupportVectorMachine).")
  }
}
