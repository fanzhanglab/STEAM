library(caret)
library(pROC)
library(ggplot2)

model.train <- function(matrix, labels, model, seed, k = 5, rf_grid = NULL, svm_grid = NULL) {
  set.seed(seed)

  # Prepare data frame
  df <- data.frame(labels = factor(labels), t(matrix))
  df[-1] <- lapply(df[-1], as.numeric)


  # Cross-validation setup with ROC metrics
  train_control <- trainControl(
    method = "cv",
    number = k,
    search = "grid",
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )

  # Create "train" directory if it doesn't exist
  if (!dir.exists("train")) {
    dir.create("train")
  }

  # Function to save ROC plot
  save_roc_plot <- function(roc_obj, model_name) {
    roc_plot <- ggroc(roc_obj) + ggtitle(paste(model_name, "ROC Curve")) +
      xlab("False Positive Rate") + ylab("True Positive Rate")
    ggsave(filename = paste0("train/", model_name, "_ROC_Curve.png"), plot = roc_plot)
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
      metric = "ROC"
    )

    message("Finished RandomForest model building with cross-validation and grid search")

    # Calculate ROC and AUC for Random Forest
    roc_obj <- roc(rf_model$pred$obs, rf_model$pred$Class1, levels = rev(levels(df$labels)))
    rf_auc <- auc(roc_obj)

    # Save model, ROC plot, and metrics
    saveRDS(rf_model, file = "train/rf_model.rds")
    save_roc_plot(roc_obj, "RandomForest")
    write.csv(rf_model$results, file = "train/rf_metrics.csv", row.names = FALSE)

    return(list(model = rf_model, auc = rf_auc, metrics = rf_model$results))

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
      metric = "ROC"
    )

    message("Finished SVM model building with cross-validation and grid search")

    # Calculate ROC and AUC for SVM
    roc_obj <- roc(svm_model$pred$obs, svm_model$pred$Class1, levels = rev(levels(df$labels)))
    svm_auc <- auc(roc_obj)

    # Save model, ROC plot, and metrics
    saveRDS(svm_model, file = "train/svm_model.rds")
    save_roc_plot(roc_obj, "SVM")
    write.csv(svm_model$results, file = "train/svm_metrics.csv", row.names = FALSE)

    return(list(model = svm_model, auc = svm_auc, metrics = svm_model$results))

  } else {
    stop("Model type not supported. Please choose 'rf' (for RandomForest) or 'svm' (for SupportVectorMachine).")
  }
}
