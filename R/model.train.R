library(randomForest)
library(e1071)
library(caret)

model.train <- function(matrix, labels, model, seed, n_tree = NULL, kernel = NULL, scale = TRUE, cv_folds = 5, folder_name = 'train.out') {
  set.seed(seed)

  # Prepare the data frame
  df <- data.frame(labels = labels, t(matrix))
  df$labels <- factor(df$labels)
  for (i in 2:ncol(df)) {
    df[, i] <- as.numeric(df[, i])
  }

  train_index <- createDataPartition(df$labels, p = 0.8, list = FALSE)
  train_data <- df[train_index, ]
  valid_data <- df[-train_index, ]

  train_control <- trainControl(method = "cv", number = cv_folds)

  if (model == "rf") {
    # Random Forest model with cross-validation
    train_model <- train(labels ~ ., data = train_data, method = "rf", trControl = train_control, ntree = n_tree)
    message("Finished RandomForest model building with cross-validation")

  } else if (model == "svm") {
    # SVM model with cross-validation
    train_model <- train(labels ~ ., data = train_data, method = "svmRadial", trControl = train_control,
                       preProcess = c("center", "scale"))
    message("Finished SVM model building with cross-validation")


  } else {
    stop("Model type is not supported. Please choose 'rf' or 'svm'.")
  }
  train_pred <- predict(train_model, train_data, type = "prob")
  valid_pred <- predict(train_model, valid_data, type = "prob")
  train_roc <- multiclass.roc(train_data$labels, train_pred)
  valid_roc <- multiclass.roc(valid_data$labels, valid_pred)

  # Create output folder if it doesn't exist
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }

  # # Generate ROC Plot for Training and Validation
  # plot_roc <- function(roc_obj, title) {
  #   aucs <- sapply(roc_obj$rocs, function(r) auc(r))
  #   avg_auc <- mean(aucs)
  #   roc_plot <- ggplot() +
  #     labs(title = paste0(title, "\nMean AUC: ", round(avg_auc, 3))) +
  #     theme_minimal()
  #   for (r in roc_obj$rocs) {
  #     roc_df <- data.frame(Sensitivity = r$sensitivities, Specificity = 1 - r$specificities)
  #     roc_plot <- roc_plot + geom_line(data = roc_df, aes(x = Specificity, y = Sensitivity), alpha = 0.5)
  #   }
  #   return(roc_plot)
  # }
  #
  # # Save Plots to PDF
  # pdf(file.path(folder_name, "train_roc_plots.pdf"), width = 10, height = 5)
  # train_roc_plot <- plot_roc(train_roc, "Training ROC Curves")
  # valid_roc_plot <- plot_roc(valid_roc, "Validation ROC Curves")
  # grid.arrange(train_roc_plot, valid_roc_plot, ncol = 2)
  # dev.off()




  # Save Model and Metrics
  train_acc <- mean(predict(train_model, train_data) == train_data$labels)
  valid_acc <- mean(predict(train_model, valid_data) == valid_data$labels)
  print(paste("Training Accuracy:", train_acc))
  print(paste("Validation Accuracy:", valid_acc))

  saveRDS(train_model, file = file.path(folder_name, "train.model.rds"))

  # Return Model and Metrics
  return(train_model)

}
