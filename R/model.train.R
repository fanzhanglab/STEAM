library(randomForest)
library(e1071)
library(caret)

model.train <- function(matrix, labels, model, seed = 123, n_tree = 200, kernel = NULL, scale = TRUE, cv_folds = 10, folder_name = 'train.out') {
  set.seed(seed)

  # Prepare the data frame
  df <- data.frame(labels = labels, t(matrix))
  df$labels <- factor(df$labels)
  for (i in 2:ncol(df)) {
    df[, i] <- as.numeric(df[, i])
  }

  train_control <- trainControl(method = "repeatedcv", number = cv_folds, repeats = 3, savePredictions = "all", classProbs = TRUE, p = 0.8, search = 'grid', allowParallel = TRUE)

  if (model == "rf") {

    train_model <- train(labels ~ ., data = df, method = "rf", trControl = train_control, ntree = n_tree, metric = 'Kappa')
    message("Finished RandomForest model building with cross-validation")

  } else if (model == "svm") {
    # SVM model with cross-validation
    train_model <- train(labels ~ ., data = df, method = "svmRadial", trControl = train_control, metric = 'Kappa')
    message("Finished SVM model building with cross-validation")


  } else {
    stop("Model type is not supported. Please choose 'rf' or 'svm'.")
  }


  # Create output folder if it doesn't exist
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }


  saveRDS(train_model, file = file.path(folder_name, "train.model.rds"))

  # Return Model and Metrics
  return(train_model)

}
