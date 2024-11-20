library(randomForest)
library(e1071)
library(caret)
library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)

model.train <- function(matrix, labels, model, seed, n_tree = NULL, kernel = NULL, scale = TRUE, cv_folds = 10, folder_name = 'train.out') {
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

  train_control <- trainControl(method = "cv", number = cv_folds, savePredictions = "all", classProbs = TRUE, p = 0.8, search = 'grid', allowParallel = TRUE)

  if (model == "rf") {
    # Random Forest model with cross-validation
    tune_grid <- expand.grid(mtry = var_seq(p = sqrt(ncol(train_data)), classification = T, len=10))
    train_model <- train(labels ~ ., data = train_data, method = "rf", trControl = train_control, ntree = 200, tunegrid = tune_grid, metric = 'Kappa')
    message("Finished RandomForest model building with cross-validation")

  } else if (model == "svm") {
    # SVM model with cross-validation
    train_model <- train(labels ~ ., data = train_data, method = "svmRadial", trControl = train_control, metric = 'Kappa')
    message("Finished SVM model building with cross-validation")


  } else {
    stop("Model type is not supported. Please choose 'rf' or 'svm'.")
  }

  stopCluster(cl)

  # Create output folder if it doesn't exist
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }


  train_acc <- mean(predict(train_model, train_data) == train_data$labels)
  print(paste("Validation Accuracy:", train_acc))

  valid_acc <- mean(predict(train_model, valid_data) == valid_data$labels)
  print(paste("Validation Accuracy:", valid_acc))

  saveRDS(train_model, file = file.path(folder_name, "train.model.rds"))

  # Return Model and Metrics
  return(train_model)

}
