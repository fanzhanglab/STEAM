model.train <- function(STEAM.obj, model, n.tree, kernel, cv.folds, cv.repeats, trainval.ratio, train.folder.name, allowParallel) {

  # Prepare the data frame
  df <- data.frame(labels = STEAM.obj$train$train.data.labels, t(STEAM.obj$train$avg.matrix))
  df$labels <- factor(df$labels)
  for (i in 2:ncol(df)) {
    df[, i] <- as.numeric(df[, i])
  }

  train_control <- trainControl(method = "repeatedcv", number = cv.folds, repeats = cv.repeats, savePredictions = "all", classProbs = TRUE, p = trainval.ratio, search = 'grid', allowParallel = allowParallel)

  if (model == "rf") {

    train_model <- train(labels ~ ., data = df, method = "rf", trControl = train_control, ntree = n.tree, metric = 'Kappa')
    message("Finished RandomForest model building with cross-validation")

  } else if (model == "svm" & kernel == 'radial') {
    # SVM model with cross-validation
    train_model <- train(labels ~ ., data = df, method = "svmRadial", trControl = train_control, metric = 'Kappa')
    message("Finished SVM model building with cross-validation")


  } else if (model == "svm" & kernel == 'radial') {
    # SVM model with cross-validation
    train_model <- train(labels ~ ., data = df, method = "svmLinear", trControl = train_control, metric = 'Kappa')
    message("Finished SVM model building with cross-validation")


  } else if (model == "xgb") {
    # SVM model with cross-validation
    train_model <- train(labels ~ ., data = df, method = "xgbtree", trControl = train_control, metric = 'Kappa')
    message("Finished SVM model building with cross-validation")


  } else if (model == "multinom") {
    # SVM model with cross-validation
    train_model <- train(labels ~ ., data = df, method = "multinom", trControl = train_control, metric = 'Kappa')
    message("Finished SVM model building with cross-validation")


  } else {
    stop("Model type is not supported. Please choose 'rf', 'svm'.")
  }


  # Create output folder if it doesn't exist
  if (!dir.exists(train.folder.name)) {
    dir.create(train.folder.name)
  }


  saveRDS(train_model, file = file.path(train.folder.name, "train.model.rds"))
  STEAM.obj$train$model <- train_model

  # Return Model and Metrics
  return(STEAM.obj)

}
