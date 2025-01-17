model.train <- function(STEAM.obj, model, n.tree, kernel, cv.folds, cv.repeats, trainval.ratio, train.folder.name, allowParallel, maxnweights = 5000) {
  
  # Prepare the data frame
  df <- data.frame(labels = STEAM.obj$train$train.data.labels, t(STEAM.obj$train$avg.matrix))
  df$labels <- factor(df$labels, levels = unique(df$labels), labels = make.names(unique(df$labels)))
  
  for (i in 2:ncol(df)) {
    df[, i] <- as.numeric(df[, i])
  }
  
  train_control <- trainControl(
    method = "repeatedcv",
    number = cv.folds,
    repeats = cv.repeats,
    savePredictions = "all",
    classProbs = TRUE,
    p = trainval.ratio,
    search = 'grid',
    allowParallel = allowParallel
  )
  
  if (model == "rf") {
    train_model <- train(labels ~ ., data = df, method = "rf", trControl = train_control, ntree = n.tree, metric = 'Kappa')
    message("Finished RandomForest model building with cross-validation")
    
  } else if (model == "svm" & kernel == 'radial') {
    train_model <- train(labels ~ ., data = df, method = "svmRadial", trControl = train_control, metric = 'Kappa')
    message("Finished SVM-Radial model building with cross-validation")
    
  } else if (model == "svm" & kernel == 'linear') {
    train_model <- train(labels ~ ., data = df, method = "svmLinear", trControl = train_control, metric = 'Kappa')
    message("Finished SVM-Linear model building with cross-validation")
    
  } else if (model == "xgb") {
    train_model <- train(labels ~ ., data = df, method = "xgbTree", trControl = train_control, metric = 'Kappa')
    message("Finished Xgboost model building with cross-validation")
    
  } else if (model == "multinom") {
    train_model <- train(labels ~ ., data = df, method = "multinom", MaxNWts = maxnweights, trControl = train_control, metric = 'Kappa')
    message("Finished Multinomial model building with cross-validation")
    
  } else {
    stop("Model type is not supported. Please choose 'rf', 'svm', 'xgb', or 'multinom'.")
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