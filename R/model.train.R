model.train <- function(STEAM.obj, model, n.tree = 50, kernel = NULL,
                        cv.folds = 2, cv.repeats = 1, trainval.ratio = 0.8,
                        train.folder.name = "./model", allowParallel = TRUE,
                        maxnweights = 5000, tune.grid = NULL) {

  # Prepare data frame
  df <- data.frame(labels = STEAM.obj$train$train.data.labels, t(STEAM.obj$train$avg.matrix))
  df$labels <- factor(df$labels, levels = unique(df$labels), labels = make.names(unique(df$labels)))
  for (i in 2:ncol(df)) df[, i] <- as.numeric(df[, i])
  
  train_control <- trainControl(
    method = "repeatedcv",
    number = cv.folds,
    repeats = cv.repeats,
    savePredictions = "all",
    classProbs = TRUE,
    search = "grid",
    allowParallel = allowParallel
  )
  
  # Train model with appropriate defaults or user grid
  if (model == "rf") {
    grid <- if (is.null(tune.grid)) expand.grid(mtry = c(5, 10, 20)) else tune.grid
    train_model <- train(labels ~ ., data = df, method = "rf",
                         trControl = train_control, tuneGrid = grid,
                         ntree = n.tree, metric = "Kappa")
    message("Finished RandomForest model building with cross-validation")
    
  } else if (model == "svm" && kernel == "radial") {
    grid <- if (is.null(tune.grid)) expand.grid(C = c(0.1, 1, 10), sigma = c(0.01, 0.05, 0.1)) else tune.grid
    train_model <- train(labels ~ ., data = df, method = "svmRadial",
                         trControl = train_control, tuneGrid = grid,
                         metric = "Kappa")
    message("Finished SVM-Radial model building with cross-validation")
    
  } else if (model == "svm" && kernel == "linear") {
    grid <- if (is.null(tune.grid)) expand.grid(C = c(0.1, 1, 10)) else tune.grid
    train_model <- train(labels ~ ., data = df, method = "svmLinear",
                         trControl = train_control, tuneGrid = grid,
                         metric = "Kappa")
    message("Finished SVM-Linear model building with cross-validation")
    
  } else if (model == "xgb") {
    grid <- if (is.null(tune.grid)) expand.grid(
      nrounds = c(100, 200),
      max_depth = c(4, 6),
      eta = c(0.01, 0.1),
      gamma = 0,
      colsample_bytree = 0.8,
      min_child_weight = 1,
      subsample = 0.8
    ) else tune.grid
    train_model <- train(labels ~ ., data = df, method = "xgbTree",
                         trControl = train_control, tuneGrid = grid,
                         metric = "Kappa")
    message("Finished XGBoost model building with cross-validation")
    
  } else if (model == "multinom") {
    grid <- if (is.null(tune.grid)) expand.grid(decay = c(0, 0.01, 0.1)) else tune.grid
    train_model <- train(labels ~ ., data = df, method = "multinom",
                         MaxNWts = maxnweights, trControl = train_control,
                         tuneGrid = grid, metric = "Kappa")
    message("Finished Multinomial model building with cross-validation")
    
  } else {
    stop("Model type is not supported. Please choose 'rf', 'svm', 'xgb', or 'multinom'.")
  }
  
  # Save model + best hyperparameters
  if (!dir.exists(train.folder.name)) dir.create(train.folder.name)
  saveRDS(train_model, file = file.path(train.folder.name, "train.model.rds"))
  
  STEAM.obj$train$model <- train_model
  STEAM.obj$train$best.params <- train_model$bestTune  # <- best hyperparameter logging
  
  return(STEAM.obj)
}

