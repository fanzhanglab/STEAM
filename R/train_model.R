train_models <- function(matrix, labels, model, seed, n_tree, kernel, scale) {
  set.seed(seed)
  
  df <- data.frame(labels = labels, t(matrix))
  df$labels <- factor(df$labels)
  
  for (i in 2:ncol(df)) {
    df[, i] <- as.numeric(df[,i])
  }
  
  if (model == "rf") {
    rf_model <- randomForest(labels ~ ., data = df, proximity = FALSE, ntree = n_tree)
    message("Finished RandomForest model building")
    return(rf_model)
    
  } else if (model == "svm") {
    svm_model <- svm(labels ~ ., data = df, kernel = kernel, scale = scale)
    message("Finished SVM model building")
    return(svm_model)
    
  } else {
    stop("Model type is not supported. Please choose 'rf' or 'svm'.")
  }

}
