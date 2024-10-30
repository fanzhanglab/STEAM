train_models <- function(train_data) {
  train_df <- data.frame(cbind(train_data$labs, t(train_data$mat)))
  train_df$labs <- factor(train_df$labs)
  rf_model <- randomForest(labs ~ ., data = train_df, proximity = FALSE, ntree = 100)
  svm_model <- svm(labs ~ ., data = train_df, kernel = "radial", scale = TRUE)
  return(list(rf = rf_model, svm = svm_model))
}
