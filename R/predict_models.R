predict_models <- function(models, test_df) {
  test_df <- data.frame(cbind(test_data$labs, t(test_data$mat)))
  test_df$labs <- factor(test_df$labs)
  rf_pred <- predict(models$rf, test_df)
  svm_pred <- predict(models$svm, test_df)
  return(list(rf_pred = rf_pred, svm_pred = svm_pred))
}
