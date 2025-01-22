#' Predict Labels and Evaluate Model Performance
#'
#' This function takes a STEAM object, predicts labels for the test dataset using the trained model, and calculates performance metrics, including accuracy, precision, recall, F1-score, mean balanced accuracy, Adjusted Rand Index (ARI), and Kappa.
#'
#' @param Object of class STEAM.Object
#' @importFrom mclust adjustedRandIndex
#' @importFrom caret confusionMatrix
#' @export
model.predict <- function(STEAM.obj) {
  # Extract test data and model from STEAM.obj
  test_matrix <- STEAM.obj$test$avg.matrix
  test_labels <- STEAM.obj$test$test.data.labels
  model <- STEAM.obj$train$model
  
  # Prepare the data frame for prediction
  df <- data.frame(labels = test_labels, t(test_matrix))
  df$labels <- factor(df$labels)
  
  for (i in 2:ncol(df)) {
    df[, i] <- as.numeric(df[, i])
  }
  
  # Make predictions
  p <- predict(model, df)
  p <- factor(p, levels = levels(df$labels))
  
  # Construct confusion matrix
  confusion_matrix <- table(Predicted = p, Actual = df$labels)
  
  # Calculate accuracy
  true_positive <- sum(diag(confusion_matrix))
  total <- sum(confusion_matrix)
  accuracy <- true_positive / total
  
  # Calculate precision, recall, and F1-score
  precision <- diag(confusion_matrix) / rowSums(confusion_matrix)
  recall <- diag(confusion_matrix) / colSums(confusion_matrix)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  f1_score[is.nan(f1_score)] <- 0 # Replace NaN values with zero, if any
  
  # Calculate specificity for each class
  specificity <- sapply(1:nrow(confusion_matrix), function(i) {
    true_negatives <- sum(confusion_matrix[-i, -i])
    false_positives <- sum(confusion_matrix[i, -i])
    true_negatives / (true_negatives + false_positives)
  })
  
  # Calculate mean balanced accuracy
  balanced_accuracy <- (recall + specificity) / 2
  mean_balanced_accuracy <- mean(balanced_accuracy, na.rm = TRUE)
  
  # Calculate Adjusted Rand Index (ARI)
  ari <- mclust::adjustedRandIndex(as.character(test_labels), as.character(p))
  
  # Calculate Kappa using confusionMatrix from caret package
  kappa <- caret::confusionMatrix(p, df$labels)$overall["Kappa"]
  
  # Compile all metrics into a list
  metrics <- list(
    accuracy = accuracy,
    mean_balanced_accuracy = mean_balanced_accuracy,
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    confusion_matrix = confusion_matrix,
    ARI = ari,
    kappa = kappa
  )
  
  # Update STEAM.obj with predictions and metrics
  STEAM.obj$test$predictions <- p
  STEAM.obj$test$true.labels <- df$labels
  STEAM.obj$test$metrics <- metrics
  
  return(STEAM.obj)
}