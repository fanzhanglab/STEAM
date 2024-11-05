model.predict <- function(model, test_matrix, test_labels) {
  df <- data.frame(labels = test_labels, t(test_matrix))
  df$labels <- factor(df$labels)

  for (i in 2:ncol(df)) {
    df[, i] <- as.numeric(df[, i])
  }

  p <- predict(model, df)
  p <- factor(p, levels = levels(df$labels))

  confusion_matrix <- table(Predicted = p, Actual = df$labels)
  true_positive <- sum(diag(confusion_matrix))
  total <- sum(confusion_matrix)
  accuracy <- true_positive / total

  precision <- diag(confusion_matrix) / rowSums(confusion_matrix)
  recall <- diag(confusion_matrix) / colSums(confusion_matrix)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  f1_score[is.nan(f1_score)] <- 0 # Replace NaN values with zero, if any

  metrics <- list(
    accuracy = accuracy,
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    confusion_matrix = confusion_matrix
  )
  return(list(predictions = p, true_labels = df$labels, metrics = metrics))
}
