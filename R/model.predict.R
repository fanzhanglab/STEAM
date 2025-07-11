model.predict <- function(STEAM.obj, saved.folder = "./") {
  # Extract test data and model
  test_matrix <- STEAM.obj$test$avg.matrix
  test_labels <- STEAM.obj$test$test.data.labels
  model <- STEAM.obj$train$model

  # Prepare data frame for prediction
  df <- data.frame(labels = test_labels, t(test_matrix))
  df$labels <- factor(df$labels)
  for (i in 2:ncol(df)) df[, i] <- as.numeric(df[, i])

  # Predict
  p <- predict(model, df)
  p <- factor(p, levels = levels(df$labels))
  true_labels <- as.character(test_labels)
  pred_labels <- as.character(p)

  # Confusion matrix
  confusion_matrix <- table(Predicted = p, Actual = df$labels)
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)

  # Precision, Recall, F1
  precision <- diag(confusion_matrix) / rowSums(confusion_matrix)
  recall <- diag(confusion_matrix) / colSums(confusion_matrix)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  f1_score[is.nan(f1_score)] <- 0

  # Specificity
  specificity <- sapply(1:nrow(confusion_matrix), function(i) {
    tn <- sum(confusion_matrix[-i, -i])
    fp <- sum(confusion_matrix[i, -i])
    tn / (tn + fp)
  })

  # Balanced Accuracy
  balanced_accuracy <- (recall + specificity) / 2
  mean_balanced_accuracy <- mean(balanced_accuracy, na.rm = TRUE)

  # External metrics
  ari <- mclust::adjustedRandIndex(true_labels, pred_labels)
  kappa <- caret::confusionMatrix(p, df$labels)$overall["Kappa"]
  nmi <- aricode::NMI(true_labels, pred_labels)
  ami <- aricode::AMI(true_labels, pred_labels)

  # Compute PAS
  compute_PAS <- function(coords, labels, k = 10, threshold = 6) {
    coords <- as.matrix(coords)
    labels <- as.character(labels)
    knn_result <- get.knn(coords, k = k)
    neighbors_idx <- knn_result$nn.index
    abnormal_flags <- sapply(seq_len(nrow(coords)), function(i) {
      neighbor_labels <- labels[neighbors_idx[i, ]]
      sum(neighbor_labels != labels[i]) >= threshold
    })
    list(score = round(mean(abnormal_flags), 3), flags = abnormal_flags)
  }
  pas_result <- compute_PAS(STEAM.obj$test$test.data.coords, pred_labels)

  # Compile results
  metrics <- list(
    accuracy = round(accuracy, 3),
    mean_balanced_accuracy = round(mean_balanced_accuracy, 3),
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    confusion_matrix = confusion_matrix,
    ARI = round(ari, 3),
    NMI = round(nmi, 3),
    AMI = round(ami, 3),
    kappa = round(as.numeric(kappa), 3),
    PAS = pas_result$score
  )

  # Final storage
  STEAM.obj$test$predictions <- p
  STEAM.obj$test$true.labels <- df$labels
  STEAM.obj$test$metrics <- metrics
  STEAM.obj$test$PAS_flags <- pas_result$flags

  return(STEAM.obj)
}
