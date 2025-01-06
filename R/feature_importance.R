feature_importance <- function(STEAM.obj, top_n = 10, title = "Top Features by Importance") {

  # Validate input
  if (!is.list(STEAM.obj) || !"train" %in% names(STEAM.obj)) {
    stop("Invalid STEAM.obj structure. It must contain 'train'.")
  }
  if (!is.numeric(top_n) || top_n <= 0) stop("'top_n' must be a positive numeric value.")

  # Extract model and method
  model <- STEAM.obj$train$model
  method <- model$method  # Extract method (e.g., "svmLinear", "xgbTree")
  importance <- NULL

  # Compute Feature Importance based on model method
  if (method == "rf") {
    # Random Forest
    importance <- varImp(model)$importance$Overall
    names(importance) <- rownames(varImp(model)$importance)

  } else if (method == "xgbTree") {
    # XGBoost
    importance <- varImp(model)$importance$Overall
    names(importance) <- rownames(varImp(model)$importance)

  } else if (method == "svmLinear") {
    # SVM Linear (kernlab)
    svm_model <- model$finalModel

    # Compute weights for linear SVM
    coefs <- colSums(svm_model@coef[[1]] * svm_model@xmatrix[[1]])  # Compute weights
    importance <- abs(coefs)  # Absolute values for importance
    names(importance) <- colnames(svm_model@xmatrix[[1]])  # Assign feature names

  } else if (method == "multinom") {
    # Multinomial Logistic Regression
    coef_matrix <- coef(model$finalModel)  # Coefficients per class
    importance <- apply(abs(coef_matrix), 2, sum)  # Sum across classes
    importance <- importance[!grepl("(Intercept)", names(importance))]

  } else {
    stop(paste("Unsupported model method:", method))
  }

  # Sort importance and get top N features
  importance <- sort(importance, decreasing = TRUE)
  top_features <- head(importance, top_n)

  # Plot Feature Importance
  library(ggplot2)
  p <- ggplot(data.frame(Feature = names(top_features), Importance = top_features),
              aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(
      title = title,
      x = "Feature",
      y = "Importance"
    ) +
    theme_minimal()

  print(p)

  return(top_features)
}
