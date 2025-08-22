#' Feature Importance
#'
#' @param STEAM.obj STEAM Object
#' @param top_n Top n
#' @param title title
#' @param nested_fold If train$model is missing (nested CV), which outer fold's model to use (integer or "best")
#' @export
feature_importance <- function(STEAM.obj, top_n = 10, title = "Top Features by Importance", nested_fold = "best") {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(caret)
  })
  
  if (!is.list(STEAM.obj)) stop("Invalid STEAM.obj.")
  if (!is.numeric(top_n) || top_n <= 0) stop("'top_n' must be a positive numeric value.")
  
  mod <- STEAM.obj$train$model
  # Fallback: nested CV model if simple-model absent
  if (is.null(mod)) {
    ncv <- try(STEAM.obj$nested$ncv, silent = TRUE)
    if (!inherits(ncv, "try-error") && !is.null(ncv) && length(ncv$outer_result) > 0) {
      if (identical(nested_fold, "best")) {
        # choose fold with max Kappa (fallback to Accuracy)
        scores <- vapply(ncv$outer_result, function(of) {
          cm_try <- try(caret::confusionMatrix(of$preds$predy, of$preds$testy), silent = TRUE)
          if (inherits(cm_try, "try-error")) return(NA_real_)
          as.numeric(cm_try$overall["Kappa"])
        }, numeric(1))
        if (all(is.na(scores))) {
          scores <- vapply(ncv$outer_result, function(of) {
            cm <- table(of$preds$predy, of$preds$testy); sum(diag(cm))/sum(cm)
          }, numeric(1))
        }
        best_i <- which.max(scores)
        mod <- ncv$outer_result[[best_i]]$fit
      } else if (is.numeric(nested_fold) && nested_fold >= 1 && nested_fold <= length(ncv$outer_result)) {
        mod <- ncv$outer_result[[nested_fold]]$fit
      }
    }
  }
  if (is.null(mod)) stop("No model found (train$model is NULL and no suitable nested model).")
  
  method <- mod$method
  vi <- NULL
  
  vi_try <- try(caret::varImp(mod, scale = FALSE)$importance, silent = TRUE)
  if (!inherits(vi_try, "try-error") && is.data.frame(vi_try) && nrow(vi_try) > 0) {
    # If multiple columns (e.g., per-class), reduce to mean absolute importance
    if (!"Overall" %in% colnames(vi_try)) {
      vi <- rowMeans(as.matrix(vi_try), na.rm = TRUE)
      names(vi) <- rownames(vi_try)
    } else {
      vi <- vi_try$Overall
      names(vi) <- rownames(vi_try)
    }
  }
  
  if (is.null(vi)) {
    if (method == "svmLinear") {
      # kernlab SVM (binary linear): compute abs weights
      svm_model <- mod$finalModel
      coefs <- try(colSums(svm_model@coef[[1]] * svm_model@xmatrix[[1]]), silent = TRUE)
      if (!inherits(coefs, "try-error")) {
        vi <- abs(coefs)
        names(vi) <- colnames(svm_model@xmatrix[[1]])
      }
    } else if (method == "multinom") {
      coef_matrix <- coef(mod$finalModel)
      vi <- apply(abs(coef_matrix), 2, sum)
      vi <- vi[!grepl("(Intercept)", names(vi), fixed = TRUE)]
    }
  }
  
  if (is.null(vi) || length(vi) == 0) {
    stop(sprintf("Variable importance not available for method '%s'. Consider permutation importance.", method))
  }
  
  vi <- sort(vi, decreasing = TRUE)
  top_features <- head(vi, top_n)
  
  p <- ggplot(data.frame(Feature = names(top_features), Importance = as.numeric(top_features)),
              aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black") +
    coord_flip() +
    labs(title = title, x = "Feature", y = "Importance") +
    theme_classic()
  
  print(p)
  return(top_features)
}
