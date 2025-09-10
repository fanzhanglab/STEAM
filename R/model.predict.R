#' Internal: Predict and Evaluate STEAM Models
#'
#' This function evaluates trained models from a STEAM object, either using a
#' simple train/test split or nested cross-validation. Not intended for direct use.
#'
#' @param STEAM.obj A STEAM object with trained models or nested CV results.
#' @param mode Either "simple" (train/test split) or "nested" (nested CV results).
#' @return The input STEAM object with added predictions and evaluation metrics.
#' @keywords internal
model.predict <- function(STEAM.obj, mode = c("simple","nested")) {
  suppressPackageStartupMessages({
    library(caret)
  })
  mode <- match.arg(mode)

  # ---------------- SIMPLE (unchanged) ----------------
  if (mode == "simple") {
    test_matrix <- STEAM.obj$test$avg.matrix
    test_labels <- STEAM.obj$test$test.data.labels
    model <- STEAM.obj$train$model

    df <- data.frame(labels = test_labels, t(test_matrix))
    df$labels <- factor(df$labels)
    for (i in 2:ncol(df)) df[, i] <- as.numeric(df[, i])

    p <- predict(model, df)
    p <- factor(p, levels = levels(df$labels))
    true_labels <- as.character(test_labels)
    pred_labels <- as.character(p)

    cm <- caret::confusionMatrix(p, df$labels)
    confusion_matrix <- cm$table
    accuracy <- unname(cm$overall["Accuracy"])
    kappa    <- unname(cm$overall["Kappa"])

    precision <- diag(confusion_matrix) / rowSums(confusion_matrix)
    recall    <- diag(confusion_matrix) / colSums(confusion_matrix)
    f1_score  <- 2 * (precision * recall) / (precision + recall)
    f1_score[is.nan(f1_score)] <- 0
    specificity <- sapply(1:nrow(confusion_matrix), function(i) {
      tn <- sum(confusion_matrix[-i, -i]); fp <- sum(confusion_matrix[i, -i])
      tn / (tn + fp)
    })
    balanced_accuracy <- (recall + specificity) / 2
    mean_balanced_accuracy <- mean(balanced_accuracy, na.rm = TRUE)

    ari <- tryCatch(mclust::adjustedRandIndex(true_labels, pred_labels), error = function(e) NA_real_)
    nmi <- tryCatch(aricode::NMI(true_labels, pred_labels), error = function(e) NA_real_)
    ami <- tryCatch(aricode::AMI(true_labels, pred_labels), error = function(e) NA_real_)

    compute_PAS <- function(coords, labels, k = 10, threshold = 6) {
      if (!requireNamespace("FNN", quietly = TRUE)) return(NA_real_)
      coords <- as.matrix(coords); labels <- as.character(labels)
      knn_result <- FNN::get.knn(coords, k = k)
      neighbors_idx <- knn_result$nn.index
      abnormal_flags <- vapply(seq_len(nrow(coords)), function(i) {
        sum(labels[neighbors_idx[i, ]] != labels[i]) >= threshold
      }, logical(1))
      round(mean(abnormal_flags), 3)
    }
    pas <- tryCatch(compute_PAS(STEAM.obj$test$test.data.coords, pred_labels), error = function(e) NA_real_)

    STEAM.obj$test$predictions <- p
    STEAM.obj$test$true.labels <- df$labels
    STEAM.obj$test$metrics <- list(
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
      PAS = pas
    )
    return(STEAM.obj)
  }

  # ---------------- NESTED CV ----------------
  ncv <- if (!is.null(STEAM.obj$train$ncv)) STEAM.obj$train$ncv else {
    if (!is.null(STEAM.obj$nested)) STEAM.obj$nested else stop("No nestedCV results found.")
  }

  # Gather all outer-fold predictions
  all_preds <- tryCatch({
    if (!is.null(ncv$output)) ncv$output else ncv$outer_result[[1]]$preds[0,]
  }, error = function(e) NULL)
  if (is.null(all_preds) && !is.null(ncv$outer_result)) {
    all_preds <- do.call(rbind, lapply(ncv$outer_result, function(of) of$preds))
  }
  if (is.null(all_preds) || nrow(all_preds) == 0)
    stop("NestedCV predictions not found.")

  obs  <- factor(all_preds$testy)
  pred <- factor(all_preds$predy, levels = levels(obs))

  cm <- caret::confusionMatrix(pred, obs)
  confusion_matrix <- cm$table
  accuracy <- unname(cm$overall["Accuracy"])
  kappa    <- unname(cm$overall["Kappa"])
  precision <- diag(confusion_matrix) / rowSums(confusion_matrix)
  recall    <- diag(confusion_matrix) / colSums(confusion_matrix)
  f1_score  <- 2 * (precision * recall) / (precision + recall)
  f1_score[is.nan(f1_score)] <- 0

  byc <- cm$byClass
  if (is.null(dim(byc))) {
    ba <- mean(c(byc["Sensitivity"], byc["Specificity"]), na.rm = TRUE)
  } else {
    sens <- suppressWarnings(as.numeric(byc[,"Sensitivity"]))
    spec <- suppressWarnings(as.numeric(byc[,"Specificity"]))
    ba <- mean((sens + spec)/2, na.rm = TRUE)
  }

  true_labels <- as.character(obs)
  pred_labels <- as.character(pred)
  ari <- tryCatch(mclust::adjustedRandIndex(true_labels, pred_labels), error = function(e) NA_real_)
  nmi <- tryCatch(aricode::NMI(true_labels, pred_labels), error = function(e) NA_real_)
  ami <- tryCatch(aricode::AMI(true_labels, pred_labels), error = function(e) NA_real_)

  pas <- NA_real_
  if (!is.null(STEAM.obj$spatial) && !is.null(rownames(STEAM.obj$spatial)) && !is.null(rownames(all_preds))) {
    if (requireNamespace("FNN", quietly = TRUE)) {
      idx <- intersect(rownames(STEAM.obj$spatial), rownames(all_preds))
      if (length(idx) > 0) {
        coords <- as.matrix(STEAM.obj$spatial[idx, , drop = FALSE])
        labs   <- pred_labels[match(idx, rownames(all_preds))]
        kn <- FNN::get.knn(coords, k = 10)$nn.index
        flags <- vapply(seq_len(nrow(coords)), function(i) {
          sum(labs[kn[i,]] != labs[i]) >= 6
        }, logical(1))
        pas <- round(mean(flags), 3)
      }
    }
  }

  STEAM.obj$nested <- list(
    metrics = list(
      Test_Accuracy = round(as.numeric(accuracy), 3),
      Mean_Balanced_Accuracy = round(as.numeric(ba), 3),
      Kappa = round(as.numeric(kappa), 3),
      ARI = round(as.numeric(ari), 3),
      NMI = round(as.numeric(nmi), 3),
      AMI = round(as.numeric(ami), 3),
      PAS = pas,
      precision = precision,
      recall = recall,
      f1_score = f1_score
    ),
    confusion = confusion_matrix,
    ncv = ncv
  )

  STEAM.obj$nested$per_fold <- lapply(seq_along(ncv$outer_result), function(i) {
    of <- ncv$outer_result[[i]]
    cm_i <- table(Predicted = of$preds$predy, Actual = of$preds$testy)
    acc_i <- sum(diag(cm_i))/sum(cm_i)
    kap_i <- tryCatch(unname(caret::confusionMatrix(cm_i)$overall["Kappa"]), error = function(e) NA_real_)
    bt <- tryCatch({
      bt_df <- of$fit$bestTune
      paste(names(bt_df), "=", as.character(unlist(bt_df)), collapse = ", ")
    }, error = function(e) NA_character_)

    prec_i <- diag(cm_i) / rowSums(cm_i)
    rec_i  <- diag(cm_i) / colSums(cm_i)
    f1_i   <- 2 * (prec_i * rec_i) / (prec_i + rec_i)
    f1_i[is.nan(f1_i)] <- 0

    list(
      fold = i,
      Accuracy = round(acc_i,3),
      Kappa = round(as.numeric(kap_i),3),
      BestTune = bt,
      confusion = cm_i,
      precision = prec_i,
      recall = rec_i,
      f1_score = f1_i
    )
  })

  return(STEAM.obj)
}
