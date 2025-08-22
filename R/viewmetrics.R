ViewMetrics <- function(STEAM.obj, fold = NULL, view = c("both","overall","folds"), return_data = FALSE) {
  suppressPackageStartupMessages({
    library(ggplot2); library(reshape2); library(grid); library(gridExtra); library(viridis); library(caret)
  })
  view <- match.arg(view)
  
  # ---------- helpers ----------
  vir_cols <- viridis::viridis(100)
  .text_col <- function(v, vmax) {
    idx <- max(1, min(100, round((v / vmax) * 99) + 1))
    rgbv <- grDevices::col2rgb(vir_cols[idx])
    br <- sum(rgbv * c(0.299, 0.587, 0.114))
    if (br > 128) "black" else "white"
  }
  .scalar_num <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA_real_)
    y <- suppressWarnings(as.numeric(x[1]))
    if (is.na(y)) NA_real_ else y
  }
  .align_levels <- function(pred, obs) {
    pr <- as.character(pred); tr <- as.character(obs)
    keep <- !(is.na(pr) | is.na(tr))
    pr <- pr[keep]; tr <- tr[keep]
    lev <- sort(unique(c(pr, tr)))
    list(pred = factor(pr, levels = lev), obs = factor(tr, levels = lev))
  }
  
  # ---------- mode detection ----------
  is_simple <- !is.null(STEAM.obj$test) && !is.null(STEAM.obj$test$metrics)
  is_nested <- !is.null(STEAM.obj$nested)
  if (!is_simple && !is_nested)
    stop("Metrics not found. Run model.predict(mode='simple' or 'nested') first.")
  
  # ====================== SIMPLE (overall only + confusion) ======================
  if (is_simple) {
    get_metric <- function(name) if (!is.null(STEAM.obj$test$metrics[[name]])) STEAM.obj$test$metrics[[name]] else NA
    metrics_data <- data.frame(
      Test_Accuracy = get_metric("accuracy"),
      Mean_Balanced_Accuracy = get_metric("mean_balanced_accuracy"),
      ARI = get_metric("ARI"),
      Kappa = get_metric("kappa"),
      NMI = get_metric("NMI"),
      AMI = get_metric("AMI"),
      PAS = get_metric("PAS")
    )
    prec <- get_metric("precision"); rec <- get_metric("recall"); f1 <- get_metric("f1_score")
    metrics_data2 <- data.frame(
      Class = names(prec),
      Precision = as.numeric(prec),
      Recall = as.numeric(rec),
      F1_Score = as.numeric(f1)
    )
    conf_matrix <- as.matrix(get_metric("confusion_matrix"))
    conf_df <- reshape2::melt(conf_matrix)
    colnames(conf_df) <- c("Predicted_Label", "Actual_Label", "Value")
    vmax <- max(conf_df$Value); conf_df$TextColor <- vapply(conf_df$Value, .text_col, character(1), vmax)
    
    p_conf <- ggplot(conf_df, aes(x = Actual_Label, y = Predicted_Label, fill = Value)) +
      geom_tile() + geom_text(aes(label = Value, color = TextColor), size = 6) +
      scale_fill_viridis_c() + scale_color_identity() +
      labs(title = "Confusion Matrix", x = "Actual", y = "Predicted") +
      theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    if (return_data) return(list(overall = metrics_data, per_fold = NULL))
    
    grid::grid.newpage()
    gridExtra::grid.arrange(
      gridExtra::tableGrob(t(metrics_data)),
      gridExtra::tableGrob(metrics_data2, rows = NULL),
      ggplot2::ggplotGrob(p_conf),
      layout_matrix = rbind(c(1,2),
                            c(3,3)),
      heights = c(1.6, 3)
    )
    return(invisible(NULL))
  }
  
  # ====================== NESTED: overall + per-fold ======================
  md  <- STEAM.obj$nested$metrics
  ncv <- STEAM.obj$nested$ncv
  if (is.null(md)) md <- list()
  if (is.null(ncv) || is.null(ncv$outer_result) || length(ncv$outer_result) == 0)
    stop("NestedCV object present but outer_result missing. Did model.predict(mode='nested') run?")
  
  need_fallback <- any(sapply(c("Test_Accuracy","Mean_Balanced_Accuracy","Kappa"),
                              function(k) is.null(md[[k]]) || length(md[[k]]) == 0))
  if (need_fallback) {
    all_preds <- tryCatch({
      if (!is.null(ncv$output)) ncv$output else ncv$outer_result[[1]]$preds[0,]
    }, error = function(e) NULL)
    if (is.null(all_preds) && !is.null(ncv$outer_result)) {
      all_preds <- do.call(rbind, lapply(ncv$outer_result, function(of) of$preds))
    }
    if (!is.null(all_preds) && nrow(all_preds) > 0) {
      al <- .align_levels(all_preds$predy, all_preds$testy)
      cm <- caret::confusionMatrix(al$pred, al$obs)
      byc <- cm$byClass
      if (is.null(dim(byc))) {
        ba <- mean(c(byc["Sensitivity"], byc["Specificity"]), na.rm = TRUE)
      } else {
        sens <- suppressWarnings(as.numeric(byc[,"Sensitivity"]))
        spec <- suppressWarnings(as.numeric(byc[,"Specificity"]))
        ba <- mean((sens + spec)/2, na.rm = TRUE)
      }
      md$Test_Accuracy <- unname(cm$overall["Accuracy"])
      md$Kappa         <- unname(cm$overall["Kappa"])
      md$Mean_Balanced_Accuracy <- ba
    }
  }
  
  # ---- overall (pooled) table ----
  metrics_data <- data.frame(
    Test_Accuracy  = round(.scalar_num(md$Test_Accuracy), 3),
    Mean_Balanced_Accuracy = round(.scalar_num(md$Mean_Balanced_Accuracy), 3),
    ARI   = round(.scalar_num(md$ARI), 3),
    Kappa = round(.scalar_num(md$Kappa), 3),
    NMI   = round(.scalar_num(md$NMI), 3),
    AMI   = round(.scalar_num(md$AMI), 3),
    PAS   = .scalar_num(md$PAS)
  )
  has_prf <- !is.null(md$precision) && !is.null(md$recall) && !is.null(md$f1_score)
  if (has_prf) {
    prf_df <- data.frame(
      Class    = names(md$precision),
      Precision= as.numeric(md$precision),
      Recall   = as.numeric(md$recall),
      F1_Score = as.numeric(md$f1_score),
      row.names = NULL
    )
  }
  
  # ---- per-fold table (square-safe) ----
  per_fold_list <- lapply(seq_along(ncv$outer_result), function(i) {
    of <- ncv$outer_result[[i]]
    al <- .align_levels(of$preds$predy, of$preds$testy)
    cm <- table(Predicted = al$pred, Actual = al$obs)
    acc_i <- sum(diag(cm)) / sum(cm)
    kap_i <- tryCatch(unname(caret::confusionMatrix(al$pred, al$obs)$overall["Kappa"]), error = function(e) NA_real_)
    bt <- tryCatch({
      bt_df <- of$fit$bestTune
      paste(names(bt_df), "=", as.character(unlist(bt_df)), collapse = ", ")
    }, error = function(e) NA_character_)
    list(fold = i, cm = cm, acc = acc_i, kappa = as.numeric(kap_i), bestTune = bt)
  })
  per_fold_df <- do.call(rbind, lapply(per_fold_list, function(z)
    data.frame(Fold = z$fold,
               Accuracy = round(z$acc, 3),
               Kappa    = round(z$kappa, 3),
               BestTune = z$bestTune,
               stringsAsFactors = FALSE)))
  
  
  chosen_fold <- if (is.null(fold)) {
    if (all(is.na(per_fold_df$Kappa))) per_fold_df$Fold[which.max(per_fold_df$Accuracy)]
    else per_fold_df$Fold[which.max(per_fold_df$Kappa)]
  } else fold
  if (!chosen_fold %in% per_fold_df$Fold)
    stop("Requested fold not found in nestedCV outer_result.")
  
  cm <- per_fold_list[[chosen_fold]]$cm
  conf_df <- reshape2::melt(as.matrix(cm))
  colnames(conf_df) <- c("Predicted_Label", "Actual_Label", "Value")
  vmax <- max(conf_df$Value)
  conf_df$TextColor <- vapply(conf_df$Value, .text_col, character(1), vmax)
  
  p_conf <- ggplot(conf_df, aes(x = Actual_Label, y = Predicted_Label, fill = Value)) +
    geom_tile() + geom_text(aes(label = Value, color = TextColor), size = 6) +
    scale_fill_viridis_c() + scale_color_identity() +
    labs(title = paste0("Confusion Matrix (Outer fold ", chosen_fold, ")"),
         x = "Actual", y = "Predicted") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (return_data) return(list(overall = metrics_data, per_fold = per_fold_df))
  
  grid::grid.newpage()
  if (view == "both") {
    if (has_prf) {
      gridExtra::grid.arrange(
        gridExtra::tableGrob(t(metrics_data)),
        gridExtra::tableGrob(per_fold_df, rows = NULL),
        gridExtra::tableGrob(prf_df, rows = NULL),
        ggplot2::ggplotGrob(p_conf),
        layout_matrix = rbind(c(1, 2),
                              c(3, 4)),
        heights = c(1.2, 3)
      )
    } else {
      gridExtra::grid.arrange(
        gridExtra::tableGrob(t(metrics_data)),
        gridExtra::tableGrob(per_fold_df, rows = NULL),
        ggplot2::ggplotGrob(p_conf),
        layout_matrix = rbind(c(1, 2),
                              c(3, 3)),
        heights = c(1.6, 3)
      )
    }
  } else if (view == "overall") {
    grobs <- list(gridExtra::tableGrob(t(metrics_data)))
    if (has_prf) grobs <- c(grobs, list(gridExtra::tableGrob(prf_df, rows = NULL)))
    grobs <- c(grobs, list(ggplot2::ggplotGrob(p_conf)))
    gridExtra::grid.arrange(grobs = grobs, ncol = 2)
  } else { # view == "folds"
    gridExtra::grid.arrange(
      gridExtra::tableGrob(per_fold_df, rows = NULL),
      ggplot2::ggplotGrob(p_conf),
      ncol = 1, heights = c(1.2, 2.5)
    )
  }
  
  invisible(NULL)
}
