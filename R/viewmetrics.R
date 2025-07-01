ViewMetrics <- function(STEAM.obj) {
  if (is.null(STEAM.obj$test$metrics)) stop("Metrics data is missing in STEAM.obj")
  
  get_metric <- function(name) {
    if (!is.null(STEAM.obj$test$metrics[[name]])) {
      return(STEAM.obj$test$metrics[[name]])
    } else {
      return(NA)
    }
  }
  
  # Pull training accuracy if available
  train_acc <- if (!is.null(STEAM.obj$train$model$results$Accuracy)) round(STEAM.obj$train$model$results$Accuracy, 3) else NA
  
  metrics_data <- data.frame(
    Train_Accuracy = train_acc,
    Test_Accuracy = get_metric("accuracy"),
    Mean_Balanced_Accuracy = get_metric("mean_balanced_accuracy"),
    ARI = get_metric("ARI"),
    Kappa = get_metric("kappa"),
    NMI = get_metric("NMI"),
    AMI = get_metric("AMI"),
    PAS = get_metric("PAS")
  )
  
  metrics_data2 <- data.frame(
    Layer = names(get_metric("precision")),
    Precision = get_metric("precision"),
    Recall = get_metric("recall"),
    F1_Score = get_metric("f1_score")
  )
  
  # Confusion Matrix Plot
  conf_matrix <- as.matrix(get_metric("confusion_matrix"))
  conf_df <- reshape2::melt(conf_matrix)
  colnames(conf_df) <- c("Predicted_Label", "Actual_Label", "Value")
  
  viridis_colors <- viridis::viridis_pal()(100)
  get_text_color <- function(value, max_value) {
    norm_value <- round((value / max_value) * 99) + 1
    brightness <- sum(grDevices::col2rgb(viridis_colors[norm_value]) * c(0.299, 0.587, 0.114))
    if (brightness > 128) "black" else "white"
  }
  max_value <- max(conf_df$Value)
  conf_df$TextColor <- sapply(conf_df$Value, get_text_color, max_value = max_value)
  
  p_conf <- ggplot2::ggplot(conf_df, ggplot2::aes(x = Actual_Label, y = Predicted_Label, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = Value, color = TextColor), size = 6) +
    ggplot2::scale_fill_viridis_c(option = "viridis") +
    ggplot2::scale_color_identity() +
    ggplot2::labs(title = "Confusion Matrix", x = "Actual Label", y = "Predicted Label") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14),
      text = ggplot2::element_text(size = 14)
    )
  
  # Final grid without silhouette plot
  grid::grid.newpage()
  gridExtra::grid.arrange(
    gridExtra::tableGrob(t(metrics_data)),
    gridExtra::tableGrob(metrics_data2, rows = NULL),
    ggplot2::ggplotGrob(p_conf),
    layout_matrix = rbind(
      c(1, 2),
      c(3, 3)
    ),
    heights = c(2, 3)
  )
}

