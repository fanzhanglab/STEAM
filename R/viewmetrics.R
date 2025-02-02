#' View Metrics visualizations
#'
#' @param STEAM.obj STEAM Object
#'
#' @importFrom ggplot2 ggplot ggplotGrob
#' @importFrom gridExtra grid.arrange tableGrob
#' @importFrom grid grid.draw grid.newpage
#' @import reshape2
#' @import viridis
#'
#' @export
ViewMetrics <- function(STEAM.obj) {

  # Check if required metrics are available
  if (is.null(STEAM.obj$test$metrics)) stop("Metrics data is missing in STEAM.obj")

  # Create metrics data frames
  metrics_data <- data.frame(
    Accuracy = STEAM.obj$test$metrics$accuracy,
    "Mean Balanced Accuracy" = STEAM.obj$test$metrics$mean_balanced_accuracy,
    ARI = STEAM.obj$test$metrics$ARI,
    Kappa = STEAM.obj$test$metrics$kappa
  )

  metrics_data2 <- data.frame(
    Layer = names(STEAM.obj$test$metric$precision),
    Precision = STEAM.obj$test$metric$precision,
    Recall = STEAM.obj$test$metric$recall,
    F1_Score = STEAM.obj$test$metric$f1_score
  )

  # Create confusion matrix
  conf_matrix <- as.matrix(STEAM.obj$test$metrics$confusion_matrix)
  conf_df <- melt(conf_matrix)
  colnames(conf_df) <- c("Predicted_Label", "Actual_Label", "Value")

  # Compute adaptive text color for tiles
  viridis_colors <- viridis_pal()(100)
  get_text_color <- function(value, max_value) {
    norm_value <- round((value / max_value) * 99) + 1
    brightness <- sum(col2rgb(viridis_colors[norm_value]) * c(0.299, 0.587, 0.114))
    if (brightness > 128) "black" else "white"
  }
  max_value <- max(conf_df$Value)
  conf_df$TextColor <- sapply(conf_df$Value, get_text_color, max_value = max_value)

  # Create confusion matrix plot
  p <- ggplot(conf_df, aes(x = Actual_Label, y = Predicted_Label, fill = Value)) +
    geom_tile() +
    geom_text(aes(label = Value, color = TextColor), size = 6) +
    scale_fill_viridis(option = "viridis") +
    scale_color_identity() +
    labs(
      title = "Confusion Matrix",
      x = "Actual Label",
      y = "Predicted Label"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      text = element_text(size = 14)
    )

  # Arrange tables and plot
  grid.newpage()
  grid.arrange(
    tableGrob(t(metrics_data)),
    tableGrob(metrics_data2, rows = NULL),
    ggplotGrob(p),
    layout_matrix = rbind(
      c(1, 2),  # First row: Two tables side by side
      c(3, 3)   # Second row: Confusion matrix spanning both columns
    ),
    heights = c(2, 3)  # Adjust heights (2 for tables, 3 for plot)
  )
}
