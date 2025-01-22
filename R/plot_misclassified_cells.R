#' Plotting Misclassified labels
#'
#' @param steam_obj STEAM Object
#' @param coordinates coordinates
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_minimal theme element_blank
#' @importFrom scales hue_pal
#' @export
plot_misclassified_cells <- function(steam_obj, coordinates = NULL) {
  # Predefined 7-color palette
  default_colors <- c(
    "#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7"
  )

  # assign colors for misclassified and base layers
  cols <- c("Misclassified" = "black",
            "Layer_1" = default_colors[1],
            "Layer_2" = default_colors[2],
            "Layer_3" = default_colors[3],
            "Layer_4" = default_colors[4],
            "Layer_5" = default_colors[5],
            "Layer_6" = default_colors[6],
            "WM"      = default_colors[7])

  # broader palette if more than 7 unique layers exist
  broader_palette <- hue_pal()(30)

  if (is.null(coordinates)) {
    coordinates <- steam_obj$spatial
    colnames(coordinates) <- c("cell_id", "col", "row")
  } else {
    colnames(coordinates) <- c("cell_id", "col", "row")
  }

  # full data for plotting
  full_data <- data.frame(
    Row = -coordinates$row,
    Col = coordinates$col,
    Labels = as.character(steam_obj$labels)  # True labels
  )

  # test coordinates, labels, and predictions
  test_coords <- steam_obj$test$test.data.coords
  test_labels <- as.character(steam_obj$test$test.data.labels)
  predictions <- as.character(steam_obj$test$predictions)

  # Identify misclassified cells
  misclassified <- test_labels != predictions
  test_indices <- rownames(test_coords)
  matching_indices <- match(test_indices, rownames(coordinates))
  full_data$Labels[matching_indices[misclassified]] <- "Misclassified"

  unique_labels <- unique(full_data$Labels)
  if (length(unique_labels) > 8) {

    missing_labels <- setdiff(unique_labels, names(cols))
    extra_colors <- setNames(broader_palette[1:length(missing_labels)], missing_labels)
    cols <- c(cols, extra_colors)
  } else {

    missing_labels <- setdiff(unique_labels, names(cols))
    for (lbl in missing_labels) {
      cols[lbl] <- broader_palette[length(cols) + 1]
    }
  }

  ggplot(full_data, aes(x = Col, y = Row, color = Labels)) +
    geom_point(size = 3) +
    scale_color_manual(values = cols) +
    labs(
      title = "",
      color = "Labels"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
}
