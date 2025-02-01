#' Plotting Misclassified labels
#'
#' @param steam_obj STEAM Object
#' @param colors list of colors for clusters
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_minimal theme element_blank
#' @importFrom scales hue_pal
#' @export
plot_misclassified_cells <- function(steam_obj, colors = NULL) {

  assign_layer_colors <- function(labels, colors = NULL) {
    unique_labels <- unique(labels)
    num_layers <- length(unique_labels)

    if (!is.null(colors)) {
      layer_colors <- colors
    } else {
      default_colors <- rainbow(num_layers)
      layer_colors <- setNames(default_colors, unique_labels)
    }

    layer_colors["Misclassified"] <- "black"

    return(layer_colors)
  }

  cols <- assign_layer_colors(steam_obj$labels, colors)
  coordinates <- steam_obj$spatial
  if (is.null(coordinates)) {
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
      color = "Layer/Label"
    ) +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),   # Removes axis lines
      panel.grid = element_blank(),  # Removes grid
      panel.border = element_blank() # Removes panel border
    )
}
