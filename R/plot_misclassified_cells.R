#' Plotting Misclassified labels
#'
#' @param label_colors list of colors for plotting labels
#' @param steam_obj STEAM Object
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_minimal theme element_blank
#' @importFrom scales hue_pal
#' @export
plot_misclassified_cells <- function(steam_obj, label_colors = NULL) {

  # Load spatial coordinates
  coordinates <- steam_obj$spatial
  if (is.null(coordinates)) stop("No spatial coordinates found in STEAM object")

  coordinates <- as.data.frame(coordinates)

  # Handle different naming conventions
  colnames_lower <- tolower(colnames(coordinates))
  x_name <- colnames(coordinates)[which(colnames_lower %in% c("x", "col", "column"))[1]]
  y_name <- colnames(coordinates)[which(colnames_lower %in% c("y", "row"))[1]]

  if (is.na(x_name) || is.na(y_name)) stop("Failed to detect spatial coordinate columns.")

  # Add cell ID if not present
  if (!"cell_id" %in% colnames(coordinates)) {
    coordinates$cell_id <- rownames(coordinates)
  }

  # Reorder columns
  coordinates <- coordinates[, c("cell_id", x_name, y_name)]
  colnames(coordinates) <- c("cell_id", "Col", "Row")

  # Prepare plotting data
  full_data <- data.frame(
    Row = -coordinates$Row,
    Col = coordinates$Col,
    Labels = as.character(steam_obj$labels),
    row.names = coordinates$cell_id
  )

  # Extract test data
  test_coords <- steam_obj$test$test.data.coords
  test_labels <- as.character(steam_obj$test$test.data.labels)
  predictions <- as.character(steam_obj$test$predictions)

  # Identify misclassified cells
  misclassified <- test_labels != predictions
  test_indices <- rownames(test_coords)
  matching_indices <- match(test_indices, rownames(full_data))
  full_data$Labels[matching_indices[misclassified]] <- "Misclassified"

  # Generate color palette
  unique_labels <- unique(full_data$Labels)
  base_labels <- setdiff(unique_labels, "Misclassified")

  if (is.null(label_colors)) {
    base_colors <- scales::hue_pal()(length(base_labels))
    label_colors <- setNames(base_colors, base_labels)
    if ("Misclassified" %in% unique_labels) {
      label_colors <- c(label_colors, Misclassified = "black")
    }
  }

  # Plot
  ggplot(full_data, aes(x = Col, y = Row, color = Labels)) +
    geom_point(size = 1) +
    scale_color_manual(values = label_colors) +
    labs(title = "", color = "Layer/Label") +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
}
