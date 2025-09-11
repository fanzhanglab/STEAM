#' Create Before vs After Comparison (Internal)
#'
#' Internal helper to generate side-by-side spatial plots showing misclassified
#' cells before correction and updated labels after \code{\link{STEAMCorrection}}.
#'
#' @param STEAM.obj A STEAM object containing spatial coordinates and correction results.
#' @param fold Optional fold number to visualize. If NULL, combines all folds.
#'
#' @return A patchwork object combining two ggplot2 plots.
#' @keywords internal
createBeforeAfterComparison <- function(STEAM.obj, fold = NULL) {

  # Validate inputs and load dependencies
  loadPatchwork()
  if (is.null(STEAM.obj$nested$ncv$outer_result)) {
    stop("No nested CV predictions found. Run model.predict with nested CV first.")
  }

  # Get spatial coordinates and align with labels
  coordinates <- extractSpatialCoordinates(STEAM.obj)
  alignment_result <- alignLabelsWithCoordinates(STEAM.obj$labels, coordinates)

  # Create base data frame
  base_data <- data.frame(
    cell_id = alignment_result$aligned_coordinates$cell_id,
    Col = alignment_result$aligned_coordinates$Col,
    Row = -alignment_result$aligned_coordinates$Row,
    Labels = as.character(alignment_result$aligned_labels),
    stringsAsFactors = FALSE,
    row.names = alignment_result$aligned_coordinates$cell_id
  )

  # Extract misclassified cells for specified fold(s)
  misclassified_cells <- getMisclassifiedCells(STEAM.obj, fold)

  # Create BEFORE plot data
  before_data <- base_data
  if (length(misclassified_cells) > 0) {
    before_data$Labels[match(misclassified_cells, before_data$cell_id)] <- "Misclassified"
  }

  # Get corrected cells
  corrected_cells <- getCorrectedCells(STEAM.obj, fold)

  # Create color palette
  master_colors <- createMasterColorPalette(STEAM.obj$labels)
  before_colors <- createColorPalette(before_data$Labels, master_colors)

  # Create BEFORE plot
  before_title <- sprintf("BEFORE: Misclassified Cells (%s)",
                          if (!is.null(fold)) paste("Fold", fold) else "all folds")

  before_plot <- createSpatialPlot(before_data, before_colors, before_title)

  # Create AFTER plot data
  after_data <- createAfterData(STEAM.obj, base_data, fold, misclassified_cells, corrected_cells)
  after_colors <- createColorPalette(after_data$Labels, master_colors)

  after_title <- sprintf("AFTER: Corrections Applied (%s)",
                         if (!is.null(fold)) paste("Fold", fold) else "All Folds")

  after_plot <- createSpatialPlot(after_data, after_colors, after_title)

  # Combine plots
  main_title <- sprintf("Before vs After: %s Analysis Impact",
                        if (!is.null(fold)) paste("Fold", fold) else "Spatial Anchor (All Folds)")

  subtitle_text <- sprintf("BEFORE: %d misclassified | AFTER: %d corrections (black = remaining misclassified)",
                           length(misclassified_cells), length(corrected_cells))

  comparison_plot <- before_plot + after_plot +
    patchwork::plot_annotation(
      title = main_title,
      subtitle = subtitle_text,
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray60")
      )
    )

  return(comparison_plot)
}

plotBeforeAfter <- function(STEAM.obj, fold = NULL) {
  return(createBeforeAfterComparison(STEAM.obj, fold = fold))
}

#' ===== VISUALIZATION FUNCTIONS =====



#' Iterative Correction Plot
#'
#' Shows misclassification corrections across iterations for a specific fold.
#'
#' @param STEAM.obj A STEAM object with correction results.
#' @param fold Numeric; specific fold to visualize (required).
#' @param iterations Optional numeric vector of iteration numbers to plot.
#'
#' @return A patchwork or ggplot2 object combining iteration plots.
#' @export
IterativePlot <- function(STEAM.obj, fold = NULL, iterations = NULL) {

  # Validate fold parameter is provided
  if (is.null(fold)) {
    stop("fold parameter is required for IterativePlot. Specify which fold to visualize.")
  }

  validateFoldParameter(STEAM.obj, fold)

  coordinates <- extractSpatialCoordinates(STEAM.obj)

  if (is.null(STEAM.obj$nested$ncv$outer_result[[fold]])) {
    stop(sprintf("Fold %d data not found in STEAM.obj$nested$ncv$outer_result", fold))
  }

  fold_data <- STEAM.obj$nested$ncv$outer_result[[fold]]

  if (is.null(fold_data$preds)) {
    stop(sprintf("No prediction data found for fold %d", fold))
  }

  fold_preds <- fold_data$preds
  if (!"cell_id" %in% colnames(fold_preds)) {
    fold_preds$cell_id <- rownames(fold_preds)
  }

  ground_truth <- setNames(as.character(fold_preds$testy), fold_preds$cell_id)
  original_predictions <- setNames(as.character(fold_preds$predy), fold_preds$cell_id)

  misclassified_cells <- names(original_predictions)[original_predictions != ground_truth]

  correction_data <- NULL

  # Get correction data from current correction_results structure
  if (!is.null(STEAM.obj$spatial_anchor_analysis$correction_results)) {
    fold_key <- paste0("fold_", fold)
    if (fold_key %in% names(STEAM.obj$spatial_anchor_analysis$correction_results)) {
      correction_data <- STEAM.obj$spatial_anchor_analysis$correction_results[[fold_key]]
    }
  }

  if (is.null(correction_data)) {
    stop(sprintf("No correction data found for fold %d. Run STEAMCorrection() first.", fold))
  }

  iteration_history <- correction_data$correction_history
  if (is.null(iteration_history) || length(iteration_history) == 0) {
    stop(sprintf("No iteration history found for fold %d", fold))
  }

  # Filter iterations if specified
  if (!is.null(iterations)) {
    if (!is.numeric(iterations)) {
      stop("iterations parameter must be numeric")
    }

    # Get available iteration numbers from the data
    available_iterations <- sapply(iteration_history, function(x) x$iteration)

    # Validate requested iterations exist
    invalid_iterations <- iterations[!iterations %in% available_iterations]
    if (length(invalid_iterations) > 0) {
      stop(paste("Invalid iteration numbers:", paste(invalid_iterations, collapse=", "),
                 "- available iterations are:", paste(available_iterations, collapse=", ")))
    }

    # Filter iteration history to only include requested iterations
    iteration_history <- iteration_history[sapply(iteration_history, function(x) x$iteration %in% iterations)]

    cat(paste("Showing", length(iterations), "selected iterations:", paste(iterations, collapse=", "), "\n"))
  } else {
    cat(paste("Showing all", length(iteration_history), "iterations\n"))
  }

  master_colors <- createMasterColorPalette(STEAM.obj$labels)

  iteration_plots <- list()

  for (i in seq_along(iteration_history)) {
    iter_data <- iteration_history[[i]]
    iteration_num <- iter_data$iteration

    current_preds <- iter_data$predictions_snapshot

    alignment_result <- alignLabelsWithCoordinates(STEAM.obj$labels, coordinates)
    labels_aligned <- alignment_result$aligned_labels
    coordinates_aligned <- alignment_result$aligned_coordinates

    plot_data <- data.frame(
      cell_id = coordinates_aligned$cell_id,
      Col = -coordinates_aligned$Row,
      Row = -coordinates_aligned$Col,
      Labels = as.character(labels_aligned),
      stringsAsFactors = FALSE
    )

    fold_cells <- intersect(names(current_preds), plot_data$cell_id)
    if (length(fold_cells) > 0) {
      fold_idx <- which(plot_data$cell_id %in% fold_cells)
      plot_data$Labels[fold_idx] <- as.character(current_preds[plot_data$cell_id[fold_idx]])
    }

    if (iteration_num > 0) {
      current_misclassified <- names(current_preds)[current_preds != ground_truth[names(current_preds)]]
      if (length(current_misclassified) > 0) {
        mis_idx <- which(plot_data$cell_id %in% current_misclassified)
        plot_data$Labels[mis_idx] <- "Still Misclassified"
      }
    } else {
      if (length(misclassified_cells) > 0) {
        mis_idx <- which(plot_data$cell_id %in% misclassified_cells)
        plot_data$Labels[mis_idx] <- "Misclassified"
      }
    }

    plot_labels <- unique(plot_data$Labels)
    plot_colors <- master_colors[names(master_colors) %in% plot_labels]

    if ("Misclassified" %in% plot_labels) {
      plot_colors["Misclassified"] <- "black"
    }
    if ("Still Misclassified" %in% plot_labels) {
      plot_colors["Still Misclassified"] <- "black"
    }

    if (iteration_num == 0) {
      plot_title <- sprintf("Iteration %d: Original\n(%d spots misclassified)",
                            iteration_num, iter_data$remaining_misclassified)
    } else {
      plot_title <- sprintf("Iteration %d: %d corrections\n(%d spots still misclassified)",
                            iteration_num, iter_data$corrections_applied,
                            iter_data$remaining_misclassified)
    }

    # Calculate dynamic sizing based on total number of iterations
    n_iterations <- length(iteration_history)

    # Dynamic point size: larger for fewer plots
    point_size <- if (n_iterations <= 3) 1.5 else if (n_iterations <= 8) 1.0 else if (n_iterations <= 15) 0.8 else 0.6

    # Dynamic title size: larger for fewer plots
    title_size <- if (n_iterations <= 3) 12 else if (n_iterations <= 8) 10 else if (n_iterations <= 15) 8 else 7

    # Dynamic margin: more space for fewer plots
    margin_size <- if (n_iterations <= 8) 8 else 5

    p <- ggplot(plot_data, aes(x = Col, y = Row, color = Labels)) +
      geom_point(size = point_size, alpha = 0.8) +
      scale_color_manual(values = plot_colors) +
      labs(title = plot_title) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = title_size),
        legend.position = "none",
        panel.border = element_rect(color = "gray80", fill = NA, size = 0.2),
        plot.margin = margin(margin_size, margin_size, margin_size, margin_size, unit = "pt")
      )

    iteration_plots[[paste0("iter_", iteration_num)]] <- p
  }

  # Calculate optimal grid layout based on number of iterations
  n_plots <- length(iteration_plots)
  if (n_plots <= 4) {
    # For 1-4 plots: single row
    grid_ncol <- n_plots
    grid_nrow <- 1
  } else if (n_plots <= 8) {
    # For 5-8 plots: 2 rows
    grid_ncol <- ceiling(n_plots / 2)
    grid_nrow <- 2
  } else if (n_plots <= 15) {
    # For 9-15 plots: 3 rows
    grid_ncol <- ceiling(n_plots / 3)
    grid_nrow <- 3
  } else if (n_plots <= 24) {
    # For 16-24 plots: 4 rows
    grid_ncol <- ceiling(n_plots / 4)
    grid_nrow <- 4
  } else {
    # For many plots: aim for roughly square layout
    grid_nrow <- ceiling(sqrt(n_plots))
    grid_ncol <- ceiling(n_plots / grid_nrow)
  }

  combined <- wrap_plots(iteration_plots, nrow = grid_nrow, ncol = grid_ncol)

  final_accuracy <- tail(iteration_history, 1)[[1]]$accuracy
  initial_accuracy <- iteration_history[[1]]$accuracy
  total_corrections <- tail(iteration_history, 1)[[1]]$total_corrections

  combined <- combined + plot_annotation(
    title = sprintf("STEAM Iterative Correction of Misclassifications"),
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray60")
    )
  )

  legend_data <- data.frame(
    x = 1, y = 1,
    Labels = names(master_colors)[names(master_colors) %in% unique(unlist(lapply(iteration_plots, function(p) unique(p$data$Labels))))]
  )

  legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = Labels)) +
    geom_point() +
    scale_color_manual(values = master_colors, name = "Cell Type") +
    theme_void() +
    theme(legend.position = "bottom")

  legend <- cowplot::get_legend(legend_plot)

  if (requireNamespace("cowplot", quietly = TRUE)) {
    final_plot <- cowplot::plot_grid(combined, legend, ncol = 1, rel_heights = c(1, 0.1))
    return(final_plot)
  } else {
    message("Install cowplot for better legend positioning: install.packages('cowplot')")
    return(combined)
  }
}



#' Accuracy Before vs After Correction
#'
#' Produces a paired bar chart showing accuracy before and after
#' correction for each fold.
#'
#' @param STEAM.obj A STEAM object with correction results.
#' @param show_improvements Logical; annotate improvements (default: TRUE).
#' @param show_corrections Logical; annotate correction counts (default: TRUE).
#'
#' @return A ggplot2 object.
#' @export
AccuracyPairedBarChart <- function(STEAM.obj, show_improvements = TRUE, show_corrections = TRUE) {

  correction_results <- STEAM.obj$spatial_anchor_analysis$correction_results

  if (is.null(correction_results)) {
    stop("No correction results found. Run STEAMCorrection() first.")
  }

  fold_data_list <- list()

  for (fold_name in names(correction_results)) {
    fold_data <- correction_results[[fold_name]]

    fold_data_list[[fold_name]] <- data.frame(
      fold_num = fold_data$fold_number,
      initial_acc = fold_data$initial_accuracy * 100,
      final_acc = fold_data$final_accuracy * 100,
      corrections = fold_data$total_corrections,
      stringsAsFactors = FALSE
    )
  }

  combined_data <- do.call(rbind, fold_data_list)
  combined_data$improvement <- combined_data$final_acc - combined_data$initial_acc

  before_data <- data.frame(
    fold = paste0("Fold ", combined_data$fold_num),
    accuracy = combined_data$initial_acc,
    timepoint = "Before",
    corrections = combined_data$corrections,
    improvement = combined_data$improvement,
    stringsAsFactors = FALSE
  )

  after_data <- data.frame(
    fold = paste0("Fold ", combined_data$fold_num),
    accuracy = combined_data$final_acc,
    timepoint = "After",
    corrections = combined_data$corrections,
    improvement = combined_data$improvement,
    stringsAsFactors = FALSE
  )

  plot_data <- rbind(before_data, after_data)
  plot_data$timepoint <- factor(plot_data$timepoint, levels = c("Before", "After"))
  plot_data$fold <- factor(plot_data$fold, levels = paste0("Fold ", sort(combined_data$fold_num)))

  if (length(unique(plot_data$timepoint)) != 2) {
    warning("Timepoint factor issues detected")
  }

  color_values <- c("Before" = "#E74C3C", "After" = "#27AE60")

  p <- ggplot(plot_data, aes(x = fold, y = accuracy, fill = timepoint)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.8, alpha = 0.8) +
    scale_fill_manual(
      values = color_values,
      name = "",
      breaks = c("Before", "After"),
      labels = c("Before", "After")
    ) +
    labs(
      title = "Accuracy Before vs After Correction",
      x = "Cross-Validation Fold",
      y = "Accuracy (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 35, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 32, color = "gray60"),
      axis.title.x = element_text(size = 28, face = "bold"),
      axis.title.y = element_text(size = 28, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 26),
      axis.text.y = element_text(size = 26),
      legend.text = element_text(size = 22),
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  if (show_improvements) {
    improvement_data <- combined_data
    improvement_data$fold_name <- paste0("Fold ", improvement_data$fold_num)
    improvement_data$y_pos <- improvement_data$final_acc + 1

    p <- p + geom_text(
      data = improvement_data,
      aes(x = fold_name, y = y_pos, label = sprintf("+%.1f%%", improvement)),
      inherit.aes = FALSE,
      color = "#2C3E50",
      fontface = "bold",
      size = 10,
      vjust = 0
    )
  }

  if (show_corrections) {
    correction_data <- combined_data
    correction_data$fold_name <- paste0("Fold ", correction_data$fold_num)
    correction_data$y_pos <- (correction_data$initial_acc + correction_data$final_acc) / 2

    p <- p + geom_text(
      data = correction_data,
      aes(x = fold_name, y = y_pos, label = paste0(corrections, "\ncorr.")),
      inherit.aes = FALSE,
      color = "black",
      fontface = "bold",
      size = 10,
      lineheight = 0.8
    )
  }

  mean_improvement <- mean(combined_data$improvement)
  total_corrections <- sum(combined_data$corrections)
  n_folds <- nrow(combined_data)

  p <- p + labs(
    subtitle = sprintf("Mean improvement: +%.2f%% | %d total corrections | %d folds",
                       mean_improvement, total_corrections, n_folds)
  )

  if (show_improvements) {
    y_max <- max(combined_data$final_acc) + 3
    y_min <- min(combined_data$initial_acc) - 1
    p <- p + coord_cartesian(ylim = c(y_min, y_max))
  }

  return(p)
}


#' Neighborhood Homogeneity Analysis
#'
#' Visualizes the relationship between local tissue homogeneity
#' and correction success rate.
#'
#' @param STEAM.obj A STEAM object with correction results.
#' @param k Number of neighbors for homogeneity (default: auto-detected).
#'
#' @return A ggplot2 object showing correction success rate by homogeneity bins.
#' @export
NeighborhoodHomogeneityAnalysis <- function(STEAM.obj, k = NULL) {

  safeRequirePackage("RANN", "Install with: install.packages('RANN')")

  if (is.null(k)) {
    if (!is.null(STEAM.obj$spatial_anchor_analysis$parameters$k)) {
      k <- STEAM.obj$spatial_anchor_analysis$parameters$k
      message(sprintf("Auto-detected k=%d from stored parameters", k))
    } else if (!is.null(STEAM.obj$correction_parameters$k)) {
      k <- STEAM.obj$correction_parameters$k
      message(sprintf("Auto-detected k=%d from correction parameters", k))
    }
  }

  coordinates <- extractSpatialCoordinates(STEAM.obj)
  coords_matrix <- as.matrix(coordinates[, c("Col", "Row")])
  rownames(coords_matrix) <- coordinates$cell_id

  knn_result <- nn2(coords_matrix, k = k + 1) # +1 to include self

  correction_results <- STEAM.obj$spatial_anchor_analysis$correction_results
  if (is.null(correction_results)) {
    stop("No correction results found. Run STEAMCorrection() first.")
  }

  homogeneity_data <- data.frame()

  for (fold_name in names(correction_results)) {
    fold_data <- correction_results[[fold_name]]
    fold_num <- fold_data$fold_number

    cv_result <- STEAM.obj$nested$ncv$outer_result[[fold_num]]
    if (is.null(cv_result$preds)) next

    ground_truth <- setNames(as.character(cv_result$preds$testy), rownames(cv_result$preds))

    correction_history <- fold_data$correction_history
    if (length(correction_history) < 2) next

    initial_preds <- correction_history[[1]]$predictions_snapshot
    final_preds <- tail(correction_history, 1)[[1]]$predictions_snapshot

    corrected_cells <- names(initial_preds)[initial_preds != final_preds]

    for (cell_id in corrected_cells) {
      if (!cell_id %in% rownames(coords_matrix)) next

      cell_idx <- which(rownames(coords_matrix) == cell_id)
      if (length(cell_idx) == 0) next

      neighbor_indices <- knn_result$nn.idx[cell_idx, -1]  # Remove first neighbor (self)
      neighbor_cell_ids <- rownames(coords_matrix)[neighbor_indices]

      neighbor_true_labels <- ground_truth[neighbor_cell_ids]
      neighbor_true_labels <- neighbor_true_labels[!is.na(neighbor_true_labels)]

      if (length(neighbor_true_labels) > 0) {
        # Homogeneity = proportion of neighbors with same true label as this cell
        cell_true_label <- ground_truth[cell_id]
        homogeneity <- sum(neighbor_true_labels == cell_true_label) / length(neighbor_true_labels)

        was_originally_wrong <- initial_preds[cell_id] != ground_truth[cell_id]
        correction_successful <- final_preds[cell_id] == ground_truth[cell_id] && was_originally_wrong

        homogeneity_data <- rbind(homogeneity_data, data.frame(
          fold = fold_num,
          cell_id = cell_id,
          homogeneity = homogeneity,
          correction_successful = correction_successful,
          original_pred = initial_preds[cell_id],
          corrected_pred = final_preds[cell_id],
          true_label = cell_true_label,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  if (nrow(homogeneity_data) == 0) {
    stop("No homogeneity data available for analysis")
  }

  homogeneity_data$homogeneity_bin <- cut(
    homogeneity_data$homogeneity,
    breaks = seq(0, 1, 0.2),
    labels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%"),
    include.lowest = TRUE
  )

  bin_summary <- homogeneity_data %>%
    group_by(homogeneity_bin) %>%
    summarise(
      total_corrections = n(),
      successful_corrections = sum(correction_successful),
      success_rate = successful_corrections / total_corrections * 100,
      .groups = 'drop'
    )

  p <- ggplot(bin_summary, aes(x = homogeneity_bin, y = success_rate)) +
    geom_col(aes(fill = success_rate), alpha = 0.8, width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", success_rate, total_corrections)),
              vjust = -0.1, size = 12, fontface = "bold") +
    scale_fill_gradient2(
      low = "#E74C3C", mid = "#F39C12", high = "#27AE60",
      midpoint = 50, name = "Success\nRate (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 35, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 32, color = "gray60"),
      axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 20)),
      axis.title.y = element_text(size = 30, face = "bold"),
      axis.text.x = element_text(size = 28),
      axis.text.y = element_text(size = 28),
      legend.title = element_text(size = 28, face = "bold", margin = margin(b = 15)),
      legend.text = element_text(size = 22),
      legend.key.width = unit(2, "cm"),
      legend.key.height = unit(1.5, "cm"),
      legend.title.align = 0.5,
      legend.spacing = unit(0.5, "cm"),
      legend.margin = margin(t = 10, b = 30, l = 15, r = 15),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = "Correction Success vs Neighborhood Homogeneity",
      subtitle = sprintf("Analysis of %d corrections across %d folds (k=%d neighbors)",
                         nrow(homogeneity_data), length(unique(homogeneity_data$fold)), k),
      x = "Neighborhood Homogeneity (% neighbors with same true tissue type)",
      y = "Correction Success Rate (%)"
    ) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.1)))

  return(p)
}


#' Error Analysis of Corrections
#'
#' Heatmap-style plot showing which tissue types were corrected successfully
#' vs unsuccessfully, with per-tissue success rates.
#'
#' @param STEAM.obj A STEAM object with correction results.
#'
#' @return A patchwork or list of ggplot2 objects (if patchwork not available).
#' @export
ErrorAnalysisPlot <- function(STEAM.obj) {

  correction_results <- STEAM.obj$spatial_anchor_analysis$correction_results
  if (is.null(correction_results)) {
    stop("No correction results found. Run STEAMCorrection() first.")
  }

  correction_outcomes <- data.frame()

  for (fold_name in names(correction_results)) {
    fold_data <- correction_results[[fold_name]]
    fold_num <- fold_data$fold_number

    cv_result <- STEAM.obj$nested$ncv$outer_result[[fold_num]]
    if (is.null(cv_result$preds)) next

    ground_truth <- setNames(as.character(cv_result$preds$testy), rownames(cv_result$preds))

    correction_history <- fold_data$correction_history
    if (length(correction_history) < 2) next

    initial_preds <- correction_history[[1]]$predictions_snapshot
    final_preds <- tail(correction_history, 1)[[1]]$predictions_snapshot

    corrected_cells <- names(initial_preds)[initial_preds != final_preds]

    if (length(corrected_cells) > 0) {
      for (cell_id in corrected_cells) {
        original_pred <- initial_preds[cell_id]
        corrected_pred <- final_preds[cell_id]
        true_label <- ground_truth[cell_id]

        was_originally_wrong <- original_pred != true_label
        correction_fixed_it <- corrected_pred == true_label
        correction_successful <- was_originally_wrong && correction_fixed_it

        correction_outcomes <- rbind(correction_outcomes, data.frame(
          fold = fold_num,
          cell_id = cell_id,
          true_tissue = true_label,
          original_prediction = original_pred,
          corrected_prediction = corrected_pred,
          correction_successful = correction_successful,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  if (nrow(correction_outcomes) == 0) {
    stop("No correction outcomes to analyze")
  }

  confusion_data <- correction_outcomes %>%
    group_by(true_tissue, corrected_prediction) %>%
    summarise(
      count = n(),
      .groups = 'drop'
    )

  all_tissues <- sort(unique(c(correction_outcomes$true_tissue, correction_outcomes$corrected_prediction)))

  complete_confusion <- expand.grid(
    true_tissue = all_tissues,
    corrected_prediction = all_tissues,
    stringsAsFactors = FALSE
  ) %>%
    left_join(confusion_data, by = c("true_tissue", "corrected_prediction")) %>%
    mutate(
      count = ifelse(is.na(count), 0, count)
    )

  p1 <- ggplot(complete_confusion, aes(x = corrected_prediction, y = factor(true_tissue, levels = rev(all_tissues)))) +
    geom_tile(aes(fill = count), color = "white", size = 0.5) +
    geom_text(aes(label = ifelse(count > 0, count, "")),
              color = "white", fontface = "bold", size = 8) +
    scale_fill_gradient(low = "#F8F9FA", high = "#2C3E50", name = "Number of\nCorrections") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 35, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 32, color = "gray60"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 28),
      axis.text.y = element_text(size = 28),
      axis.title.x = element_text(size = 30),
      axis.title = element_text(size = 30),
      legend.title = element_text(size = 30, face = "bold", margin = margin(b = 20)),
      legend.text = element_text(size = 28),
      legend.key.width = unit(2, "cm"),
      legend.key.height = unit(1.5, "cm"),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = "Correction Confusion Matrix",
      subtitle = "Rows: True tissue type | Columns: Corrected prediction",
      x = "Corrected Prediction",
      y = "True Tissue Type (Actual)"
    )

  tissue_success <- correction_outcomes %>%
    group_by(true_tissue) %>%
    summarise(
      total_corrections = n(),
      successful_corrections = sum(correction_successful),
      success_rate = successful_corrections / total_corrections * 100,
      .groups = 'drop'
    ) %>%
    arrange(desc(success_rate))

  p2 <- ggplot(tissue_success, aes(x = reorder(true_tissue, success_rate), y = success_rate)) +
    geom_col(aes(fill = success_rate), alpha = 0.8, width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", success_rate, total_corrections)),
              hjust = -0.1, size = 8, fontface = "bold") +
    scale_fill_gradient2(
      low = "#E74C3C", mid = "#F39C12", high = "#27AE60",
      midpoint = 50, name = "Success\nRate (%)"
    ) +
    coord_flip() +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 28),
      legend.title = element_text(size = 30, face = "bold", margin = margin(b = 15)),
      legend.text = element_text(size = 28),
      legend.key.width = unit(2, "cm"),
      legend.key.height = unit(1.5, "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    labs(
      title = "Success Rate by True Tissue Type",
      x = "True Tissue Type",
      y = "Correction Success Rate (%)"
    ) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.1)))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- p1 / p2 + plot_annotation(
      title = "Correction Error Analysis",
      subtitle = sprintf("Analysis of %d corrections across %d tissue types",
                         nrow(correction_outcomes), length(all_tissues)),
      theme = theme(
        plot.title = element_text(hjust = 0, size = 32, face = "bold"),
        plot.subtitle = element_text(hjust = 0, size = 35, color = "gray60")
      )
    )
    return(combined)
  } else {
    return(list(confusion_matrix = p1, success_by_tissue = p2))
  }
}
