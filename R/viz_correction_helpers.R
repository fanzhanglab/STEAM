#' STEAM Visualization Helper Functions
#'
#' Helper functions used internally by STEAM visualization routines. These
#' functions provide reusable components for extracting data, applying
#' corrections, building color palettes, and generating standardized plots.
#'
#' @keywords internal

#' Get Misclassified Cells
#'
#' Extracts misclassified cell IDs from nested CV predictions in a STEAM object.
#'
#' @param STEAM.obj A STEAM object with nested CV results.
#' @param fold Optional numeric fold index. If \code{NULL}, returns across all folds.
#'
#' @return A character vector of misclassified cell IDs.
#' @keywords internal
getMisclassifiedCells <- function(STEAM.obj, fold = NULL) {
  fold_indices <- if (!is.null(fold)) fold else seq_along(STEAM.obj$nested$ncv$outer_result)

  all_preds <- do.call(rbind, lapply(fold_indices, function(i) {
    if (i > length(STEAM.obj$nested$ncv$outer_result)) return(NULL)
    p <- STEAM.obj$nested$ncv$outer_result[[i]]$preds
    if (is.null(p)) return(NULL)
    p$outer_fold <- i
    p
  }))

  if (is.null(all_preds) || nrow(all_preds) == 0) {
    stop("No nested CV predictions found for the specified fold(s).")
  }

  if (!"cell_id" %in% colnames(all_preds)) {
    if (!is.null(rownames(all_preds))) {
      all_preds$cell_id <- rownames(all_preds)
    } else {
      stop("Predictions don't have rownames or 'cell_id' to match spatial coordinates.")
    }
  }

  misclassified_mask <- all_preds$testy != all_preds$predy
  if (!any(misclassified_mask)) return(character(0))

  misclassified_cells <- all_preds$cell_id[misclassified_mask]

  if (!is.null(fold)) {
    fold_mask <- all_preds$outer_fold[misclassified_mask] == fold
    misclassified_cells <- misclassified_cells[fold_mask]
  }

  return(misclassified_cells)
}

#' Get Corrected Cells
#'
#' Extracts cell IDs that were corrected by STEAMCorrection.
#'
#' @param STEAM.obj A STEAM object containing correction results.
#' @param fold Optional numeric fold index. If \code{NULL}, returns across all folds.
#'
#' @return A character vector of corrected cell IDs.
#' @keywords internal
getCorrectedCells <- function(STEAM.obj, fold = NULL) {
  correction_results <- STEAM.obj$spatial_anchor_analysis$correction_results
  if (is.null(correction_results)) return(character(0))

  if (!is.null(fold)) {
    fold_key <- paste0("fold_", fold)
    if (!fold_key %in% names(correction_results)) return(character(0))

    fold_data <- correction_results[[fold_key]]
    if (is.null(fold_data$correction_history) || length(fold_data$correction_history) <= 1) {
      return(character(0))
    }

    initial_preds <- fold_data$correction_history[[1]]$predictions_snapshot
    final_preds <- tail(fold_data$correction_history, 1)[[1]]$predictions_snapshot

    return(names(initial_preds)[initial_preds != final_preds])
  } else {
    all_corrected <- character(0)
    for (fold_name in names(correction_results)) {
      fold_data <- correction_results[[fold_name]]
      if (!is.null(fold_data$correction_history) && length(fold_data$correction_history) > 1) {
        initial_preds <- fold_data$correction_history[[1]]$predictions_snapshot
        final_preds <- tail(fold_data$correction_history, 1)[[1]]$predictions_snapshot
        fold_corrected <- names(initial_preds)[initial_preds != final_preds]
        all_corrected <- c(all_corrected, fold_corrected)
      }
    }
    return(unique(all_corrected))
  }
}

#' Create Color Palette
#'
#' Generates a color palette for spatial visualization, ensuring
#' misclassified cells are colored black.
#'
#' @param labels A vector of label names.
#' @param master_colors A named vector of colors for tissue layers.
#'
#' @return A named vector of colors for the given labels.
#' @keywords internal
createColorPalette <- function(labels, master_colors) {
  unique_labels <- unique(labels)
  colors <- master_colors[names(master_colors) %in% unique_labels]

  if ("Misclassified" %in% unique_labels) colors["Misclassified"] <- "black"
  if ("Still Misclassified" %in% unique_labels) colors["Still Misclassified"] <- "black"

  missing_labels <- setdiff(unique_labels, names(colors))
  if (length(missing_labels) > 0) {
    missing_colors <- setNames(scales::hue_pal()(length(missing_labels)), missing_labels)
    colors <- c(colors, missing_colors)
  }

  return(colors)
}

#' Create Spatial Plot
#'
#' Generates a standardized ggplot2 scatter plot of spatial cells.
#'
#' @param data A data.frame with \code{Col}, \code{Row}, and \code{Labels}.
#' @param colors A named vector of colors for labels.
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_classic theme
#' @keywords internal
createSpatialPlot <- function(data, colors, title) {
  ggplot(data, aes(x = Col, y = Row, color = Labels)) +
    geom_point(size = 3) +
    scale_color_manual(values = colors) +
    labs(title = title, color = "Layer / Status") +
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

#' Create Post-Correction Data
#'
#' Updates a plotting dataset with corrected labels and marks cells
#' that remain misclassified.
#'
#' @param STEAM.obj A STEAM object containing correction results.
#' @param base_data Base plotting data.frame (with \code{cell_id} and \code{Labels}).
#' @param fold Fold index (numeric).
#' @param misclassified_cells Character vector of misclassified cells.
#' @param corrected_cells Character vector of corrected cells.
#'
#' @return A data.frame with updated \code{Labels}.
#' @keywords internal
createAfterData <- function(STEAM.obj, base_data, fold, misclassified_cells, corrected_cells) {
  after_data <- base_data

  correction_results <- STEAM.obj$spatial_anchor_analysis$correction_results
  if (!is.null(correction_results) && !is.null(fold)) {
    fold_key <- paste0("fold_", fold)
    if (fold_key %in% names(correction_results)) {
      fold_data <- correction_results[[fold_key]]
      if (!is.null(fold_data$corrected_labels)) {
        corrected_labels <- fold_data$corrected_labels
        fold_cells <- intersect(names(corrected_labels), after_data$cell_id)
        if (length(fold_cells) > 0) {
          cell_idx <- match(fold_cells, after_data$cell_id)
          after_data$Labels[cell_idx] <- as.character(corrected_labels[fold_cells])
        }
      }
    }
  }

  if (length(misclassified_cells) > 0) {
    uncorrected_cells <- setdiff(misclassified_cells, corrected_cells)
    if (length(uncorrected_cells) > 0) {
      cell_idx <- match(uncorrected_cells, after_data$cell_id)
      cell_idx <- cell_idx[!is.na(cell_idx)]
      after_data$Labels[cell_idx] <- "Still Misclassified"
    }
  }

  return(after_data)
}
