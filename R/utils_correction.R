#' STEAM Shared Utility Functions
#'
#' This file contains shared utility functions used across STEAM analysis components.
#' These functions handle common operations like coordinate extraction, validation, and error handling.
#'
#' @name utils_correction
#' @keywords internal
NULL
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 package is required. Install with: install.packages('ggplot2')")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("scales package is required. Install with: install.packages('scales')")
}



#' Extract and Validate Spatial Coordinates
#'
#' Extracts spatial coordinates from STEAM object and validates column structure
#'
#' @param STEAM.obj STEAM object containing spatial coordinates
#' @param require_cell_id Whether to require cell_id column (default: TRUE)
#' @return data.frame with standardized coordinate columns (cell_id, Col, Row)
#' @export
extractSpatialCoordinates <- function(STEAM.obj, require_cell_id = TRUE) {
  # Validate STEAM object has spatial coordinates
  if (is.null(STEAM.obj$spatial)) {
    stop("No spatial coordinates found in STEAM object")
  }

  coordinates <- as.data.frame(STEAM.obj$spatial)

  # Detect coordinate columns using flexible naming
  cl <- tolower(colnames(coordinates))
  x_name <- colnames(coordinates)[which(cl %in% c("x", "col", "column"))[1]]
  y_name <- colnames(coordinates)[which(cl %in% c("y", "row"))[1]]

  if (is.na(x_name) || is.na(y_name)) {
    stop("Failed to detect spatial coordinate columns. Expected columns named 'x'/'Col'/'Column' and 'y'/'Row'.")
  }

  # Ensure cell_id column exists if required
  if (require_cell_id && !"cell_id" %in% colnames(coordinates)) {
    if (!is.null(rownames(coordinates)) && !all(rownames(coordinates) == as.character(seq_len(nrow(coordinates))))) {
      coordinates$cell_id <- rownames(coordinates)
    } else {
      coordinates$cell_id <- as.character(seq_len(nrow(coordinates)))
    }
  }

  # Standardize column names
  result <- data.frame(
    cell_id = if (require_cell_id) coordinates$cell_id else rownames(coordinates),
    Col = coordinates[[x_name]],
    Row = coordinates[[y_name]],
    stringsAsFactors = FALSE
  )

  if (!require_cell_id) {
    result$cell_id <- NULL
  }

  return(result)
}

#' Validate Fold Parameter
#'
#' Validates fold parameter against available folds in STEAM object
#'
#' @param STEAM.obj STEAM object containing fold data
#' @param fold Fold number to validate
#' @return TRUE if valid, stops execution with error message if invalid
#' @export
validateFoldParameter <- function(STEAM.obj, fold) {
  if (is.null(fold)) {
    return(TRUE)  # NULL fold is always valid (means use all folds)
  }

  if (!is.numeric(fold) || length(fold) != 1) {
    stop("fold must be a single numeric value")
  }

  # Check nested CV structure
  if (!is.null(STEAM.obj$nested) && !is.null(STEAM.obj$nested$ncv$outer_result)) {
    if (fold > length(STEAM.obj$nested$ncv$outer_result)) {
      stop(sprintf("fold %d not found. Available folds: 1-%d", fold, length(STEAM.obj$nested$ncv$outer_result)))
    }
    return(TRUE)
  }

  # Check corrections data for available folds
  if (!is.null(STEAM.obj$spatial_anchor_analysis$corrections)) {
    available_folds <- unique(STEAM.obj$spatial_anchor_analysis$corrections$fold)
    if (!fold %in% available_folds) {
      stop(sprintf("fold %d not found. Available folds: %s", fold, paste(sort(available_folds), collapse = ", ")))
    }
    return(TRUE)
  }

  # Check fold_summaries
  if (!is.null(STEAM.obj$spatial_anchor_analysis$fold_summaries)) {
    n_folds <- length(STEAM.obj$spatial_anchor_analysis$fold_summaries)
    if (fold > n_folds) {
      stop(sprintf("fold %d not found. Available folds: 1-%d", fold, n_folds))
    }
    return(TRUE)
  }

  stop("No fold information found in STEAM object")
}

#' Create Master Color Palette
#'
#' Creates a consistent color palette for tissue layers across all visualizations
#'
#' @param labels Vector of all possible label names
#' @param include_misclassified Whether to include black color for "Misclassified" (default: TRUE)
#' @return Named vector of colors
#' @export
createMasterColorPalette <- function(labels, include_misclassified = TRUE) {
  # Get unique tissue layer names, excluding special labels
  tissue_layers <- unique(as.character(labels))
  tissue_layers <- tissue_layers[!tissue_layers %in% c("Misclassified", "Still Misclassified")]
  tissue_layers <- sort(tissue_layers)  # Sort for consistent order

  # Create master color palette
  master_colors <- setNames(scales::hue_pal()(length(tissue_layers)), tissue_layers)

  # Add black for misclassified cells if requested
  if (include_misclassified) {
    master_colors["Misclassified"] <- "black"
    master_colors["Still Misclassified"] <- "black"
  }

  return(master_colors)
}

#' Align Labels with Coordinates
#'
#' Aligns label vector with coordinate data frame, handling size mismatches
#'
#' @param labels Vector of labels
#' @param coordinates Data frame with coordinates
#' @return List with aligned_labels and aligned_coordinates
#' @export
alignLabelsWithCoordinates <- function(labels, coordinates) {
  if (length(labels) >= nrow(coordinates)) {
    # If we have more labels than coordinates, truncate labels to match coordinates
    aligned_labels <- labels[seq_len(nrow(coordinates))]
    aligned_coordinates <- coordinates
  } else {
    # If we have fewer labels than coordinates, truncate coordinates to match labels
    aligned_coordinates <- coordinates[seq_along(labels), , drop = FALSE]
    aligned_labels <- labels
  }

  return(list(
    aligned_labels = aligned_labels,
    aligned_coordinates = aligned_coordinates
  ))
}

#' Get Misclassified Cells from Nested CV
#'
#' Extracts misclassified cells from nested cross-validation results
#'
#' @param STEAM.obj STEAM object with nested CV results
#' @param fold Optional specific fold to filter by
#' @return Character vector of misclassified cell IDs
#' @export
getMisclassifiedCells <- function(STEAM.obj, fold = NULL) {
  if (is.null(STEAM.obj$nested) || is.null(STEAM.obj$nested$ncv$outer_result)) {
    return(character(0))
  }

  # Filter by fold if specified
  fold_indices <- if (!is.null(fold)) fold else seq_along(STEAM.obj$nested$ncv$outer_result)

  all_preds <- do.call(rbind, lapply(fold_indices, function(i) {
    if (i > length(STEAM.obj$nested$ncv$outer_result)) return(NULL)
    p <- STEAM.obj$nested$ncv$outer_result[[i]]$preds
    if (is.null(p)) return(NULL)
    p$outer_fold <- i
    p
  }))

  if (is.null(all_preds) || nrow(all_preds) == 0) {
    return(character(0))
  }

  if (!"cell_id" %in% colnames(all_preds)) {
    if (!is.null(rownames(all_preds))) {
      all_preds$cell_id <- rownames(all_preds)
    } else {
      all_preds$cell_id <- as.character(seq_len(nrow(all_preds)))
    }
  }

  # Identify misclassified cells
  mis <- all_preds$testy != all_preds$predy
  if (!any(mis)) {
    return(character(0))
  }

  mis_ids <- all_preds$cell_id[mis]
  mis_folds <- all_preds$outer_fold[mis]

  if (!is.null(fold)) {
    # Filter to only show misclassified cells from the current fold
    current_fold_mask <- mis_folds == fold
    return(mis_ids[current_fold_mask])
  }

  return(mis_ids)
}

#' Create Standard Plot Theme
#'
#' Creates a consistent theme for spatial plots
#'
#' @return ggplot2 theme object
#' @export
createSpatialPlotTheme <- function() {
  theme_void() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
}

#' Safe Require Package
#'
#' Safely requires a package with informative error message
#'
#' @param package_name Name of the package
#' @param install_message Custom installation message
#' @return TRUE if package is available
#' @export
safeRequirePackage <- function(package_name, install_message = NULL) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    if (is.null(install_message)) {
      install_message <- sprintf("Install with: install.packages('%s')", package_name)
    }
    stop(sprintf("%s package is required. %s", package_name, install_message))
  }
  return(TRUE)
}

#' Load Patchwork Package Safely
#'
#' Centralized function to load patchwork for combined plots
#'
#' @return TRUE if loaded successfully
#' @export
loadPatchwork <- function() {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package is required for combined plots. Install with: install.packages('patchwork')")
  }
  suppressPackageStartupMessages(library(patchwork))
  return(TRUE)
}

#' Extract Ground Truth from STEAM Object
#'
#' Centralized function to extract ground truth from various STEAM object locations
#'
#' @param STEAM.obj STEAM object
#' @param fold Optional fold number to get fold-specific ground truth
#' @return Vector of ground truth labels or NULL if not found
#' @export
extractGroundTruth <- function(STEAM.obj, fold = NULL) {
  # Check direct ground truth field
  if (!is.null(STEAM.obj$ground_truth)) {
    return(STEAM.obj$ground_truth)
  }

  # Check in spatial_anchor_analysis
  if (!is.null(STEAM.obj$spatial_anchor_analysis$iterative_results$ground_truth)) {
    return(STEAM.obj$spatial_anchor_analysis$iterative_results$ground_truth)
  }

  # Check nested CV for fold-specific ground truth
  if (!is.null(fold) && !is.null(STEAM.obj$nested$ncv$outer_result)) {
    if (fold <= length(STEAM.obj$nested$ncv$outer_result)) {
      fold_result <- STEAM.obj$nested$ncv$outer_result[[fold]]
      if (!is.null(fold_result$preds$testy)) {
        return(setNames(fold_result$preds$testy, rownames(fold_result$preds)))
      }
    }
  }

  return(NULL)
}
