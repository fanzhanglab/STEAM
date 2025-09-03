# Enhanced STEAM visualization tools for spatial analysis
# Supports both single-pass and iterative anchor analysis results

# Load required libraries
if (!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)

# Helper operator for null coalescing
`%||%` <- function(lhs, rhs) {
  if (!is.null(lhs) && length(lhs) > 0) lhs else rhs
}

#' Anchor summary panel (flips + Δ accuracy)
#'
#' Creates a two-panel summary: (A) flips per outer fold (correct vs wrong),
#' and (B) Δ accuracy per fold.
#'
#' @param STEAM.obj A STEAM object after [STEAM_anchor()].
#' @return A patchwork object (printable).
#' @examples
#' \dontrun{
#'   STEAM.obj <- STEAM_anchor(STEAM.obj)
#'   p <- AnchorPanel(STEAM.obj)
#'   print(p)
#' }
#' @importFrom ggplot2 ggplot aes geom_col geom_text labs theme_classic theme element_text
#' @importFrom ggplot2 position_stack scale_fill_manual geom_hline scale_y_continuous
#' @importFrom tidyr pivot_longer
#' @importFrom patchwork plot_layout
#' @importFrom scales percent_format
#' @export
AnchorPanel <- function(STEAM.obj) {
  # Basic checks
  if (!is.list(STEAM.obj) || is.null(STEAM.obj$spatial_anchor_analysis)) {
    stop("AnchorPanel: run STEAM_anchor() first; no $spatial_anchor_analysis found.", call. = FALSE)
  }
  saa <- STEAM.obj$spatial_anchor_analysis
  
  # Compute flip_summary if absent but we have final_by_fold
  flip_summary <- saa$flip_summary
  if ((is.null(flip_summary) || !nrow(flip_summary)) &&
      !is.null(STEAM.obj$nested$ncv) && !is.null(saa$final_by_fold)) {
    if (!exists("FlipStats", mode = "function")) {
      stop("AnchorPanel: FlipStats() not found to derive flip_summary.", call. = FALSE)
    }
    flip_summary <- FlipStats(STEAM.obj$nested$ncv, saa$final_by_fold)
  }
  
  # Validate required columns
  req <- c("Outer_Fold","Flips_To_Correct","Flips_To_Wrong","Delta_Accuracy")
  nm  <- if (is.null(flip_summary)) character() else names(flip_summary)
  if (is.null(flip_summary) || !all(req %in% nm)) {
    missing <- paste(setdiff(req, nm), collapse = ", ")
    stop("AnchorPanel: flip_summary missing required columns: ", missing, call. = FALSE)
  }
  if (!nrow(flip_summary)) {
    stop("AnchorPanel: flip_summary is empty; no flips to display.", call. = FALSE)
  }
  
  # Order folds as encountered
  flip_summary$Outer_Fold <- factor(flip_summary$Outer_Fold,
                                    levels = unique(flip_summary$Outer_Fold))
  
  # Long form for flips
  flips_long <- tidyr::pivot_longer(
    flip_summary[, c("Outer_Fold","Flips_To_Correct","Flips_To_Wrong")],
    names_to = "Category", values_to = "Count",
    cols = c("Flips_To_Correct","Flips_To_Wrong")
  )
  map_levels <- c(Flips_To_Correct = "Flipped to Correct",
                  Flips_To_Wrong   = "Flipped to Wrong")
  flips_long$Category <- factor(map_levels[as.character(flips_long$Category)],
                                levels = c("Flipped to Correct","Flipped to Wrong"))
  
  p_flips <- ggplot2::ggplot(flips_long, ggplot2::aes(x = Outer_Fold, y = Count, fill = Category)) +
    ggplot2::geom_col(position = "stack", width = 0.75) +
    ggplot2::geom_text(ggplot2::aes(label = Count),
                       position = ggplot2::position_stack(vjust = 0.5), size = 3, color = "white") +
    ggplot2::scale_fill_manual(values = c("Flipped to Correct" = "forestgreen",
                                          "Flipped to Wrong"   = "firebrick")) +
    ggplot2::labs(x = "Outer Fold", y = "# of Flips", fill = "Type",
                  title = "Flips per Outer Fold (Correct vs Wrong)") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  p_delta <- ggplot2::ggplot(flip_summary, ggplot2::aes(x = Outer_Fold, y = Delta_Accuracy)) +
    ggplot2::geom_col(width = 0.75, fill = "steelblue") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(x = "Outer Fold", y = "Δ Accuracy",
                  title = "Accuracy Gain from Spatial Corrections") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))
  
  # Compose without requiring patchwork to be attached
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(p_flips, p_delta, ncol = 2))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p_flips, p_delta, ncol = 2)
    return(invisible(NULL))
  } else {
    warning("Neither 'patchwork' nor 'gridExtra' is installed. Returning a list of plots.")
    return(list(flips = p_flips, delta = p_delta))
  }
}








#' Spatial flip map for a given outer fold
#'
#' Highlights cells flipped by STEAM_anchor() while showing all other cells in
#' a softly faded version of the layer palette.
#'
#' @param STEAM.obj STEAM object after [STEAM_anchor()].
#' @param fold Outer fold index to visualize.
#' @param layer_colors optional named vector to override the layer palette.
#' @param fade Amount (0..1) to fade non-flipped layers toward white.
#' @param point_size Scatter point size.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#'   p <- AnchorMap(STEAM.obj, fold = 1)
#'   print(p)
#' }
#' @export
AnchorMap <- function(steam_obj,
                      final_by_fold = NULL,   # <- now optional
                      fold = 1,
                      layer_colors = NULL,
                      fade = 0.35,
                      point_size = 3) {
  
  suppressPackageStartupMessages({
    library(ggplot2); library(scales)
  })
  
  # pull final_by_fold from STEAM.obj if not passed
  if (is.null(final_by_fold)) {
    faa <- steam_obj$spatial_anchor_analysis$final_by_fold
    if (is.null(faa) || !length(faa)) {
      stop("AnchorMap: no final_by_fold found. Run STEAM_anchor() first.", call. = FALSE)
    }
    final_by_fold <- faa
  }
  
  # ---- coords ----
  coordinates <- as.data.frame(steam_obj$spatial)
  if (!nrow(coordinates)) stop("steam_obj$spatial is empty.")
  cl <- tolower(colnames(coordinates))
  x_name <- colnames(coordinates)[which(cl %in% c("x","col","column"))[1]]
  y_name <- colnames(coordinates)[which(cl %in% c("y","row"))[1]]
  if (is.na(x_name) || is.na(y_name)) stop("Could not detect X/Y columns in steam_obj$spatial.")
  if (!"cell_id" %in% colnames(coordinates)) coordinates$cell_id <- rownames(coordinates)
  coordinates <- coordinates[, c("cell_id", x_name, y_name)]
  colnames(coordinates) <- c("cell_id","Col","Row")
  
  # ---- labels aligned to coords ----
  lab_vec <- steam_obj$labels
  if (is.null(names(lab_vec)) && !is.null(rownames(steam_obj$spatial))) {
    names(lab_vec) <- rownames(steam_obj$spatial)
  }
  base_lab <- lab_vec[coordinates$cell_id]
  
  full_data <- data.frame(
    cell_id = coordinates$cell_id,
    Col     = coordinates$Col,
    Row     = -coordinates$Row,
    Labels  = as.character(base_lab),
    stringsAsFactors = FALSE,
    row.names = coordinates$cell_id
  )
  
  # ---- fold preds to mark flips ----
  of <- steam_obj$nested$ncv$outer_result[[fold]]
  pd <- as.data.frame(of$preds)
  if (is.null(rownames(pd))) rownames(pd) <- as.character(of$test_id %||% of$test_ix)
  ids_pred <- rownames(pd)
  
  truth_test  <- setNames(as.character(pd$testy), ids_pred)
  before_test <- setNames(as.character(pd$predy), ids_pred)
  after_test  <- final_by_fold[[fold]]
  if (is.null(names(after_test))) names(after_test) <- ids_pred
  
  common_ids <- intersect(ids_pred, full_data$cell_id)
  truth  <- truth_test[common_ids]
  before <- before_test[common_ids]
  after  <- after_test[common_ids]
  
  flipped <- before != after
  if (any(flipped)) {
    ids_flip   <- common_ids[flipped]
    to_correct <- after[flipped] == truth[flipped]
    lab_vec2 <- ifelse(to_correct, "Flipped to Correct", "Flipped to Wrong")
    full_data$Labels[match(ids_flip, full_data$cell_id)] <- lab_vec2
  }
  
  # palette
  pal <- layer_colors %||% steam_obj$viz$label_colors %||% steam_obj$plot$label_colors
  if (is.null(pal)) {
    all_layers <- sort(unique(as.character(steam_obj$labels)))
    pal <- setNames(scales::hue_pal()(length(all_layers)), all_layers)
  }
  
  layer_levels <- sort(unique(setdiff(full_data$Labels,
                                      c("Flipped to Correct","Flipped to Wrong", NA_character_))))
  base_cols <- pal[intersect(names(pal), layer_levels)]
  missing <- setdiff(layer_levels, names(base_cols))
  if (length(missing)) {
    base_cols <- c(base_cols, setNames(scales::hue_pal()(length(missing)), missing))
  }
  base_cols <- base_cols[layer_levels]
  
  fade_color <- function(hex, f = fade) {
    oc <- col2rgb(hex)/255; rc <- (1-f)*oc + f*1
    rgb(rc[1], rc[2], rc[3])
  }
  faded_cols <- vapply(base_cols, fade_color, character(1))
  
  status_cols <- c("Flipped to Correct" = "forestgreen",
                   "Flipped to Wrong"   = "firebrick")
  
  labs_all <- unique(full_data$Labels)
  cols <- c(status_cols[names(status_cols) %in% labs_all],
            faded_cols[names(faded_cols) %in% labs_all])
  
  full_data$Labels <- factor(full_data$Labels,
                             levels = c(names(status_cols), names(pal)))
  
  ggplot2::ggplot(full_data, ggplot2::aes(x = Col, y = Row, color = Labels)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_manual(values = cols, drop = FALSE, na.value = "#D9D9D9") +
    ggplot2::labs(title = paste0("Flips vs True Labels (Fold ", fold, ")"),
                  color = "Layer / Status") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank())
}

# ========================================================================
# ITERATIVE STEAM ANCHOR VISUALIZATION FUNCTIONS
# ========================================================================

#' Plot Iterative Progress Overview
#'
#' @description
#' Shows the progress of corrections across iterations
#'
#' @param STEAM.obj STEAM object with iterative_anchor_analysis results
#' @return ggplot object
#' @export
plotIterativeProgress <- function(STEAM.obj) {
  
  if (is.null(STEAM.obj$iterative_anchor_analysis)) {
    stop("No iterative_anchor_analysis results found. Run STEAM_anchor with iterative=TRUE first.")
  }
  
  iter_analysis <- STEAM.obj$iterative_anchor_analysis
  
  # Extract iteration data
  iteration_data <- data.frame()
  
  if (iter_analysis$method == "iterative_with_ground_truth") {
    # Ground truth method
    for (i in seq_along(iter_analysis$iteration_results)) {
      corrections <- if (!is.null(iter_analysis$iteration_results[[i]]$corrections)) {
        nrow(iter_analysis$iteration_results[[i]]$corrections)
      } else {
        iter_analysis$iteration_results[[i]]$total_corrections
      }
      
      iteration_data <- rbind(iteration_data, data.frame(
        iteration = i,
        corrections = corrections,
        cumulative = sum(iteration_data$corrections, corrections),
        method = "Ground Truth"
      ))
    }
  } else {
    # No ground truth method
    for (i in seq_along(iter_analysis$correction_history)) {
      corrections <- if (!is.null(iter_analysis$correction_history[[i]])) {
        nrow(iter_analysis$correction_history[[i]])
      } else {
        0
      }
      
      iteration_data <- rbind(iteration_data, data.frame(
        iteration = i,
        corrections = corrections,
        cumulative = sum(iteration_data$corrections, corrections),
        method = "Alternative Validation"
      ))
    }
  }
  
  # Create plot
  p <- ggplot(iteration_data, aes(x = iteration)) +
    geom_col(aes(y = corrections), fill = "steelblue", alpha = 0.7, width = 0.6) +
    geom_line(aes(y = cumulative), color = "red", linewidth = 1.2, group = 1) +
    geom_point(aes(y = cumulative), color = "red", size = 3) +
    geom_text(aes(y = corrections, label = corrections), vjust = -0.5, size = 3) +
    geom_text(aes(y = cumulative, label = cumulative), vjust = -0.5, color = "red", size = 3) +
    theme_minimal() +
    labs(
      title = sprintf("Iterative STEAM Progress (%s)", iter_analysis$method),
      subtitle = sprintf("Total: %d corrections across %d iterations", 
                        iter_analysis$summary$total_corrections, 
                        iter_analysis$iterations_completed),
      x = "Iteration",
      y = "Number of Corrections"
    ) +
    scale_x_continuous(breaks = iteration_data$iteration) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(p)
}

#' Plot Spatial Corrections by Iteration
#'
#' @description
#' Shows spatial distribution of corrections for each iteration
#'
#' @param STEAM.obj STEAM object with iterative results
#' @param max_iterations Maximum number of iterations to plot (default: 4)
#' @param point_size Size of points (default: 1.5)
#' @return List of ggplot objects (one per iteration)
#' @export
plotCorrectionsByIteration <- function(STEAM.obj, max_iterations = 4, point_size = 1.5) {
  
  if (is.null(STEAM.obj$iterative_anchor_analysis)) {
    stop("No iterative_anchor_analysis results found. Run STEAM_anchor with iterative=TRUE first.")
  }
  
  iter_analysis <- STEAM.obj$iterative_anchor_analysis
  max_iter <- min(max_iterations, iter_analysis$iterations_completed)
  
  plots <- list()
  
  for (iter in 1:max_iter) {
    # Skip iteration 1 if it's just baseline
    if (iter == 1 && iter_analysis$method == "iterative_with_ground_truth") {
      next
    }
    
    plot_title <- sprintf("Iteration %d Corrections", iter)
    p <- visualizeSpatialFlips(STEAM.obj, iteration = iter, 
                                point_size = point_size, 
                                show_iteration_info = FALSE) +
      labs(title = plot_title)
    
    plots[[paste0("iteration_", iter)]] <- p
  }
  
  return(plots)
}

#' Plot Iterative Improvement Metrics
#'
#' @description
#' Shows improvement metrics across iterations (for ground truth method)
#'
#' @param STEAM.obj STEAM object with iterative results
#' @return ggplot object
#' @export
plotIterativeImprovement <- function(STEAM.obj) {
  
  if (is.null(STEAM.obj$iterative_anchor_analysis)) {
    stop("No iterative_anchor_analysis results found.")
  }
  
  iter_analysis <- STEAM.obj$iterative_anchor_analysis
  
  if (iter_analysis$method == "iterative_with_ground_truth") {
    # Extract improvement data from ground truth method
    improvement_data <- data.frame()
    
    for (i in 2:length(iter_analysis$iteration_results)) {
      if (!is.null(iter_analysis$iteration_results[[i]])) {
        # Try to extract improvement metrics if available
        corrections <- iter_analysis$iteration_results[[i]]$total_corrections %||% 0
        
        improvement_data <- rbind(improvement_data, data.frame(
          iteration = i,
          corrections = corrections,
          method = "Ground Truth"
        ))
      }
    }
    
    title_text <- "Iterative Improvement (Ground Truth Validation)"
    
  } else {
    # No ground truth method - use spatial coherence metrics
    improvement_data <- data.frame()
    
    for (i in seq_along(iter_analysis$iteration_metrics)) {
      metrics <- iter_analysis$iteration_metrics[[i]]
      
      improvement_data <- rbind(improvement_data, data.frame(
        iteration = i,
        spatial_coherence = metrics$spatial_coherence,
        expression_coherence = metrics$expression_coherence,
        method = "Alternative Validation"
      ))
    }
    
    title_text <- "Iterative Improvement (Alternative Validation)"
  }
  
  # Create appropriate plot based on available data
  if (iter_analysis$method == "iterative_with_ground_truth") {
    p <- ggplot(improvement_data, aes(x = iteration, y = corrections)) +
      geom_line(color = "blue", linewidth = 1.2) +
      geom_point(color = "blue", size = 3) +
      geom_text(aes(label = corrections), vjust = -0.5) +
      theme_minimal() +
      labs(
        title = title_text,
        x = "Iteration",
        y = "Number of Corrections"
      )
  } else {
    # Reshape data for multiple metrics
    plot_data <- tidyr::pivot_longer(
      improvement_data, 
      cols = c("spatial_coherence", "expression_coherence"),
      names_to = "metric",
      values_to = "value"
    )
    
    p <- ggplot(plot_data, aes(x = iteration, y = value, color = metric)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      theme_minimal() +
      labs(
        title = title_text,
        x = "Iteration",
        y = "Coherence Score",
        color = "Metric"
      ) +
      scale_color_manual(
        values = c("spatial_coherence" = "blue", "expression_coherence" = "red"),
        labels = c("Spatial Coherence", "Expression Coherence")
      )
  }
  
  return(p)
}

#' Create Comprehensive Iterative Analysis Dashboard
#'
#' @description
#' Creates a multi-panel dashboard showing iterative analysis results
#'
#' @param STEAM.obj STEAM object with iterative results
#' @param max_spatial_plots Maximum number of spatial plots to include (default: 3)
#' @return Combined plot object (requires patchwork)
#' @export
createIterativeDashboard <- function(STEAM.obj, max_spatial_plots = 3) {
  
  if (is.null(STEAM.obj$iterative_anchor_analysis)) {
    stop("No iterative_anchor_analysis results found.")
  }
  
  # Create individual plots
  p_progress <- plotIterativeProgress(STEAM.obj)
  p_improvement <- plotIterativeImprovement(STEAM.obj)
  p_overall <- visualizeSpatialFlips(STEAM.obj, iteration = NULL, 
                                      show_iteration_info = TRUE)
  
  # Get spatial plots for individual iterations
  spatial_plots <- plotCorrectionsByIteration(STEAM.obj, 
                                                 max_iterations = max_spatial_plots,
                                                 point_size = 1.2)
  
  # Combine plots
  if (requireNamespace("patchwork", quietly = TRUE)) {
    # Top row: progress and improvement
    top_row <- patchwork::wrap_plots(p_progress, p_improvement, ncol = 2)
    
    # Middle row: overall spatial view
    middle_row <- p_overall
    
    # Bottom row: individual iteration spatial plots
    if (length(spatial_plots) > 0) {
      bottom_row <- patchwork::wrap_plots(plotlist = spatial_plots, ncol = min(3, length(spatial_plots)))
      
      # Combine all
      combined <- top_row / middle_row / bottom_row
    } else {
      combined <- top_row / middle_row
    }
    
    combined <- combined + patchwork::plot_annotation(
      title = "STEAM Iterative Anchor Analysis Dashboard",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
    
    return(combined)
    
  } else {
    warning("patchwork package not available. Returning list of plots.")
    return(list(
      progress = p_progress,
      improvement = p_improvement,
      overall_spatial = p_overall,
      spatial_by_iteration = spatial_plots
    ))
  }
}

# ========================================================================
# ENHANCED STEAM PLOTTING FUNCTIONS
# Comprehensive visualization tools for enhanced STEAM analysis
# ========================================================================

#' Quick Visualization for Iterative STEAM Results
#'
#' @description
#' One-function approach to visualize all iterative STEAM anchor results
#'
#' @param STEAM.obj STEAM object with iterative results
#' @param plot_type Type of visualization: "dashboard", "progress", "spatial", "improvement", "corrected", "comparison"
#' @param save_plots Whether to save plots to files (default: FALSE)
#' @param output_dir Directory to save plots (default: current directory)
#' @param fold Which fold to show for corrected plots ("all" or specific fold number)
#' @param iteration Which iteration to show (NULL for all, 0 for single-pass)
#' @return Plot object or list of plots
#' @export
quickIterativeViz <- function(STEAM.obj, 
                                plot_type = "dashboard", 
                                save_plots = FALSE,
                                output_dir = ".",
                                fold = "all",
                                iteration = NULL) {
  
  if (is.null(STEAM.obj$iterative_anchor_analysis) && !plot_type %in% c("corrected", "comparison")) {
    stop("No iterative_anchor_analysis results found. Run STEAM_anchor with iterative=TRUE first.")
  }
  
  result <- switch(plot_type,
    "dashboard" = createIterativeDashboard(STEAM.obj),
    "progress" = plotIterativeProgress(STEAM.obj),
    "spatial" = visualizeSpatialFlips(STEAM.obj, iteration = NULL),
    "improvement" = plotIterativeImprovement(STEAM.obj),
    "corrected" = plotCorrectedCells(STEAM.obj, fold = fold, iteration = iteration),
    "comparison" = plotBeforeAfterComparison(STEAM.obj, fold = fold, iteration = iteration),
    "all" = list(
      dashboard = createIterativeDashboard(STEAM.obj),
      progress = plotIterativeProgress(STEAM.obj),
      spatial = visualizeSpatialFlips(STEAM.obj, iteration = NULL),
      improvement = plotIterativeImprovement(STEAM.obj),
      corrected = plotCorrectedCells(STEAM.obj, fold = fold, iteration = iteration),
      comparison = plotBeforeAfterComparison(STEAM.obj, fold = fold, iteration = iteration)
    ),
    stop("plot_type must be one of: 'dashboard', 'progress', 'spatial', 'improvement', 'corrected', 'comparison', 'all'")
  )
  
  # Save plots if requested
  if (save_plots && requireNamespace("ggplot2", quietly = TRUE)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    if (plot_type == "all") {
      for (name in names(result)) {
        if (inherits(result[[name]], "ggplot") || inherits(result[[name]], "patchwork")) {
          filename <- file.path(output_dir, paste0("STEAM_iterative_", name, ".png"))
          ggplot2::ggsave(filename, result[[name]], width = 12, height = 8, dpi = 300)
          message("Saved: ", filename)
        }
      }
    } else {
      filename <- file.path(output_dir, paste0("STEAM_iterative_", plot_type, ".png"))
      ggplot2::ggsave(filename, result, width = 12, height = 8, dpi = 300)
      message("Saved: ", filename)
    }
  }
  
  return(result)
}

# ========================================================================

# Required packages check for enhanced plotting
checkPlotPackages <- function() {
  required_pkgs <- c("ggplot2", "viridis")
  missing <- setdiff(required_pkgs, rownames(installed.packages()))
  if (length(missing) > 0) {
    warning(paste("Enhanced plotting packages missing:", paste(missing, collapse = ", ")))
  }
}

#' Extract anchor corrections from analysis results
#' @param spatial_anchor_analysis Results from STEAM_anchor
#' @param STEAM.obj STEAM object
#' @return Data frame with corrections
extractAnchorCorrections <- function(spatial_anchor_analysis, STEAM.obj) {
  
  # Get final predictions by fold
  final_by_fold <- spatial_anchor_analysis$final_by_fold
  if (is.null(final_by_fold)) {
    return(data.frame())
  }
  
  # Get original predictions from nested CV
  ncv <- STEAM.obj$nested$ncv
  all_corrections <- data.frame()
  
  for (f in seq_along(ncv$outer_result)) {
    preds_df <- ncv$outer_result[[f]]$preds
    if (!is.null(preds_df) && f <= length(final_by_fold)) {
      # Get cell IDs for this fold
      cell_ids <- rownames(preds_df)
      
      # Get original predictions
      original_pred <- preds_df$predy
      true_label <- preds_df$testy
      
      # Get final predictions for this fold
      final_pred_fold <- final_by_fold[[f]]
      
      # Align predictions by cell ID
      if (!is.null(names(final_pred_fold))) {
        # Reorder final predictions to match the order of cell_ids in preds_df
        final_pred <- final_pred_fold[cell_ids]
        
        # Check for any missing cells (NA values)
        missing_cells <- is.na(final_pred)
        if (any(missing_cells)) {
          # For missing cells, use original predictions (no correction)
          final_pred[missing_cells] <- original_pred[missing_cells]
        }
      } else {
        # If final_pred_fold doesn't have names, assume same order as preds_df
        if (length(final_pred_fold) != length(original_pred)) {
          warning(sprintf("Fold %d: Length mismatch between original (%d) and final (%d) predictions. Using original predictions.", 
                         f, length(original_pred), length(final_pred_fold)))
          final_pred <- original_pred  # No corrections for this fold
        } else {
          final_pred <- final_pred_fold
        }
      }
      
      # Ensure vectors are same length (should be guaranteed now)
      if (length(original_pred) != length(final_pred) || length(original_pred) != length(true_label)) {
        warning(sprintf("Fold %d: Vector length mismatch after alignment. orig=%d, final=%d, true=%d. Skipping.", 
                       f, length(original_pred), length(final_pred), length(true_label)))
        next
      }
      
      # Find corrections (where final != original)
      corrected_idx <- original_pred != final_pred
      
      if (any(corrected_idx)) {
        fold_corrections <- data.frame(
          cell_id = cell_ids[corrected_idx],
          original_pred = original_pred[corrected_idx],
          suggested_correction = final_pred[corrected_idx],
          true_label = true_label[corrected_idx],
          fold = f,
          stringsAsFactors = FALSE
        )
        
        # Check if correction is to correct label
        fold_corrections$correction_correct <- fold_corrections$suggested_correction == fold_corrections$true_label
        
        all_corrections <- rbind(all_corrections, fold_corrections)
      }
    }
  }
  
  return(all_corrections)
}

#' Extract summary statistics from anchor analysis
#' @param spatial_anchor_analysis Results from STEAM_anchor
#' @param STEAM.obj STEAM object
#' @return List with summary statistics
extractAnchorSummaryStats <- function(spatial_anchor_analysis, STEAM.obj) {
  
  # Get nested CV results
  ncv <- STEAM.obj$nested$ncv
  
  # Count total cells and misclassifications
  total_cells <- 0
  total_misclassified <- 0
  
  for (f in seq_along(ncv$outer_result)) {
    preds_df <- ncv$outer_result[[f]]$preds
    if (!is.null(preds_df)) {
      fold_cells <- nrow(preds_df)
      fold_misclassified <- sum(preds_df$predy != preds_df$testy)
      
      total_cells <- total_cells + fold_cells
      total_misclassified <- total_misclassified + fold_misclassified
    }
  }
  
  # Get corrections
  corrections <- extractAnchorCorrections(spatial_anchor_analysis, STEAM.obj)
  total_corrections <- nrow(corrections)
  
  corrections_to_correct <- 0
  corrections_to_wrong <- 0
  if (total_corrections > 0) {
    corrections_to_correct <- sum(corrections$correction_correct)
    corrections_to_wrong <- sum(!corrections$correction_correct)
  }
  
  # Calculate rates
  misclassification_rate <- if (total_cells > 0) total_misclassified / total_cells else 0
  correction_rate <- if (total_misclassified > 0) total_corrections / total_misclassified else 0
  
  return(list(
    total_cells = total_cells,
    total_misclassified = total_misclassified,
    total_corrections = total_corrections,
    corrections_to_correct = corrections_to_correct,
    corrections_to_wrong = corrections_to_wrong,
    misclassification_rate = misclassification_rate,
    correction_rate = correction_rate
  ))
}

# ========================================================================
# STEAM ANCHOR VISUALIZATION FUNCTIONS
# ========================================================================

#' Visualize Spatial Distribution of Anchor Corrections (Iterative Compatible)
#'
#' @description
#' Shows where corrections happened spatially, supports both single-pass and iterative results
#'
#' @param STEAM.obj STEAM object with spatial_anchor_analysis results
#' @param iteration Which iteration to visualize (NULL for all corrections, 0 for single-pass)
#' @param point_size Size of points (default: 2)
#' @param highlight_corrections Whether to highlight corrected cells (default: TRUE)
#' @param show_iteration_info Whether to show iteration information in title (default: TRUE)
#' 
#' @return ggplot object
#' @export
visualizeSpatialFlips <- function(STEAM.obj, iteration = NULL, point_size = 2, 
                                   highlight_corrections = TRUE, show_iteration_info = TRUE) {
  
  
  # Check for iterative analysis results
  has_iterative <- !is.null(STEAM.obj$iterative_anchor_analysis)
  has_single_pass <- !is.null(STEAM.obj$spatial_anchor_analysis)
  
  if (!has_iterative && !has_single_pass) {
    stop("No spatial_anchor_analysis or iterative_anchor_analysis results found. Run STEAM_anchor first.")
  }
  
  # Extract corrections based on analysis type and iteration parameter
  corrections <- data.frame()
  analysis_info <- ""
  
  if (has_iterative) {
    # Iterative analysis
    iter_analysis <- STEAM.obj$iterative_anchor_analysis
    
    if (is.null(iteration)) {
      # Show all corrections from all iterations
      if (iter_analysis$method == "iterative_with_ground_truth") {
        for (i in seq_along(iter_analysis$iteration_results)) {
          if (!is.null(iter_analysis$iteration_results[[i]]$corrections)) {
            iter_corr <- iter_analysis$iteration_results[[i]]$corrections
            iter_corr$iteration <- i
            corrections <- rbind(corrections, iter_corr)
          }
        }
        analysis_info <- sprintf("All iterations (%d total)", iter_analysis$iterations_completed)
      } else {
        # No ground truth - use correction_history
        for (i in seq_along(iter_analysis$correction_history)) {
          if (!is.null(iter_analysis$correction_history[[i]])) {
            iter_corr <- iter_analysis$correction_history[[i]]
            iter_corr$iteration <- i
            corrections <- rbind(corrections, iter_corr)
          }
        }
        analysis_info <- sprintf("All iterations (%d total)", iter_analysis$iterations_completed)
      }
    } else if (iteration == 0) {
      # Show single-pass results only
      if (has_single_pass) {
        corrections <- extractAnchorCorrections(STEAM.obj$spatial_anchor_analysis, STEAM.obj)
        analysis_info <- "Single-pass only"
      } else {
        stop("No single-pass results available. Use iteration = NULL for all iterative results.")
      }
    } else {
      # Show specific iteration
      if (iteration > iter_analysis$iterations_completed) {
        stop(sprintf("Iteration %d not available. Maximum iteration: %d", 
                    iteration, iter_analysis$iterations_completed))
      }
      
      if (iter_analysis$method == "iterative_with_ground_truth") {
        if (!is.null(iter_analysis$iteration_results[[iteration]]$corrections)) {
          corrections <- iter_analysis$iteration_results[[iteration]]$corrections
          corrections$iteration <- iteration
        }
      } else {
        if (!is.null(iter_analysis$correction_history[[iteration]])) {
          corrections <- iter_analysis$correction_history[[iteration]]
          corrections$iteration <- iteration
        }
      }
      analysis_info <- sprintf("Iteration %d", iteration)
    }
  } else {
    # Single-pass only
    corrections <- extractAnchorCorrections(STEAM.obj$spatial_anchor_analysis, STEAM.obj)
    analysis_info <- "Single-pass"
  }
  
  # Extract summary statistics
  if (has_iterative && is.null(iteration)) {
    # Use iterative summary
    total_corrections <- iter_analysis$summary$total_corrections
    if (iter_analysis$method == "iterative_with_ground_truth") {
      total_misclassified <- sum(sapply(STEAM.obj$nested$ncv$outer_result, function(x) {
        if (!is.null(x$preds)) sum(x$preds$testy != x$preds$predy) else 0
      }))
    } else {
      total_misclassified <- total_corrections  # Conservative estimate for no ground truth
    }
  } else {
    # Use single-pass or specific iteration stats
    stats <- extractAnchorSummaryStats(STEAM.obj$spatial_anchor_analysis, STEAM.obj)
    total_corrections <- ifelse(nrow(corrections) > 0, nrow(corrections), stats$total_corrections)
    total_misclassified <- stats$total_misclassified
  }
  
  # Get spatial coordinates
  coords <- as.data.frame(STEAM.obj$spatial)
  colnames(coords) <- c("x", "y")
  coords$cell_id <- rownames(coords)
  
  # Get all predictions
  ncv <- STEAM.obj$nested$ncv
  all_preds <- data.frame()
  
  for (f in seq_along(ncv$outer_result)) {
    preds_df <- ncv$outer_result[[f]]$preds
    if (!is.null(preds_df) && all(c("predy", "testy") %in% colnames(preds_df))) {
      fold_preds <- data.frame(
        cell_id = rownames(preds_df),
        original_pred = preds_df$predy,
        true_label = preds_df$testy,
        fold = f,
        stringsAsFactors = FALSE
      )
      all_preds <- rbind(all_preds, fold_preds)
    }
  }
  
  # Remove duplicates (keep first occurrence)
  all_preds <- all_preds[!duplicated(all_preds$cell_id), ]
  
  # Merge with coordinates
  plot_data <- merge(coords, all_preds, by = "cell_id", all.x = TRUE)
  
  # Add correction status
  plot_data$status <- "Correct"
  plot_data$status[plot_data$original_pred != plot_data$true_label] <- "Misclassified"
  
  if (nrow(corrections) > 0) {
    # Mark corrected cells
    plot_data$status[plot_data$cell_id %in% corrections$cell_id] <- "Corrected"
    
    # Update predictions for corrected cells
    for (i in seq_len(nrow(corrections))) {
      idx <- plot_data$cell_id == corrections$cell_id[i]
      plot_data$display_pred <- plot_data$original_pred
      plot_data$display_pred[idx] <- corrections$suggested_correction[i]
    }
  } else {
    plot_data$display_pred <- plot_data$original_pred
  }
  
  plot_data$status <- factor(plot_data$status, 
                            levels = c("Correct", "Misclassified", "Corrected"))
  
  # Create color palette
  status_colors <- c("Correct" = "gray80", 
                    "Misclassified" = "red3", 
                    "Corrected" = "blue3")
  
  # Create title text
  if (show_iteration_info) {
    title_text <- sprintf("Spatial Corrections (%s)\n%d corrections, %d total misclassified", 
                         analysis_info, total_corrections, total_misclassified)
  } else {
    title_text <- sprintf("Spatial Corrections (%s)", analysis_info)
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = x, y = y)) +
    geom_point(aes(color = status), size = point_size, alpha = 0.8) +
    scale_color_manual(values = status_colors, name = "Status") +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = title_text)
  
  if (highlight_corrections && nrow(corrections) > 0) {
    # Add highlights for corrected cells
    corrected_data <- plot_data[plot_data$status == "Corrected", ]
    p <- p + 
      geom_point(data = corrected_data, 
                aes(x = x, y = y), 
                color = "blue3", size = point_size * 1.5, 
                shape = 21, stroke = 1.5, fill = NA)
  }
  
  return(p)
}

#' Show Anchor Correction Summary
#'
#' @description
#' Print summary of corrections
#'
#' @param STEAM.obj STEAM object with spatial_anchor_analysis results
#' @export
showFlipSummary <- function(STEAM.obj) {
  
  if (is.null(STEAM.obj$spatial_anchor_analysis)) {
    stop("No spatial_anchor_analysis results found. Run STEAM_anchor first.")
  }
  
  # Extract statistics
  stats <- extractAnchorSummaryStats(STEAM.obj$spatial_anchor_analysis, STEAM.obj)
  corrections <- extractAnchorCorrections(STEAM.obj$spatial_anchor_analysis, STEAM.obj)
  
  cat("\n=== STEAM Anchor Flip Summary ===\n")
  cat(sprintf("Total cells analyzed: %d\n", stats$total_cells))
  cat(sprintf("Misclassified cells: %d (%.1f%%)\n", 
              stats$total_misclassified, 
              stats$misclassification_rate * 100))
  cat(sprintf("Correct flips applied: %d (%.1f%% of misclassified)\n", 
              stats$total_corrections,
              stats$correction_rate * 100))
  
  if (stats$total_corrections > 0) {
    cat(sprintf("  - Corrections to correct label: %d (%.1f%%)\n",
                stats$corrections_to_correct,
                stats$corrections_to_correct / stats$total_corrections * 100))
    cat(sprintf("  - Corrections to wrong label: %d (%.1f%%)\n",
                stats$corrections_to_wrong,
                stats$corrections_to_wrong / stats$total_corrections * 100))
    
    # Show correction breakdown by cell type
    cat("\nCorrections by original cell type:\n")
    correction_summary <- corrections %>%
      group_by(original_pred) %>%
      summarise(
        n_corrections = n(),
        n_correct = sum(correction_correct),
        .groups = 'drop'
      ) %>%
      arrange(desc(n_corrections))
    
    for (i in seq_len(nrow(correction_summary))) {
      cat(sprintf("  %s: %d corrections (%d correct)\n",
                  correction_summary$original_pred[i],
                  correction_summary$n_corrections[i],
                  correction_summary$n_correct[i]))
    }
  } else {
    cat("\nNo corrections were made.\n")
  }
  
  cat("\n")
}

#' Visualize Anchor Correction Details
#'
#' @description
#' Show detailed view of corrections with before/after labels
#'
#' @param STEAM.obj STEAM object with spatial_anchor_analysis results
#' @param max_cells Maximum number of cells to show (default: 20)
#' @export
visualizeCorrectionDetails <- function(STEAM.obj, max_cells = 20) {
  
  if (is.null(STEAM.obj$spatial_anchor_analysis)) {
    stop("No spatial_anchor_analysis results found. Run STEAM_anchor first.")
  }
  
  # Extract corrections
  corrections <- extractAnchorCorrections(STEAM.obj$spatial_anchor_analysis, STEAM.obj)
  
  if (nrow(corrections) == 0) {
    cat("No corrections were made.\n")
    return(invisible(NULL))
  }
  
  # Limit to max_cells
  if (nrow(corrections) > max_cells) {
    corrections <- corrections[seq_len(max_cells), ]
    cat(sprintf("Showing first %d corrections out of %d total\n\n", 
                max_cells, nrow(corrections)))
  }
  
  # Create visualization data
  viz_data <- corrections %>%
    select(cell_id, original_pred, suggested_correction, true_label) %>%
    mutate(
      correction_type = ifelse(suggested_correction == true_label, 
                              "Correct", "Incorrect"),
      change = paste(original_pred, "→", suggested_correction)
    )
  
  # Create plot
  p <- ggplot(viz_data, aes(y = reorder(cell_id, seq_len(nrow(viz_data))), x = 1)) +
    geom_tile(aes(fill = correction_type), width = 0.8, height = 0.8) +
    geom_text(aes(label = change), size = 3) +
    scale_fill_manual(values = c("Correct" = "lightgreen", "Incorrect" = "lightcoral"),
                     name = "Correction") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(title = "Correction Details (Original → Corrected)",
         y = "Cell ID")
  
  return(p)
}

#' Plot spatial overview of tissue with cell types
#' @param STEAM.obj STEAM object containing spatial coordinates and labels
#' @param color_by Column name to color points by (default: "labels")
#' @param colors Named vector of colors for cell types (optional)
#' @param point_size Size of points (default: 0.8)
#' @param alpha Transparency (default: 0.7)
#' @param title Plot title
#' @return ggplot object
#' @export
plotSpatialOverview <- function(STEAM.obj, color_by = "labels", colors = NULL, 
                                 point_size = 0.8, alpha = 0.7, title = "Spatial Overview") {
  
  # Extract spatial coordinates
  coords <- STEAM.obj$spatial
  if (is.null(coords)) stop("No spatial coordinates found in STEAM.obj$spatial")
  
  # Get color values
  if (color_by == "labels") {
    color_values <- STEAM.obj$labels
  } else {
    color_values <- STEAM.obj[[color_by]]
  }
  
  if (is.null(color_values)) stop(paste("Column", color_by, "not found in STEAM.obj"))
  
  # Create plot data
  plot_data <- data.frame(
    x = coords[, 1],
    y = coords[, 2], 
    color_var = color_values
  )
  
  # Create base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = color_var)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14)
    ) +
    ggplot2::labs(title = title, x = "", y = "", color = tools::toTitleCase(color_by)) +
    ggplot2::coord_fixed()
  
  # Apply custom colors if provided
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors)
  } else if (requireNamespace("viridis", quietly = TRUE)) {
    if (is.numeric(color_values)) {
      p <- p + ggplot2::scale_color_viridis_c()
    } else {
      p <- p + ggplot2::scale_color_viridis_d()
    }
  }
  
  return(p)
}

#' Plot spatial confidence scores from enhanced STEAM analysis
#' @param enhanced_results STEAM_anchor_enhanced results
#' @param confidence_type Which confidence metric to plot
#' @param color_palette Color palette name (default: "viridis")
#' @param point_size Point size (default: 1.2)
#' @param title Plot title (auto-generated if NULL)
#' @return ggplot object
#' @export
plotConfidenceSpatial <- function(enhanced_results, 
                                   confidence_type = "multimodal_confidence",
                                   color_palette = "viridis", point_size = 1.2, 
                                   title = NULL) {
  
  # Extract confidence data
  confidence_data <- enhanced_results$spatial_anchor_analysis$multimodal_confidence
  if (is.null(confidence_data)) {
    stop("No multimodal confidence data found in results")
  }
  
  if (!(confidence_type %in% colnames(confidence_data))) {
    stop(paste("Confidence type", confidence_type, "not found in data"))
  }
  
  # Get spatial coordinates from enhanced results or global environment
  cell_ids <- confidence_data$cell_id
  
  # Try to get coordinates from STEAM.obj in global environment
  if (exists("STEAM.obj", envir = .GlobalEnv)) {
    coords <- get("STEAM.obj", envir = .GlobalEnv)$spatial[cell_ids, ]
  } else {
    stop("STEAM.obj not found in global environment - needed for spatial coordinates")
  }
  
  confidence_values <- confidence_data[[confidence_type]]
  
  # Create plot data
  plot_data <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    confidence = confidence_values
  )
  
  # Auto-generate title if not provided
  if (is.null(title)) {
    title <- paste("Spatial", tools::toTitleCase(gsub("_", " ", confidence_type)))
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = confidence)) +
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      legend.title = ggplot2::element_text(size = 12)
    ) +
    ggplot2::labs(title = title, x = "", y = "", 
         color = tools::toTitleCase(gsub("_", " ", confidence_type))) +
    ggplot2::coord_fixed()
  
  # Apply color palette
  if (requireNamespace("viridis", quietly = TRUE) && color_palette == "viridis") {
    p <- p + ggplot2::scale_color_viridis_c()
  } else {
    p <- p + ggplot2::scale_color_gradient(low = "blue", high = "red")
  }
  
  return(p)
}

#' Plot expression gradients spatially
#' @param enhanced_results STEAM_anchor_enhanced results
#' @param gradient_type Type of gradient to plot ("magnitude")
#' @param color_palette Color palette (default: "plasma")
#' @param title Plot title
#' @return ggplot object
#' @export
plotExpressionGradients <- function(enhanced_results, gradient_type = "magnitude",
                                     color_palette = "plasma", title = "Expression Gradients") {
  
  # Extract gradient data
  gradients <- enhanced_results$spatial_anchor_analysis$expression_gradients
  if (is.null(gradients)) {
    stop("No expression gradients found in results")
  }
  
  # Get spatial coordinates
  if (exists("STEAM.obj", envir = .GlobalEnv)) {
    coords <- get("STEAM.obj", envir = .GlobalEnv)$spatial[names(gradients), ]
  } else {
    stop("STEAM.obj not found in global environment")
  }
  
  # Create plot data
  plot_data <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    gradient = as.numeric(gradients)
  )
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = gradient)) +
    ggplot2::geom_point(size = 1.0, alpha = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14)
    ) +
    ggplot2::labs(title = title, x = "", y = "", color = "Gradient\nMagnitude") +
    ggplot2::coord_fixed()
  
  # Apply color palette
  if (requireNamespace("viridis", quietly = TRUE) && color_palette == "plasma") {
    p <- p + ggplot2::scale_color_viridis_c(option = "plasma")
  } else {
    p <- p + ggplot2::scale_color_gradient(low = "blue", high = "yellow")
  }
  
  return(p)
}

#' Plot expression boundaries
#' @param enhanced_results STEAM_anchor_enhanced results
#' @param boundary_threshold Threshold for highlighting boundaries (default: 0.8)
#' @param highlight_color Color for boundary cells (default: "red")
#' @param title Plot title
#' @return ggplot object
#' @export
plotExpressionBoundaries <- function(enhanced_results, boundary_threshold = 0.8,
                                      highlight_color = "red", title = "Expression Boundaries") {
  
  # Extract boundary data
  boundaries <- enhanced_results$spatial_anchor_analysis$expression_boundaries
  if (is.null(boundaries)) {
    stop("No expression boundaries found in results")
  }
  
  # Get spatial coordinates
  if (exists("STEAM.obj", envir = .GlobalEnv)) {
    coords <- get("STEAM.obj", envir = .GlobalEnv)$spatial[boundaries$cell_id, ]
  } else {
    stop("STEAM.obj not found in global environment")
  }
  
  # Create plot data
  plot_data <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    boundary_score = boundaries$boundary_score,
    is_boundary = boundaries$boundary_score > boundary_threshold
  )
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(color = boundary_score), size = 0.8, alpha = 0.6) +
    ggplot2::geom_point(data = plot_data[plot_data$is_boundary, ], 
               ggplot2::aes(x = x, y = y), color = highlight_color, size = 1.2, 
               alpha = 0.8, shape = 1, stroke = 1) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14)
    ) +
    ggplot2::labs(title = title, x = "", y = "", color = "Boundary\nScore") +
    ggplot2::coord_fixed()
  
  # Apply color palette
  if (requireNamespace("viridis", quietly = TRUE)) {
    p <- p + ggplot2::scale_color_viridis_c()
  }
  
  return(p)
}

#' Plot gene expression spatially
#' @param STEAM.obj STEAM object with expression data
#' @param gene Gene name to plot
#' @param color_palette Color palette (default: "Blues")
#' @param title Plot title
#' @return ggplot object
#' @export
plotGeneExpressionSpatial <- function(STEAM.obj, gene, color_palette = "Blues", title = NULL) {
  
  # Extract expression data
  expr_matrix <- STEAM.obj$count_exp
  if (is.null(expr_matrix)) {
    stop("No expression data found in STEAM.obj$count_exp")
  }
  
  # Check gene exists
  if (gene %in% rownames(expr_matrix)) {
    expr_values <- expr_matrix[gene, ]
  } else if (gene %in% colnames(expr_matrix)) {
    expr_values <- expr_matrix[, gene]
  } else {
    stop(paste("Gene", gene, "not found in expression matrix"))
  }
  
  # Get spatial coordinates
  coords <- STEAM.obj$spatial[names(expr_values), ]
  
  # Create plot data
  plot_data <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    expression = as.numeric(expr_values)
  )
  
  # Auto-generate title
  if (is.null(title)) {
    title <- paste(gene, "Expression")
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = expression)) +
    ggplot2::geom_point(size = 1.0, alpha = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14)
    ) +
    ggplot2::labs(title = title, x = "", y = "", color = "Expression\nLevel") +
    ggplot2::coord_fixed()
  
  # Apply color palette
  p <- p + ggplot2::scale_color_gradient(low = "white", high = "darkblue")
  
  return(p)
}

#' Plot corrected cells with same color scheme as plot_misclassified_cells
#'
#' @description
#' Shows corrected cells in spatial context with consistent color scheme
#' for easy comparison with misclassified cells
#'
#' @param steam_obj STEAM object with spatial_anchor_analysis or iterative results
#' @param label_colors Named vector of colors for cell types (same as plot_misclassified_cells)
#' @param fold Which fold to show ("all" for all folds, or specific fold number)
#' @param iteration Which iteration to show (NULL for all, 0 for single-pass, specific number for iteration)
#' @return ggplot object
#' @export
plotCorrectedCells <- function(steam_obj, label_colors = NULL, fold = "all", iteration = NULL) {
  
  suppressPackageStartupMessages({
    library(ggplot2); library(scales)
  })
  
  # ---- coords ----
  coordinates <- steam_obj$spatial
  if (is.null(coordinates)) stop("No spatial coordinates found in STEAM object")
  coordinates <- as.data.frame(coordinates)
  
  # detect coordinate columns
  cl <- tolower(colnames(coordinates))
  x_name <- colnames(coordinates)[which(cl %in% c("x","col","column"))[1]]
  y_name <- colnames(coordinates)[which(cl %in% c("y","row"))[1]]
  if (is.na(x_name) || is.na(y_name)) stop("Failed to detect spatial coordinate columns.")
  
  # ensure cell_id
  if (!"cell_id" %in% colnames(coordinates)) coordinates$cell_id <- rownames(coordinates)
  coordinates <- coordinates[, c("cell_id", x_name, y_name)]
  colnames(coordinates) <- c("cell_id","Col","Row")
  
  # base plotting df (flip y for image-like orientation)
  full_data <- data.frame(
    cell_id = coordinates$cell_id,
    Col     = coordinates$Col,
    Row     = -coordinates$Row,
    Labels  = as.character(steam_obj$labels),
    stringsAsFactors = FALSE,
    row.names = coordinates$cell_id
  )
  
  # Get corrections based on analysis type
  has_iterative <- !is.null(steam_obj$iterative_anchor_analysis)
  has_single_pass <- !is.null(steam_obj$spatial_anchor_analysis)
  
  if (has_iterative && (is.null(iteration) || iteration > 0)) {
    # Iterative analysis corrections
    iter_analysis <- steam_obj$iterative_anchor_analysis
    
    if (is.null(iteration)) {
      # Show all corrections from all iterations
      if (iter_analysis$method == "iterative_with_ground_truth") {
        all_corrections <- data.frame()
        for (i in seq_along(iter_analysis$iteration_results)) {
          if (!is.null(iter_analysis$iteration_results[[i]]$corrections)) {
            iter_corr <- iter_analysis$iteration_results[[i]]$corrections
            iter_corr$iteration <- i
            all_corrections <- rbind(all_corrections, iter_corr)
          }
        }
      } else {
        all_corrections <- data.frame()
        for (i in seq_along(iter_analysis$correction_history)) {
          if (!is.null(iter_analysis$correction_history[[i]])) {
            iter_corr <- iter_analysis$correction_history[[i]]
            iter_corr$iteration <- i
            all_corrections <- rbind(all_corrections, iter_corr)
          }
        }
      }
      title_suffix <- "(All Iterations)"
    } else {
      # Specific iteration
      if (iter_analysis$method == "iterative_with_ground_truth") {
        if (iteration <= length(iter_analysis$iteration_results) && 
            !is.null(iter_analysis$iteration_results[[iteration]]$corrections)) {
          all_corrections <- iter_analysis$iteration_results[[iteration]]$corrections
          all_corrections$iteration <- iteration
        } else {
          all_corrections <- data.frame()
        }
      } else {
        if (iteration <= length(iter_analysis$correction_history) && 
            !is.null(iter_analysis$correction_history[[iteration]])) {
          all_corrections <- iter_analysis$correction_history[[iteration]]
          all_corrections$iteration <- iteration
        } else {
          all_corrections <- data.frame()
        }
      }
      title_suffix <- paste0("(Iteration ", iteration, ")")
    }
    
  } else if (has_single_pass || iteration == 0) {
    # Single-pass analysis corrections
    if (!has_single_pass) {
      stop("No single-pass results available. Use iteration = NULL or > 0 for iterative results.")
    }
    
    # Extract corrections using helper function
    all_corrections <- extractAnchorCorrections(steam_obj$spatial_anchor_analysis, steam_obj)
    title_suffix <- "(Single-pass)"
    
  } else {
    stop("No analysis results found. Run STEAM_anchor first.")
  }
  
  # Mark corrected cells
  if (nrow(all_corrections) > 0) {
    
    if (identical(fold, "all")) {
      # Show all folds together
      if ("fold" %in% colnames(all_corrections)) {
        # Mark with fold information
        for (f in unique(all_corrections$fold)) {
          fold_corrections <- all_corrections[all_corrections$fold == f, ]
          corr_ids <- fold_corrections$cell_id
          corr_tags <- paste0("Corrected (Fold ", f, ")")
          full_data$Labels[match(corr_ids, full_data$cell_id)] <- corr_tags
        }
      } else {
        # No fold info, just mark as corrected
        corr_ids <- all_corrections$cell_id
        full_data$Labels[match(corr_ids, full_data$cell_id)] <- "Corrected"
      }
      
    } else {
      # Specific fold
      if (!is.numeric(fold)) {
        stop("Invalid 'fold'. Use 'all' or a valid fold number.")
      }
      
      if ("fold" %in% colnames(all_corrections)) {
        fold_corrections <- all_corrections[all_corrections$fold == fold, ]
      } else {
        # If no fold column, assume all corrections are for the requested fold
        fold_corrections <- all_corrections
      }
      
      if (nrow(fold_corrections) > 0) {
        corr_ids <- fold_corrections$cell_id
        full_data$Labels[match(corr_ids, full_data$cell_id)] <- paste0("Corrected (Fold ", fold, ")")
      }
    }
  }
  
  # ---- colors (matching plot_misclassified_cells style) ----
  labs_all <- unique(full_data$Labels)
  base_labs <- setdiff(labs_all, grep("^Corrected", labs_all, value = TRUE))
  corr_labs <- grep("^Corrected", labs_all, value = TRUE)
  
  if (is.null(label_colors)) {
    base_cols <- setNames(scales::hue_pal()(length(base_labs)), base_labs)
  } else {
    base_cols <- label_colors
    missing_base <- setdiff(base_labs, names(base_cols))
    if (length(missing_base)) {
      base_cols <- c(base_cols, setNames(scales::hue_pal()(length(missing_base)), missing_base))
    }
  }
  
  # corrected fold colors (distinct but complementary to misclassified colors)
  if (length(corr_labs)) {
    corr_cols <- setNames(scales::hue_pal(h = c(120, 240), l = 45, c = 90)(length(corr_labs)), corr_labs)
    cols <- c(base_cols, corr_cols)
  } else {
    cols <- base_cols
  }
  
  # ---- plot ----
  title_text <- if (identical(fold, "all")) {
    paste("Corrected Cells", title_suffix, "(colored by fold)")
  } else {
    paste("Corrected Cells", title_suffix, paste0("(Fold ", fold, ")"))
  }
  
  ggplot(full_data, aes(x = Col, y = Row, color = Labels)) +
    geom_point(size = 3) +
    scale_color_manual(values = cols) +
    labs(
      title = title_text,
      color = "Layer / Status"
    ) +
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

#' Compare misclassified vs corrected cells side by side
#'
#' @description
#' Creates a side-by-side comparison of misclassified and corrected cells
#'
#' @param steam_obj STEAM object with analysis results
#' @param label_colors Named vector of colors for consistent coloring
#' @param fold Which fold to show ("all" or specific fold number)
#' @param iteration Which iteration to show for corrections (NULL for all)
#' @return Combined plot object (requires patchwork)
#' @export
plotBeforeAfterComparison <- function(steam_obj, label_colors = NULL, fold = "all", iteration = NULL) {
  
  # Create misclassified plot
  p_misclassified <- plotMisclassifiedCells(steam_obj, label_colors = label_colors, fold = fold)
  
  # Create corrected plot
  p_corrected <- plotCorrectedCells(steam_obj, label_colors = label_colors, fold = fold, iteration = iteration)
  
  # Combine plots
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- p_misclassified | p_corrected
    combined <- combined + patchwork::plot_annotation(
      title = "Before vs After STEAM Correction",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    )
    return(combined)
  } else {
    warning("patchwork package not available. Returning list of plots.")
    return(list(misclassified = p_misclassified, corrected = p_corrected))
  }
}

#' Plot misclassified cells (copied from user's code for consistency)
#'
#' @description
#' Shows misclassified cells in spatial context
#'
#' @param steam_obj STEAM object with predictions
#' @param label_colors Named vector of colors for cell types
#' @param fold Which fold to show ("all" or specific fold number)
#' @return ggplot object
#' @export
plotMisclassifiedCells <- function(steam_obj, label_colors = NULL, fold = "all") {
  suppressPackageStartupMessages({
    library(ggplot2); library(scales)
  })
  
  # ---- coords ----
  coordinates <- steam_obj$spatial
  if (is.null(coordinates)) stop("No spatial coordinates found in STEAM object")
  coordinates <- as.data.frame(coordinates)
  
  # detect coordinate columns
  cl <- tolower(colnames(coordinates))
  x_name <- colnames(coordinates)[which(cl %in% c("x","col","column"))[1]]
  y_name <- colnames(coordinates)[which(cl %in% c("y","row"))[1]]
  if (is.na(x_name) || is.na(y_name)) stop("Failed to detect spatial coordinate columns.")
  
  # ensure cell_id
  if (!"cell_id" %in% colnames(coordinates)) coordinates$cell_id <- rownames(coordinates)
  coordinates <- coordinates[, c("cell_id", x_name, y_name)]
  colnames(coordinates) <- c("cell_id","Col","Row")
  
  # base plotting df (flip y for image-like orientation)
  full_data <- data.frame(
    cell_id = coordinates$cell_id,
    Col     = coordinates$Col,
    Row     = -coordinates$Row,
    Labels  = as.character(steam_obj$labels),
    stringsAsFactors = FALSE,
    row.names = coordinates$cell_id
  )
  
  # ---- SIMPLE CV path ----
  if (!is.null(steam_obj$test) && !is.null(steam_obj$test$predictions)) {
    test_coords <- steam_obj$test$test.data.coords
    test_lbls  <- as.character(steam_obj$test$test.data.labels)
    preds      <- as.character(steam_obj$test$predictions)
    
    mis <- test_lbls != preds
    test_ids <- rownames(as.data.frame(test_coords))
    idx <- match(test_ids, full_data$cell_id)
    # mark misclassified
    full_data$Labels[idx[mis]] <- "Misclassified"
    
    # ---- NESTED CV path ----
  } else if (!is.null(steam_obj$nested) && !is.null(steam_obj$nested$ncv$outer_result)) {
    
    if (identical(fold, "all")) {
      all_preds <- do.call(rbind, lapply(seq_along(steam_obj$nested$ncv$outer_result), function(i) {
        p <- steam_obj$nested$ncv$outer_result[[i]]$preds
        if (is.null(p)) return(NULL)
        p$outer_fold <- i
        p
      }))
      if (is.null(all_preds) || nrow(all_preds) == 0) stop("No nested CV predictions found.")
      
      
      if (!"cell_id" %in% colnames(all_preds)) {
        if (!is.null(rownames(all_preds))) {
          all_preds$cell_id <- rownames(all_preds)
        } else {
          stop("Predictions don't have rownames or 'cell_id' to match spatial coordinates.")
        }
      }
      
      # build per-fold mislabel strings, e.g. "Misclassified (Fold 3)"
      mis <- all_preds$testy != all_preds$predy
      if (any(mis)) {
        mis_ids   <- all_preds$cell_id[mis]
        mis_folds <- all_preds$outer_fold[mis]
        mis_tags  <- paste0("Misclassified (Fold ", mis_folds, ")")
        full_data$Labels[match(mis_ids, full_data$cell_id)] <- mis_tags
      }
      
    } else {
      # Specific outer fold
      if (!is.numeric(fold) || fold < 1 || fold > length(steam_obj$nested$ncv$outer_result)) {
        stop("Invalid 'fold'. Use 'all' or a valid outer fold index.")
      }
      p <- steam_obj$nested$ncv$outer_result[[fold]]$preds
      if (is.null(p) || nrow(p) == 0) stop(sprintf("No predictions for outer fold %s.", fold))
      if (!"cell_id" %in% colnames(p)) {
        if (!is.null(rownames(p))) p$cell_id <- rownames(p) else {
          stop("Fold predictions lack rownames and 'cell_id'; cannot align to spatial coords.")
        }
      }
      mis_ids <- p$cell_id[p$testy != p$predy]
      if (length(mis_ids)) {
        full_data$Labels[match(mis_ids, full_data$cell_id)] <- paste0("Misclassified (Fold ", fold, ")")
      }
    }
    
  } else {
    stop("No predictions found. Run model.predict(mode='simple' or 'nested') first.")
  }
  
  # ---- colors ----
  # Base labels keep their own palette; misclassified folds get distinct colors.
  labs_all <- unique(full_data$Labels)
  base_labs <- setdiff(labs_all, grep("^Misclassified", labs_all, value = TRUE))
  mis_labs  <- grep("^Misclassified", labs_all, value = TRUE)
  
  if (is.null(label_colors)) {
    base_cols <- setNames(scales::hue_pal()(length(base_labs)), base_labs)
  } else {
    base_cols <- label_colors
    missing_base <- setdiff(base_labs, names(base_cols))
    if (length(missing_base)) {
      base_cols <- c(base_cols, setNames(scales::hue_pal()(length(missing_base)), missing_base))
    }
  }
  
  # misclassified fold colors (distinct hues to see fold source)
  if (length(mis_labs)) {
    mis_cols <- setNames(scales::hue_pal(h = c(0, 360), l = 35, c = 100)(length(mis_labs)), mis_labs)
    cols <- c(base_cols, mis_cols)
  } else {
    cols <- base_cols
  }
  
  # ---- plot ----
  ggplot(full_data, aes(x = Col, y = Row, color = Labels)) +
    geom_point(size = 3) +
    scale_color_manual(values = cols) +
    labs(
      title = if (!is.null(steam_obj$nested) && identical(fold, "all"))
        "Misclassified Cells (colored by outer fold)"
      else if (!is.null(steam_obj$nested) && is.numeric(fold))
        paste0("Misclassified Cells (Fold ", fold, ")")
      else "Misclassified Cells",
      color = "Layer / Status"
    ) +
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

#' Plot accuracy improvement across folds
#' @param results STEAM_anchor results
#' @return ggplot object
#' @export
plotAccuracyImprovement <- function(results) {
  
  eval_data <- results$eval
  if (is.null(eval_data)) {
    stop("No evaluation data found in results")
  }
  
  # Create plot
  p <- ggplot2::ggplot(eval_data, ggplot2::aes(x = factor(Outer_Fold))) +
    ggplot2::geom_col(ggplot2::aes(y = Acc_Before), alpha = 0.5, fill = "lightblue", width = 0.4, 
                      position = ggplot2::position_nudge(x = -0.2)) +
    ggplot2::geom_col(ggplot2::aes(y = Acc_After), alpha = 0.7, fill = "darkblue", width = 0.4, 
                      position = ggplot2::position_nudge(x = 0.2)) +
    ggplot2::geom_text(ggplot2::aes(y = Acc_After + 0.01, label = paste("+", round(Delta, 3))), 
                       vjust = 0, size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Accuracy Improvement by Fold",
      x = "CV Fold",
      y = "Accuracy"
    ) +
    ggplot2::ylim(0, max(eval_data$Acc_After) + 0.05)
  
  return(p)
}

#' Plot confidence score distributions
#' @param enhanced_results STEAM_anchor_enhanced results
#' @return ggplot object  
#' @export
plotConfidenceDistributions <- function(enhanced_results) {
  
  confidence_data <- enhanced_results$spatial_anchor_analysis$multimodal_confidence
  if (is.null(confidence_data)) {
    stop("No confidence data found")
  }
  
  # Get confidence columns
  conf_cols <- c("multimodal_confidence", "spatial_confidence", "expression_consistency")
  conf_cols <- conf_cols[conf_cols %in% colnames(confidence_data)]
  
  if (length(conf_cols) == 0) {
    stop("No confidence columns found in data")
  }
  
  # Create long format data
  plot_data <- data.frame()
  for (col in conf_cols) {
    plot_data <- rbind(plot_data, data.frame(
      metric = col,
      score = confidence_data[[col]]
    ))
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = score, fill = metric)) +
    ggplot2::geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
    ggplot2::facet_wrap(~ metric, scales = "free_y", 
                        labeller = ggplot2::labeller(metric = function(x) tools::toTitleCase(gsub("_", " ", x)))) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Confidence Score Distributions",
      x = "Score",
      y = "Frequency"
    ) +
    ggplot2::theme(legend.position = "none")
  
  return(p)
}

#' Plot performance comparison summary
#' @param results STEAM_anchor results
#' @return ggplot object
#' @export
plotPerformanceComparison <- function(results) {
  
  eval_data <- results$eval
  if (is.null(eval_data)) {
    stop("No evaluation data found")
  }
  
  # Create summary metrics
  summary_data <- data.frame(
    Metric = c("Mean Accuracy Before", "Mean Accuracy After", "Mean Improvement"),
    Value = c(mean(eval_data$Acc_Before), mean(eval_data$Acc_After), mean(eval_data$Delta))
  )
  
  p <- ggplot2::ggplot(summary_data, ggplot2::aes(x = Metric, y = Value)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = round(Value, 3)), vjust = -0.5) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = "STEAM Performance Summary",
      x = "",
      y = "Value"
    )
  
  return(p)
}

# Print message when enhanced plotting functions are loaded
message("\n=== Enhanced STEAM Plotting Functions Loaded ===")
message("New functions available:")
message("- plotSpatialOverview() - Tissue architecture overview")
message("- plotConfidenceSpatial() - Confidence scores in spatial context")
message("- plotExpressionGradients() - Expression gradient visualization")
message("- plotExpressionBoundaries() - Expression boundary detection")
message("- plotGeneExpressionSpatial() - Gene-specific expression patterns") 
message("- plotAccuracyImprovement() - Performance across folds")
message("- plotConfidenceDistributions() - Score distributions")
message("- plotPerformanceComparison() - Overall performance metrics")
message("==================================================")
