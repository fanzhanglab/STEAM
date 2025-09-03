#!/usr/bin/env Rscript

#' STEAM Iterative Analysis Dashboard
#' 
#' Creates a cleaner, less confusing dashboard with better visual hierarchy
#' and clearer explanations of what each plot shows.

library(ggplot2)
library(dplyr)

#' Create Clearer Progress Plot
#' 
#' Shows progress with separate plots for per-iteration and cumulative
createClearProgressPlot <- function(STEAM.obj) {
  
  if (is.null(STEAM.obj$iterative_anchor_analysis)) {
    stop("No iterative_anchor_analysis results found.")
  }
  
  iter_analysis <- STEAM.obj$iterative_anchor_analysis
  
  # Extract iteration data
  iteration_data <- data.frame()
  
  if (iter_analysis$method == "iterative_with_ground_truth") {
    for (i in seq_along(iter_analysis$iteration_results)) {
      corrections <- if (!is.null(iter_analysis$iteration_results[[i]]$corrections)) {
        nrow(iter_analysis$iteration_results[[i]]$corrections)
      } else {
        iter_analysis$iteration_results[[i]]$total_corrections %||% 0
      }
      
      iteration_data <- rbind(iteration_data, data.frame(
        iteration = i,
        corrections = corrections
      ))
    }
  } else {
    for (i in seq_along(iter_analysis$correction_history)) {
      corrections <- if (!is.null(iter_analysis$correction_history[[i]])) {
        nrow(iter_analysis$correction_history[[i]])
      } else {
        0
      }
      
      iteration_data <- rbind(iteration_data, data.frame(
        iteration = i,
        corrections = corrections
      ))
    }
  }
  
  # Add cumulative column
  iteration_data$cumulative <- cumsum(iteration_data$corrections)
  
  # Create cleaner plot - just show corrections per iteration
  p <- ggplot(iteration_data, aes(x = factor(iteration), y = corrections)) +
    geom_col(fill = "#3498DB", alpha = 0.8, width = 0.7) +
    geom_text(aes(label = corrections), vjust = -0.3, size = 4, fontface = "bold") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray60"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = "Corrections Made Per Iteration",
      subtitle = sprintf("Total: %d corrections across %d iterations", 
                        max(iteration_data$cumulative), 
                        nrow(iteration_data)),
      x = "Iteration",
      y = "Number of Corrections"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  return(p)
}

#' Create Accuracy Summary Plot (Simple Before/After)
#' 
#' Shows just the key result - accuracy improvement
createAccuracySummaryPlot <- function(STEAM.obj) {
  
  if (is.null(STEAM.obj$iterative_anchor_analysis)) {
    return(NULL)
  }
  
  iter_analysis <- STEAM.obj$iterative_anchor_analysis
  
  # Only works with ground truth method
  if (iter_analysis$method != "iterative_with_ground_truth") {
    return(createCoherenceSummaryPlot(STEAM.obj))
  }
  
  # Get number of folds
  n_folds <- length(STEAM.obj$nested$ncv$outer_result)
  
  # Calculate baseline and final accuracy
  baseline_accuracy <- numeric(n_folds)
  final_accuracy <- numeric(n_folds)
  
  for (f in 1:n_folds) {
    fold_result <- STEAM.obj$nested$ncv$outer_result[[f]]
    if (!is.null(fold_result$preds)) {
      preds_df <- fold_result$preds
      baseline_accuracy[f] <- mean(preds_df$testy == preds_df$predy, na.rm = TRUE)
      
      # Apply all corrections to get final accuracy
      current_preds <- setNames(preds_df$predy, rownames(preds_df))
      true_labels <- setNames(preds_df$testy, rownames(preds_df))
      
      # Get all corrections
      all_corrections <- NULL
      if (!is.null(STEAM.obj$spatial_anchor_analysis$corrections)) {
        all_corrections <- STEAM.obj$spatial_anchor_analysis$corrections
      }
      
      if (!is.null(all_corrections)) {
        fold_corrections <- all_corrections[all_corrections$fold == f, ]
        if (nrow(fold_corrections) > 0) {
          current_preds[fold_corrections$cell_id] <- fold_corrections$suggested_correction
        }
      }
      
      final_accuracy[f] <- mean(current_preds == true_labels, na.rm = TRUE)
    }
  }
  
  # Create simple before/after plot
  plot_data <- data.frame(
    Accuracy = c(baseline_accuracy, final_accuracy),
    Stage = factor(rep(c("Before", "After"), each = n_folds),
                   levels = c("Before", "After")),
    Fold = rep(1:n_folds, 2)
  )
  
  improvement <- mean(final_accuracy - baseline_accuracy)
  
  p <- ggplot(plot_data, aes(x = Stage, y = Accuracy)) +
    geom_boxplot(aes(fill = Stage), alpha = 0.7, outlier.shape = NA) +
    geom_point(aes(color = factor(Fold)), size = 2.5, 
               position = position_jitter(width = 0.2)) +
    scale_fill_manual(values = c("Before" = "#E74C3C", "After" = "#27AE60")) +
    scale_color_discrete(name = "Fold") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray60"),
      axis.title = element_text(size = 12),
      legend.position = "none"
    ) +
    labs(
      title = "Overall Accuracy Improvement",
      subtitle = sprintf("Mean improvement: %+.4f (%+.2f%%)", 
                        improvement, improvement * 100),
      x = "",
      y = "Accuracy"
    )
  
  return(p)
}

#' Create Coherence Summary Plot (for no ground truth)
createCoherenceSummaryPlot <- function(STEAM.obj) {
  
  iter_analysis <- STEAM.obj$iterative_anchor_analysis
  
  if (length(iter_analysis$iteration_metrics) < 2) {
    return(NULL)
  }
  
  # Extract coherence data
  coherence_data <- data.frame()
  for (i in seq_along(iter_analysis$iteration_metrics)) {
    metrics <- iter_analysis$iteration_metrics[[i]]
    coherence_data <- rbind(coherence_data, data.frame(
      iteration = i,
      spatial_coherence = metrics$spatial_coherence,
      expression_coherence = metrics$expression_coherence %||% 0.5
    ))
  }
  
  # Show improvement from first to last iteration
  first_spatial <- coherence_data$spatial_coherence[1]
  last_spatial <- coherence_data$spatial_coherence[nrow(coherence_data)]
  spatial_improvement <- last_spatial - first_spatial
  
  plot_data <- data.frame(
    Coherence = c(first_spatial, last_spatial),
    Stage = factor(c("Initial", "Final"), levels = c("Initial", "Final")),
    Type = "Spatial Coherence"
  )
  
  p <- ggplot(plot_data, aes(x = Stage, y = Coherence)) +
    geom_col(aes(fill = Stage), alpha = 0.8, width = 0.6) +
    geom_text(aes(label = sprintf("%.3f", Coherence)), vjust = -0.3, 
              size = 4, fontface = "bold") +
    scale_fill_manual(values = c("Initial" = "#E74C3C", "Final" = "#27AE60")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray60"),
      axis.title = element_text(size = 12),
      legend.position = "none"
    ) +
    labs(
      title = "Spatial Coherence Improvement",
      subtitle = sprintf("Improvement: %+.4f", spatial_improvement),
      x = "",
      y = "Spatial Coherence"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  return(p)
}

#' Create Summary Stats Text Plot
createSummaryStatsPlot <- function(STEAM.obj) {
  
  if (is.null(STEAM.obj$iterative_anchor_analysis)) {
    return(NULL)
  }
  
  iter_analysis <- STEAM.obj$iterative_anchor_analysis
  
  # Calculate key statistics
  total_iterations <- iter_analysis$iterations_completed
  total_corrections <- iter_analysis$summary$total_corrections
  
  # Create text summary
  stats_text <- sprintf(
    "Summary Statistics\n\n• Method: %s\n• Iterations completed: %d\n• Total corrections: %d\n• Average per iteration: %.1f",
    ifelse(iter_analysis$method == "iterative_with_ground_truth", 
           "Ground Truth Validation", "Alternative Validation"),
    total_iterations,
    total_corrections,
    total_corrections / total_iterations
  )
  
  # Create a text plot
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = stats_text, 
             hjust = 0.5, vjust = 0.5, size = 4.5, lineheight = 1.2) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "#F8F9FA", color = "#DEE2E6"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  
  return(p)
}

#' Create Simplified Spatial Overview
createSpatialOverview <- function(STEAM.obj, max_points = 2000) {
  
  # Check for any corrections
  all_corrections <- NULL
  if (!is.null(STEAM.obj$spatial_anchor_analysis$corrections)) {
    all_corrections <- STEAM.obj$spatial_anchor_analysis$corrections
  }
  
  if (is.null(all_corrections) || nrow(all_corrections) == 0) {
    return(NULL)
  }
  
  # Get spatial coordinates
  coords <- STEAM.obj$spatial
  if (is.null(coords)) {
    return(NULL)
  }
  
  # Sample points if too many
  if (nrow(coords) > max_points) {
    sample_indices <- sample(nrow(coords), max_points)
    coords <- coords[sample_indices, ]
  }
  
  # Mark corrected cells
  corrected_cells <- unique(all_corrections$cell_id)
  coords$corrected <- rownames(coords) %in% corrected_cells
  
  # Create simple spatial plot
  p <- ggplot(coords, aes(x = x, y = y)) +
    geom_point(aes(color = corrected), size = 0.8, alpha = 0.7) +
    scale_color_manual(
      values = c("FALSE" = "#BDC3C7", "TRUE" = "#E74C3C"),
      labels = c("Unchanged", "Corrected"),
      name = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      title = "Spatial Distribution of Corrections",
      subtitle = sprintf("%d cells corrected", length(corrected_cells)),
      x = "", y = ""
    ) +
    coord_fixed()
  
  return(p)
}

#' Create Improved Dashboard
#' 
#' @param STEAM.obj STEAM object with iterative results
#' @return Combined plot or list of plots
createImprovedDashboard <- function(STEAM.obj) {
  
  if (is.null(STEAM.obj$iterative_anchor_analysis)) {
    stop("No iterative_anchor_analysis results found.")
  }
  
  # Create individual plots
  p_progress <- createClearProgressPlot(STEAM.obj)
  p_accuracy <- createAccuracySummaryPlot(STEAM.obj)
  p_stats <- createSummaryStatsPlot(STEAM.obj)
  p_spatial <- createSpatialOverview(STEAM.obj)
  
  # Filter out NULL plots
  plots <- list(
    progress = p_progress,
    accuracy = p_accuracy,
    stats = p_stats,
    spatial = p_spatial
  )
  plots <- plots[!sapply(plots, is.null)]
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    # Arrange plots in a clean layout
    if (length(plots) >= 4) {
      # 2x2 layout
      combined <- patchwork::wrap_plots(plotlist = plots, ncol = 2)
    } else if (length(plots) == 3) {
      # Top row: 2 plots, bottom row: 1 plot
      top_row <- patchwork::wrap_plots(plots[1:2], ncol = 2)
      combined <- top_row / plots[[3]]
    } else {
      # Simple horizontal layout
      combined <- patchwork::wrap_plots(plotlist = plots, ncol = length(plots))
    }
    
    combined <- combined + patchwork::plot_annotation(
      title = "STEAM Iterative Analysis - Clear Summary",
      subtitle = "Simplified view of key results and improvements",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, color = "gray60")
      )
    )
    
    return(combined)
  } else {
    warning("patchwork package not available. Returning list of plots.")
    return(plots)
  }
}

#' Create Before/After Comparison Plot
#' 
#' Shows misclassifications before anchor analysis vs same plot with corrections overlaid
#' Uses exact same logic as plot_misclassified_cells function
createBeforeAfterComparison <- function(STEAM.obj) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    cat("patchwork package needed for comparison plots. Install with: install.packages('patchwork')\n")
    return(NULL)
  }
  
  suppressPackageStartupMessages({
    library(ggplot2); library(scales)
  })
  
  # ---- coords (EXACT same logic as plot_misclassified_cells) ----
  coordinates <- STEAM.obj$spatial
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
  
  # base plotting df (flip y for image-like orientation) - EXACT same as original
  full_data <- data.frame(
    cell_id = coordinates$cell_id,
    Col     = coordinates$Col,
    Row     = -coordinates$Row,
    Labels  = as.character(STEAM.obj$labels),
    stringsAsFactors = FALSE,
    row.names = coordinates$cell_id
  )
  
  # ---- Get misclassifications from nested CV (EXACT same logic) ----
  if (!is.null(STEAM.obj$nested) && !is.null(STEAM.obj$nested$ncv$outer_result)) {
    
    all_preds <- do.call(rbind, lapply(seq_along(STEAM.obj$nested$ncv$outer_result), function(i) {
      p <- STEAM.obj$nested$ncv$outer_result[[i]]$preds
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
    misclassified_cells <- c()
    if (any(mis)) {
      mis_ids   <- all_preds$cell_id[mis]
      mis_folds <- all_preds$outer_fold[mis]
      mis_tags  <- paste0("Misclassified (Fold ", mis_folds, ")")
      full_data$Labels[match(mis_ids, full_data$cell_id)] <- mis_tags
      misclassified_cells <- mis_ids
    }
    
  } else {
    stop("No nested CV predictions found. Run model.predict with nested CV first.")
  }
  
  # ---- colors (EXACT same logic as plot_misclassified_cells) ----
  labs_all <- unique(full_data$Labels)
  base_labs <- setdiff(labs_all, grep("^Misclassified", labs_all, value = TRUE))
  mis_labs  <- grep("^Misclassified", labs_all, value = TRUE)
  
  base_cols <- setNames(scales::hue_pal()(length(base_labs)), base_labs)
  
  # misclassified fold colors (distinct hues to see fold source)
  if (length(mis_labs)) {
    mis_cols <- setNames(scales::hue_pal(h = c(0, 360), l = 35, c = 100)(length(mis_labs)), mis_labs)
    cols_before <- c(base_cols, mis_cols)
  } else {
    cols_before <- base_cols
  }
  
  # ---- BEFORE plot (EXACT same as plot_misclassified_cells) ----
  before_plot <- ggplot(full_data, aes(x = Col, y = Row, color = Labels)) +
    geom_point(size = 3) +
    scale_color_manual(values = cols_before) +
    labs(
      title = "BEFORE: Misclassified Cells (colored by outer fold)",
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
  
  # ---- AFTER plot - same structure but add "Corrected" layer ----
  
  # Get corrections
  corrected_cells <- c()
  if (!is.null(STEAM.obj$spatial_anchor_analysis$corrections)) {
    corrections_df <- STEAM.obj$spatial_anchor_analysis$corrections
    corrected_cells <- corrections_df$cell_id
  }
  
  # Create after data - start fresh with original tissue labels (no misclassifications)
  after_data <- data.frame(
    cell_id = coordinates$cell_id,
    Col     = coordinates$Col,
    Row     = -coordinates$Row,
    Labels  = as.character(STEAM.obj$labels),
    stringsAsFactors = FALSE,
    row.names = coordinates$cell_id
  )
  
  # Mark corrected cells as "Corrected" to highlight them
  if (length(corrected_cells) > 0) {
    for (cell_id in corrected_cells) {
      cell_idx <- which(after_data$cell_id == cell_id)
      if (length(cell_idx) > 0) {
        # Mark as corrected to highlight the spatial anchor corrections
        after_data$Labels[cell_idx] <- "Corrected"
      }
    }
  }
  
  # ---- colors for AFTER plot ----
  labs_all_after <- unique(after_data$Labels)
  base_labs_after <- setdiff(labs_all_after, grep("^Corrected", labs_all_after, value = TRUE))
  corrected_labs <- grep("^Corrected", labs_all_after, value = TRUE)
  
  # Use same base colors as BEFORE plot for tissue types
  base_cols_after <- setNames(scales::hue_pal()(length(base_labs_after)), base_labs_after)
  
  # Start with base tissue colors
  cols_after <- base_cols_after
  
  # Add corrected color (bright green)
  if (length(corrected_labs)) {
    corrected_cols <- setNames(c("#00FF00"), corrected_labs)
    cols_after <- c(cols_after, corrected_cols)
  }
  
  # ---- AFTER plot (same style as original function) ----
  after_plot <- ggplot(after_data, aes(x = Col, y = Row, color = Labels)) +
    geom_point(size = 3) +
    scale_color_manual(values = cols_after) +
    labs(
      title = "AFTER: Tissue + Corrections Applied",
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
  
  # Combine plots
  comparison_plot <- before_plot + after_plot + 
    patchwork::plot_annotation(
      title = "Before vs After: Spatial Anchor Analysis Impact",
      subtitle = sprintf("BEFORE: %d misclassified cells | AFTER: %d corrections applied (green)", 
                        length(misclassified_cells), length(corrected_cells)),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray60")
      )
    )
  
  return(comparison_plot)
}

#' Quick Iterative Visualization
#' 
#' @param STEAM.obj STEAM object
#' @param plot_type Type of plot to create ("dashboard", "progress", "accuracy", "stats", "spatial")
#' @return Plot object or dashboard
quickIterativeViz <- function(STEAM.obj, plot_type = "dashboard") {
  
  if (is.null(STEAM.obj$iterative_anchor_analysis)) {
    stop("No iterative analysis results found. Run STEAM_anchor with iterative=TRUE first.")
  }
  
  # Validate plot_type
  valid_types <- c("dashboard", "progress", "accuracy", "stats", "spatial", "comparison")
  if (!plot_type %in% valid_types) {
    stop(sprintf("plot_type must be one of: %s", paste(valid_types, collapse = ", ")))
  }
  
  # Create individual plots based on type
  switch(plot_type,
    "progress" = createClearProgressPlot(STEAM.obj),
    "accuracy" = createAccuracySummaryPlot(STEAM.obj),
    "stats" = createSummaryStatsPlot(STEAM.obj),
    "spatial" = createSpatialOverview(STEAM.obj),
    "comparison" = createBeforeAfterComparison(STEAM.obj),
    "dashboard" = {
      # Show dashboard with instructions
      cat("=== STEAM Iterative Analysis Dashboard ===\n\n")
      
      iter_analysis <- STEAM.obj$iterative_anchor_analysis
      
      cat("Analysis Summary:\n")
      cat(sprintf("• Method: %s\n", 
                 ifelse(iter_analysis$method == "iterative_with_ground_truth",
                       "Ground Truth Validation", "Alternative Validation")))
      cat(sprintf("• Iterations: %d\n", iter_analysis$iterations_completed))
      cat(sprintf("• Total corrections: %d\n", iter_analysis$summary$total_corrections))
      cat("\nGenerating improved dashboard...\n\n")
      
      # Create and return the dashboard
      dashboard <- createImprovedDashboard(STEAM.obj)
      
      cat("Dashboard components:\n")
      cat("1. Progress Plot: Shows corrections made per iteration\n")
      cat("2. Accuracy Plot: Before vs after comparison (if ground truth available)\n")
      cat("3. Statistics Panel: Key numbers and metrics\n")
      cat("4. Spatial Overview: Where corrections were made\n\n")
      
      return(dashboard)
    }
  )
}

#' Legacy function name for backwards compatibility
showClearDashboard <- function(STEAM.obj) {
  quickIterativeViz(STEAM.obj, plot_type = "dashboard")
}
