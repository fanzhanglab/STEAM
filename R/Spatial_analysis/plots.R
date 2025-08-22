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
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  
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