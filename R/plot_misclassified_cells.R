#' Plot Misclassified Cells
#'
#' Visualize spatial positions of misclassified cells from STEAM predictions,
#' either from a simple train/test split or nested cross-validation.
#'
#' @param steam_obj A STEAM object with spatial coordinates and predictions.
#' @param label_colors Optional named vector of colors for labels.
#' @param fold Either "all" (show misclassifications from all outer folds)
#'   or an integer specifying one outer fold.
#'
#' @return A ggplot2 plot object.
#' @examples
#' \dontrun{
#'   plot_misclassified_cells(STEAM.obj)
#'   plot_misclassified_cells(STEAM.obj, fold = 1)
#' }
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_classic theme element_blank
#' @importFrom scales hue_pal
#' @export
plot_misclassified_cells <- function(steam_obj, label_colors = NULL, fold = "all") {

  # ---- coords ----
  coordinates <- steam_obj$spatial
  if (is.null(coordinates)) stop("No spatial coordinates found in STEAM object")
  coordinates <- as.data.frame(coordinates)

  # detect coordinate columns
  cl <- tolower(colnames(coordinates))
  x_name <- colnames(coordinates)[which(cl %in% c("x", "col", "column"))[1]]
  y_name <- colnames(coordinates)[which(cl %in% c("y", "row"))[1]]
  if (is.na(x_name) || is.na(y_name)) stop("Failed to detect spatial coordinate columns.")

  # ensure cell_id
  if (!"cell_id" %in% colnames(coordinates)) coordinates$cell_id <- rownames(coordinates)
  coordinates <- coordinates[, c("cell_id", x_name, y_name)]
  colnames(coordinates) <- c("cell_id", "Col", "Row")

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
    test_lbls   <- as.character(steam_obj$test$test.data.labels)
    preds       <- as.character(steam_obj$test$predictions)

    mis <- test_lbls != preds
    test_ids <- rownames(as.data.frame(test_coords))
    idx <- match(test_ids, full_data$cell_id)
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

      mis <- all_preds$testy != all_preds$predy
      if (any(mis)) {
        mis_ids   <- all_preds$cell_id[mis]
        mis_folds <- all_preds$outer_fold[mis]
        mis_tags  <- paste0("Misclassified (Fold ", mis_folds, ")")
        full_data$Labels[match(mis_ids, full_data$cell_id)] <- mis_tags
      }

    } else {
      if (!is.numeric(fold) || fold < 1 || fold > length(steam_obj$nested$ncv$outer_result)) {
        stop("Invalid 'fold'. Use 'all' or a valid outer fold index.")
      }
      p <- steam_obj$nested$ncv$outer_result[[fold]]$preds
      if (is.null(p) || nrow(p) == 0) stop(sprintf("No predictions for outer fold %s.", fold))
      if (!"cell_id" %in% colnames(p)) {
        if (!is.null(rownames(p))) {
          p$cell_id <- rownames(p)
        } else {
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
  labs_all <- unique(full_data$Labels)
  base_labs <- setdiff(labs_all, grep("^Misclassified", labs_all, value = TRUE))
  mis_labs  <- grep("^Misclassified", labs_all, value = TRUE)

  if (is.null(label_colors)) {
    all_possible_labels <- sort(unique(as.character(steam_obj$labels)))
    all_label_colors <- setNames(scales::hue_pal()(length(all_possible_labels)), all_possible_labels)
    base_cols <- all_label_colors[names(all_label_colors) %in% base_labs]
  } else {
    base_cols <- label_colors
    missing_base <- setdiff(base_labs, names(base_cols))
    if (length(missing_base)) {
      all_possible_labels <- sort(unique(as.character(steam_obj$labels)))
      all_label_colors <- setNames(scales::hue_pal()(length(all_possible_labels)), all_possible_labels)
      base_cols <- c(base_cols, all_label_colors[missing_base])
    }
  }

  if (length(mis_labs)) {
    if (identical(fold, "all")) {
      mis_cols <- setNames(scales::hue_pal(h = c(0, 360), l = 35, c = 100)(length(mis_labs)), mis_labs)
    } else {
      mis_cols <- setNames(rep("black", length(mis_labs)), mis_labs)
    }
    cols <- c(base_cols, mis_cols)
  } else {
    cols <- base_cols
  }

  # ---- plot ----
  p <- ggplot2::ggplot(full_data, ggplot2::aes(x = Col, y = Row, color = Labels)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(
      title = if (!is.null(steam_obj$nested) && identical(fold, "all"))
        "Misclassified Cells (colored by outer fold) - Spatial CV"
      else if (!is.null(steam_obj$nested) && is.numeric(fold))
        paste0("Misclassified Cells (Fold ", fold, ") - Spatial CV")
      else if (!is.null(steam_obj$test))
        "Misclassified Cells - Spatial CV"
      else "Misclassified Cells",
      color = "Layer / Status"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )

  return(p)
}
