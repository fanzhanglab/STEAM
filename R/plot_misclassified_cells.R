#' Plotting Misclassified labels
#'
#' @param steam_obj STEAM Object
#' @param colors list of colors for clusters
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_minimal theme element_blank
#' @importFrom scales hue_pal
#' @export
plot_misclassified_cells <- function(steam_obj, label_colors = NULL, fold = "all") {
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
