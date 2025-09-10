#' Internal: Neighborhood Averaging for STEAM Objects
#'
#' Applies local spatial neighborhood averaging to training or test data
#' within a STEAM object. For each cell, expression profiles are replaced
#' with the mean of all neighbors within an \code{n.size} window in X/Y
#' coordinates (if at least two neighbors are found).
#'
#' @param STEAM.obj A STEAM object containing training and/or test data.
#' @param n.size Numeric, neighborhood radius (coordinate units).
#' @param is_train Logical, whether to apply to training (\code{TRUE}) or
#'   test (\code{FALSE}) data.
#'
#' @return The updated STEAM object with an added \code{avg.matrix} in
#'   \code{$train} or \code{$test}.
#' @keywords internal
neighborhood.avg <- function(STEAM.obj, n.size, is_train) {

  # Helper function to find coordinate columns
  get_coord_columns <- function(wc) {
    colnames_lower <- tolower(colnames(wc))

    x_candidates <- c("x", "col", "column")
    y_candidates <- c("y", "row")

    x_col <- colnames(wc)[which(colnames_lower %in% x_candidates)]
    y_col <- colnames(wc)[which(colnames_lower %in% y_candidates)]

    if (length(x_col) != 1 || length(y_col) != 1) {
      stop("Could not uniquely identify x and y coordinate columns.")
    }

    return(list(x = x_col, y = y_col))
  }

  if (is_train) {
    avg_mat <- STEAM.obj$train$train.data.matrix
    rn <- unique(as.character(STEAM.obj$train$train.data.labels))

    if (length(rn) < 2) {
      stop("Training labels contain only one class, classification not possible.")
    }

    for (i in seq_along(rn)) {
      mat <- STEAM.obj$train$train.data.matrix[, STEAM.obj$train$train.data.labels == rn[i], drop = FALSE]
      wc <- as.data.frame(STEAM.obj$train$train.data.coords[STEAM.obj$train$train.data.labels == rn[i], ])

      coords <- get_coord_columns(wc)
      x_col <- coords$x
      y_col <- coords$y

      wc[[x_col]] <- as.numeric(as.character(wc[[x_col]]))
      wc[[y_col]] <- as.numeric(as.character(wc[[y_col]]))

      for (j in seq_len(ncol(mat))) {
        roi <- wc[j, x_col]
        coi <- wc[j, y_col]

        neighs <- which(abs(wc[[x_col]] - roi) <= n.size &
                          abs(wc[[y_col]] - coi) <= n.size)

        if (length(neighs) < 2) next

        avg_mat[, colnames(mat)[j]] <- rowMeans(mat[, neighs, drop = FALSE])
      }
    }
    STEAM.obj$train$avg.matrix <- avg_mat

  } else {
    avg_mat <- STEAM.obj$test$test.data.matrix
    wc <- as.data.frame(STEAM.obj$test$test.data.coords)

    coords <- get_coord_columns(wc)
    x_col <- coords$x
    y_col <- coords$y

    wc[[x_col]] <- as.numeric(as.character(wc[[x_col]]))
    wc[[y_col]] <- as.numeric(as.character(wc[[y_col]]))

    for (j in seq_len(ncol(avg_mat))) {
      roi <- wc[j, x_col]
      coi <- wc[j, y_col]

      neighs <- which(abs(wc[[x_col]] - roi) <= n.size &
                        abs(wc[[y_col]] - coi) <= n.size)

      if (length(neighs) < 2) next

      avg_mat[, j] <- rowMeans(STEAM.obj$test$test.data.matrix[, neighs, drop = FALSE])
    }
    STEAM.obj$test$avg.matrix <- avg_mat
  }

  return(STEAM.obj)
}
