#' Perform Neighborhood Averaging
#'
#' This function performs neighborhood averaging on the expression data in a STEAM object. 
#' For each cell, it computes the average expression of neighboring cells within a specified spatial distance (\code{n.size}).
#'
#' @param STEAM.obj 
#' @param n.size An integer specifying the neighborhood size (distance) for averaging.
#' @param is_train A logical value indicating whether to apply neighborhood averaging on the training dataset (\code{TRUE}) or the test dataset (\code{FALSE}).
#' @export
neighborhood.avg <- function(STEAM.obj, n.size, is_train) {

  if (is_train) {
    avg_mat <- STEAM.obj$train$train.data.matrix
    rn <- levels(factor(STEAM.obj$train$train.data.labels))

    for (i in seq_along(rn)) {
      mat <- STEAM.obj$train$train.data.matrix[, STEAM.obj$train$train.data.labels == rn[i]]
      wc <- STEAM.obj$train$train.data.coords[STEAM.obj$train$train.data.labels == rn[i], ]

      for (j in seq_len(ncol(mat))) {
        roi <- wc[j, 2]
        coi <- wc[j, 3]

        neighs <- which((wc[, 2] %in% (roi - n.size):(roi + n.size)) &
                          (wc[, 3] %in% (coi - n.size):(coi + n.size)))

        if (length(neighs) < 2) next

        avg_mat[, colnames(mat)[j]] <- rowMeans(mat[, neighs])
      }
    }
    STEAM.obj$train$avg.matrix <- avg_mat
  } else {
    avg_mat <- STEAM.obj$test$test.data.matrix
    for (j in seq_len(ncol(STEAM.obj$test$test.data.matrix))) {
      roi <- STEAM.obj$test$test.data.coords [j, 2]
      coi <- STEAM.obj$test$test.data.coords [j, 3]

      neighs <- which((STEAM.obj$test$test.data.coords [, 2] %in% (roi - n.size):(roi + n.size)) &
                        (STEAM.obj$test$test.data.coords [, 3] %in% (coi - n.size):(coi + n.size)))

      if (length(neighs) < 2) next

      avg_mat[, j] <- rowMeans(STEAM.obj$test$test.data.matrix[, neighs])
    }
    STEAM.obj$test$avg.matrix <- avg_mat
  }


  return(STEAM.obj)
}
