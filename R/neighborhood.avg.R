neighborhood.avg() <- function(matrix, coordinates, labels = NULL, n_size, seed, is_train = TRUE) {
  set.seed(seed)
  avg_mat <- matrix

  if (is_train) {
    rn <- levels(factor(labels))

    for (i in seq_along(rn)) {
      mat <- matrix[, labels == rn[i]]
      wc <- coordinates[labels == rn[i], ]

      for (j in seq_len(ncol(mat))) {
        roi <- wc[j, 2]
        coi <- wc[j, 3]

        neighs <- which((wc[, 2] %in% (roi - n_size):(roi + n_size)) &
                          (wc[, 3] %in% (coi - n_size):(coi + n_size)))

        if (length(neighs) < 2) next

        avg_mat[, colnames(mat)[j]] <- rowMeans(mat[, neighs])
      }
    }
  } else {
    for (j in seq_len(ncol(matrix))) {
      roi <- coordinates[j, 2]
      coi <- coordinates[j, 3]

      neighs <- which((coordinates[, 2] %in% (roi - n_size):(roi + n_size)) &
                        (coordinates[, 3] %in% (coi - n_size):(coi + n_size)))

      if (length(neighs) < 2) next

      avg_mat[, j] <- rowMeans(matrix[, neighs])
    }
  }

  return(avg_mat)
}
