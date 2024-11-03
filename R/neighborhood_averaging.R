neighborhood_averaging <- function(matrix, coordinates, labels = NULL , n_size, seed, is_train = TRUE ) {
  set.seed(seed)
  avg_mat <- matrix
  
  if (is_train) {
    rn <- levels(factor(labels))
  
  for (i in seq_along(rn)) {
    mat <- matrix[, which(labels == rn[i])]
    wc <- coordinates[which(labels == rn[i]),]
    
    for (j in 1:ncol(mat)) {
      roi <- wc[j, 2]
      coi <- wc[j, 3]
      
      allrows <- wc[, 2]
      allcols <- wc[, 3]
      
      neighs <- which((allrows %in% c((roi - n_size):(roi + n_size))) & 
                        (allcols %in% c((coi - n_size):(coi + n_size))))
      
      if (length(neighs) < 2) next
      
      avg_mat[,colnames(mat)[j]] <- rowMeans(mat[, neighs])
    }
  }
} else {
  for (j in 1:ncol(matrix)) {
    roi <- coordinates[j, 2]
    coi <- coordinates[j, 3]
    
    allrows <- coordinates[, 2]
    allcols <- coordinates[, 3]
    
    neighs <- which((allrows %in% c((roi - n_size):(roi + n_size))) & 
                      (allcols %in% c((coi - n_size):(coi + n_size))))
    
    if (length(neighs) < 2) next
    
    avg_mat[, j] <- rowMeans(matrix[, neighs])
    
  }
}

  return(avg_mat)
}
