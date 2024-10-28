neighborhood_averaging <- function(matrix, coordinates, labels, n_size) {
  
  labs <- levels(factor(labels))
  avg_mat <- matrix
  
  for (i in 1:length(labs)) {
    mat <- matrix[, which(labels == labs[i])]
    wc <- coordinates[which(labels == labs[i]),]
    
    for (j in 1:ncol(mat)) {
      roi <- wc[j, 2]
      coi <- wc[j, 3]
      
      allrows <- wc[, 2]
      allcols <- wc[, 3]
      
      neighs <- which((allrows %in% c((roi - n_size):(roi + n_size))) & 
                        (allcols %in% c((coi - n_size):(coi + n_size))))
      
      if (length(neighs) < 2) {
        next
      }
      
      newj <- rowMeans(mat[, neighs])

      avg_mat[, colnames(mat)[j]] <- newj
    }
  }
  
  message("Finished neighborhood averaging")
  
  return(avg_mat)
}