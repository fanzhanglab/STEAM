data.split <- function(matrix, coordinates, labels, train_ratio, seed){
  set.seed(seed)
  ind <- createDataPartition(labels, p = train_ratio, list = FALSE)

  train_matrix <- matrix[, ind]
  train_coords <- coordinates[ind, ]
  train_labels <- labels[ind]

  test_matrix <- matrix[, -ind]
  test_coords <- coordinates[-ind,]
  test_labels <- labels[-ind]

  return(list(train = list(matrix = train_matrix, coordinates = train_coords, labels = train_labels),
              test = list(matrix = test_matrix, coordinates = test_coords, labels = test_labels)))

}
