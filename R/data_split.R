data_split <- function(matrix, coordinates, labels, train_ratio, seed){
  set.seed(seed)
  ind <- sample(2, ncol(matrix), replace = TRUE, prob = c(train_ratio, 1 - train_ratio))
  
  train_matrix <- matrix[, ind == 1]
  train_coords <- coordinates[ind == 1, ]
  train_labels <- labels[ind == 1]
  
  test_matrix <- matrix[, ind == 2]
  test_coords <- coordinates[ind == 2,]
  test_labels <- labels[ind == 2]
  
  return(list(train = list(matrix = train_matrix, coordinates = train_coords, labels = train_labels),
              test = list(matrix = test_matrix, coordinates = test_coords, labels = test_labels)))
  
}