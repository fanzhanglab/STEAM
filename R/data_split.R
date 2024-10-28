data_split(matrix, coordinates, labels, split_ratio = 0.7){
  
  ind <- sample(2, ncol(matrix), replace = TRUE, prob = c(split_ratio, 1 - split_ratio))
  
  train_matrix <- matrix[,ind == 1]
  train_coords <- coordinates[ind == 1,]
  train_labels <- labels[ind == 1]
  
  test_matrix <- matrix[,ind == 2]
  test_coords <- coordinates[ind == 2,]
  test_labels <- labels[ind == 2]
  
  return(list(train = list(matrix = train_matrix, coords = train_coords, labels = train_labels),
              test = list(matrix = test_matrix, coords = test_coords, labels = test_labels)))
  
}