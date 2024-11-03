RunSTEAM <- function(matrix, coordinates, labels, train_ratio, n_size, seed, n_tree = 100, model = "rf", kernel = "radial") {
  data <- data_split(matrix, coordinates, labels, train_ratio, seed)
  
  train_matrix <- data$train$matrix
  test_matrix <- data$test$matrix
  train_labels <- data$train$labels
  test_labels <- data$test$labels
  train_coords <- data$train$coordinates
  test_coords <- data$test$coordinates
  
  avg_train_matrix <- neighborhood_averaging(train_matrix, train_coords, train_labels, n_size, seed)
  avg_test_matrix <- neighborhood_averaging(test_matrix, test_coords, test_labels, n_size, seed)
  
  message("Finished neighborhood averaging")
  
  model <- train_model(avg_train_matrix, train_labels, model = model, seed = seed, n_tree = n_tree, kernel = kernel)
  
  results <- predict_model(model, avg_test_matrix, test_labels)
  
  return(results)
}


#results <- RunSTEAM(matrix, coordinates, labels, train_ratio, n_size, seed, n_tree = 100, model = "rf", kernel = "radial")