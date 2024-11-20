#' Training and Evaluating Clustering performance with RunSTEAM()
#'
#' @param matrix gene expression matrix
#' @param coordinates spatial/image coordinates
#' @param labels cell annotations
#' @param train_ratio Fraction to set for training (if 0.8 is given, 80% of data will be split for training. Default is 0.8)
#' @param n_size Neighbourhood
#' @param seed
#' @param n_tree
#' @param model
#' @param kernel
#'
#' @return
#' @export
#'
#' @examples
RunSTEAM <- function(matrix, coordinates, labels, train_ratio, n_size, seed, n_tree = 100, model = "rf", kernel = "radial") {
  data <- STEAM:::data.split(matrix, coordinates, labels, train_ratio, seed)

  train_matrix <- data$train$matrix
  test_matrix <- data$test$matrix
  train_labels <- data$train$labels
  test_labels <- data$test$labels
  train_coords <- data$train$coordinates
  test_coords <- data$test$coordinates

  avg_train_matrix <- STEAM:::neighborhood.avg(train_matrix, train_coords, train_labels, n_size, seed, is_train = TRUE)
  avg_test_matrix <- STEAM:::neighborhood.avg(test_matrix, test_coords, n_size = n_size, seed = seed, is_train = FALSE)

  message("Finished neighborhood averaging")

  model <- model.train(avg_train_matrix, train_labels, model, seed = seed, n_tree = n_tree, kernel = kernel)

  results <- model.predict(model, avg_test_matrix, test_labels)

  return(results)
}

