RunSTEAM <- function(matrix, coordinates, labels, train_ratio, n_size, seed, n_tree = 100, model = "rf", kernel = "radial") {
  data <- data.split(matrix, coordinates, labels, train_ratio, seed)

  train_matrix <- data$train$matrix
  test_matrix <- data$test$matrix
  train_labels <- data$train$labels
  test_labels <- data$test$labels
  train_coords <- data$train$coordinates
  test_coords <- data$test$coordinates

  avg_train_matrix <- neighborhood.avg(train_matrix, train_coords, train_labels, n_size, seed, is_train = TRUE)
  avg_test_matrix <- neighborhood.avg(test_matrix, test_coords, n_size = n_size, seed = seed, is_train = FALSE)

  message("Finished neighborhood averaging")

  model <- model.train(avg_train_matrix, train_labels, model = model, seed = seed, n_tree = n_tree, kernel = kernel)

  results <- model.predict(model, avg_test_matrix, test_labels)

  return(results)
}

# library(caret)
# library(randomForest)
# library(e1071)
#
# DLPFCs <- readRDS("DLPFC.RDS")
# gl <- readRDS("DLPFCGeneList.RDS")
# DLPFC <- DLPFCs$`151669`
# labels <- DLPFC$Labs
# coordinates <- DLPFC@images$slice1@coordinates[,c(1:3)]
# matrix <- DLPFC@assays$SCT@scale.data
# seed <- 22
# train_ratio <- 0.7
# n_size <- 3
#
# results <- RunSTEAM(matrix, coordinates, labels, train_ratio, n_size, seed, n_tree = 100, model = "rf", kernel = "radial")
