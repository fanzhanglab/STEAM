#' Data Splitting for Training and Testing Sets
#'
#' This function splits the data in a STEAM object into training and testing sets based on a specified ratio. 
#'
#' @param STEAM.obj Object of class STEAM.Object
#' @param train.ratio A numeric value between 0 and 1 specifying the proportion of the data to include in the training set.
#' @importFrom caret createDataPartition
#' @export
data.split <- function(STEAM.obj, train.ratio){

  ind <- createDataPartition(STEAM.obj$labels, p = train.ratio, list = FALSE)

  STEAM.obj$train$train.data.matrix <- STEAM.obj$count_exp[, ind]
  STEAM.obj$train$train.data.coords <- STEAM.obj$spatial[ind, ]
  STEAM.obj$train$train.data.labels <- STEAM.obj$labels[ind]

  STEAM.obj$test$test.data.matrix <- STEAM.obj$count_exp[, -ind]
  STEAM.obj$test$test.data.coords <- STEAM.obj$spatial[-ind,]
  STEAM.obj$test$test.data.labels <- STEAM.obj$labels[-ind]

  return(STEAM.obj)

}
