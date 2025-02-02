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
