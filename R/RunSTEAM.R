#' Training and Evaluating Clustering performance with RunSTEAM()
#'
#'
#' @param STEAM.obj Object of class STEAM.Object from LoadSTEAM()
#' @param seed Seed Value. Default is 123
#' @param model Currently available methods - rf (for random forest) , svm (for Support Vector Machine), xgb (for Extreme Gradient Boosting), multinom (for Penalized Multinomial Regression)
#' @param kernel Default is Linear (argument valid only for svm)
#' @param train.ratio training set ratio. Default is 0.8 (Splits 80% data for training and 20% for testing)
#' @param n.size Neighborhood size for neighborhood averaging
#' @param cv.folds #folds for Cross-validation
#' @param cv.repeats #repeats for repeated cross validation
#' @param trainval.ratio training set ratio for cross-validation. Default is 0.8 (Splits 80% data for training and 20% for validation)
#' @param n.tree #trees (argument valid only for Random Forest). Deafult is 500
#' @param train.folder.name Training Output folder name
#' @param allowParallel Boolean to perform parallel processing, if resources are available. Default is FALSE
#'
#' @return S4 Object of class 'STEAM.Object' with all outputs
#'
#'@import caret
#'randomForest
#'e1071
#'
#' @export
RunSTEAM <- function(STEAM.obj, train.ratio = 0.8, n.size = 5, seed = 123, cv.folds = 10, cv.repeats = 3, trainval.ratio = 0.8, model = "rf", n.tree = 500, kernel = 'linear', train.folder.name = 'train.out', allowParallel = FALSE) {
  set.seed(seed)
  STEAM.obj <- data.split(STEAM.obj, train.ratio = train.ratio)
  message("Finished Data Splitting")

  STEAM.obj <- neighborhood.avg(STEAM.obj, n.size = n.size, is_train = TRUE)
  STEAM.obj <- neighborhood.avg(STEAM.obj, n.size = n.size, is_train = FALSE)
  message("Finished neighborhood averaging")


  STEAM.obj <- model.train(STEAM.obj, model, n.tree = n.tree, kernel = kernel, cv.folds = cv.folds, cv.repeats = cv.repeats, trainval.ratio = trainval.ratio, train.folder.name = train.folder.name, allowParallel = allowParallel)
  message("Finished Model training")

  STEAM.obj <- model.predict(STEAM.obj)
  message("Finished Evaluation")

  return(STEAM.obj)
}

