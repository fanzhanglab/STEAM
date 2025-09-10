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
RunSTEAM <- function(
    STEAM.obj,
    mode = c("simple","nested"),
    train.ratio = 0.8,
    n.size = 5,
    seed = 123,
    # simple mode
    cv.folds = 10, cv.repeats = 3, trainval.ratio = 0.8,
    # nested mode
    n_outer_folds = 5, n_inner_folds = 3,
    # model
    model = "rf", n.tree = 500, kernel = "linear",
    train.folder.name = "train.out",
    allowParallel = FALSE,              # inner caret parallel
    metric = "Kappa",
    tune.grid = NULL,
    outer_folds = NULL, inner_folds = NULL,
    verbose = TRUE,
    cv.cores = 1,
    parallel_mode = NULL,
    maxnweights = 10000
) {
  set.seed(seed)
  mode <- match.arg(mode)
  
  if (mode == "simple") {
    STEAM.obj <- data.split(STEAM.obj, train.ratio = train.ratio, seed = seed)
    message("Finished Data Splitting")
    
    STEAM.obj <- neighborhood.avg(STEAM.obj, n.size = n.size, is_train = TRUE)
    STEAM.obj <- neighborhood.avg(STEAM.obj, n.size = n.size, is_train = FALSE)
    message("Finished neighborhood averaging")
    
    STEAM.obj <- model.train(
      STEAM.obj, mode = "simple",
      model = model, kernel = kernel,
      cv.folds = cv.folds, cv.repeats = cv.repeats,
      metric = metric, tune.grid = tune.grid,
      allowParallel = allowParallel, seed = seed, n.size = n.size,
      maxnweights = maxnweights 
    )
    message("Finished Model training")
    
    STEAM.obj <- model.predict(STEAM.obj, mode = "simple")
    message("Finished Evaluation")
    return(STEAM.obj)
  }
  
  # --- nested mode ---
  STEAM.obj <- model.train(
    STEAM.obj, mode = "nested",
    model = model, kernel = kernel,
    n_outer_folds = n_outer_folds, n_inner_folds = n_inner_folds,
    metric = metric, tune.grid = tune.grid,
    allowParallel = allowParallel,   # inner caret parallel (will be disabled if cv.cores>1)
    seed = seed, n.size = n.size,
    outer_folds = outer_folds, inner_folds = inner_folds,
    verbose = verbose,
    cv.cores = cv.cores,
    parallel_mode = parallel_mode,
    maxnweights = maxnweights
  )
  message("Finished NestedCV training (outer+inner with leakage-safe averaging)")
  
  STEAM.obj <- model.predict(STEAM.obj, mode = "nested")
  message("Collected NestedCV results")
  STEAM.obj
}

