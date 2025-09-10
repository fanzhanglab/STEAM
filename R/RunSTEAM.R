#' Run STEAM: Training and Evaluation
#'
#' Main entry point to train and evaluate STEAM models. Supports both
#' simple train/test splitting and leakage-safe nested cross-validation.
#'
#' @param STEAM.obj A STEAM object created by \code{\link{LoadSTEAM}}.
#' @param mode Either "simple" (train/test split) or "nested" (nested CV).
#' @param train.ratio Proportion of data for training (simple mode only).
#' @param n.size Neighborhood size for neighborhood averaging (spatial smoothing).
#' @param seed Random seed (default 123).
#' @param cv.folds Number of folds for CV (simple mode).
#' @param cv.repeats Number of repeats for repeated CV (simple mode).
#' @param n_outer_folds Outer folds for nested CV.
#' @param n_inner_folds Inner folds for nested CV.
#' @param model Model type: "rf", "svm", "xgb", or "multinom".
#' @param kernel Kernel type for SVM ("linear" or "radial").
#' @param metric Metric to optimize during training (default "Kappa").
#' @param tune.grid Optional custom tuning grid.
#' @param allowParallel Allow caret-level parallelism (default FALSE).
#' @param outer_folds Optional precomputed outer folds (nested mode).
#' @param inner_folds Optional precomputed inner folds (nested mode).
#' @param verbose Logical, print messages (default TRUE).
#' @param cv.cores Number of cores for nestedcv parallelism.
#' @param parallel_mode Optional parallel backend ("multisession", etc.).
#' @param maxnweights Maximum weights (for multinom models).
#' @param trainval.ratio Train set ratio for train-val split
#' @param n.tree tree depth
#' @param train.folder.name folder name to save
#'
#' @return A STEAM object with trained models, predictions, and evaluation metrics.
#' @export
#' @import caret
#' @import randomForest
#' @import e1071
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


  return(STEAM.obj)
}
