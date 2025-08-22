model.train <- function(
    STEAM.obj,
    mode = c("simple","nested"),
    model = c("rf","svm","xgb","multinom"),
    kernel = "linear",
    cv.folds = 10, cv.repeats = 1,
    n_outer_folds = 5, n_inner_folds = 3,
    metric = "Kappa",
    tune.grid = NULL,
    allowParallel = FALSE,
    seed = 42,
    n.size = 5,
    outer_folds = NULL, inner_folds = NULL,
    verbose = TRUE,
    cv.cores = 1,
    parallel_mode = NULL,
    maxnweights = 10000
) {
  set.seed(seed)
  mode  <- match.arg(mode)
  model <- match.arg(model)
  
  .global_svm_levels <- levels(factor(STEAM.obj$train$train.data.labels))
  
  # ---- mapper ----
  map_model <- function(model, kernel, p, tg) {
    if (model == "rf") {
      grid <- if (is.null(tg)) {
        data.frame(mtry = unique(pmax(1L, c(floor(sqrt(p)), floor(p/3), pmax(2L, floor(p/10))))))
      } else tg
      list(method = "rf", grid = grid)
      
    } else if (model == "svm") {
      m <- if (kernel == "radial") "svmRadial" else "svmLinear"
      grid <- if (is.null(tg)) {
        if (m == "svmLinear") {
          expand.grid(C = 2^(-3:5))
        } else {
          expand.grid(C = 2^(-3:5), sigma = 2^(-8:-2))
        }
      } else tg
      list(method = m, grid = grid)
      
    } else if (model == "xgb") {
      grid <- if (is.null(tg)) {
        expand.grid(
          nrounds = c(100, 200),
          max_depth = c(2, 3, 4),
          eta = c(0.05, 0.1, 0.2),
          gamma = c(0, 1),
          colsample_bytree = c(0.6, 0.8, 1.0),
          min_child_weight = c(1, 3, 5),
          subsample = c(0.7, 1.0)
        )
      } else tg
      list(method = "xgbTree", grid = grid)
      
    } else { # multinom
      grid <- if (is.null(tg)) expand.grid(decay = c(0, 0.01, 0.1)) else tg
      list(method = "multinom", grid = grid)
    }
  }
  
  # ================= SIMPLE mode =================
  if (mode == "simple") {
    stopifnot(!is.null(STEAM.obj$train$avg.matrix))
    X <- t(as.matrix(STEAM.obj$train$avg.matrix))
    y <- factor(STEAM.obj$train$train.data.labels)
    
    mm <- map_model(model, kernel, ncol(X), tune.grid)
    
    uses_probs <- TRUE  # Always use probabilities for consistency
    
    if (metric %in% c("ROC","PR","logLoss","mnLogLoss")) {
      sumfun <- if (nlevels(y) > 2) caret::multiClassSummary else caret::twoClassSummary
    } else {
      sumfun <- caret::multiClassSummary  # Use multiClassSummary for Kappa, Accuracy, etc.
    }
    
    preproc <- if (mm$method == "svmLinear") c("center","scale") else NULL
    
    ctrl <- caret::trainControl(
      method = if (cv.repeats > 1) "repeatedcv" else "cv",
      number = cv.folds, repeats = cv.repeats,
      classProbs = uses_probs,      # Always TRUE now
      summaryFunction = sumfun,
      allowParallel = allowParallel,
      savePredictions = "final"
    )
    
    fit <- caret::train(
      x = X, y = y,
      method = mm$method,
      tuneGrid = mm$grid,
      trControl = ctrl,
      metric = metric,
      preProcess = preproc,
      MaxNWts = if (mm$method == "multinom") maxnweights else NULL
    )
    
    STEAM.obj$train$model <- fit
    return(STEAM.obj)
  }
  
  # ================= NESTED mode =================
  clean_labels <- function(labels) {
    labels <- as.character(labels)
    labels <- gsub("/", ".", labels); labels <- gsub("-", ".", labels)
    labels <- gsub(" ", "_", labels); labels <- gsub("[()]", "", labels)
    labels <- gsub(",", "", labels);  labels <- gsub("\\.\\.", ".", labels)
    labels <- trimws(labels); labels <- iconv(labels, to = "ASCII//TRANSLIT"); labels
  }
  y <- factor(clean_labels(STEAM.obj$labels))
  X <- t(as.matrix(STEAM.obj$count_exp))
  coords_all <- as.data.frame(STEAM.obj$spatial)
  if (is.null(rownames(X))) stop("count_exp must have column names matching spatial rownames.")
  if (is.null(rownames(coords_all))) rownames(coords_all) <- colnames(STEAM.obj$count_exp)
  
  if (is.null(outer_folds)) {
    outer_folds <- caret::createFolds(y, k = n_outer_folds)  # stratified
  }
  if (is.null(inner_folds)) {
    inner_folds <- lapply(outer_folds, function(te) {
      ytr <- y[-te]
      caret::createFolds(ytr, k = n_inner_folds)
    })
  }
  
  mm <- map_model(model, kernel, ncol(X), tune.grid)
  
  uses_probs <- TRUE  # Always use probabilities for consistency
  
  if (metric %in% c("ROC","PR","logLoss","mnLogLoss")) {
    sumfun <- if (nlevels(y) > 2) caret::multiClassSummary else caret::twoClassSummary
  } else {
    sumfun <- caret::multiClassSummary  # Use multiClassSummary for Kappa, Accuracy, etc.
  }
  
  preproc <- if (mm$method == "svmLinear") c("center","scale") else NULL
  
  tr_ctrl_inner <- caret::trainControl(
    method = "cv", number = n_inner_folds,
    classProbs = uses_probs,
    summaryFunction = sumfun,
    savePredictions = "final",
    allowParallel = allowParallel
  )
  
  if (cv.cores > 1 && isTRUE(allowParallel)) {
    message("Note: cv.cores > 1 and allowParallel = TRUE → disabling inner caret parallel to avoid nested parallelism.")
    tr_ctrl_inner$allowParallel <- FALSE
  }
  
  ncv_args <- list(
    y = y, x = X,
    method = mm$method,
    tuneGrid = mm$grid,
    metric = metric,
    trControl = tr_ctrl_inner,
    outer_folds = outer_folds,
    inner_folds = inner_folds,
    modifyX = "neighavg_fit",
    modifyX_useY = TRUE,
    modifyX_options = list(coords_all = coords_all, n.size = as.integer(n.size)),
    filterFUN = NULL,
    filter_options = list(),
    finalCV = NA,
    cv.cores = cv.cores,
    parallel_mode = parallel_mode,
    verbose = as.logical(verbose)
  )
  
  if (!is.null(preproc)) {
    ncv_args$preProcess <- preproc
  }
  
  if (mm$method == "multinom") {
    ncv_args$MaxNWts <- maxnweights
  }
  
  ncv <- do.call(nestedcv::nestcv.train, ncv_args)
  
  STEAM.obj$train$ncv <- ncv
  STEAM.obj
}
