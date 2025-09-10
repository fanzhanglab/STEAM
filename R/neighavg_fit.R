# --- 0) helper to guess coordinate column names ---
.guess_coord_names <- function(df) {
  cn <- colnames(df)
  x_candidates <- c("x", "X", "Xcoord", "Xcorr", "xcoord", "xcorr", "pos_x", "px", "X_centroid", "x_centroid")
  y_candidates <- c("y", "Y", "Ycoord", "Ycorr", "ycoord", "ycorr", "pos_y", "py", "Y_centroid", "y_centroid")
  
  x_name <- intersect(x_candidates, cn)
  y_name <- intersect(y_candidates, cn)
  if (length(x_name) == 0L || length(y_name) == 0L)
    stop("Could not infer x/y column names in coords_all. Please ensure columns like 'x'/'y' or 'Xcorr'/'Ycorr' are present.")
  setNames(c(x_name[1], y_name[1]), c("x_name","y_name"))
}

# --- 1) fitter (stores what predict() needs) ---
neighavg_fit <- function(x, y = NULL, coords_all, n.size = 5L) {
  if (is.null(rownames(x))) stop("Training matrix x must have rownames (cell ids).")
  if (is.null(colnames(x))) stop("Training matrix x must have colnames (features).")
  if (is.null(rownames(coords_all))) stop("coords_all must have rownames (cell ids).")
  cn <- .guess_coord_names(coords_all)
  
  ids_train <- rownames(x)
  miss <- setdiff(ids_train, rownames(coords_all))
  if (length(miss)) {
    stop(sprintf("Missing coordinates for %d TRAIN samples. First few: %s",
                 length(miss), paste(utils::head(miss, 5), collapse = ", ")))
  }
  
  coords_train <- coords_all[ids_train, c(cn["x_name"], cn["y_name"]), drop = FALSE]
  coords_train[[cn["x_name"]]] <- as.numeric(coords_train[[cn["x_name"]]])
  coords_train[[cn["y_name"]]] <- as.numeric(coords_train[[cn["y_name"]]])
  
  X_train <- as.matrix(x)
  storage.mode(X_train) <- "double"
  
  structure(
    list(
      ids_train    = ids_train,
      X_train      = X_train,                  # samples x features
      coords_all   = coords_all,               # data.frame, rownames = ALL cells
      coords_train = coords_train,             # subset aligned to ids_train
      x_name       = unname(cn["x_name"]),
      y_name       = unname(cn["y_name"]),
      n.size       = as.integer(n.size)
    ),
    class = "neighavg_fit"
  )
}

# --- 2) predictor (the hardened version I sent earlier) ---
predict.neighavg_fit <- function(object, newdata, ...) {
  stopifnot(!is.null(newdata))
  if (is.null(rownames(newdata))) stop("NEW data must have rownames (cell IDs).")
  ids_new <- rownames(newdata)
  if (anyNA(ids_new)) stop("NEW rownames contain NA.")
  if (any(duplicated(ids_new))) {
    d <- unique(ids_new[duplicated(ids_new)])
    stop(sprintf("NEW rownames contain duplicates. First few: %s",
                 paste(head(d, 5), collapse = ", ")))
  }
  
  req <- c("ids_train","X_train","coords_all","n.size","x_name","y_name")
  miss <- setdiff(req, names(object))
  if (length(miss)) stop("neighavg_fit object missing: ", paste(miss, collapse=", "))
  
  ids_train  <- object$ids_train
  X_train    <- object$X_train
  coords_all <- object$coords_all
  n.size     <- as.numeric(object$n.size)
  x_name     <- object$x_name
  y_name     <- object$y_name
  
  ctrain <- object$coords_train
  if (is.null(ctrain)) {
    if (is.null(rownames(coords_all))) stop("coords_all must have rownames (cell IDs).")
    missing_tr <- setdiff(ids_train, rownames(coords_all))
    if (length(missing_tr))
      stop(sprintf("Missing coordinates for %d TRAIN samples (e.g., %s).",
                   length(missing_tr), paste(head(missing_tr, 5), collapse=", ")))
    ctrain <- coords_all[ids_train, , drop = FALSE]
  }
  
  need_cols <- c(x_name, y_name)
  if (!all(need_cols %in% colnames(ctrain)))
    stop("coords_train lacks required columns: ", paste(setdiff(need_cols, colnames(ctrain)), collapse=", "))
  if (!all(need_cols %in% colnames(coords_all)))
    stop("coords_all lacks required columns: ", paste(setdiff(need_cols, colnames(coords_all)), collapse=", "))
  
  cn_new <- colnames(newdata)
  cn_tr  <- colnames(X_train)
  if (is.null(cn_new) || is.null(cn_tr)) stop("Both newdata and X_train must have column names (features).")
  common <- intersect(cn_new, cn_tr)
  if (length(common) == 0L) stop("No common features between newdata and X_train.")
  common <- sort(common)
  newdata <- as.matrix(newdata[, common, drop = FALSE])
  X_train <- as.matrix(X_train[, common, drop = FALSE])
  storage.mode(newdata) <- "double"
  storage.mode(X_train) <- "double"
  
  if (is.null(rownames(coords_all)))
    stop("coords_all must have rownames (cell IDs).")
  miss_new <- setdiff(ids_new, rownames(coords_all))
  if (length(miss_new))
    stop(sprintf("Missing coordinates for %d NEW samples (e.g., %s).",
                 length(miss_new), paste(head(miss_new, 5), collapse=", ")))
  
  cnew <- coords_all[ids_new, need_cols, drop = FALSE]
  ctrain[[x_name]] <- as.numeric(ctrain[[x_name]])
  ctrain[[y_name]] <- as.numeric(ctrain[[y_name]])
  cnew[[x_name]]   <- as.numeric(cnew[[x_name]])
  cnew[[y_name]]   <- as.numeric(cnew[[y_name]])
  
  if (is.null(rownames(ctrain)) || anyNA(rownames(ctrain)))
    stop("coords_train must have non-NA rownames matching ids_train.")
  ctrain <- ctrain[ids_train, , drop = FALSE]
  
  res <- matrix(NA_real_, nrow = nrow(newdata), ncol = ncol(newdata),
                dimnames = list(ids_new, colnames(newdata)))
  
  xtr <- ctrain[[x_name]]; ytr <- ctrain[[y_name]]
  for (i in seq_len(nrow(cnew))) {
    dx <- abs(xtr - cnew[[x_name]][i])
    dy <- abs(ytr - cnew[[y_name]][i])
    neigh_idx <- which(dx <= n.size & dy <= n.size)
    if (length(neigh_idx) <= 1L) res[i, ] <- newdata[i, ] else
      res[i, ] <- colMeans(X_train[neigh_idx, , drop = FALSE], na.rm = TRUE)
  }
  
  if (!identical(rownames(res), ids_new))
    stop("Internal invariant violated: output rownames != input rownames")
  res
}