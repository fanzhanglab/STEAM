neighavg_fit <- function(x, y = NULL, coords_all, n.size = 5L) {
  if (is.null(rownames(x))) stop("Training matrix x must have rownames (cell ids).")
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
  
  structure(
    list(
      ids_train   = ids_train,
      X_train     = as.matrix(x),
      coords_all  = coords_all,
      coords_train= coords_train,
      x_name      = unname(cn["x_name"]),
      y_name      = unname(cn["y_name"]),
      n.size      = as.integer(n.size)
    ),
    class = "neighavg_fit"   # IMPORTANT: class must match predict.neighavg_fit
  )
}
