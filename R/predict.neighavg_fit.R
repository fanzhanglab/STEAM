predict.neighavg_fit <- function(object, newdata, ...) {
  stopifnot(!is.null(rownames(newdata)))
  ids_new <- rownames(newdata)
  
  if (anyNA(ids_new)) stop("NEW rownames contain NA.")
  if (any(duplicated(ids_new))) {
    d <- unique(ids_new[duplicated(ids_new)])
    stop(sprintf("NEW rownames contain duplicates. First few: %s",
                 paste(head(d, 5), collapse = ", ")))
  }
  
  ids_train  <- object$ids_train
  X_train    <- object$X_train
  coords_all <- object$coords_all
  ctrain     <- object$coords_train
  x_name     <- object$x_name
  y_name     <- object$y_name
  n.size     <- object$n.size
  
  # align columns (defensive): ensure X_train columns match newdata columns
  if (!identical(colnames(newdata), colnames(X_train))) {
    common <- intersect(colnames(newdata), colnames(X_train))
    if (length(common) == 0L)
      stop("No common features between newdata and X_train.")
    # reorder X_train to newdata's order, and drop missing cols from both sides consistently
    X_train  <- X_train[, common, drop = FALSE]
    newdata  <- as.matrix(newdata[, common, drop = FALSE])
  } else {
    newdata <- as.matrix(newdata)
  }
  
  # check coordinates for new ids
  missing_coords <- setdiff(ids_new, rownames(coords_all))
  if (length(missing_coords)) {
    stop(sprintf("Missing coordinates for %d NEW samples. First few: %s",
                 length(missing_coords), paste(head(missing_coords, 5), collapse = ", ")))
  }
  
  cnew <- coords_all[ids_new, c(x_name, y_name), drop = FALSE]
  cnew[[x_name]] <- as.numeric(cnew[[x_name]])
  cnew[[y_name]] <- as.numeric(cnew[[y_name]])
  
  # output matrix with same dims as (possibly reduced) newdata
  res <- matrix(NA_real_, nrow = nrow(newdata), ncol = ncol(newdata),
                dimnames = list(ids_new, colnames(newdata)))
  
  # simple rectangular neighborhood (Chebyshev-like via axis thresholds)
  for (i in seq_len(nrow(cnew))) {
    dx <- abs(ctrain[[x_name]] - cnew[[x_name]][i])
    dy <- abs(ctrain[[y_name]] - cnew[[y_name]][i])
    neigh_idx <- which(dx <= n.size & dy <= n.size)
    
    if (length(neigh_idx) <= 1L) {
      # keep original row if too few neighbors
      res[i, ] <- as.numeric(newdata[i, ])
    } else {
      res[i, ] <- colMeans(X_train[neigh_idx, , drop = FALSE], na.rm = TRUE)
    }
  }
  
  # sanity check: rownames must be preserved
  if (!identical(rownames(res), ids_new))
    stop("Internal invariant violated: output rownames != input rownames")
  res
}
