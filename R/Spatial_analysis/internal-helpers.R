# ======================================================================
# internal-helpers.R
#
# Internal helper functions for STEAM
# ----------------------------------------------------------------------
# Building blocks for spatial anchor analysis and nested CV evaluation.
# 
# Responsibilities include:
#   • Fold handling (inferring train/test IDs for outer CV folds)
#   • Spatial utilities (kNN search, Gaussian weighting, consensus)
#   • Confidence scoring (agreement, purity, model-based confidence)
#   • Anchor-based corrections (rules for flipping predictions)
# ======================================================================


inferfolds <- function(ncv, coords_all, strict = TRUE) {
  if (is.null(ncv$outer_result) || !length(ncv$outer_result))
    stop("ncv$outer_result is missing or empty.")
  all_ids <- rownames(coords_all)
  if (is.null(all_ids))
    stop("coords_all must have rownames (cell IDs).")
  if (anyNA(all_ids)) stop("coords_all rownames contain NA.")
  if (any(duplicated(all_ids))) {
    d <- unique(all_ids[duplicated(all_ids)])
    stop("coords_all rownames contain duplicates (e.g., ", paste(head(d, 5), collapse=", "), ").")
  }
  
  nfolds <- length(ncv$outer_result)
  out <- vector("list", nfolds)
  
  pick_from_preds <- function(of) {
    prn <- rownames(of$preds)
    if (is.null(prn)) return(NULL)
    # keep preds' order
    keep <- prn %in% all_ids
    if (!all(keep)) {
      msg <- sprintf("Only %d/%d preds rownames match coords_all.", sum(keep), length(prn))
      if (strict) stop(msg) else warning(msg, call. = FALSE)
      if (!any(keep)) return(NULL)
    }
    test_id <- prn[keep]
    list(test_id = test_id, test_ix = match(test_id, all_ids))
  }
  
  pick_from_ix <- function(ix) {
    if (is.null(ix)) return(NULL)
    ix <- as.integer(ix)
    if (anyNA(ix) || any(ix < 1) || any(ix > length(all_ids)))
      stop("test_ix out of bounds or contains NA.")
    list(test_ix = ix, test_id = all_ids[ix])
  }
  
  pick_from_id <- function(id) {
    if (is.null(id)) return(NULL)
    id <- as.character(id)
    if (anyNA(id)) stop("test_id contains NA.")
    ix <- match(id, all_ids)
    if (anyNA(ix)) {
      bad <- id[is.na(ix)]
      msg <- paste0("Some test_id not in coords_all (e.g., ", paste(head(bad, 5), collapse=", "), ").")
      if (strict) stop(msg) else warning(msg, call. = FALSE)
      id  <- id[!is.na(ix)]; ix <- ix[!is.na(ix)]
      if (!length(ix)) return(NULL)
    }
    list(test_ix = ix, test_id = id)
  }
  
  for (i in seq_len(nfolds)) {
    of <- ncv$outer_result[[i]]
    
    got <- NULL
    if (!is.null(of$preds)) got <- pick_from_preds(of)
    if (is.null(got) && !is.null(of$test_ix)) got <- pick_from_ix(of$test_ix)
    if (is.null(got) && !is.null(of$test_id)) got <- pick_from_id(of$test_id)
    if (is.null(got) && !is.null(ncv$outer_folds) && length(ncv$outer_folds) >= i)
      got <- pick_from_ix(ncv$outer_folds[[i]])
    
    if (is.null(got))
      stop(sprintf("Outer fold %d: cannot infer train/test IDs.", i))
    
    # finalize train set = complement
    test_ix  <- got$test_ix
    test_id  <- got$test_id
    train_ix <- setdiff(seq_along(all_ids), test_ix)
    out[[i]] <- list(
      test_ix  = test_ix,
      train_ix = train_ix,
      test_id  = test_id,
      train_id = all_ids[train_ix]
    )
  }
  
  out
}



tuneanchors <- function(outer_result, coords_all, 
                        k_grid = c(4,6,8,12),
                        anchor_cut_grid = seq(0.70, 0.90, 0.05),
                        consensus_grid  = seq(0.55, 0.80, 0.05),
                        class_precision = NULL,
                        min_class_precision = 0.70,
                        train_ids) {
  fit <- outer_result$fit
  if (is.null(fit)) stop("Outer fit missing.")
  if (is.null(train_ids)) stop("train_ids required (rowIndex -> cell_id).")
  
  inner <- .get_inner_cv_preds(fit, train_ids)
  if (!"pred" %in% names(inner) && "predy" %in% names(inner)) names(inner)[names(inner)=="predy"] <- "pred"
  if (!"obs"  %in% names(inner) && "testy" %in% names(inner)) names(inner)[names(inner)=="testy"] <- "obs"
  
  best <- list(score = -Inf, k = NA_integer_, anchor_cut = NA_real_, consensus = NA_real_)
  for (k in k_grid) for (ac in anchor_cut_grid) for (cc in consensus_grid) {
    s <- tryCatch(
      .score_setting_inner(inner, coords_all, k, ac, cc,
                           class_precision = class_precision,
                           min_class_precision = min_class_precision),
      error = function(e) { message(sprintf("[tuneanchors] k=%d ac=%.2f cc=%.2f ERROR: %s", k, ac, cc, e$message)); NA_real_ }
    )
    if (!is.finite(s)) {
      message(sprintf("[tuneanchors] k=%d ac=%.2f cc=%.2f -> non-finite score", k, ac, cc))
      next
    }
    if (s >= best$score) best <- list(score = s, k = k, anchor_cut = ac, consensus = cc)
  }
  if (!is.finite(best$score)) message("[tuneanchors] No finite score found across grid.")
  best
}




AnchorKNN <- function(coordinates, k = 8) {
  stopifnot(requireNamespace("FNN", quietly = TRUE))
  coords <- as.matrix(coordinates)
  if (nrow(coords) < 2) stop("Need at least 2 points for kNN.")
  k_eff <- min(k, nrow(coords) - 1L)
  kn <- FNN::get.knn(coords, k = k_eff)
  list(
    nn.index = kn$nn.index,
    nn.dist  = kn$nn.dist,
    idx_list = split(kn$nn.index, row(kn$nn.index)),
    dist_list= split(kn$nn.dist,  row(kn$nn.dist))
  )
}


# helpers (internal) ----
conf_features_from_probs <- function(probs) {
  if (is.null(probs)) return(list(conf=NA_real_, margin=NA_real_, entropy=NA_real_))
  stopifnot(is.matrix(probs))
  pmax <- apply(probs, 1, max)
  p2   <- apply(probs, 1, function(v) sort(v, TRUE)[2])
  ent  <- -rowSums(probs * log(pmax(probs, .Machine$double.eps)))
  list(conf=pmax, margin=pmax - p2, entropy=ent)
}

gaussian_weights <- function(d, sigma) {
  if (!is.finite(sigma) || sigma <= 0) sigma <- 1
  exp(-(d^2)/(2*sigma^2))
}

agreement_purity <- function(neigh_preds, self_pred, w) {
  if (!length(neigh_preds)) return(list(agree=0.5, purity=1))
  sw <- sum(w); if (sw <= 0 || !is.finite(sw)) w <- rep(1/length(w), length(w)) else w <- w/sum(w)
  agree <- sum(w * (neigh_preds == self_pred))
  vote  <- tapply(w, neigh_preds, sum)
  purity <- max(vote, na.rm = TRUE)
  list(agree=ifelse(is.finite(agree), agree, 0), purity=ifelse(is.finite(purity), purity, 0))
}

combine_confidence <- function(agree, purity, conf, w_ag=0.55, w_pu=0.25, w_cf=0.20) {
  if (is.na(conf)) { s <- w_ag + w_pu; w_ag <- w_ag/s; w_pu <- w_pu/s; w_cf <- 0 }
  w_ag*agree + w_pu*purity + w_cf*conf
}

# main helper (internal) ----
confscore <- function(predictions, knn, probs = NULL, sigma = NULL,
                      w_ag = 0.55, w_pu = 0.25, w_cf = 0.20) {
  n <- length(predictions)
  stopifnot(length(knn$idx_list) == n, length(knn$dist_list) == n)
  
  # model confidence features
  cf <- if (!is.null(probs)) conf_features_from_probs(probs) else list(conf=rep(NA_real_, n),
                                                                       margin=rep(NA_real_, n),
                                                                       entropy=rep(NA_real_, n))
  # sigma
  if (is.null(sigma)) {
    sigma <- stats::median(as.numeric(knn$nn.dist), na.rm = TRUE)
    if (!is.finite(sigma) || sigma <= 0) sigma <- 1
  }
  
  out_score <- numeric(n)
  out_agree <- numeric(n)
  out_purity<- numeric(n)
  
  for (i in seq_len(n)) {
    idx <- knn$idx_list[[i]]; dst <- knn$dist_list[[i]]
    if (!length(idx)) { out_score[i] <- combine_confidence(0.5,1, cf$conf[i], w_ag,w_pu,w_cf); next }
    w <- gaussian_weights(dst, sigma)
    ap <- agreement_purity(predictions[idx], predictions[i], w)
    out_agree[i] <- ap$agree; out_purity[i] <- ap$purity
    # rescale conf to [0,1] per-vector (avoid NA issues)
    conf_i <- cf$conf[i]
    if (!all(is.na(cf$conf))) {
      rng <- range(cf$conf, na.rm = TRUE); conf_i <- if (diff(rng)>0) (conf_i - rng[1])/diff(rng) else 0
    }
    out_score[i] <- combine_confidence(ap$agree, ap$purity, conf_i, w_ag, w_pu, w_cf)
  }
  
  data.frame(
    cell_id            = names(predictions),
    spatial_confidence = out_score,
    neighbor_agreement = out_agree,
    neighborhood_purity= out_purity,
    model_confidence   = cf$conf,
    model_margin       = cf$margin,
    model_entropy      = cf$entropy,
    stringsAsFactors = FALSE
  )
}


# --- helpers (internal) -------------------------------------------------------

# 1) Which neighbors are reliable anchors
.anchor_mask <- function(spatial_anchors, cut) {
  stopifnot(is.data.frame(spatial_anchors),
            "spatial_confidence" %in% names(spatial_anchors))
  spatial_anchors$spatial_confidence >= cut
}

# 2) Weighted class consensus among a set of neighbor indices
.weighted_consensus <- function(neigh_idx, neigh_dist, cls_vec, sigma) {
  if (!length(neigh_idx)) return(list(winner = NA_character_, consensus = 0, n = 0))
  w <- exp(-(neigh_dist^2) / (2 * sigma^2))
  if (!is.finite(sum(w)) || sum(w) <= 0) return(list(winner = NA_character_, consensus = 0, n = 0))
  vote <- tapply(w, cls_vec[neigh_idx], sum)
  winner <- names(which.max(vote))
  consensus <- as.numeric(max(vote)) / sum(w)
  list(winner = winner, consensus = consensus, n = length(neigh_idx))
}

# 3) Guard rails; return TRUE if a flip is allowed
.passes_guards <- function(candidate_class,
                           consensus,
                           orig_class,
                           threshold,
                           class_precision = NULL,
                           min_class_precision = 0.7,
                           margin_value = NA_real_,
                           max_allowed_margin = 0.15) {
  if (is.na(candidate_class) || candidate_class == orig_class) return(FALSE)
  if (consensus < threshold) return(FALSE)
  if (!is.null(class_precision)) {
    cp <- unname(class_precision[candidate_class])
    if (is.na(cp) || cp < min_class_precision) return(FALSE)
  }
  if (!is.na(margin_value) && margin_value > max_allowed_margin) return(FALSE)
  TRUE
}

# 4) Row builder (avoids repeated rbind)
.make_flip_row <- function(id, orig, sug, cons, n_anchor) {
  data.frame(
    cell_id = id,
    original_prediction = orig,
    suggested_correction = sug,
    correction_confidence = min(0.95, 0.10 + cons),
    n_anchor_neighbors = n_anchor,
    anchor_consensus = cons,
    correction_rationale = sprintf("%s consensus %.2f via %d anchors", sug, cons, n_anchor),
    stringsAsFactors = FALSE
  )
}

# --- thin wrapper (internal) --------------------------------------------------

FlipAnchors <- function(predictions, knn, spatial_anchors,
                        correction_threshold = 0.6, anchor_cut = 0.8,
                        class_precision = NULL, min_class_precision = 0.7,
                        model_margin = NULL, max_allowed_margin = 0.15,
                        sigma = NULL) {
  
  # Basic validation
  ids <- names(predictions)
  if (is.null(ids)) stop("predictions must be a named vector (cell IDs as names).")
  stopifnot(is.list(knn), length(knn$idx_list) == length(predictions),
            length(knn$dist_list) == length(predictions))
  if (is.null(sigma)) {
    sigma <- stats::median(as.numeric(knn$nn.dist), na.rm = TRUE)
    if (!is.finite(sigma) || sigma <= 0) sigma <- 1.0
  }
  
  # Precompute masks & lookups
  reliable <- .anchor_mask(spatial_anchors, anchor_cut)
  margin_lookup <- if (!is.null(model_margin)) setNames(as.numeric(model_margin), names(model_margin)) else NULL
  pred_vec <- as.character(predictions)  # ensure character
  
  # Collect rows in a list (faster than rbind in a loop)
  out_rows <- vector("list", length(pred_vec))
  keep <- 0L
  
  for (i in seq_along(pred_vec)) {
    idx <- knn$idx_list[[i]]; dst <- knn$dist_list[[i]]
    if (!length(idx)) next
    
    a_idx <- idx[reliable[idx]]; a_dst <- dst[reliable[idx]]
    if (!length(a_idx)) next
    
    wc <- .weighted_consensus(neigh_idx = a_idx,
                              neigh_dist = a_dst,
                              cls_vec = pred_vec,
                              sigma = sigma)
    
    # Guard rails
    mg <- if (!is.null(margin_lookup)) margin_lookup[ids[i]] else NA_real_
    if (!.passes_guards(candidate_class = wc$winner,
                        consensus = wc$consensus,
                        orig_class = pred_vec[i],
                        threshold = correction_threshold,
                        class_precision = class_precision,
                        min_class_precision = min_class_precision,
                        margin_value = mg,
                        max_allowed_margin = max_allowed_margin)) next
    
    keep <- keep + 1L
    out_rows[[keep]] <- .make_flip_row(ids[i], pred_vec[i], wc$winner, wc$consensus, wc$n)
  }
  
  if (keep == 0L) {
    return(data.frame(
      cell_id = character(0),
      original_prediction = character(0),
      suggested_correction = character(0),
      correction_confidence = numeric(0),
      n_anchor_neighbors = integer(0),
      anchor_consensus = numeric(0),
      correction_rationale = character(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, out_rows[seq_len(keep)])
}