#' Spatial Anchor Analysis (tuning + corrections over outer folds)
#'
#' @description
#' Runs STEAM's spatial anchor analysis on nested-CV results:
#' per-fold hyperparameter tuning (k, anchor_cut, consensus), applies
#' leak-free spatial flips on outer test sets, and aggregates evaluation.
#'
#' @param STEAM.obj A STEAM object produced by \code{RunSTEAM(mode = "nested")}.
#'                  Must contain \code{$nested$ncv} and \code{$spatial}.
#' @param k_grid Integer vector of k-NN neighborhood sizes to try.
#' @param anchor_cut_grid Numeric vector in (0,1]; minimum spatial confidence
#'        for a neighbor to be considered an anchor.
#' @param consensus_grid Numeric vector in (0,1]; minimum weighted neighbor
#'        consensus required to flip a label.
#' @param class_precision Optional named numeric vector of per-class precision.
#' @param min_class_precision Numeric in [0,1]. Minimum allowed class precision.
#' @param probs_by_fold Optional list of class-probability matrices per fold.
#' @param parallel One of \code{"none"}, \code{"multicore"}, \code{"multisession"}.
#' @param workers Integer number of workers if \code{parallel!="none"}.
#' @param progress Logical; show progress bar when not using futures.
#' @param verbose Logical; emit per-fold messages.
#'
#' @return The input \code{STEAM.obj} augmented at
#'   \code{$spatial_anchor_analysis}.
#' @export
STEAM_anchor <- function(
    STEAM.obj,
    k_grid = c(4,6,8,12),
    anchor_cut_grid = seq(0.70, 0.90, 0.05),
    consensus_grid  = seq(0.55, 0.80, 0.05),
    class_precision = NULL,
    min_class_precision = 0.70,
    probs_by_fold = NULL,
    parallel = c("none","multicore","multisession"),
    workers  = max(1, parallel::detectCores() - 1),
    progress = TRUE,
    verbose  = TRUE
) {
  t0 <- Sys.time()
  parallel <- match.arg(parallel)
  
  .stop <- function(...) stop(sprintf(...), call. = FALSE)
  if (is.null(STEAM.obj$nested) || is.null(STEAM.obj$nested$ncv)) {
    .stop("STEAM.obj lacks nested CV results. Run RunSTEAM(mode='nested') first.")
  }
  if (is.null(STEAM.obj$spatial)) .stop("STEAM.obj has no spatial coordinates.")
  ncv <- STEAM.obj$nested$ncv
  coords_all <- STEAM.obj$spatial
  
  # infer train/test IDs per fold
  fold_ids <- inferfolds(ncv, coords_all)
  nfolds <- length(ncv$outer_result)
  
  tuned_list     <- vector("list", nfolds)
  final_by_fold  <- vector("list", nfolds)
  corr_list      <- vector("list", nfolds)
  eval_rows      <- vector("list", nfolds)
  
  .process_fold <- function(i) {
    of   <- ncv$outer_result[[i]]
    ids  <- fold_ids[[i]]
    pd   <- as.data.frame(of$preds)
    if (is.null(rownames(pd))) rownames(pd) <- ids$test_id
    test_ids <- rownames(pd)
    
    preds <- setNames(as.character(pd$predy), test_ids)
    truth <- setNames(as.character(pd$testy), test_ids)
    
    # tune with inner-CV
    best <- tuneanchors(
      outer_result        = of,
      coords_all          = coords_all,
      k_grid              = k_grid,
      anchor_cut_grid     = anchor_cut_grid,
      consensus_grid      = consensus_grid,
      class_precision     = class_precision,
      min_class_precision = min_class_precision,
      train_ids           = ids$train_id
    )
    
    if (!is.finite(best$score)) {
      return(list(
        tuned = best, corrections = NULL, final = preds,
        eval  = data.frame(Outer_Fold = i,
                           Acc_Before = mean(preds == truth),
                           Acc_After  = mean(preds == truth),
                           Delta      = 0, Flips = 0,
                           stringsAsFactors = FALSE)
      ))
    }
    
    # apply on test
    Ctest <- coords_all[test_ids, , drop = FALSE]
    knn   <- AnchorKNN(Ctest, k = best$k)
    probs <- if (!is.null(probs_by_fold)) probs_by_fold[[i]] else NULL
    sc    <- confscore(preds, knn, probs = probs)
    
    corr  <- FlipAnchors(
      predictions          = preds,
      knn                  = knn,
      spatial_anchors      = sc,
      correction_threshold = best$consensus,
      anchor_cut           = best$anchor_cut,
      class_precision      = class_precision,
      min_class_precision  = min_class_precision
    )
    
    final <- preds
    if (!is.null(corr) && nrow(corr)) final[corr$cell_id] <- corr$suggested_correction
    
    acc_b <- mean(preds == truth)
    acc_a <- mean(final == truth)
    if (!is.null(corr) && nrow(corr)) {
      corr$true_label       <- truth[corr$cell_id]
      corr$would_be_correct <- corr$suggested_correction == corr$true_label
      corr$Outer_Fold       <- i
    }
    
    list(
      tuned = best,
      corrections = corr,
      final = final,
      eval = data.frame(Outer_Fold = i,
                        Acc_Before = acc_b,
                        Acc_After  = acc_a,
                        Delta      = acc_a - acc_b,
                        Flips      = sum(preds != final),
                        stringsAsFactors = FALSE)
    )
  }
  
  # run across folds (serial or parallel)
  results <- lapply(seq_len(nfolds), .process_fold)
  
  for (i in seq_len(nfolds)) {
    r <- results[[i]]
    if (is.null(r)) next
    tuned_list[[i]]     <- list(tuned = r$tuned, corrections = r$corrections,
                                delta_accuracy = r$eval$Delta)
    final_by_fold[[i]]  <- r$final
    corr_list[[i]]      <- r$corrections
    eval_rows[[i]]      <- r$eval
  }
  
  corr_table <- do.call(rbind, Filter(Negate(is.null), corr_list))
  eval_tbl   <- do.call(rbind, Filter(Negate(is.null), eval_rows))
  flip_sum   <- if (length(Filter(Negate(is.null), final_by_fold)) == nfolds) {
    FlipStats(ncv, final_by_fold)
  } else NULL
  
  STEAM.obj$spatial_anchor_analysis <- list(
    per_fold      = tuned_list,
    final_by_fold = final_by_fold,
    eval          = eval_tbl,
    flip_summary  = flip_sum,
    corr_table    = corr_table,
    grids = list(
      k_grid = k_grid,
      anchor_cut_grid = anchor_cut_grid,
      consensus_grid  = consensus_grid
    ),
    runtime_info = list(
      started  = as.character(t0),
      finished = as.character(Sys.time()),
      parallel = parallel,
      workers  = if (parallel == "none") 1L else workers
    )
  )
  class(STEAM.obj) <- unique(c("STEAM_SpatialAnchors", class(STEAM.obj)))
  STEAM.obj
}