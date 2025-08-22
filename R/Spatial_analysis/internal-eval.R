# ======================================================================
# internal-eval.R
#
# Internal evaluation helpers for STEAM spatial anchor analysis
# ----------------------------------------------------------------------
# Functions here compute fold-level and cell-level summaries:
#   • EvalAnchor – per-fold accuracy before/after flips
#   • FlipStats  – aggregate flip statistics (counts, rates, Δaccuracy)
#   • FlipLong   – long-format data for plotting flips
#
# summary/plotting functions (e.g., STEAM_anchor_summary, plot_flips_*).
# ======================================================================


EvalAnchor <- function(ncv, final_by_fold) {
  stopifnot(length(ncv$outer_result) == length(final_by_fold))
  rows <- lapply(seq_along(ncv$outer_result), function(i) {
    of <- ncv$outer_result[[i]]
    pd <- as.data.frame(of$preds)
    rownames(pd) <- rownames(pd) %||% as.character(of$test_id %||% of$test_ix)
    ids <- rownames(pd)
    truth  <- setNames(as.character(pd$testy), ids)
    before <- setNames(as.character(pd$predy), ids)
    after  <- final_by_fold[[i]][ids]
    
    acc_b <- mean(before == truth)
    acc_a <- mean(after  == truth)
    data.frame(Outer_Fold = i, Acc_Before = acc_b, Acc_After = acc_a,
               Delta = acc_a - acc_b, Flips = sum(before != after), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}


FlipStats <- function(ncv, final_by_fold) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  stopifnot(length(ncv$outer_result) == length(final_by_fold))
  
  rows <- lapply(seq_along(ncv$outer_result), function(i) {
    of <- ncv$outer_result[[i]]
    pd <- as.data.frame(of$preds)
    # robust rownames: prefer rownames(preds); else test_id; else test_ix
    rownames(pd) <- rownames(pd) %||% as.character(of$test_id %||% of$test_ix)
    
    ids    <- rownames(pd)
    truth  <- setNames(as.character(pd$testy), ids)
    before <- setNames(as.character(pd$predy), ids)
    
    after  <- final_by_fold[[i]]
    if (is.null(names(after))) names(after) <- ids
    
    common <- intersect(ids, names(after))
    truth  <- truth[common]
    before <- before[common]
    after  <- after[common]
    
    # classify each cell
    status <- ifelse(before != after,
                     ifelse(after == truth, "Flipped to Correct", "Flipped to Wrong"),
                     "Unchanged")
    
    n_all       <- length(common)
    n_flip      <- sum(status != "Unchanged")
    n_correct   <- sum(status == "Flipped to Correct")
    n_wrong     <- sum(status == "Flipped to Wrong")
    n_unchanged <- sum(status == "Unchanged")
    
    acc_b <- mean(before == truth)
    acc_a <- mean(after  == truth)
    
    data.frame(
      Outer_Fold         = i,
      N_Test             = n_all,
      Flips_Total        = n_flip,
      Flips_To_Correct   = n_correct,
      Flips_To_Wrong     = n_wrong,
      Unchanged          = n_unchanged,
      Flip_Rate          = n_flip / n_all,
      Correct_Flip_Rate  = if (n_flip > 0) n_correct / n_flip else 0,
      Wrong_Flip_Rate    = if (n_flip > 0) n_wrong   / n_flip else 0,
      Acc_Before         = acc_b,
      Acc_After          = acc_a,
      Delta_Accuracy     = acc_a - acc_b,
      stringsAsFactors   = FALSE
    )
  })
  
  do.call(rbind, rows)
}


FlipLong <- function(ncv, final_by_fold) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  stopifnot(length(ncv$outer_result) == length(final_by_fold))
  
  rows <- lapply(seq_along(ncv$outer_result), function(i) {
    of <- ncv$outer_result[[i]]
    pd <- as.data.frame(of$preds)
    rownames(pd) <- rownames(pd) %||% as.character(of$test_id %||% of$test_ix)
    
    ids    <- rownames(pd)
    truth  <- setNames(as.character(pd$testy), ids)
    before <- setNames(as.character(pd$predy), ids)
    
    after  <- final_by_fold[[i]]
    if (is.null(names(after))) names(after) <- ids
    
    common <- intersect(ids, names(after))
    truth  <- truth[common]
    before <- before[common]
    after  <- after[common]
    
    status <- ifelse(before != after,
                     ifelse(after == truth, "Flipped to Correct", "Flipped to Wrong"),
                     "Unchanged")
    
    data.frame(
      Outer_Fold = i,
      cell_id    = common,
      Truth      = truth,
      Before     = before,
      After      = after,
      Flip_Status= status,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, rows)
}
