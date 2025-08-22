#' Summarize and visualize STEAM spatial-anchor corrections
#'
#' Produces a per-fold before/after accuracy table, a flip summary panel
#' (flipped-to-correct / flipped-to-wrong vs ΔAccuracy), and an optional
#' spatial plot highlighting flips for a chosen fold. If no anchor analysis
#' is present (or \code{force = TRUE}), it will run \code{STEAM_anchor()}
#' first, passing any extra tuning arguments via \code{...}.
#'
#' @param STEAM.obj A STEAM object with nested CV results (from \code{RunSTEAM(..., mode="nested")}).
#' @param force Logical; if TRUE, re-runs \code{STEAM_anchor()} even if results exist.
#' @param plot Logical; if TRUE, prints the flip panel and spatial flip plot.
#' @param fold_plot Integer scalar; which outer fold to visualize spatially.
#' @param layer_colors Optional named character vector of hex colors for layers
#'   (names must match label names, e.g., \code{c(Layer_1="#...", ..., WM="#...")}).
#' @param fade Numeric in [0,1]; fade factor for base-layer points in spatial plot (higher = lighter).
#' @param verbose Logical; print progress/info messages.
#'
#' @return (Invisibly) the updated \code{STEAM.obj} with:
#' \itemize{
#'   \item \code{$spatial_anchor_analysis$eval}: per-fold before/after accuracy.
#'   \item \code{$spatial_anchor_analysis$flip_summary}: flip counts & Δaccuracy per fold.
#'   \item \code{$spatial_anchor_analysis$final_by_fold}: post-correction predictions (list).
#'   \item \code{$viz$label_colors}: persisted palette (if provided).
#'   \item \code{$spatial_anchor_analysis$summary_quick}: a tiny summary list.
#' }
#'
#' @examples
#' # After RunSTEAM(..., mode="nested"):
#' # STEAM.obj <- STEAM_anchor_summary(STEAM.obj, plot = TRUE, fold_plot = 1)
#'
#' @export
STEAM_anchor_summary <- function(
    STEAM.obj,
    force = FALSE,
    plot  = TRUE,
    fold_plot = 1L,
    layer_colors = NULL,
    fade = 0.30,
    verbose = TRUE,
    ...
) {
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  .stop <- function(...) stop(sprintf(...), call. = FALSE)
  .msg  <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # ---- Core validation -------------------------------------------------------
  if (!is.list(STEAM.obj)) .stop("`STEAM.obj` must be a list-like STEAM object.")
  if (is.null(STEAM.obj$nested) || is.null(STEAM.obj$nested$ncv))
    .stop("Nested CV not found. Run RunSTEAM(mode = 'nested') first.")
  if (is.null(STEAM.obj$spatial))
    .stop("`STEAM.obj$spatial` is missing.")
  if (is.null(rownames(STEAM.obj$spatial)))
    .stop("`STEAM.obj$spatial` must have rownames (cell IDs).")
  
  if (!is.logical(force) || length(force) != 1L)
    .stop("`force` must be TRUE/FALSE.")
  if (!is.logical(plot) || length(plot) != 1L)
    .stop("`plot` must be TRUE/FALSE.")
  
  if (!is.numeric(fold_plot) || length(fold_plot) != 1L || is.na(fold_plot) || fold_plot < 1)
    .stop("`fold_plot` must be a positive integer fold index.")
  
  if (!is.numeric(fade) || length(fade) != 1L || is.na(fade) || fade < 0 || fade > 1)
    .stop("`fade` must be a number in [0, 1].")
  
  # Optional palette validation
  if (!is.null(layer_colors)) {
    if (!is.character(layer_colors) || is.null(names(layer_colors)))
      .stop("`layer_colors` must be a *named* character vector of hex colors.")
    # quick hex check (loose)
    bad_hex <- which(!grepl("^#?[0-9A-Fa-f]{6}$", layer_colors))
    if (length(bad_hex)) .stop("`layer_colors` contains invalid hex at positions: %s",
                               paste(bad_hex, collapse=", "))
  }
  
  # ---- Run analysis if needed -----------------------------------------------
  if (is.null(STEAM.obj$spatial_anchor_analysis) || isTRUE(force)) {
    .msg("Running STEAM_anchor() ...")
    STEAM.obj <- STEAM_anchor(STEAM.obj, verbose = verbose, ...)
  }
  
  # Persist palette for future plots
  if (!is.null(layer_colors)) {
    STEAM.obj$viz <- STEAM.obj$viz %||% list()
    STEAM.obj$viz$label_colors <- layer_colors
  }
  
  saa  <- STEAM.obj$spatial_anchor_analysis %||% list()
  eval <- saa$eval %||% data.frame()
  
  # ---- Print summary table ---------------------------------------------------
  cat("=== STEAM Spatial Anchor Summary ===\n")
  nfolds <- length(STEAM.obj$nested$ncv$outer_result)
  cat(sprintf("- Outer folds: %d\n", nfolds))
  
  if (nrow(eval)) {
    mean_delta <- mean(eval$Delta, na.rm = TRUE)
    total_flips <- sum(eval$Flips, na.rm = TRUE)
    cat(sprintf("- Mean ΔAccuracy: %+0.2f%%\n", 100*mean_delta))
    cat(sprintf("- Total flips applied: %d\n\n", total_flips))
    
    eval_print <- transform(eval,
                            Acc_Before = round(Acc_Before, 3),
                            Acc_After  = round(Acc_After, 3),
                            Delta      = round(Delta, 3))
    print(eval_print, row.names = FALSE)
  } else {
    cat("- No evaluation rows found.\n")
  }
  
  # ---- Corrections overview --------------------------------------------------
  if (!is.null(saa$corr_table) && nrow(saa$corr_table)) {
    cat("\nTop correction transitions (original → suggested):\n")
    df <- saa$corr_table
    tab <- as.data.frame(table(df$original_prediction, df$suggested_correction))
    names(tab) <- c("Original","Suggested","n")
    tab <- tab[order(tab$n, decreasing = TRUE), ]
    print(utils::head(tab, 10), row.names = FALSE)
  } else {
    cat("\nNo correction rows.\n")
  }
  
  # ---- Save quick summary back on object ------------------------------------
  STEAM.obj$spatial_anchor_analysis$summary_quick <- list(
    mean_delta_accuracy = if (nrow(eval)) mean(eval$Delta, na.rm = TRUE) else NA_real_,
    total_flips         = if (nrow(eval)) sum(eval$Flips, na.rm = TRUE) else 0L
  )
  
  # ---- Plotting (safe) -------------------------------------------------------
  if (isTRUE(plot)) {
    # 1) Panel plot
    if (!is.null(saa$flip_summary) && nrow(saa$flip_summary)) {
      if (exists("plot_flips_panel", mode = "function")) {
        tryCatch(
          print(plot_flips_panel(saa$flip_summary)),
          error = function(e) .msg("plot_flips_panel() failed: %s", e$message)
        )
      } else {
        .msg("plot_flips_panel() not found; skipping panel plot.")
      }
    } else {
      .msg("flip_summary is empty; skipping panel plot.")
    }
    
    # 2) Spatial flip plot for `fold_plot`
    if (fold_plot > nfolds) {
      .msg("Requested fold_plot=%d exceeds available folds (%d); skipping spatial plot.", fold_plot, nfolds)
    } else if (!is.null(saa$final_by_fold) && length(saa$final_by_fold) >= fold_plot) {
      if (exists("plot_flips_fold_like_misclassified", mode = "function")) {
        tryCatch(
          print(plot_flips_fold_like_misclassified(
            STEAM.obj,
            saa$final_by_fold,
            fold = fold_plot,
            layer_colors = layer_colors %||% STEAM.obj$viz$label_colors,
            fade = fade
          )),
          error = function(e) .msg("plot_flips_fold_like_misclassified() failed: %s", e$message)
        )
      } else {
        .msg("plot_flips_fold_like_misclassified() not found; skipping spatial plot.")
      }
    } else {
      .msg("final_by_fold missing or too short; skipping spatial plot.")
    }
  }
  
  invisible(STEAM.obj)
}