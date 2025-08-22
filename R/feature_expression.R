#' Feature Expression
#'
#' @param steam_obj STEAM Obj
#' @param feature_name Feature Name
#' @param title title
#'
#' @export
feature_expression <- function(steam_obj, feature_name,
                               view = c("pooled","facet"), title = "Expression Across Layers") {
  suppressPackageStartupMessages(library(ggplot2))
  view <- match.arg(view)
  
  # -------- Simple CV path (same as before, shows test split) --------
  if (!is.null(steam_obj$test) && !is.null(steam_obj$test$avg.matrix)) {
    mat  <- steam_obj$test$avg.matrix                  # features x test cells
    labs <- steam_obj$test$test.data.labels            # factor/char len = ncol(mat)
    
    if (!feature_name %in% rownames(mat))
      stop(sprintf("Feature '%s' not found in test$avg.matrix.", feature_name))
    
    expr <- as.numeric(mat[feature_name, ])
    df <- data.frame(Layer = factor(as.character(labs)), Expression = expr)
    
    p <- ggplot(df, aes(x = Layer, y = Expression, fill = Layer)) +
      geom_boxplot(outlier.shape = NA, width = 0.7) +
      geom_jitter(width = 0.15, alpha = 0.35, size = 0.8, color = "black") +
      labs(title = title, x = "Layer / Label", y = paste("Expression:", feature_name)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                              legend.position = "none")
    print(p); return(invisible(p))
  }
  
  # -------- Nested CV path: use ONLY held-out cells from each outer fold --------
  if (is.null(steam_obj$nested) || is.null(steam_obj$nested$ncv) ||
      is.null(steam_obj$count_exp) || is.null(steam_obj$labels))
    stop("For nested mode, need nested$ncv, count_exp, and labels in the STEAM object.")
  
  ncv <- steam_obj$nested$ncv
  if (is.null(ncv$outer_result) || !length(ncv$outer_result))
    stop("No outer folds found (nested$ncv$outer_result missing).")
  
  mat_all  <- steam_obj$count_exp
  labs_all <- steam_obj$labels
  if (is.null(rownames(mat_all))) stop("count_exp must have rownames (features).")
  if (is.null(colnames(mat_all))) stop("count_exp must have colnames (cell IDs).")
  
  hit <- which(rownames(mat_all) == feature_name)
  if (!length(hit)) {
    ci <- which(tolower(rownames(mat_all)) == tolower(feature_name))
    if (length(ci)) hit <- ci[1]
  }
  if (!length(hit)) {
    m <- grep(feature_name, rownames(mat_all), ignore.case = TRUE, value = TRUE)
    msg <- if (length(m)) paste0(" Did you mean: ", paste(utils::head(m, 8), collapse = ", "), "?") else ""
    stop(sprintf("Feature '%s' not found in count_exp.%s", feature_name, msg))
  }
  
  .fold_df <- function(of_idx) {
    p <- ncv$outer_result[[of_idx]]$preds
    if (is.null(p) || !nrow(p)) return(NULL)
    
    if (!"cell_id" %in% colnames(p)) {
      if (is.null(rownames(p))) stop("Preds lack rownames and 'cell_id'; cannot align.")
      p$cell_id <- rownames(p)
    }
    test_ids <- as.character(p$cell_id)
    
    j <- match(test_ids, colnames(mat_all))
    keep <- which(!is.na(j))
    if (!length(keep)) return(NULL)
    
    expr <- as.numeric(mat_all[hit, j[keep], drop = TRUE])
    
    if ("testy" %in% colnames(p)) {
      labs <- as.character(p$testy[keep])
    } else {
      if (!is.null(names(labs_all))) {
        labs <- as.character(labs_all[test_ids[keep]])
      } else {
        labs <- as.character(labs_all[j[keep]])
      }
    }
    
    data.frame(
      OuterFold = factor(paste0("Fold ", of_idx)),
      Layer     = factor(labs),
      Expression= expr,
      stringsAsFactors = FALSE
    )
  }
  
  folds <- seq_along(ncv$outer_result)
  df_list <- lapply(folds, .fold_df)
  df_all  <- do.call(rbind, df_list)
  if (is.null(df_all) || !nrow(df_all)) stop("Could not assemble held-out cells for any fold.")
  
  base <- ggplot(df_all, aes(x = Layer, y = Expression, fill = Layer)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.35, size = 0.8, color = "black", show.legend = FALSE) +
    labs(title = title, x = "Layer / Label", y = paste("Expression:", feature_name)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  p <- if (view == "facet") {
    base + facet_wrap(~ OuterFold, nrow = 1)
  } else {
    base
  }
  
  print(p)
  invisible(p)
}
