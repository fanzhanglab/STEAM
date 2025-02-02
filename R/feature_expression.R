#' Feature Expression
#'
#' @param steam_obj STEAM Obj
#' @param feature_name Feature Name
#' @param title title
#'
#' @export
feature_expression <- function(steam_obj, feature_name, title = "Expression Across Layers") {
  feature_expression <- data.frame(
    Layer = steam_obj$test$test.data.labels,
    Expression = steam_obj$test$avg.matrix[feature_name, ]
  )

  ggplot(feature_expression, aes(x = Layer, y = Expression, fill = Layer)) +
    geom_boxplot() +
    labs(
      title = title,
      x = "Layer/Label",
      y = "Expression Level",
      fill = "Layer/Label"
    ) +
    theme_classic()
}
