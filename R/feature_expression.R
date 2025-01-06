feature_expression <- function(steam_obj, feature_name, title = "Expression Across Layers") {
  feature_expression <- data.frame(
    Layer = steam_obj$test$test.data.labels,
    Expression = steam_obj$test$avg.matrix[feature_name, ]
  )

  ggplot(feature_expression, aes(x = Layer, y = Expression, fill = Layer)) +
    geom_boxplot() +
    labs(
      title = title,
      x = "Layer",
      y = "Expression Level"
    ) +
    theme_minimal()
}

feature_expression(STEAM.obj, "MBP", title = "Expression of MBP Across Layers")
