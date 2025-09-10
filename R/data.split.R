#' Split STEAM object into training and testing sets
#'
#' @param STEAM.obj A STEAM object with expression, labels, and spatial slots
#' @param train.ratio Numeric (0–1), proportion of data used for training
#' @return Modified STEAM object with train/test subslots
#' @keywords internal
data.split <- function(STEAM.obj, train.ratio) {

  # Universal labels
  clean_labels <- function(labels) {
    labels <- as.character(labels)
    labels <- gsub("/", ".", labels)
    labels <- gsub("-", ".", labels)
    labels <- gsub(" ", "_", labels)
    labels <- gsub("[()]", "", labels)
    labels <- gsub(",", "", labels)
    labels <- gsub("\\.\\.", ".", labels)
    labels <- trimws(labels)
    labels <- iconv(labels, to = "ASCII//TRANSLIT")
    return(labels)
  }

  STEAM.obj$labels <- clean_labels(STEAM.obj$labels)

  # Create train/test split
  ind <- createDataPartition(STEAM.obj$labels, p = train.ratio, list = FALSE)

  STEAM.obj$train$train.data.matrix <- STEAM.obj$count_exp[, ind]
  STEAM.obj$train$train.data.coords <- STEAM.obj$spatial[ind, ]
  STEAM.obj$train$train.data.labels <- STEAM.obj$labels[ind]

  STEAM.obj$test$test.data.matrix <- STEAM.obj$count_exp[, -ind]
  STEAM.obj$test$test.data.coords <- STEAM.obj$spatial[-ind, ]
  STEAM.obj$test$test.data.labels <- STEAM.obj$labels[-ind]

  return(STEAM.obj)
}
