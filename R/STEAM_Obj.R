#' STEAM object - Loading data for STEAM
#'
#' @param count_exp Count matrix/ Gene expression data
#' @param spatial Spatial Coordinates
#' @param labels Cluster info/annotations
#' @param Seurat.obj SEURAT Obj
#'
#' @return Object of class 'STEAM.Object'
#' @export
LoadSTEAM <- function(count_exp = NULL, spatial = NULL, labels = NULL, Seurat.obj = NULL) {
  # Initialize STEAM object as a list
  STEAM.obj <- list()

  # Check if a Seurat object is provided
  if (is.null(Seurat.obj)) {
    # Use provided inputs directly
    STEAM.obj$count_exp <- count_exp
    STEAM.obj$spatial <- spatial
    STEAM.obj$labels <- labels
  } else {
    # Extract data from the Seurat object
    STEAM.obj$count_exp <- Seurat.obj@assays$SCT@scale.data
    STEAM.obj$spatial <- Seurat.obj@images$slice1@coordinates
    STEAM.obj$labels <- Seurat.obj$seurat_clusters
  }

  # Assign class to the STEAM object
  class(STEAM.obj) <- 'STEAM.Object'

  # Return the STEAM object
  return(STEAM.obj)
}

