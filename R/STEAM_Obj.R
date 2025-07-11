#' STEAM object - Loading data for STEAM
#'
#' @param count_exp Count matrix/ Gene expression data
#' @param spatial Spatial Coordinates
#' @param labels Cluster info/annotations
#' @param Seurat.obj SEURAT Obj
#' @param pca PC embeddings (if using Seurat.obj)
#' @param umap UMAP embeddings (if using Seurat.obj)
#' @param clusters cluster labels (if using Seurat.obj)
#' @param label.column cluster labels
#' @param assay assay name (if using Seurat.obj)
#'
#' @return Object of class 'STEAM.Object'
#' @export
LoadSTEAM <- function(count_exp = NULL, spatial = NULL, labels = NULL,
                      pca = NULL, umap = NULL, clusters = NULL,
                      Seurat.obj = NULL, label.column = NULL, assay = "SCT") {

  STEAM.obj <- list()

  # ── Use Seurat object if provided ───────────────────────────────
  if (!is.null(Seurat.obj)) {
    # Extract scaled expression from specified assay
    if (!(assay %in% names(Seurat.obj@assays))) {
      stop(paste("Assay", assay, "not found in Seurat object."))
    }
    STEAM.obj$count_exp <- Seurat.obj@assays[[assay]]@scale.data

    # Extract spatial coordinates (supports Seurat v5+ and earlier)
    coords <- tryCatch({
      Seurat.obj@images$slice1@coordinates
    }, error = function(e) {
      NULL
    })

    if (!is.null(coords)) {
      STEAM.obj$spatial <- coords
    } else {
      message("Using GetTissueCoordinates() to extract spatial data.")
      coords <- Seurat::GetTissueCoordinates(Seurat.obj)
      rownames(coords) <- coords$cell
      coords$cell <- NULL
      STEAM.obj$spatial <- coords
    }

    # Label column is required
    if (is.null(label.column)) {
      stop("Please specify the label column from Seurat metadata using `label.column`.")
    }
    if (!(label.column %in% colnames(Seurat.obj@meta.data))) {
      stop(paste("Label column", label.column, "not found in Seurat metadata."))
    }
    STEAM.obj$labels <- as.character(Seurat.obj@meta.data[[label.column]])

    # Add reductions if available
    if ("pca" %in% names(Seurat.obj@reductions)) {
      STEAM.obj$pca <- Seurat.obj@reductions[["pca"]]@cell.embeddings
    }
    if ("umap" %in% names(Seurat.obj@reductions)) {
      STEAM.obj$umap <- Seurat.obj@reductions[["umap"]]@cell.embeddings
    }

    # Use labels as clusters if none are specified
    STEAM.obj$clusters <- as.character(STEAM.obj$labels)

    # ── Manual input mode ───────────────────────────────────────────
  } else {
    if (is.null(count_exp) || is.null(spatial) || is.null(labels)) {
      stop("Must provide `count_exp`, `spatial`, and `labels` if not using Seurat.obj.")
    }

    STEAM.obj$count_exp <- count_exp
    STEAM.obj$spatial <- spatial
    STEAM.obj$labels <- as.character(labels)

    if (!is.null(pca))    STEAM.obj$pca <- pca
    if (!is.null(umap))   STEAM.obj$umap <- umap
    if (!is.null(clusters)) {
      STEAM.obj$clusters <- as.character(clusters)
    } else {
      STEAM.obj$clusters <- as.character(labels)
    }
  }

  # Final check for NA in cluster assignments
  if (anyNA(STEAM.obj$clusters)) {
    stop("STEAM.obj$clusters contains NA values. Please check your label or cluster input.")
  }

  class(STEAM.obj) <- "STEAM.Object"
  return(STEAM.obj)
}
