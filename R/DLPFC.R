#' LIBD Human DLPFC Spatial Transcriptomics Data (Sample 151673)
#'
#' A sample dataset containing human dorsolateral prefrontal cortex (DLPFC)
#' spatial transcriptomics data from the LIBD project, generated with the
#' 10x Genomics Visium platform (sample ID 151673).
#'
#' This dataset is included as an example to demonstrate STEAM workflows.
#' It is not required for package functionality.
#'
#' @docType data
#' @name DLPFC
#' @aliases dlpfc
#' @usage data(DLPFC)
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{spatial}{A data frame of cell spatial coordinates.}
#'   \item{count_exp}{A gene expression count matrix, with genes as rows and cells as columns.}
#'   \item{labels}{A factor or character vector of ground-truth labels.}
#' }
#'
#' @source LIBD: \url{http://research.libd.org/spatialLIBD/}
#'
"DLPFC"
