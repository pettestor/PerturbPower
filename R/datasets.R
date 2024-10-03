#' Example gRNA Counts Data
#'
#' This dataset contains example gRNA counts used in the analysis.
#' Each row corresponds to a gRNA, and each column corresponds to counts.
#'
#' @format A data frame with X rows and Y variables:
#' \describe{
#'   \item{target_gene_symbol}{A unique identifier for each gene}
#'   \item{count}{The count of each gRNA in the sample.}
#'   \item{X}{Any additional variables relevant to your dataset.}
#'   \item{ID}{A unique identifier for each gRNA.}
#' }
#' @source Simulated data, or provide a real data source.
"gRNA_counts"
#' Example Seurat Metadata
#'
#' This dataset contains metadata for a Seurat object used in the analysis.
#' Each row corresponds to a cell, and the columns provide information about various features.
#'
#' @format A data frame with X rows and the following variables:
#' \describe{
#'   \item{CC.Difference}{Numeric, representing the difference in cell cycle phases.}
#'   \item{Cycling}{Logical, indicating whether the cell is cycling.}
#'   \item{G2M.Score}{Numeric, the G2M phase score.}
#'   \item{G2M.Score.log}{Numeric, log-transformed G2M phase score.}
#'   \item{NamedClusters}{Factor, named clusters assigned to cells.}
#'   \item{Phase}{Factor, representing the cell cycle phase (e.g., G1, S, G2M).}
#'   \item{RNA_snn_res.0.07}{Numeric, RNA single-nearest neighbor resolution.}
#'   \item{S.Score}{Numeric, S phase score.}
#'   \item{S.Score.log}{Numeric, log-transformed S phase score.}
#'   \item{X}{Numeric, another feature (please describe what this represents).}
#'   \item{cell_id}{Character, unique identifier for each cell.}
#'   \item{groups}{Factor, groups or categories assigned to the cells.}
#'   \item{nCount_RNA}{Numeric, RNA molecule count per cell.}
#'   \item{nFeature_RNA}{Numeric, number of features (genes) detected per cell.}
#'   \item{orig.ident}{Factor, original identity or sample identifier.}
#'   \item{seurat_clusters}{Factor, Seurat-assigned cluster identities.}
#' }
#' @source Simulated data, or provide a real data source.
"seurat_obj.metadata"
