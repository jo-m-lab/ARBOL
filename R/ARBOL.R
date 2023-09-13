#' @name ARBOL
#' @title Wraps GenTieredClusters in a larger script to output a cell-level
#'    dataframe of cluster identity through the tree.
#' 
#' @description Performs ARBOL tiered hierarchical clustering on a seurat
#'    object. Outputs a dataframe of tier membership per cell
#' 
#' @param srobj \code{SeuratObject} single cell dataset.
#' @param PreProcess_fun function to use for preproccessing defaults to 
#'   PreProcess_sctransform. PreProcess_sctransform_harmony now available.
#' @param min_cluster_size integer, minimum number of cells to allow further
#'   clustering. Defaults to 100.
#' @param max_tiers integer, maximum number of tiers to allow for recursive
#'   clustering. Defaults to 10.
#' @param ChooseOptimalClustering_fun function that returns srobj with clusters
#'   in `srobj$Best.Clusters` after choosing optimal clustering resolution.
#' @param res_scan_step integer, number of resolutions of decreasing score
#'   before an early silhouette analysis resolution scan stop. Defaults to 5.
#' @param res_scan_min integer, smallest resolution to scan. Defaults to 0.01.
#' @param res_scan_max integer, highest resolution to scan. Defaults to 1.
#' @param res_scan_n integer, number of resolutions to scan between
#'   \code{res_scan_min} and \code{res_scan_max} using a logarithmic sequence.
#' @param harmony_var character string, variable over which to run harmony
#'   integration.
#' @param DownsampleNum integer, number of cells to downsample to when running
#'  the silhouette scan.
#' @return dataframe with tierN cluster membership per cell
#' 
#' @details
#' ARBOL aims to obtain tiered hierarchical clustering of all the cells in the
#'   dataset. To do so it iteratively subclusters the data, renormalizes, 
#'   finds variable features, scales the data computes PC and finds clusters at
#'   various resolutions. Then it assesses the robustness of the different
#'   cluster resolutions by using the silhouette score. Once that resolution is
#'   chosen each cluster is processed as previouslt described. Stopping criteria
#'   when it gets 1. below the minimum number of cells specified 2. if the best 
#'   resolution returns 1 cluster.
#' 
#' @examples 
#' endclustDF <- ARBOL(srobj = srobj,
#'                            saveSROBJdir = "~/tieredoutput/srobjs",
#'                            figdir = "~/tieredoutput/figdir",
#'                            SaveEndNamesDir = "~/tieredoutput/endclusts")
#' @rdname ARBOL
#' @export
#' 
ARBOL <- function(
        srobj,
        .normalize = .normalize_log1p,
        min_cluster_size = 100,
        max_tiers = 10,
        ChooseOptimalClustering_fun = ChooseOptimalClustering_default,
        res_scan_step = 5,
        res_scan_min = 0.01,
        res_scan_max = 1,
        res_scan_n = 40,
        harmony_var = NULL,
        DownsampleNum = 7500) {
    
    # Some input quality checks
    stopifnot(
        is(srobj, "Seurat"),
        is(.normalize, "function"),
        is(ChooseOptimalClustering_fun, "function"),
        is.numeric(min_cluster_size) & length(min_cluster_size) == 1,
        is.numeric(max_tiers) & length(max_tiers) == 1,
        is.numeric(res_scan_step) & length(res_scan_step) == 1,
        is.numeric(res_scan_min) & length(res_scan_min) == 1,
        is.numeric(res_scan_max) & length(res_scan_max) == 1,
        is.numeric(res_scan_n) & length(res_scan_n) == 1,
        is.numeric(DownsampleNum) & length(DownsampleNum) == 1,
        is.character(harmony_var) | is.null(harmony_var),
        is.character(assay.use))
    
    tieredsrobjs <- GenTieredClusters(
        srobj = srobj,
        PreProcess_fun = .preprocess,
        min_cluster_size = min_cluster_size,
        max_tiers = max_tiers,
        ChooseOptimalClustering_fun = ChooseOptimalClustering_fun,
        res_scan_step = res_scan_step,
        res_scan_min = res_scan_min,
        res_scan_max = res_scan_max,
        res_scan_n = res_scan_n,
        harmony_var = harmony_var,
        DownsampleNum = DownsampleNum)

  srobjslist <- unlist(tieredsrobjs)
  tryCatch({
    tierNidents <- lapply(srobjslist, function(tsrobj) {
      tryCatch({
        idents <- data.frame(
            CellID = row.names(tsrobj@meta.data),
            tierNident = tsrobj@meta.data$tierNident)
        return(idents)
      })
    })
  },
    error = function(e) {message('unnesting srobjs failed, probably because clustering failed. reason:'); print(e)}
  )

  endclustDF <- bind_rows(tierNidents)
  endclustDF$tierNident <- sub('_T0C0_', '',endclustDF$tierNident)
  endclustDF$tierNident <- gsub('_','.',endclustDF$tierNident)

  endclustDF <- endclustDF %>% 
      tidyr::separate(
          tierNident,
          into = paste0('tier', 1:max_tiers),
          sep = '\\.',
          remove = FALSE)

  return(endclustDF)
}
