#############################################################
## ARBOL supported function to choose the optimal cluster  ##
#############################################################

# wrapper for chooseResolution_SilhouetteAnalysisParameterScan
#
#
# examples 
# nPCs <- ChoosePCs_default(srobj, figure_dir=fig_dir)
#' 
#' @param srobj \code{Seurat} single cell object
#' @param DownsampleNum integer, number of cells to downsample to when running
#'   the silhouette scan.
#' @param .normalize_data function that return a normalized seurat object. By 
#'   default \code(.normalize_log1p) but can also choose \code(.normalize_sct)
#'   or pass your own.
#' @param harmony_var character string, variable over which to run harmony
#'   integration.
#' @return vectors of PCs chosen
#' @importFrom Seurat RunPCA JackStraw ScoreJackStraw
#' @importFrom SeuratObject DefaultAssay Reductions
#' @noRd
.choose_optimal_clustering_default <- function(
        srobj,
        downsample_num = Inf,
        PreProcess_fun = .normalize_log1p,
        harmony_var = harmony_var,
        ...) {
    # get resolution from paramsweep of heavily downsampled dataset
    
    # downsample to downsample_num
    if (downsample_num <  ncol(srobj) & ncol(srobj) >= 100) {
        srobj <- srobj[, sample(colnames(srobj), downsample_num)]
        srobj <- .preprocess(srobj, harmony_var = harmony_var)
    }
    #setting resolution choice to arbitrary number in case of error
    ndrsave <- srobj@misc$nDRs # TODO is this needed?
    srobj@misc$resolution.choice <- 1
    srobj@misc <- list("resolution.choice" = resolution.choice)
    
    tryCatch( { srobj <- .choose_resolution_silhouette_scan(srobj, ...)
    }, error = function(e) {
        message('Harmony Silhouette Analysis failed. resolution set to 1.',
                ' WARNING: Double check your result.')
        message(paste("Harmony Silhouette Analysis error: ", e))
    })
    

    
    
    return(srobj@misc$resolution.choice)
    
}
