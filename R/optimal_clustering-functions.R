#############################################################
## ARBOL supported function to choose the optimal cluster  ##
#############################################################

# wrapper for chooseResolution_SilhouetteAnalysisParameterScan
#' @param n.drs number of PCs to include in clustering
#' @inheritParams ARBOL

#' @return vectors of PCs chosen
#' @importFrom Seurat RunPCA JackStraw ScoreJackStraw
#' @importFrom SeuratObject DefaultAssay Reductions
#' @noRd
.choose_optimal_clustering_default <- function(
        srobj,
        DownsampleNum = Inf,
        .normalize_data = .normalize_log1p,
        harmony_var = harmony_var,
        n.drs = srobj@misc$nDRs,
        res_scan_min = res_scan_min,
        res_scan_max = res_scan_max,
        res_scan_n = res_scan_n,
        res_scan_step = res_scan_step) {
    # get resolution from paramsweep of heavily downsampled dataset
    
    # downsample to downsample_num
    if (DownsampleNum <  ncol(srobj) & ncol(srobj) >= 100) {
        srobj <- srobj[, sample(colnames(srobj), DownsampleNum)]
        # Reprocess the data
        srobj <- .preprocess(srobj, harmony_var = harmony_var)
    }
    #setting resolution choice to arbitrary number in case of error
    resolution.choice <- 1
    
    tryCatch( { resolution.choice <- .choose_resolution_silhouette_scan(
        srobj,
        n.drs = srobj@misc$nDRs,
        res_scan_min = res_scan_min,
        res_scan_max = res_scan_max,
        res_scan_n = res_scan_n,
        res_scan_step = res_scan_step,
        harmony_var = harmony_var)
    
    }, error = function(e) {
        message('Silhouette Analysis failed. resolution set to 1.',
                ' WARNING: Double check your result.')
        message(paste("Silhouette Analysis error: ", e))
    })
    

    
    
    return(resolution.choice)
    
}
