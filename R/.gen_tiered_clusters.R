#' @name .gen_tiered_clusters
#' @title Performs iterative clustering (ARBOL) on a v4 Seurat object
#'
#' @description Iteratively clusters seurat object from single cell datasets,
#'   choosing optimum resolution parameters at each stage of clustering.
#'#' @param tier integer, defines at which iteration tier we are. Necessary for 
#'   the recursive approach to make sure we go through all the tiers.
#'   Starts from 0.
#' @param clustN integer, defines which cluster in a tier we are on.
#'    Necessary for the recursive approach to make sure we go through all the
#'    cluster in a tier. Starts from 0.
#' @inheritParams ARBOL
#' @return list of lists with all seurat objects 
#' 
#' @examples
#' srobj <- readRDS("/path/to/seurat_object.rds")
#' tiers <- GenTieredclusters(srobj)
#'
#' @rdname .gen_tiered_clusters
#' @export
#' 
.gen_tiered_clusters <- function(
        srobj,
        PreProcess_fun = PreProcess_sctransform,
        tier = 0,
        clustN = 0,
        min_cluster_size = 100,
        max_tiers = 10,
        ChooseOptimalClustering_fun = ChooseOptimalClustering_default,
        res_scan_step = 5,
        res_scan_min = 0.01,
        res_scan_max = 1,
        res_scan_n = 40,
        harmony_var = NULL,
        DownsampleNum = 7500) {
    
    # keep track of position in tree for logging
    message("Number of cells: ", paste0(ncol(srobj), collapse = "\n"))
    message(paste0("Starting tier: ", tier, ", with ", ncol(srobj), " cells" ))
    
    ######################################################################################################
    # check if too few cells or past tier 10. if so, return end-node without processing
    ######################################################################################################
    
    if (min_cluster_size < 5) {
        message(
        'WARNING: setting minimum cluster size to < 5 will throw errors in processing functions. \n 
        the full run should finish, but all processing of srobj < 5 cells is unreliable. \n
        Search "failure" in logs for information')
    }
    
    if ( (ncol(srobj) < min_cluster_size) | (srobj@misc$tier > max_tiers) ) {
        message(cbind(
            "found end-node below min number of cells or above max tier. num cells: ",
            ncol(srobj),' < ', min_cluster_size))
        #add tierNident to metadata and stop recursion
        srobj@meta.data$tierNident <- SaveEndFileName %>% str_replace_all('/','.')
        return(srobj)
    }
    
    # Preprocessing: Run SCtransform or log transform, PCA, choose num PCs for
    #  downstream analysis, find neighbors.
    tryCatch({srobj <- PreProcess_fun(srobj, regressVar = harmony_var)
    },
    error = function(e) {
        message('Pre-processing failure')
        print(paste("Pre-processing error: ", e))
    })
    
    # Check if SCT ran successfully, if not, default back to log1p
    if(!(cluster_assay %in% names(srobj@assays))){
        cluster_assay = "RNA"
    }
    
    # Get optimum cluster resolution and cluster
    tryCatch({res <- ChooseOptimalClustering_fun(
        srobj,
        PreProcess_fun = PreProcess_fun,
        res.low = res_scan_min,
        res.high = res_scan_max,
        res.n = res_scan_n,
        res.step = res_scan_step,
        harmony_var = harmony_var,
        downsample_num = DownsampleNum,
        )
    },
    error = function(e) {
        message('Resolution choice failure')
        print(paste("Resolution choice error: ", e))
    })
    tryCatch({srobj <- FindClusters(srobj, resolution = res)
    },
    error = function(e) {
        message('FindClusters failure');
        print(paste("Clustering error: ", e))
    })
    
    ######################################################################################################
    # check if too few clusters to call end-node and end recursion
    ######################################################################################################
    
    # if one cluster:
    if ( !(length(unique(Idents(srobj))) > 1) ) {
        message("found end-node with one cluster")
        srobj@meta.data$tierNident <- SaveEndFileName %>% str_replace_all('/','.')
        return(srobj)
    }
    
    ######################################################################################################
    # recurse
    ######################################################################################################
    
    # split all clusters into separate srobjs
    message("continuing to recurse")
    
    subsets <- SplitSrobjOnIdents(srobj, paste0("tier", tier))
    
    print(subsets)
    # recurse along subsets
    return(lapply(seq_along(subsets),
                  function(i) {
                      .gen_tiered_clusters(
                          # Pass one each of the new clusters as a new seurat object
                          subsets[[i]],
                          # Increase the tier by 1 level from the previous
                          tier = tier + 1,
                          # evaluated cluster
                          clustN = i - 1,
                          PreProcess_fun = PreProcess_fun,
                          ChooseOptimalClustering_fun = ChooseOptimalClustering_fun,
                          min_cluster_size = min_cluster_size,
                          max_tiers = max_tiers,
                          harmony_var = harmony_var
                      )
                  }))
}
