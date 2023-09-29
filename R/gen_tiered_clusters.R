#' @name .gen_tiered_clusters
#' @title Performs iterative clustering (ARBOL) on a v4 Seurat object
#'
#' @description Iteratively clusters seurat object from single cell datasets,
#'   choosing optimum resolution parameters at each stage of clustering.
#' @param tier integer, defines at which iteration tier we are. Necessary for 
#'   the recursive approach to make sure we go through all the tiers.
#'   Starts from 0.
#' @param clustN integer, defines which cluster in a tier we are on.
#'    Necessary for the recursive approach to make sure we go through all the
#'    cluster in a tier. Starts from 0.
#' @param id_tc 0 or character vector indicating the tier and cluster the
#'   function is iterating on. ## TODO
#' @inheritParams ARBOL
#' @return list of lists with all seurat objects 
#' 
#' @examples
#' srobj <- readRDS("/path/to/seurat_object.rds")
#' tiers <- .gen_tiered_clusters(srobj)
#'
#' @rdname .gen_tiered_clusters
#' @noRd

.gen_tiered_clusters <- function(
        srobj,
        id_tc,
        tier = 0,
        clustN = 0,
        min_cluster_size = 100,
        max_tiers = 10,
        .choose_optimal_clustering_fun = .choose_optimal_clustering_default,
        .normalize_data = .normalize_log1p,
        res_scan_step = 5,
        res_scan_min = 0.01,
        res_scan_max = 1,
        res_scan_n = 40,
        reduction = "pca",
        harmony_var = NULL,
        DownsampleNum = 7500) {
    
    # Some input quality checks
    stopifnot(
        is(srobj, "Seurat"),
        is.character(id_tc) & length(id_tc) == 1,
        is.numeric(tier) & length(tier) == 1,
        is.numeric(clustN) & length(clustN) == 1,
        "nCount_RNA" %in% colnames(srobj@meta.data),
        "nFeature_RNA" %in% colnames(srobj@meta.data),
        is.numeric(min_cluster_size) & length(min_cluster_size) == 1,
        is.numeric(max_tiers) & length(max_tiers) == 1,
        is(.choose_optimal_clustering_fun, "function"),
        is(.normalize_data, "function"),
        is.numeric(res_scan_step) & length(res_scan_step) == 1,
        is.numeric(res_scan_min) & length(res_scan_min) == 1,
        is.numeric(res_scan_max) & length(res_scan_max) == 1,
        is.numeric(res_scan_n) & length(res_scan_n) == 1,
        is.numeric(DownsampleNum) & length(DownsampleNum) == 1,
        DownsampleNum < Inf,
        is.character(harmony_var) | is.null(harmony_var)
    )
    # keep track of position in tree for logging
    message("Number of cells: ", paste0(ncol(srobj), collapse = "\n"))
    message(paste0("Starting tier: ", tier, ", cluster: ", clustN,
                   ", identN: ", id_tc, ", with ", ncol(srobj), " cells" ))
    # ID with the tier and cluster
    if (id_tc == "") {
        id_tc <- paste0(id_tc, sprintf("T%sC%s", tier, clustN))
    } else {
        id_tc <- paste0(id_tc, sprintf("_T%sC%s", tier, clustN))
    }
    srobj@misc$tier <- tier
    srobj@misc$clustN <- clustN
    
    ######################################################################################################
    # check if too few cells or past tier 10. if so, return end-node without processing
    ######################################################################################################
    
    if (min_cluster_size < 5) {
        message(
        'WARNING: setting minimum cluster size to < 5 will throw errors in',
        'processing functions.\nThe full run should finish, but all processing',
        'of srobj < 5 cells is unreliable.\n',
        'Search "failure" in logs for information')
    }
    
    if ( (ncol(srobj) < min_cluster_size) | (tier > max_tiers) ) {
        message(cbind(
            "found end-node below min number of cells or above max tier. num cells: ",
            ncol(srobj),' < ', min_cluster_size))
        #add tierNident to metadata and stop recursion
        srobj@meta.data$tierNident <- gsub('/', '.', id_tc)
    }
    
    # Preprocessing: Run SCtransform or log transform, PCA, choose num PCs for
    #  downstream analysis, find neighbors.
    tryCatch({
        srobj <- .preprocess(
            srobj,
            .normalize_data = .normalize_data,
            harmony_var = harmony_var)
    },
    error = function(e) {
        message('Pre-processing failure')
        print(paste("Pre-processing error: ", e))
    })
    
    # Check if SCT ran successfully, if not, default back to log1p
    if (!("SCT" %in% names(srobj@assays)))
        cluster_assay <- "RNA"
    
    # Get optimal cluster resolution for clustering
    tryCatch({res <- .choose_optimal_clustering_fun(
        srobj,
        DownsampleNum = DownsampleNum,
        .normalize_data = .normalize_log1p,
        harmony_var = harmony_var,
        n.drs = srobj@misc$nDRs,
        res_scan_min = res_scan_min,
        res_scan_max = res_scan_max,
        res_scan_n = res_scan_n,
        res_scan_step = res_scan_step
        )
    },
    error = function(e) {
        message('Resolution choice failure')
        print(paste("Resolution choice error: ", e))
    })
    
    # Find clusters
    tryCatch({srobj <- FindClusters(srobj, resolution = res, verbose = FALSE)
    },
    error = function(e) {
        message('FindClusters failure');
        print(paste("Clustering error: ", e))
    })
    
    ######################################################################################################
    # check if too few clusters to call end-node and end recursion
    ######################################################################################################
    
    # if one cluster:
    ## TODO Idents doesn't seem the right approach here cause its not reset
    # by FindClusters :suspicious_look
    # Not happy with my approach, but it should be robust
    res_nm <- paste0("_snn_res.", res, "$")
    idx <- grep(res, colnames(srobj@meta.data))
    
    if ( length(unique(srobj@meta.data[, idx])) == 1 ) {
        message("found end-node with one cluster")
        srobj@meta.data$tierNident <- gsub('/', '.', id_tc)
        return(srobj)
    }
    
    if ((ncol(srobj) < min_cluster_size) | (tier > max_tiers)) {
        message(cbind("found end-node below min number of cells or above",
                      " max tier. num cells: ",
                      ncol(srobj),
                      ' < ',
                      min_cluster_size))
        return(srobj)
    }
    
    ######################################################################################################
    # recurse
    ######################################################################################################
    
    # split all clusters into separate srobjs
    message("continuing to recurse")
    Idents(srobj) <- srobj@meta.data[, idx]
    subsets <- SplitObject(srobj, split.by = "ident")
    
    # recurse along subsets
    return(
        lapply(seq_along(subsets),
               function(i) {
                   .gen_tiered_clusters(
                       # Pass one each of the new clusters as a new seurat object
                       subsets[[i]],
                       id_tc = id_tc,
                       # Increase the tier by 1 level from the previous
                       tier = tier + 1,
                       # evaluated cluster, -1 because of the 0 indexing
                       clustN = i - 1,
                       .normalize_data = .normalize_data,
                       .choose_optimal_clustering_fun = .choose_optimal_clustering_fun,
                       min_cluster_size = min_cluster_size,
                       max_tiers = max_tiers,
                       harmony_var = harmony_var
                      )
                  }))
}
