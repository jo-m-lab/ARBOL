#' Performs a parameter scan based on silhouette analysis on a seurat object
#'
#' @description Silhouette Analysis, a parameter scan method for choosing the most distinct clusters, is performed across 40 possible resolution on a seurat object. In ARBOL, srobj is downsampled for input to the parameter scan when there are >7.5k cells
#'
#' Also saves plots describing the parameter scan results
#'
#' @examples srobj <- readRDS("/path/to/seurat_object.rds")
#' srobj <- chooseResolution_SilhouetteAnalysisParameterScan(srobj)
#' 
#' @param n.drs number of PCs to include in clustering
#' @param bias bias toward resolution choice, when resolution results are similar
#' @inheritParams ARBOL
#' @return seurat object with slot misc$resolution.choice
#' @noRd

#' @importFrom Seurat FindClusters

.choose_resolution_silhouette_scan <- function(
        srobj,
        n.drs = srobj@misc$nDRs,
        res_scan_min = res_scan_min,
        res_scan_max = res_scan_max,
        res_scan_n = res_scan_n,
        res_scan_step = res_scan_step,
        harmony_var = NULL,
        bias = "under") {
    
    ######## step 1: save the input seurat object as a new temporary object, 
    ########         dont want to overwrite or change the original one with all of the parameter scans
    
    # in case there have been other clusterings calculated in the metadata,
    # just cut down to simplify/avoid errors, we have previouslt checked these
    # are present in the metadata
    srobj@meta.data <- srobj@meta.data[, c("CellID", "nCount_RNA","nFeature_RNA")]
    
    
    ## TODO harmony is called but not necessary
    if (is.null(harmony_var)) {
        red <- "pca"
    } else {
        red <- "harmony"
    }
    ## TODO maybe change to dist() for euclidean - important because correlation the higher the value the better and inverse with euclidean distance
    # dist.temp <- cor(t(srobj@reductions[[red]]@cell.embeddings[, n.drs]),
    #                  method = "pearson")
    dist.temp <- dist(srobj@reductions[[red]]@cell.embeddings[, n.drs],
                     method = "euclidean")
    # dist.temp <- dist.temp[random.cells.choose, random.cells.choose]
    sil.all.matrix <- matrix(data = NA, nrow = ncol(srobj), ncol = res_scan_n)
    
    # set up empty additive data structures for silhouette scan
    set.res <- round(exp(seq(log(res_scan_min), log(res_scan_max), length.out = res_scan_n)),
                     digits = 3)
    res.that.work <- rep(TRUE, length(set.res))
    
    ######## step 2: calculate the FindClusters over a large range of resolutions
    message("Performing parameter scan over multiple resolutions...")
    
    
    for (i in seq_len(length(set.res))) {
        
        ###### Step 2.1 computes FindClusters for each resolution individually
        # tryCatch allows testing resolutions that might cause errors in algorithm
        tryCatch({
            srobj <- FindClusters(srobj, resolution = set.res[i], verbose = FALSE)
        }, error = function(e) {
            res.that.work[i] <<- FALSE
            message(
                paste("errored on", set.res[i],
                      "flagging that resolution as it doesn't work"))})
        # print(paste("          ", round(100*i/length(set.res)),
        #             "% done with parameter scan", sep=""))
        #if an error is encountered, skip this resolution and don't add it to the testing matrix
        if (res.that.work[i] == FALSE) {
            next
        }
        
        ######## step 3: calculate the silhouette width for each resolution (iterated in for loop)
        # print("Computing a silhouette width for each cell")
        # Extract which metadata column corresponds to the resolution just computed 
        res <- paste0("_snn_res.", set.res[i], "$")
        indexes <- grep(res, colnames(srobj@meta.data))
        
        # Extract cluster for each cell at that resolution
        clusters.temp <- as.numeric(as.vector(srobj@meta.data[, indexes]))
        
        # Save silhouettes for each cell computed for that resolution's cluster
        if (length(table(clusters.temp)) > 1) {
            # sil.out <- silhouette(clusters.temp, as.dist(1 - as.matrix(dist.temp))) # For correlation
            sil.out <- silhouette(clusters.temp, dist.temp) # For euclidean
            sil.all.matrix[, i] <- sil.out[, 3]
        } else if (length(table(clusters.temp)) == 1) {
            sil.all.matrix[, i] <- 0
        }
        
        ## TODO kyle please explain better what is going on here
        # only perform calculations if you've done more iterations than res.step (res.step < 2 will throw errors)
        if (i > res_scan_step) {
            sil.average <- setNames(colMeans(sil.all.matrix), set.res[1:i])
            sil.medians <- setNames(apply(sil.all.matrix, 2, median), set.res[1:i])
    
            ######## step 4: automate choosing resolution that maximizes the silhouette
            # TODO this can be simplified to the commented line below
            sil_keep <- which(sil.average == max(sil.average, na.rm = TRUE))
            sil_max <- sil.average[sil_keep]
            
            if (bias == "over") {
                resolution.choice <- sil_max[length(sil_max)]
            }
            if (bias == "under") {
                resolution.choice <- sil_max[1]
            }
    
            if (sil.average[i] < sil.average[i-1] &
                !(is.unsorted(-sil.average[(i-(res_scan_step-1)):i], FALSE))) {
                message("found early stopping point in silhouette analysis")
                break
            }
        }
    }
    
    ######## step 5: calculate summary metric to compare the silhouette distributions,
    ########         average has worked well so far... could get fancier
    # get the silhouette of the best resolution: 
    silhouette.best <- as.numeric(resolution.choice)
    # Extract clustering with best silhouette
    res <- paste0("_snn_res.", names(resolution.choice), "$")
    res_idx <- grep(res, colnames(srobj@meta.data))
    
    n.clusters <- length(unique(srobj@meta.data[, res_idx]))
    
    message(paste("Best Resolution Choice: ", names(resolution.choice),
                ", with silhouette score of: ", 
                round(silhouette.best, digits = 3),
                ", giving ", n.clusters,
                " clusters", sep = ""))
    
    # Only need to return the best resolution, since the function
    return(as.numeric(names(resolution.choice)))
}
