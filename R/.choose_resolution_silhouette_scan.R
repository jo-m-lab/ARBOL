#' Performs a parameter scan based on silhouette analysis on a seurat object
#'
#' @description Silhouette Analysis, a parameter scan method for choosing the most distinct clusters, is performed across 40 possible resolution on a seurat object. In ARBOL, srobj is downsampled for input to the parameter scan when there are >7.5k cells
#'
#' Also saves plots describing the parameter scan results
#'
#' @examples srobj <- readRDS("/path/to/seurat_object.rds")
#' srobj <- chooseResolution_SilhouetteAnalysisParameterScan(srobj)
#' 
#' @param input.srobj v4 seurat object
#' @param assay seurat object assay to use for silhouette analysis
#' @param n.pcs number of PCs to include in clustering
#' @param sample.name names resulting parameter scan results graphs
#' @param res.low lowest resolution to include in parameter scan
#' @param res.high highest resolution
#' @param res.n number of resolutions to test
#' @param res.step number of resolution steps before testing for local maximum. used to reduce runtime.
#' @param bias bias toward resolution choice, when resolution results are similar
#' @param figdir directory to save results graphs
#' @return seurat object with slot misc$resolution.choice
#' @noRd

.choose_resolution_silhouette_scan <- function(
        srobj,
        assay = cluster_assay,
        n.drs = srobj@misc$nDRs,
        res.low = res_scan_min,
        res.high = res_scan_max,
        res.n = res_scan_n,
        res.step = res_scan_step,
        bias = "under") {
    
    ######## step 1: save the input seurat object as a new temporary object, 
    ########         dont want to overwrite or change the original one with all of the parameter scans
    
    srobj.tmp <- srobj # TODO find a workaround bc we can't just double the memory requirements
    # in case there have been other things calculated in the metadata, just cut down to simplify/avoid errors
    srobj.tmp@meta.data <- srobj.tmp@meta.data[, c("nCount_RNA","nFeature_RNA")] # should just be the nUMI and nGene 
    
    ## TODO Scanning is over at max 500 cells - the dataset is already downsampled to DownsampleNum
    # if ( nrow(dist.temp) > 500) {
    #     random.cells.choose <- sample(1:nrow(dist.temp), round(nrow(dist.temp)/10, digits=0))
    # } else {
    #     random.cells.choose <- seq_len(nrow(dist.temp))
    # }
        
    dist.temp <- cor(t(srobj.tmp@reductions$harmony@cell.embeddings[, n.drs]), method = "pearson")
    # dist.temp <- dist.temp[random.cells.choose, random.cells.choose]
    sil.all.matrix <- matrix(data = NA, nrow = nrow(dist.temp), ncol = res.n)
    
    # set up empty additive data structures for silhouette scan
    set.res <- round(exp(seq(log(res.low), log(res.high), length.out = res.n)), digits=3)
    res.that.work <- rep(TRUE, length(set.res))
    
    ######## step 2: calculate the FindClusters over a large range of resolutions
    print("Performing parameter scan over multiple resolutions...")
    
    
    for (i in seq_len(length(set.res))) {
        
        ###### Step 2.1 computes FindClusters for each resolution individually
        # tryCatch allows testing resolutions that might cause errors in algorithm
        tryCatch({
            srobj.tmp <- FindClusters(srobj.tmp, resolution = set.res[i], verbose = FALSE)
        }, error = function(e) {res.that.work[i] <<- FALSE; message(paste("errored on", set.res[i], "flagging that resolution as it doesn't work"))})
        print(paste("          ", round(100*i/length(set.res)), "% done with parameter scan", sep=""))
        #if an error is encountered, skip this resolution and don't add it to the testing matrix
        if (res.that.work[i] == FALSE) {
            next
        }
        
        ######## step 3: calculate the silhouette width for each resolution (iterated in for loop)
        print("Computing a silhouette width for each cell")
        # Extract which metadata column corresponds to the resolution just computed 
        res <- paste0("_snn_res.", set.res[i], "$")
        indexes <- grep(res, colnames(srobj.tmp@meta.data))
        
        # Extract cluster for each cell at that resolution
        clusters.temp <- as.numeric(as.vector(srobj.tmp@meta.data[, indexes]))
        
        # Save silhouettes for each cell computed for that resolution's cluster
        if (length(table(clusters.temp)) > 1) {
            sil.out <- silhouette(clusters.temp, as.dist(1 - as.matrix(dist.temp)))
            sil.all.matrix[, i] <- sil.out[, 3]
        } else if (length(table(clusters.temp)) == 1) {
            sil.all.matrix[, i] <- 0
        }
    }
    
    
    # only perform calculations if you've done more iterations than res.step (res.step < 2 will throw errors)
    ## TODO don't understand code below why not just do this:
    # This returns a named vector with 1 character where the name is the
    # resolution and the value is the sil. score.
    resolution.choice <- sil.average[which.max(sil.average)]
    
    # if (i > res.step) {
    #     sil.average <- setNames(colMeans(sil.all.matrix), set.res[1:i])
    #     sil.medians <- setNames(apply(sil.all.matrix, 2, median), set.res[1:i])
    #     
    #     ######## step 4: automate choosing resolution that maximizes the silhouette 
    #     hist.out <- hist(sil.average, length(sil.average)/1.2,  plot=T)
    #     
    #     #  take the ones that fall into the top bin, 
    #     #  and the max OR MIN of those  ******* can change this to under vs over cluster
    #     ## TODO don't understand this...
    #     names(sil.average[which(sil.average > hist.out$breaks[length(hist.out$breaks)-1])])
    #     if (bias == "over") {
    #         resolution.choice <- as.numeric(max())
    #     }
    #     if (bias == "under") {
    #         resolution.choice <- as.numeric(min(
    #             names(sil.average[which(sil.average>hist.out$breaks[length(hist.out$breaks)-1])])))
    #     }
    #     
    #     if (sil.average[i] < sil.average[i-1] & 
    #         !(is.unsorted(-sil.average[(i-(res.step-1)):i], FALSE))) {
    #         message("found early stopping point in silhouette analysis")
    #         break
    #     }
    # }
    
    
    ######## step 5: calculate summary metric to compare the silhouette distributions,
    ########         average has worked well so far... could get fancier
    # get the silhouette of the best resolution: 
    silhouette.best <- as.numeric(resolution.choice)
    # Extract clustering with best silhouette
    res <- paste0("_snn_res.", names(resolution.choice), "$")
    res_idx <- grep(res, colnames(srobj.tmp@meta.data))
    
    n.clusters <- length(unique(srobj.tmp@meta.data[, res_idx]))
    
    print(paste("Best Resolution Choice: ", names(resolution.choice),
                ", with average silhouette score of: ", 
                round(silhouette.best, digits = 3),
                ", giving ", n.clusters,
                " clusters", sep = ""))
    
    ######## step 8: return the original seurat object, with the metadata containing a 
    ########         concatenated vector with the clusters defined by the best choice here,
    ########         as well as the ident set to this new vector
    # TODO we only need to return the best resolution, since the function
    # .choose_optimal_clustering_default - return(srobj@misc$resolution.choice)
    
    # Best.Clusters <- srobj.tmp@meta.data[, res_idx]
    # 
    # srobj.tmp$Best.Clusters <- Best.Clusters
    # Idents(srobj.tmp) <- srobj.tmp$Best.Clusters
    # srobj.tmp@misc <- list("resolution.choice" = resolution.choice)
    return(as.numeric(names(resolution.choice)))
}
