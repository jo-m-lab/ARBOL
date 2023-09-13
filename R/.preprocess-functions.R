#############################################
## ARBOL supported normalization functions ##
#############################################
.normalize_log1p <- function(srobj) {
    NormalizeData(object = srobj, verbose = FALSE)
    }

.normalize_sct <- function(srobj) {
    SCTransform(object = srobj, verbose = FALSE, method = "glmGamPoi")
    }

##################################
## ARBOL supported HVG choosing ##
##################################
# .choose_hvg_default <- function(srobj) {
#     FindVariableFeatures(srobj, nfeatures=2000, selection.method = "disp")
# }

#############################################
## ARBOL supported dimensionality setting  ##
#############################################

# Chooses the number of principle components to include in clustering
# Seurat's JackStraw is used when there are less than 500 cells in a
# clustering event. Otherwise, PCs are used in clustering if they explain 15%
# more variance than the following PC.
#
#
# examples 
# nPCs <- ChoosePCs_default(srobj, figure_dir=fig_dir)
#' 
#' param srobj v4 seurat object
#' param improved_diff_quantile percent variance explained of next PC to choose PC
#' param significance JackStraw significance required to choose PC
#' param reduction name of the reduction slot where to select the reduced
#'   dimensions for downstream analysis
#' return vectors of PCs chosen
#' @importFrom Seurat RunPCA JackStraw ScoreJackStraw
#' @importFrom SeuratObject DefaultAssay Reductions
#' 
.choose_dims_default <- function(
        srobj,
        improved_diff_quantile = .85,
        significance = 0.01,
        reduction = "pca") {
    
    # Some input quality checks
    stopifnot(
        is(srobj, "Seurat"),
        is.numeric(improved_diff_quantile) & length(improved_diff_quantile == 1),
        is.numeric(significance)  & length(significance == 1),
        is.character(reduction) & reduction %in% Reductions(srobj))
    
    # if num cells > 500 use elbow plot if small use jackstraw, 
    #   if none sig use first 2 pcs
    if ((reduction == "pca" & ncol(srobj) > 500) | reduction != "pca") {
        eigValues <- (srobj[[reduction]]@stdev)^2  ## EigenValues
        varExplained <- eigValues / sum(eigValues)
        nPCs <- 1:max(
            which(diff(varExplained) < 
                      quantile(diff(varExplained),
                               1 - improved_diff_quantile)) + 1)
    } else {
        # Change default assay to run PCA on RNA
        # DefaultAssay(srobj) <- 'RNA' ## TODO Is this needed? We can run JackStraw on any normalized data
        #Run PCA under RNA for Jackstraw
        srobj <- RunPCA(
            object = srobj,
            npcs = min(50, round(ncol(srobj)/2)),
            verbose = FALSE)
        suppressWarnings({srobj <- JackStraw(srobj, dims = 50)})
        srobj <- ScoreJackStraw(
            srobj,
            dims = 1:ncol(srobj@reductions$pca@cell.embeddings))
        # Extract number of PCs
        nPCs <- which(JS(srobj$pca)@overall.p.values[, "Score"] < significance)

        # DefaultAssay(srobj) <- 'SCT' ## TODO related to the comment 7 lines above
    }
    
    # Set vector for corner case 2
    if (length(nPCs) <= 2) {
        nPCs <- seq_len(2)
    }
    return(nPCs)
}


####################################
## ARBOL supported preprocessing  ##
####################################
# Helper function to process the seurat object at each step
#' @importFrom Seurat FindVariableFeatures ScaleData RunPCA FindNeighbors
.preprocess <- function(
        srobj,
        .normalize = .normalize_log1p,
        .choose_hvg = .choose_hvg_default,
        harmony_var = NULL,
        assay.use = DefaultAssay(srobj)) {
    
    # Some input quality checks
    stopifnot(
        is(srobj, "Seurat"),
        is(.normalize, "function"),
        is(.choose_hvg, "function"),
        is.character(harmony_var) | is.null(harmony_var),
        is.character(assay.use))
    
    # Standard data preprocessing
    srobj <- .normalize_data(object = srobj)
    srobj <- FindVariableFeatures(
        srobj,
        nfeatures = 2000,
        selection.method = "disp")
    srobj <- ScaleData(object = srobj, features = VariableFeatures(srobj))
    srobj <- RunPCA(
        object = srobj,
        npcs = min(50, round(ncol(srobj) / 2)),
        verbose = FALSE)
    # Specify this for the .choose_dims_default step
    reduction <- "pca"
    
    # Run Harmony integration if specified
    if (!is.null(harmony_var)) {
        srobj <- RunHarmony(
            object = srobj,
            group.by.vars = harmony_var,
            reduction = "pca")
        # Specify this for the .choose_dims_default step
        reduction <- "harmony"
    }
    
    # Choose number of reduced dimensions
    nPCs <- .choose_dims_default(srobj, reduction = reduction)
    
    # Save the PCs to use in misc
    srobj@misc$nPCs <- nPCs
    message(
        paste0("chose ", length(nPCs), "PCs, added choices to srobj@misc$nPCs"))
    
    # Compute KNN graph
    srobj <- FindNeighbors(
        object = srobj,
        dims = nPCs,
        k.param = ceiling(0.5*sqrt(ncol(srobj))),
        reduction = reduction,
    )
    
    return(srobj)
}