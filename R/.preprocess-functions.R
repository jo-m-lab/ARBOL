#############################################
## ARBOL supported normalization functions ##
#############################################
#' @noRd
#' @importFrom Seurat NormalizeData

.normalize_log1p <- function(srobj) {
    NormalizeData(object = srobj, verbose = FALSE)
    }

#' @noRd
#' @importFrom Seurat SCTransform
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
# nDRs <- ChoosePCs_default(srobj, figure_dir=fig_dir)
#' 
#' @param srobj v4 seurat object
#' @param improved_diff_quantile percent variance explained of next PC to choose PC
#' @param significance JackStraw significance required to choose PC
#' @param reduction name of the reduction slot where to select the reduced
#'   dimensions for downstream analysis
#' return vectors of PCs chosen
#' @importFrom Seurat RunPCA JackStraw ScoreJackStraw
#' @importFrom SeuratObject DefaultAssay Reductions
#' @noRd
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
        nDRs <- 1:max(which(
            diff(varExplained) <
                quantile(diff(varExplained), 1 - improved_diff_quantile)) + 1)
    } else {
        # Change default assay to run PCA on RNA
        DefaultAssay(srobj) <- "RNA" # JackStraw cannot be run on SCTransform-normalized data
        #Run PCA under RNA for Jackstraw
        npcs <- min(50, round(ncol(srobj)/2))
        srobj <- RunPCA(
            object = srobj,
            npcs = npcs,
            assay = "RNA",
            verbose = FALSE)
        suppressWarnings({srobj <- JackStraw(srobj, dims = npcs, assay = "RNA")})
        srobj <- ScoreJackStraw(
            srobj,
            dims = 1:ncol(srobj[["pca"]]@cell.embeddings))
        # Extract number of PCs
        nDRs <- which(JS(srobj$pca)@overall.p.values[, "Score"] < significance)

        # DefaultAssay(srobj) <- 'SCT' ## TODO related to the comment 7 lines above
    }
    
    # Set vector for corner case 2
    if (length(nDRs) <= 2) {
        nDRs <- seq_len(2)
    }
    return(nDRs)
}


####################################
## ARBOL supported preprocessing  ##
####################################
# Helper function to process the seurat object at each step
#' @importFrom Seurat FindVariableFeatures ScaleData RunPCA FindNeighbors
#' @importFrom harmony RunHarmony
#' @noRd
.preprocess <- function(
        srobj,
        .normalize_data = .normalize_log1p,
        harmony_var = NULL) {
    
    # Some input quality checks
    stopifnot(
        is(srobj, "Seurat"),
        is(.normalize_data, "function"),
        is.character(harmony_var) | is.null(harmony_var))
    
    # Standard data preprocessing
    srobj <- .normalize_data(srobj = srobj)
    srobj <- FindVariableFeatures(
        srobj,
        nfeatures = 2000,
        selection.method = "disp",
        verbose = FALSE)
    srobj <- ScaleData(
        object = srobj,
        features = VariableFeatures(srobj),
        verbose = FALSE)
    srobj <- RunPCA(
        object = srobj,
        npcs = min(50, round(ncol(srobj) / 2)),
        verbose = FALSE)
    # Specify this for the .choose_dims_default step
    reduction <- "pca"
    
    # Run Harmony integration if specified
    if (!is.null(harmony_var)) {
        suppressWarnings(
            srobj <- RunHarmony(
                object = srobj,
                group.by.vars = harmony_var,
                reduction = "pca")
        )
        
        # Specify this for the .choose_dims_default step
        reduction <- "harmony"
    }
    
    # Choose number of reduced dimensions
    nDRs <- .choose_dims_default(srobj, reduction = reduction)
    
    # Save the PCs to use in misc
    srobj@misc$nDRs <- nDRs
    message(
        paste0("chose ", length(nDRs), "nDRs, added choices to srobj@misc$nDRs"))
    
    # Compute KNN graph
    srobj <- FindNeighbors(
        object = srobj,
        dims = nDRs,
        k.param = ceiling(0.5*sqrt(ncol(srobj))),
        reduction = reduction,
    )
    
    return(srobj)
}
