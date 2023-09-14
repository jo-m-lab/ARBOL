set.seed(321)
# mock up some single-cell, mixture & marker data
sce <- mockSC(ng = 500, nc = 50, nt = 3)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----  Check normalization methods  -------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# .normalize_log1p with Seurat object ----
test_that(".normalize_log1p", {
    # Normalize Seurat object
    sce <- .normalize_log1p(sce)
    
    # Run checks Seurat object
    expect_is(sce[["RNA"]]@data, "dgCMatrix")
    expect_false(identical(sce[["RNA"]]@counts, sce[["RNA"]]@data))
    expect_identical(dim(sce[["RNA"]]@counts), dim(sce[["RNA"]]@data))
})

# .normalize_sct with Seurat object ----
test_that(".normalize_sct", {
    # Normalize Seurat object
    sce <- .normalize_sct(sce)
    
    # Run checks Seurat object
    expect_is(sce[["SCT"]]@data, "dgCMatrix")
    expect_false(identical(sce[["RNA"]]@counts, sce[["SCT"]]@data))
    expect_identical(dim(sce[["RNA"]]@counts), dim(sce[["SCT"]]@data))
})

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----  Check nDR selection  ---------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# .choose_dims_default PCA with Seurat object ----
test_that(".choose_dims_default PCA", {
    sce <- NormalizeData(object = sce, verbose = FALSE)
    sce <- FindVariableFeatures(sce, nfeatures = 2000, verbose = FALSE)
    sce <- ScaleData(object = sce, features = VariableFeatures(sce), verbose = FALSE)
    sce <- RunPCA(
        object = sce,
        npcs = 10,
        verbose = FALSE)
    
    # Normalize Seurat object
    nDRs <- .choose_dims_default(sce, reduction = "pca")
    
    # Run checks object
    expect_is(nDRs, "integer")
    expect_true(length(nDRs) > 1)
    expect_true(length(nDRs) <= ncol(sce[["pca"]]))
})

# .choose_dims_default not PCA with Seurat object ----
test_that(".choose_dims_default not PCA", {
    sce <- NormalizeData(object = sce, verbose = FALSE)
    sce <- FindVariableFeatures(sce, nfeatures = 2000, verbose = FALSE)
    sce <- ScaleData(object = sce, features = VariableFeatures(sce), verbose = FALSE)
    sce <- RunPCA(
        object = sce,
        npcs = 10,
        verbose = FALSE)
    sce[["harmony"]] <- sce[["pca"]]
    
    # Normalize Seurat object
    nDRs <- .choose_dims_default(sce, reduction = "harmony")
    
    # Run checks object
    expect_is(nDRs, "integer")
    expect_true(length(nDRs) > 1)
    expect_true(length(nDRs) <= ncol(sce[["pca"]]))
})

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----  Check .preprocess  -----------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# check .preprocessing() ----
test_that(".preprocessing pca", {
    sce <- NormalizeData(object = sce, verbose = FALSE)
    sce <- FindVariableFeatures(sce, nfeatures = 2000, verbose = FALSE)
    sce <- ScaleData(object = sce, features = VariableFeatures(sce), verbose = FALSE)
    sce <- RunPCA(
        object = sce,
        npcs = 10,
        verbose = FALSE)
    
    # Normalize Seurat object
    sce <- .preprocess(sce, .normalize_data = .normalize_log1p, harmony_var = NULL)
    
    # Run checks object
    nDRs <- sce@misc$nDRs
    expect_is(nDRs, "integer")
    expect_true(length(nDRs) > 1)
    expect_true(length(nDRs) <= ncol(sce[["pca"]]))
    expect_true(length(sce@graphs) > 0)
    expect_true(all(names(sce@graphs) == c("RNA_nn", "RNA_snn")))
    expect_true(identical(dim(sce@graphs$RNA_snn), c(ncol(sce), ncol(sce))))
    
})

# check .preprocessing() with harmony ----
test_that(".preprocessing harmony", {
    sce <- NormalizeData(object = sce, verbose = FALSE)
    sce <- FindVariableFeatures(sce, nfeatures = 2000, verbose = FALSE)
    sce <- ScaleData(object = sce, features = VariableFeatures(sce), verbose = FALSE)
    sce <- RunPCA(
        object = sce,
        npcs = 10,
        verbose = FALSE)
    
    # Normalize Seurat object
    sce <- .preprocess(sce, .normalize_data = .normalize_log1p, harmony_var = "type")
    
    # Run checks object
    nDRs <- sce@misc$nDRs
    expect_is(nDRs, "integer")
    expect_true(length(nDRs) > 1)
    expect_true(length(nDRs) <= ncol(sce[["pca"]]))
    expect_true(length(sce@graphs) > 0)
    expect_true(all(names(sce@graphs) == c("RNA_nn", "RNA_snn")))
    expect_true(identical(dim(sce@graphs$RNA_snn), c(ncol(sce), ncol(sce))))
    
})
