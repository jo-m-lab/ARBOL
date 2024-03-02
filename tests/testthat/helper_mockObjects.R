#! /usr/bin/Rscript

CreateMockSeuratObject <- function() {
  
  suppressMessages({
    # Simulate some data
    data <- Matrix(sample(0:500, 2000, replace = TRUE), nrow = 100, ncol = 20, sparse=TRUE)
    rownames(data) <- paste0("Gene", 1:100)
    colnames(data) <- paste0("Cell", 1:20)

    # Create a basic Seurat object
    mock_seurat <- CreateSeuratObject(counts = data)

    # Add necessary metadata
    mock_seurat <- AddMetaData(mock_seurat, metadata = data.frame(
        CellID = colnames(data),
        tierNident = sample(c("Type1", "Type2", "Type3", "Type4"),
                            20, replace = TRUE),
        condition = sample(c("Control", "Treatment"), 20, replace = TRUE),
        batch = sample(c("Batch1", "Batch2"), 20, replace = TRUE),
        sample = sample(c("Sample1", "Sample2"), 20, replace = TRUE),
        row.names = colnames(data)
    ))

    mock_seurat <- NormalizeData(mock_seurat)
    mock_seurat <- FindVariableFeatures(mock_seurat)
    mock_seurat <- ScaleData(mock_seurat)
    mock_seurat <- RunPCA(mock_seurat, npcs = 5)
    mock_seurat <- SCTransform(mock_seurat, method='glmGamPoi')
    DefaultAssay(mock_seurat) <- "RNA"


    # Placeholder for @misc slot modifications, assuming these are required for your tests
    # (Adjust as necessary based on your package's specific requirements)
    mock_seurat@misc$pvclust <- list() # Empty list, to be filled by your function as needed
  })
  # Return the mock Seurat object
  return(mock_seurat)
}

