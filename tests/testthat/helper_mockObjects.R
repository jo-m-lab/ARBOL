#! /usr/bin/Rscript

CreateMockSeuratObject <- function(num_cells = 20,
                                   num_cell_types = 4,
                                   num_genes = 100) {
  suppressMessages({
    library(Seurat)
    library(dplyr)
    library(Matrix)
    # Simulate some counts matrix data
    data <- Matrix(sample(0:500, num_genes * num_cells, replace = TRUE),
                   nrow = num_genes,
                   ncol = num_cells,
                   sparse = TRUE)
    rownames(data) <- paste0("Gene", 1:num_genes)
    colnames(data) <- paste0("Cell", 1:num_cells)

    # Create a basic Seurat object
    mock_seurat <- CreateSeuratObject(counts = data)

    # Simulate basic metadata with 3 ARBOL columns: CellID, tierNident, sample
    metadata <- data.frame(
      CellID = colnames(data),
      tierNident = sample(paste0("Type", 1:num_cell_types),
                          num_cells,
                          replace = TRUE),
      condition = sample(c("Control", "Treatment"),
                         num_cells,
                         replace = TRUE),
      batch = sample(c("Batch1", "Batch2"),
                     num_cells,
                     replace = TRUE),
      sample = sample(c("Sample1", "Sample2"),
                      num_cells,
                      replace = TRUE),
      row.names = colnames(data)
    )

    # Add metadata to seurat object
    mock_seurat <- AddMetaData(mock_seurat, metadata = metadata)

    # Simulate 2 cell type lineages (e.g. immune or stroma) for reference frames / tree plots.
    tierNidents <- unique(metadata$tierNident)
    half_index <- ceiling(length(tierNidents) / 2)
    group1 <- tierNidents[1:half_index]
    group2 <- tierNidents[(half_index + 1):length(tierNidents)] # Second half
    # Assign lineages to seurat object
    ctdf <- data.frame(tierNident = tierNidents)
    ctdf$celltype <- ifelse(ctdf$tierNident %in% group1, "T1", "T2")
    ctdf <- metadata %>% left_join(ctdf, by = "tierNident")
    row.names(ctdf) <- ctdf$CellID
    mock_seurat <- AddMetaData(mock_seurat, ctdf)

    mock_seurat <- NormalizeData(mock_seurat)
    mock_seurat <- FindVariableFeatures(mock_seurat)
    mock_seurat <- ScaleData(mock_seurat)
    mock_seurat <- RunPCA(mock_seurat, npcs = 5)
    mock_seurat <- SCTransform(mock_seurat, method = "glmGamPoi")
    DefaultAssay(mock_seurat) <- "RNA"

    # Placeholder for @misc slot modifications if needed for tests
    mock_seurat@misc$pvclust <- list()
  })
  # Return the mock Seurat object
  return(mock_seurat)
}
