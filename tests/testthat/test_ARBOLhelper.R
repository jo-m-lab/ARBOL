source("helper_mockObjects.R")
library(testthat)
# test_local here to make all ARBOL functions available to test environment.
# testthat::test_local()

test_that("ARBOLcentroidTaxonomy workflow tests\n
           centroidTaxonomy adds a pvclust obj to srobj@misc$pvclust\n
           binaryTreeToDF DF contains 2(nSpecies)-1 nodes\n
           treeAllotment adds correct metadata to binarydf", {
  mock_seurat <- CreateMockSeuratObject() # Create your mock Seurat object
  # Test centroidTaxonomy
  mock_seurat <- centroidTaxonomy(mock_seurat)
  expect_true(inherits(mock_seurat@misc$pvclust, "pvclust"))
  binarydf <- binaryTreeToDF(pvclust_tree=mock_seurat@misc$pvclust)
  expect_true(length(unique(binarydf$pathString)) ==
                2 * length(unique(mock_seurat$tierNident)) - 1)
  dataTree <- treeAllotment(srobj = mock_seurat,
              treedf = binarydf,
              categories = c("sample", "batch", "condition"),
              diversities = c("sample", "batch", "condition"),
              diversity_metric = "simpson",
              counts = c("sample", "batch", "condition"),
              totals = "nCount_RNA"
            )
  expect_true(inherits(dataTree, "Node")) # dataTree is a data.tree object
  expect_true(all(c("ids", # make sure attributes are added to dataTree
                "tierNident",
                "plotHeight",
                "batch_diversity",
                "condition_n_Treatment",
                "nCount_RNA"
              ) %in%
                dataTree$attributesAll))
  dataTree <- propagateTree(dataTree,
                            srobj = mock_seurat,
                            categorical_attributes = c("sample", "batch", "condition"),
                            diversity_attributes = c("sample", "batch", "condition"),
                            numerical_attributes = c("sample", "batch", "condition"),
                            total_attributes = "nCount_RNA")

  # make sure number of cells in object matches n cells in tree
  expect_true(length(unique(unlist(dataTree$Get("ids")))) ==
                length(Cells(mock_seurat)))
  # make sure all cell names are preserved correctly
  expect_equal(sort(unique(unlist(dataTree$Get("ids")))),
               sort(Cells(mock_seurat)))
  # make sure diversities are in end nodes
  # and make sure they"re propagated correctly
  end_node_diversities <- vegan::diversity(mock_seurat@meta.data %>%
                             dplyr::count(sample, tierNident) %>%
                             tidyr::pivot_wider(names_from = sample,
                                         values_from = n, values_fill = 0) %>%
                             column_to_rownames("tierNident"),
                             "simpson") %>%
                             data.frame %>% t %>% data.frame
  expect_equal(FindNode(dataTree,"Type4")$sample_diversity,
               end_node_diversities$Type4)
  expect_equal(dataTree$root$sample_diversity,
               diversity(mock_seurat@meta.data %>%
                           dplyr::count(sample) %>%
                           tidyr::pivot_wider(names_from = sample, values_from = n), 
                           "simpson"))
  # make sure summed totals work as expected
  expect_true(sum(mock_seurat[["RNA"]]$counts) ==
                dataTree$Get("nCount_RNA")[1])
  tidyGraph <- data.tree_to_ggraph(dataTree,
                                   categories = c("batch", "condition"),
                                   diversities = c("batch", "condition"),
                                   counts = c("batch", "condition"),
                                   totals = "nCount_RNA",
                                   "numChildren")
  expect_true(inherits(tidyGraph, "tbl_graph"))
  tidyGraph <- tidyGraph %>% activate(nodes) %>%
               mutate(tier = str_count(label, "\\."))
  tidyGraph <- tidyGraph %>% activate(nodes) %>%
               mutate(string = label,
                      name = basename(label) %>%
               str_replace_all("T0C0.", ""))
  rm(dataTree)
})

test_that("sv.cbind combines 3 int sparse vectors correctly", {
  # Create mock sparse vectors
  v1 <- Matrix::sparseVector(x = c(1,2), i = c(1,4), length = 7)
  v2 <- Matrix::sparseVector(x = c(3,4), i = c(2,5), length = 7)
  v3 <- Matrix::sparseVector(x = c(5,6), i = c(3,6), length = 7)
  # Expected combined matrix
  expected_mat <- Matrix::sparseMatrix(
    i = c(1,4,2,5,3,6),
    j = c(1,1,2,2,3,3),
    x = c(1,2,3,4,5,6),
    dims = c(7, 3)
  )
  # Run sv.cbind
  result <- sv.cbind(list(v1, v2, v3))
  # Test if result matches expected
  expect_equal(result, expected_mat)
})

test_that("diversityPerGroup returns diversities as expected", {
  mock_seurat <- CreateMockSeuratObject()
  end_node_diversities <- vegan::diversity(mock_seurat@meta.data %>%
                             dplyr::count(sample, tierNident) %>%
                             tidyr::pivot_wider(names_from = sample,
                                         values_from = n, values_fill = 0) %>%
                             column_to_rownames("tierNident"),
                             "simpson") %>%
                             data.frame %>% rownames_to_column('tierNident')
  colnames(end_node_diversities) <- c('tierNident', 'sample_diversity')
  # DiversityPerGroup tests unwrapped
  df <- mock_seurat@meta.data
  species <- "tierNident"
  group <- "sample"
  species_sym <- rlang::sym(species)
  group_sym <- rlang::sym(group)
  
  # Count groups per species directly using curly-curly
  tierNcount <- df %>%
    group_by({{species_sym}}) %>%
    count({{group_sym}}, name = "n") %>% ungroup
  
  # Pivot table to allow vegan::diversity call
  tierNwide <- tierNcount %>%
    pivot_wider(names_from = {{group_sym}}, values_from = n, values_fill = list(n = 0))

  # Use rownames from species for the diversity function, which requires a matrix or dataframe with rownames
  tierNwide_df <- as.data.frame(tierNwide)
  # Ensure species column is the first column for row names
  tierNwide_df <- tierNwide_df %>% select({{species}}, everything())
  rownames(tierNwide_df) <- tierNwide_df[, 1]
  tierNwide_df <- tierNwide_df[, -1]

  # Calculate diversity
  diversity_values <- vegan::diversity(tierNwide_df, index = diversity_metric)
  
  # Prepare the result as a dataframe
  result <- data.frame(
    species = rownames(tierNwide_df),
    diversity = diversity_values,
    row.names = NULL
  )
  # Rename the diversity column
  names(result)[1] <- species
  names(result)[2] <- sprintf('%s_diversity', group)
  wrapped <- diversityPerGroup(mock_seurat@meta.data, 'tierNident', 'sample')

  # manual calculation vs. unwrapped run check
  expect_equal(end_node_diversities, result)
  # unwrapped run vs. wrapped check
  expect_equal(result, wrapped)
})

test_that("sv.cbind combines two random vectors correctly", {
  set.seed(687) # Ensure reproducibility
  # Generate controlled random sparse vectors
  v1 <- Matrix::sparseVector(x = runif(5), i = sample(1:10, 5), length = 10)
  v2 <- Matrix::sparseVector(x = runif(5), i = sample(1:10, 5), length = 10)
  # Manually construct the expected result
  # Here, we simply combine the two vectors knowing their structure
  expected_matrix <- Matrix::sparseMatrix(
    i = c(v1@i, v2@i),
    j = c(rep(1, length(v1@x)), rep(2, length(v2@x))),
    x = c(v1@x, v2@x),
    dims = c(10, 2)
  )
  # Combine vectors using sv.cbind
  result_matrix <- sv.cbind(list(v1, v2))
  # Compare result with expected
  expect_equal(result_matrix, expected_matrix)
})

test_that("dgcMatrix.aggregate aggregates correctly", {
  # Step 1: Create a reproducible sparse matrix
  set.seed(687)  # Ensure reproducibility
  mat <- Matrix(rnorm(20), 5, 4, sparse = TRUE)
  # Step 2: Define groupings
  groupings <- c(1, 2, 1, 2, 1)
  # Step 3: Select an aggregation function
  aggregation_function <- mean # Aggregating by mean, as it"s default
  # First we make sure aggregate is producing a matrix of the correct size.
  expected_rows <- length(unique(groupings))
  expected_cols <- ncol(mat)
  # Step 4: Use the function and compare output
  aggregated_mat <- dgcMatrix.aggregate(mat, groupings, aggregation_function)
  # Check if dimensions match expected results
  expect_equal(dim(aggregated_mat), c(expected_rows, expected_cols))
})

test_that("dgcMatrix.aggregate manually calculates correctly", {
  # Define a non-random sparse matrix
  values <- c(1, 2, 3, 4, 5, 6, 7, 8)
  i <- c(1, 2, 3, 1, 2, 3, 4, 4)  # Row indices
  j <- c(1, 1, 1, 2, 2, 2, 1, 2)  # Column indices
  mat <- sparseMatrix(i = i, j = j, x = values, dims = c(4, 2))
  # Define groupings (align with row indices)
  groupings <- c(1, 1, 2, 2)  # Two groups
  # Use mean for agg fun as it"s default
  aggregation_function <- mean  # Aggregating by mean
  # Manually calculate the expected result
  expected_values <- c(mean(c(1, 2)), mean(c(4, 5)),
                       mean(c(3, 7)), mean(c(6, 8)))
  expected_i <- c(1, 1, 2, 2) 
  expected_j <- c(1, 2, 1, 2)
  expected_mat <- sparseMatrix(i = expected_i,
                               j = expected_j,
                               x = expected_values,
                               dims = c(2, 2))
  # Use the function and compare output
  aggregated_mat <- dgcMatrix.aggregate(mat, groupings, aggregation_function)
  # Check if the resulting matrix matches the expected matrix
  expect_equal(as(aggregated_mat, "CsparseMatrix"),
               as(expected_mat, "CsparseMatrix"))
})

test_that("getCentroids returns correct output", {
  mock_seurat <- CreateMockSeuratObject() # Create your mock Seurat object
  # Test getCentroids
  centroids <- getCentroids(srobj = mock_seurat,
                            tree_reduction = "pca",
                            reduction_dims = 1:5,
                            centroid_method = "mean",
                            centroid_assay = "RNA",
                            gene_list = Features(mock_seurat))                            
  # Check if centroids is a data frame and has correct dimensions
  expect_true(is.data.frame(centroids))
  expect_equal(dim(centroids)[1],
               dim(Embeddings(mock_seurat,
                                 reduction = "pca"))[2])
  expect_equal(dim(centroids)[2],
               length(unique(mock_seurat@meta.data$tierNident)))

  centroids <- getCentroids(srobj = mock_seurat,
                            tree_reduction = "centroids",
                            centroid_method = "median",
                            centroid_assay = "RNA",
                            centroid_layer = "scale.data",
                            gene_list = Features(mock_seurat))
  expect_true(is.data.frame(centroids))
})

test_that("standardNames adds one name per species to srobj", {
  # create mock seurat object with extra cells to make DE work
  set.seed(687)
  mock_seurat <- CreateMockSeuratObject(num_cells = 200,
                                        num_cell_types = 6,
                                        num_genes = 500)
  # Test getCentroids
  mock_seurat <- getStandardNames(mock_seurat,
                                  figdir = NULL,
                                  celltype_col = "celltype",
                                  standardname_col = "standardName",
                                  n_genes = 2,
                                  logfc.threshold = 0,
                                  min.pct = 0)
  standardNames <- unique(mock_seurat@meta.data$standardName)
  pattern <- "T[0-9]+\\.Gene[0-9]+"
  expect_true(all(sapply(standardNames, function(x) grepl(pattern, x))))
})
