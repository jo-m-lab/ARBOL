#! /usr/bin/Rscript
require(Seurat)
require(tidyverse)
require(devtools)
require(data.table)
require(glmGamPoi)
require(future)

#' Performs iterative clustering (ARBOL) on a v4 Seurat object
#'
#' @description Iteratively clusters seurat object from single cell datasets, choosing optimum resolution parameters at each stage of clustering. Outputs QC plots for each tier and stage.
#'
#' can return output directories figs/ srobjs/ (see examples) containing subset seurat objects and QC plots matching tiered structure of found in dataset
#'
#' can also return endclusts/ folder which contains tsv files named with the path to the endcluster (T0C0T1C1.tsv) and contains the cellnames (one per line) that comprise that cluster
#'
#' @examples srobj <- readRDS("/path/to/seurat_object.rds")
#' tiers <- GenTieredclusters(srobj,
#'                            saveSROBJdir = "~/tieredoutput/srobjs",
#'                            figdir = "~/tieredoutput/figdir",
#'                            SaveEndNamesDir = "~/tieredoutput/endclusts")
#'
#' @param srobj v4 seurat object
#' @param cluster_assay assay to use for clustering defaults to "SCT"
#' @param cells cellnames if tiered clustering should start on subset of object
#' @param tier starting level defaults to 0
#' @param clustN cluster starting from default to 0
#' @param PreProcess_fun function to use for preproccessing defaults to PreProcess_sctransform
#' @param min_cluster_size minimum number of cells to allow further clustering. defaults to 100
#' @param max_tiers maximum number of tiers to allow further clustering. defaults to 10
#' @param EnoughDiffUp minimum number of up-regulated genes to call clusters unique. Differential expression is performed when clustering finds 2 clusters. defaults to 5
#' @param EnoughDiffDown minimum number of down-regulated genes. If either up or down is not met, the 2 clusters are joined, and further clustering is stopped. defaults to 5
#' @param tierAllowedRecomb minimum tier where differential expression can be called to decide on recombination. defaults to 0. clustering may stop early when clustering finds 2 clusters with high cell numbers, as Wilcoxon effect sizes may be low in their DE.
#' @param ChooseOptimalClustering_fun function that returns srobj with clusters in `srobj$Best.Clusters` after choosing optimal clustering resolution 
#' @param saveSROBJdir where to save seurat objects for each tier and cluster, if null does not save
#' @param figdir where to save QC figures for each tier and cluster, if null does not save
#' @param saveEndNamesDir where to save directory of end clusters, if null does not save
#' @param SaveEndFileName prefix for all end cluster files
#' @return list of lists with all seurat objects (highly recommend using folder arguments for saving outputs)
#' @export

GenTieredClusters <- function(srobj, cluster_assay = "SCT", cells = NULL, tier=0, clustN = 0,
                              PreProcess_fun = PreProcess_sctransform,
                              ChooseOptimalClustering_fun = ChooseOptimalClustering_default,
                              saveSROBJdir=NULL, figdir=NULL, SaveEndNamesDir=NULL, SaveEndFileName=NULL,
                              min_cluster_size = 100, max_tiers = 10, EnoughDiffUp = 5, EnoughDiffDown = 5,
                              tierAllowedRecomb=0) {
  ######################################################################################################
  #' make sure output directories exist
  ######################################################################################################
  if (!is.null(saveSROBJdir)) dir.create(saveSROBJdir, showWarnings = T, recursive = T)
  if (!is.null(figdir)) dir.create(figdir, showWarnings = T, recursive = T)
  if (!is.null(SaveEndNamesDir) & tier==0) dir.create(SaveEndNamesDir, showWarnings = T, recursive = T)
  SaveEndFileName <- paste0(SaveEndFileName, sprintf("_T%sC%s", tier, clustN))
  
  ######################################################################################################
  #' make sure QC metadata exists
  ######################################################################################################
  if (tier == 0) {srobj$nFeature_RNA <- Matrix::colSums(srobj@assays$RNA@counts > 0)}
  if (tier == 0) {srobj$nCount_RNA <- Matrix::colSums(srobj@assays$RNA@counts)}
  if (tier == 0) {srobj$percent.mt <- PercentageFeatureSet(srobj, pattern = "^MT-")}

  ######################################################################################################
  #' basic processing
  ###################################################################################################### 
  #' subset to provided cells
  if (is.null(cells)) {
    cells <- colnames(srobj)
  }
  working_srobj <- subset(srobj, cells=cells)
  working_srobj@misc$tier <- tier

  ######################################################################################################

  #' keep track of position in tree for logging
  message("Number of cells: ",paste0(ncol(working_srobj),collapse='\n'))
  message(paste0("Starting tier: ", tier, ", cluster: ", clustN, ", with ", ncol(working_srobj), " cells" ))
  message(paste0("file name: ", sprintf("%s/%s.tsv", SaveEndNamesDir, SaveEndFileName)))

  #' Basic QC plotting of end-nodes before interruption of recursion
  tryCatch({QC_Plotting(working_srobj, fig_dir = figdir)
          },
          error = function(e) {message('QC Plotting failure'); print(paste("error: ",e))
          })
  
  ######################################################################################################
  #' check if too few cells or past tier 10. if so, return end-node without processing
  ######################################################################################################

  if (min_cluster_size < 5) {
    message('WARNING: setting minimum cluster size to < 5 will throw errors in processing functions. \n 
             the full run should finish, but all processing of srobj < 5 cells is unreliable. \n
             Search "failure" in logs for information')
  }

  if ( (ncol(working_srobj) < min_cluster_size) | (working_srobj@misc$tier > max_tiers) ) {
    message(cbind("found end-node below min number of cells or above max tier. num cells: ", ncol(working_srobj),' < ', min_cluster_size))
    #write end-node table and srobj
    EndNode_Write(working_srobj, srobj_dir = saveSROBJdir, endclust_dir = SaveEndNamesDir, filename = SaveEndFileName)
    #stop recursion
    return(working_srobj)
  }
  
  #' Preprocessing: Run SCtransform or log transform, PCA, choose num PCs for downstream analysis, find neighbors.

  tryCatch({working_srobj <- PreProcess_fun(working_srobj, fig_dir=figdir)
          },
          error = function(e) {message('Pre-processing failure'); print(paste("Pre-processing error: ", e))
        })
  
  #' Get optimum cluster resolution and cluster

  tryCatch({res <- ChooseOptimalClustering_fun(working_srobj,
                                     assay=cluster_assay,
                                     PreProcess_fun = PreProcess_fun,
                                     downsample_num = 7500,
                                     figdir=figdir)
            },
            error = function(e) {message('Resolution choice failure'); print(paste("Resolution choice error: ", e))
          })
  tryCatch({working_srobj <- FindClusters(working_srobj,
                                assay = cluster_assay,
                                resolution = res,
                                k.param=ceiling(0.5*sqrt(ncol(working_srobj))))
            },
            error = function(e) {message('FindClusters failure'); print(paste("Clustering error: ", e))
          })
  
  ######################################################################################################
  #' logging & plotting: pre-processing and cluster results
  ###################################################################################################### 

  tryCatch({working_srobj <- Plotting(working_srobj, fig_dir = figdir)
  },
  error = function(e) {message('Plotting failure'); print(paste("Plotting error: ", e))
  })

  ######################################################################################################
  #' check if too few clusters to call end-node and end recursion
  ######################################################################################################

  #' if one cluster:
  if ( !(length(unique(Idents(working_srobj))) > 1) ) {
    #write end-node table and srobj
    EndNode_Write(working_srobj, srobj_dir = saveSROBJdir, endclust_dir = SaveEndNamesDir, filename = SaveEndFileName)
    
    message("found end-node with one cluster")
    return(working_srobj)
  }

 #' if tier > 2 and
  #' if two indistinguishable clusters, defined by # genes differential expression with FDR < 0.05

  if ( (working_srobj@misc$tier > tierAllowedRecomb) ) {

    if ( (length(unique(Idents(working_srobj))) == 2) ) {
      if ( !(EnoughDiffExp(working_srobj, levels(Idents(working_srobj))[1], levels(Idents(working_srobj))[2], nUpregulated = EnoughDiffUp, nDownregulated = EnoughDiffDown)) ) {
        message("found end-node with two indistinguishable clusters")
        
        Idents(working_srobj) <- 0
        working_srobj$two_not_quite_clusters <- working_srobj$seurat_clusters
        working_srobj$seurat_clusters <- 0

        #write end-node table and srobj
        EndNode_Write(working_srobj, srobj_dir = saveSROBJdir, endclust_dir = SaveEndNamesDir, filename = SaveEndFileName)

        return(working_srobj)
      }
    }
  }

  #Write srobj of internal node
  Node_Write(working_srobj, srobj_dir = saveSROBJdir)

  ######################################################################################################
  #' recurse
  ######################################################################################################

  #' split all clusters into separate srobjs
  message("continuing to recurse")

  subsets <- SplitSrobjOnIdents(working_srobj, paste0("tier", tier))

  print(subsets)
  #' recurse along subsets
  return(lapply(seq_along(subsets),
                function(i) {
                  saveSROBJdir_nextTier=NULL
                  if (!is.null(saveSROBJdir)) {
                    saveSROBJdir_nextTier <- sprintf("%s/%s", saveSROBJdir, paste0("T", tier+1, "C", i-1))
                  }
                  figdir_nextTier=NULL
                  if (!is.null(figdir)) {
                    figdir_nextTier <- sprintf("%s/%s", figdir, paste0("T", tier+1, "C", i-1))
                  }
                  min_preserve = min_cluster_size
                  max_preserve = max_tiers
                  EnoughDiffUp_preserve = EnoughDiffUp
                  EnoughDiffDown_preserve = EnoughDiffDown
                  tierAllowedRecomb_preserve = tierAllowedRecomb
                  GenTieredClusters(srobj,
                                    cells=colnames(subsets[[i]]),
                                    tier=tier+1, clustN = i-1,
                                    saveSROBJdir = saveSROBJdir_nextTier,
                                    figdir = figdir_nextTier,
                                    SaveEndNamesDir = SaveEndNamesDir,
                                    SaveEndFileName = SaveEndFileName,
                                    min_cluster_size=min_preserve,
                                    max_tiers= max_preserve,
                                    EnoughDiffUp = EnoughDiffUp_preserve,
                                    EnoughDiffDown = EnoughDiffDown_preserve,
                                    tierAllowedRecomb = tierAllowedRecomb_preserve
                  )
                }))
}



PreProcess <- function(srobj, ChoosePCs_fun = ChoosePCs_default, ChooseHVG_fun = ChooseHVG_default, fig_dir=NULL, ...) {
  #' seurat v3 processing
  srobj <- NormalizeData(object = srobj)
  srobj <- ChooseHVG_fun(srobj)
  srobj <- ScaleData(object = srobj, features=VariableFeatures(srobj))
  srobj <- RunPCA(object = srobj, npcs = min(50, round(ncol(srobj)/2)), verbose=F)
  nPCs  <- ChoosePCs_fun(srobj, figure_dir=fig_dir)
  srobj@misc$nPCs <- nPCs
  message(paste0("chose ", length(nPCs), "PCs, added choices to srobj@misc$nPCs "))
  srobj <- FindNeighbors(object = srobj, dims=nPCs)
  return(srobj)
}

PreProcess_sctransform <- function(srobj, ChoosePCs_fun = ChoosePCs_default, fig_dir=NULL) {
  #' seurat v4 processing
  
  tryCatch({
            srobj <- SCTransform(object = srobj, verbose=TRUE, method='glmGamPoi')
           }, error=function(e) 
           {message(sprintf('SCTransform failed to run, likely due to too few cells. cell num: %s .... Defaulting back to log1p normalization',ncol(srobj)))
           })
  
  #If it fails, normalization defaults back to log1p
  srobj <- NormalizeData(srobj, assay = "RNA")
  #Add ScaleData slot for regular normalization, as it may be called in downstream functions
  srobj <- ScaleData(srobj, assay = "RNA")
  #Additionally FindVariableFeatures to allow PCA
  srobj <- FindVariableFeatures(srobj, assay="RNA")
  tryCatch({
    message(sprintf('dimensions of SCT matrix: %s', paste0(capture.output(dim(srobj[["SCT"]]@scale.data)),collapse='\n')
      ))
  }, error = function(e)
  {message('skipping SCT metrics..')
  })

  srobj <- RunPCA(object = srobj, npcs = min(50, round(ncol(srobj)/2)), verbose=F)

  nPCs  <- ChoosePCs_fun(srobj, figure_dir=fig_dir)
  srobj@misc$nPCs <- nPCs
  message(paste0("chose ", length(nPCs), "PCs, added choices to srobj@misc$nPCs "))
  srobj <- FindNeighbors(object = srobj, dims=nPCs)
  return(srobj)
}

#' Chooses the number of principle components to include in clustering
#'
#' @description Seurat's JackStraw is used when there are less than 500 cells in a clustering event. 
#' otherwise, PCs are used in clustering if they explain 15% more variance than the following PC.
#'
#' Also saves plots describing PC choice
#'
#' @examples nPCs <- ChoosePCs_default(srobj, figure_dir=fig_dir)
#' 
#' @param srobj v4 seurat object
#' @param improved_diff_quantile percent variance explained of next PC to choose PC
#' @param significance JackStraw significance required to choose PC
#' @param figure_dir directory where PC choice plots are saved
#' @return list of PCs chosen
#' @export

ChoosePCs_default <- function(srobj, improved_diff_quantile=.85, significance=0.01, figure_dir=NULL) {
  #' if num cells > 500 use elbow plot if small use jackstraw, 
  #'   if none sig use first 2 pcs

  if (ncol(srobj) > 500) {
    eigValues <- (srobj@reductions$pca@stdev)^2  ## EigenValues
    varExplained <- eigValues / sum(eigValues)
    nPCs <- 1:max(which(diff(varExplained) < quantile(diff(varExplained), 1-improved_diff_quantile)) + 1)
    if (!is.null(figure_dir)) {
      ElbowPlot(srobj, ndims=ncol(srobj@reductions$pca@cell.embeddings)) + 
        geom_vline(xintercept  = length(nPCs), color="red")
      ggsave(sprintf("%s/ChosenPCs.pdf", figure_dir))
    }
  } 
  else {
    DefaultAssay(srobj) <- 'RNA'
    #Run PCA under RNA for Jackstraw

    srobj <- RunPCA(object = srobj, npcs = min(50, round(ncol(srobj)/2)), verbose=F)
    suppressWarnings({srobj <- JackStraw(srobj, dims=50)})
    srobj <- ScoreJackStraw(srobj, dims=1:ncol(srobj@reductions$pca@cell.embeddings))
    nPCs <- which(JS(srobj$pca)@overall.p.values[,"Score"] < significance)
    if (!is.null(figure_dir)) {
      JackStrawPlot(srobj, dims = 1:ncol(srobj@reductions$pca@cell.embeddings)) + NoLegend()
      ggsave(sprintf("%s/ChosenPCs.pdf", figure_dir))
    }
    DefaultAssay(srobj) <- 'SCT'
  }
  if (length(nPCs) <= 2) {
    nPCs <- 1:2
  }
  return(nPCs)
}

ChoosePCs_JackStraw <- function(srobj, threshold=0.01, figure_dir=NULL) {
  suppressWarnings({srobj <- JackStraw(srobj, dims=50)})
  srobj <- ScoreJackStraw(srobj, dims=1:ncol(srobj@reductions$pca@cell.embeddings))
  nPCs <- which(JS(srobj$pca)@overall.p.values[,"Score"] < threshold)
  if (!is.null(figure_dir)) {
    JackStrawPlot(srobj, dims = 1:ncol(srobj@reductions$pca@cell.embeddings)) + NoLegend()
    ggsave(sprintf("%s/ChosenPCs.pdf", figure_dir))
  }
  if (length(nPCs) < 2) {
    nPCs <- 1:2
  }
  return(nPCs)
}

ChooseHVG_default <- function(srobj) {
  FindVariableFeatures(srobj, nfeatures=2000, selection.method = "disp")
}

ChooseOptimalClustering_default <- function(srobj, downsample_num = Inf,
                                            PreProcess_fun = PreProcess, ...) {
  #' get resolution from paramsweep of heavily downsampled dataset
  
  # downsample to 10% of data
  if ((downsample_num < Inf) & (downsample_num <  ncol(srobj)) & ncol(srobj) >= 100) {
    cells <- sample(colnames(srobj), downsample_num)
    srobj <- subset(srobj, cells=cells)
    srobj <- PreProcess_fun(srobj)
  }
  #setting resolution choice to arbitrary number in case of error
  srobj@misc$resolution.choice = 1
  tryCatch({srobj <- chooseResolution_SilhouetteAnalysisParameterScan(srobj, ...)
          }, error = function(e) {
            message('Silhouette Analysis failed. resolution set to 1. WARNING: Double check your result.')
            print(paste("Silhouette Analysis error: ", e))
          })
  return(srobj@misc$resolution.choice)
  
}


EnoughDiffExp <- function(srobj, ident1, ident2, nUpregulated = 5, nDownregulated = 5) {
  tryCatch({ 
    markers <- FindMarkers(srobj, ident.1=ident1, ident.2=ident2)
    fccol <- grep("FC",colnames(markers))
    sig.markers <- markers[
    which((abs(markers[,fccol]) >= 1.5) & (markers$p_val_adj < 0.05)), ]
    return((sum(sig.markers[,fccol] > 0) > nUpregulated) & 
       (sum(sig.markers[,fccol] < 0) > nDownregulated))

  }, error = function(e) {
    message("Failed to compare end clusters when exactly two. combining anyway. Error: ", e);
    return(FALSE)    
  })
}


SplitSrobjOnIdents <- function(srobj, proj_nm_suffix) {
  numbers_only <- function(x) !grepl("\\D", x)
  if (any(!numbers_only(levels(srobj)))) {
    idts <- sort(levels(Idents(srobj)))
  } else {
    idts <- sort(as.numeric(levels(Idents(srobj))))
  }
  lapply(idts,
         function(idt) {tmp <- subset(srobj, idents=idt); 
         tmp@project.name=paste0(tmp@project.name, idt, proj_nm_suffix, sep="_");
         return(tmp)})
}

EndNode_Write <- function(working_srobj, srobj_dir = NULL, endclust_dir = NULL, filename = NULL) {
  if (!is.null(srobj_dir)) {
    message(paste("saving srobj to", srobj_dir))
    saveRDS(working_srobj, sprintf("%s/subset_srobj.rds", srobj_dir))
  }
  if (!is.null(endclust_dir)) {
    write.table(as.matrix(colnames(working_srobj)),
                sprintf("%s/%s.tsv", endclust_dir, filename),
                sep="\t", row.names = F, col.names = F)
  }    
}

Node_Write <- function(working_srobj, srobj_dir = NULL) {
  if (!is.null(srobj_dir)) {
    message(paste("saving srobj to", srobj_dir))
    saveRDS(working_srobj, sprintf("%s/subset_srobj.rds", srobj_dir))
  }
}

QC_Plotting <- function(working_srobj, fig_dir=NULL) {
  Idents(working_srobj) <- rep("all.data", ncol(working_srobj))
    
  VlnPlot(working_srobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = -1)
  ggsave(sprintf("%s/QC_vln_plot.pdf", fig_dir))
    
  #' basic QC scatter plots
  FeatureScatter(working_srobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  ggsave(sprintf("%s/QC_scatter_plot_nCount_pctMito.pdf", fig_dir))
  FeatureScatter(working_srobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ggsave(sprintf("%s/QC_scatter_plot_nCount_nFeat.pdf", fig_dir))
  
}

Plotting <- function(working_srobj,fig_dir=NULL) {

  cell.num.per.clust <- table(Idents(working_srobj))
  message("Cell counts per cluster: ",paste0(capture.output(cell.num.per.clust),collapse='\n'))

  if (!is.null(fig_dir)) {
    #' save basic cell counts
    cat(sprintf("Srobj: with %s genes & %s cells", nrow(working_srobj), ncol(working_srobj)),
        file = sprintf("%s/basic_counts.txt", fig_dir), sep="\n")
    cat(sprintf("\n Choose %s PCs \n", working_srobj@misc$nPCs), file = sprintf("%s/basic_counts.txt", fig_dir), append = T)
    cat("\nTable of num cells per cluster:\n", file = sprintf("%s/basic_counts.txt", fig_dir), append = T)
    write.table(cell.num.per.clust, file = sprintf("%s/basic_counts.txt", fig_dir), row.names = F, append = T)  
     
    #' basic umap of clusters (using umap-learn, as only it is allowed on premade graphs)
    working_srobj <- RunUMAP(working_srobj,
                             dims=working_srobj@misc$nPCs, 
                             n.neighbors = min(30L, ncol(working_srobj)-1))
    DimPlot(working_srobj, reduction = "umap", label = T)
    ggsave(sprintf("%s/cluster_umap.pdf", fig_dir))
    
    if (!is.null(working_srobj@meta.data$orig.ident)) {
      DimPlot(working_srobj, group.by = "orig.ident", reduction = "umap", label = F)
      ggsave(sprintf("%s/orig_ident_umap.pdf", fig_dir))
    }

    markers <- FindAllMarkers(working_srobj, assay = "RNA", only.pos = T, max.cells.per.ident = 2000)
    
    if ( !(nrow(markers) < 1) ) {
      fccol <- grep("FC",colnames(markers))
      top10markers <- markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>%  slice_max(markers[,fccol],n=10)
      write.table(top10markers, sprintf("%s/top10markers.tsv", fig_dir), sep="\t", row.names = F)
      write.table(markers, sprintf("%s/allmarkers.tsv", fig_dir), sep="\t", row.names = F)
      
      #' output heatmap of markers
      #if (nrow(top10markers) < 5) {top10markers <- markers %>% group_by(cluster) %>% slice_max(markers[,fccol],n=10)}
      tryCatch({ 
        working_srobj <- ScaleData(working_srobj, features = top10markers$gene)
        DoHeatmap(working_srobj, features = top10markers$gene, raster = F)
        ggsave(sprintf("%s/top10markers_heatmap.pdf", fig_dir))
      }, error = function(e) {message("Failed to produce heatmap"); print(paste("HEATMAP ERRORED ON:  ", e))})
    }
  }  
  return(working_srobj)
}

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
#' @param bias bias toward resolution choice, when resolution results are similar
#' @param figdir directory to save results graphs
#' @return seurat object with slot misc$resolution.choice
#' @export

chooseResolution_SilhouetteAnalysisParameterScan <- function(
  input.srobj, assay="RNA", n.pcs = input.srobj@misc$nPCs, sample.name =  format(Sys.time(), "%a-%b-%d-%X-%Y-%Z"),
  res.low = .01, res.high=3, res.n = 40, bias = "under", figdir=NULL) {
  
  ######## step 1: save the input seurat object as a new temporary object, 
  ########         dont want to overwrite or change the original one with all of the parameter scans
  
  srobj.tmp = input.srobj 
  # in case there have been other things calculated in the metadata, just cut down to simplify/avoid errors
  srobj.tmp@meta.data = srobj.tmp@meta.data[,c("nCount_RNA","nFeature_RNA")] # should just be the nUMI and nGene  
  
  ######## step 2: calculate the FindClusters over a large range of resolutions
  print("Performing parameter scan over multiple resolutions...")

  set.res = round(exp(seq(log(res.low), log(res.high), length.out=res.n)), digits=3)
  
  srobj.tmp = FindClusters(srobj.tmp, assay=assay, dims.use = n.pcs, k.param=ceiling(0.5*sqrt(ncol(srobj.tmp))),
                           resolution=set.res[1], save.SNN=T, plot.SNN=F,
                           force.recalc=TRUE, verbose=FALSE)
  
  res.that.work <- rep(T, length(set.res))
  for(i in 2:length(set.res)){
    tryCatch({
      srobj.tmp = FindClusters(
        srobj.tmp, assay=assay, resolution=set.res[i], reuse.SNN=TRUE, plot.SNN=FALSE, verbose=FALSE, force.recalc=FALSE)
    }, error = function(e) {res.that.work[i] <<- F; message(paste("errored on", set.res[i], "Flagging as doesn't work"))})
    print(paste("          ", round(100*i/length(set.res)), "% done with parameter scan", sep=""))
  }
  set.res <- set.res[res.that.work]
  
  
  ######## step 3: output plot of how the resolution changes the number of clusters you get
  n.clusters = vector(mode="numeric", length=length(set.res))
  names(n.clusters) = set.res
  for(i in 1:length(n.clusters)){
    n.clusters[i] = length(table(as.vector(srobj.tmp@meta.data[,paste0(assay, "_snn_res.", names(n.clusters)[i])])))
  }
  
  ######## step 4: calculate the silhouette width for each resolution
  print("Computing a silhouette width for each cell, for each resolution...")
  require(cluster)
  
  dist.temp = cor(t(srobj.tmp@reductions$pca@cell.embeddings[,n.pcs]), method="pearson")
  # consider changing to sampling 10% of each cluster, runs into issues calc silhouette score if there is only 1 cell in each cluster
  random.cells.choose = if ( nrow(dist.temp) > 500 ) sample(1:nrow(dist.temp), round(nrow(dist.temp)/10, digits=0)) else 1:nrow(dist.temp)
  dist.temp.downsample = dist.temp[random.cells.choose, random.cells.choose]
  sil.all.matrix = matrix(data=NA, nrow=nrow(dist.temp.downsample), ncol=0)
  
  for(i in 1:length(set.res)){
    clusters.temp = as.numeric(as.vector(
      srobj.tmp@meta.data[random.cells.choose,paste0(assay, "_snn_res.", set.res[i])]))
    if(length(table(clusters.temp))>1){
      sil.out = silhouette(clusters.temp, as.dist(1-as.matrix(dist.temp.downsample)))
      sil.all.matrix = cbind(sil.all.matrix, sil.out[,3])
    }
    if(length(table(clusters.temp))==1){
      sil.all.matrix = cbind(sil.all.matrix, rep(0, length(clusters.temp)))
    }
    print(paste("          ", round(100*i/length(set.res)), "% done with silhouette calculation", sep=""))
    
  }
  
  
  ######## step 5: calculate summary metric to compare the silhouette distributions,
  ########         average has worked well so far... could get fancier
  
  print("Identifying a best resolution to maximize silhouette width")
  sil.average = setNames(colMeans(sil.all.matrix), set.res)
  sil.medians <- setNames(apply(sil.all.matrix, 2, median), set.res)
  
  ######## step 6: automate choosing resolution that maximizes the silhouette 
  hist.out = hist(sil.average, length(sil.average)/1.2,  plot=FALSE)
  
  #  take the ones that fall into the top bin, 
  #  and the max OR MIN of those  ******* can change this to under vs over cluster
  if(bias=="over"){
    resolution.choice = as.numeric(max(
      names(sil.average[which(sil.average>hist.out$breaks[length(hist.out$breaks)-1])])))
  }
  if(bias=="under"){
    resolution.choice = as.numeric(min(
      names(sil.average[which(sil.average>hist.out$breaks[length(hist.out$breaks)-1])])))
  }
  
  # get the silhouette of the best resolution: 
  silhouette.best = as.numeric(sil.average[paste(resolution.choice)])
  
  print(paste("Best Resolution Choice: ", resolution.choice, ", with average silhouette score of: ",
              round(silhouette.best, digits=3), ", giving ", as.numeric(n.clusters[paste(resolution.choice)]),
              " clusters", sep=""))
  
  ######### step 7: output plot and data 
  
  if (!is.null(figdir)) {
    
    print(paste0("Ouptutting summary statistics and returning seurat object... ",
                 "This will create a pdf in your output directory,",
                 " and will return your input seurat object ammended with the best choice",
                 " for clusters (found as Best.Clusters in the meta.data matrix, and set to your new ident)..."))
    
    pdf(paste(figdir, "/", sample.name, ".pdf", sep=""),
        width=10, height=4, useDingbats=FALSE)
    par(mfrow=c(1,3))
    # Resolution vs # of Clusters
    plot(set.res, n.clusters, col="black", pch=19,
         type="p", xlab="Resolution", ylab="# Clusters",
         main="Resolution vs. # Clusters")
    
    # Resolution vs Average Silhouette
    plot(set.res, sil.average, col="black", pch=19,
         type="p", xlab="Resolution", ylab="Average Silhouette",
         main="Resolution vs. Average Silhouette")
    points(set.res, sil.medians, col="red", pch=15)
    abline(h=hist.out$breaks[length(hist.out$breaks)-1], col="firebrick3", lty=2)
    abline(v=resolution.choice, col="dodgerblue2", lty=2)
    
    # N Clusters vs Average Silhouette
    plot(n.clusters, sil.average, col="black", pch=19,
         type="p", xlab="# Clusters", ylab="Average Silhouette",
         main="# Clusters vs. Average Silhouette")
    points(n.clusters, sil.medians, col="red", pch=15)
    abline(h=hist.out$breaks[length(hist.out$breaks)-1], col="firebrick3", lty=2)
    abline(v=as.numeric(n.clusters[paste(resolution.choice)]), col="dodgerblue2", lty=2)
    dev.off()
  }
  
  ######## step 8: return the original seurat object, with the metadata containing a 
  ########         concatenated vector with the clusters defined by the best choice here,
  ########         as well as the ident set to this new vector
  
  Best.Clusters = srobj.tmp@meta.data[,paste0(assay, "_snn_res.", resolution.choice)]
  
  input.srobj$Best.Clusters = Best.Clusters
  Idents(input.srobj) = input.srobj$Best.Clusters
  input.srobj@misc <- list("resolution.choice" = resolution.choice)
  return(input.srobj)
}
