% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ARBOL.R
\name{GenTieredClusters}
\alias{GenTieredClusters}
\title{Performs iterative clustering (ARBOL) on a v4 Seurat object}
\usage{
GenTieredClusters(
  srobj,
  cluster_assay = "SCT",
  cells = NULL,
  tier = 0,
  clustN = 0,
  PreProcess_fun = PreProcess_sctransform,
  ChooseOptimalClustering_fun = ChooseOptimalClustering_default,
  saveSROBJdir = NULL,
  figdir = NULL,
  SaveEndNamesDir = NULL,
  SaveEndFileName = NULL,
  min_cluster_size = 100,
  max_tiers = 10,
  EnoughDiffUp = 5,
  EnoughDiffDown = 5,
  res_scan_step = 5,
  res_scan_min = 0.01,
  res_scan_max = 3,
  res_scan_n = 40,
  tierAllowedRecomb = 0,
  harmony_var = NULL,
  DownsampleNum = 7500
)
}
\arguments{
\item{srobj}{v4 seurat object}

\item{cluster_assay}{assay to use for clustering defaults to "SCT"}

\item{cells}{cellnames if tiered clustering should start on subset of object}

\item{tier}{starting level defaults to 0}

\item{clustN}{cluster starting from default to 0}

\item{PreProcess_fun}{function to use for preproccessing defaults to PreProcess_sctransform. PreProcess_sctransform_harmony now available.}

\item{ChooseOptimalClustering_fun}{function that returns srobj with clusters in `srobj$Best.Clusters` after choosing optimal clustering resolution}

\item{saveSROBJdir}{where to save seurat objects for each tier and cluster, if null does not save}

\item{figdir}{where to save QC figures for each tier and cluster, if null does not save}

\item{SaveEndFileName}{prefix for all end cluster files}

\item{min_cluster_size}{minimum number of cells to allow further clustering. defaults to 100}

\item{max_tiers}{maximum number of tiers to allow further clustering. defaults to 10}

\item{EnoughDiffUp}{minimum number of up-regulated genes to call clusters unique. Differential expression is performed when clustering finds 2 clusters. defaults to 5}

\item{EnoughDiffDown}{minimum number of down-regulated genes. If either up or down is not met, the 2 clusters are joined, and further clustering is stopped. defaults to 5}

\item{res_scan_step}{number of resolutions of decreasing score before an early silhouette analysis resolution scan stop}

\item{res_scan_min}{the smallest resolution to scan}

\item{res_scan_max}{the highest resolution to scan (default is 3 due to legacy, we recommend lowering this)}

\item{res_scan_n}{the number of resolutions to scan. A logarithmic sequence of numbers between min and max are scanned}

\item{tierAllowedRecomb}{minimum tier where differential expression can be called to decide on recombination. defaults to 0. clustering may stop early when clustering finds 2 clusters with high cell numbers, as Wilcoxon effect sizes may be low.}

\item{harmony_var}{variable over which to iterate harmony integration}

\item{saveEndNamesDir}{where to save directory of end clusters, if null does not save}
}
\value{
list of lists with all seurat objects (highly recommend using folder arguments for saving outputs)
}
\description{
Iteratively clusters seurat object from single cell datasets, choosing optimum resolution parameters at each stage of clustering. Outputs QC plots for each tier and stage.

can return output directories figs/ srobjs/ (see examples) containing subset seurat objects and QC plots matching tiered structure of found in dataset

can also return endclusts/ folder which contains tsv files named with the path to the endcluster (T0C0T1C1.tsv) and contains the cellnames (one per line) that comprise that cluster
}
\examples{
srobj <- readRDS("/path/to/seurat_object.rds")
tiers <- GenTieredclusters(srobj,
                           saveSROBJdir = "~/tieredoutput/srobjs",
                           figdir = "~/tieredoutput/figdir",
                           SaveEndNamesDir = "~/tieredoutput/endclusts")

}
