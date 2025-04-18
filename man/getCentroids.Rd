% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ARBOLhelper.R
\name{getCentroids}
\alias{getCentroids}
\title{Calculate pvclust() tree (a binary tree of distances between end-clusters) for ARBOL results
tree based on euclidean distance between cluster centroids based on gene medians with complete linkage}
\usage{
getCentroids(
  srobj = srobj,
  tree_reduction = tree_reduction,
  reduction_dims = reduction_dims,
  centroid_method = "mean",
  centroid_assay = "RNA",
  centroid_layer = "scale.data",
  gene_list = Features(srobj)
)
}
\arguments{
\item{srobj}{a seurat object with ARBOL 'tierNident' column}

\item{tree_reduction}{cell x gene reduction space to work from. 'centroids' uses full cell x gene.}

\item{reduction_dims}{the dimensions of the reduction layer to use for centroid calculation. defaults to 1:25}

\item{centroid_method}{function to calculate gene centroids per ARBOL end cluster, see Matrix.utils::aggregate.Matrix for options}

\item{centroid_assay}{if using cell x gene data (not any srobj@reduction), the assay within which to calculate centroids}

\item{gene_list}{if using cell x gene data (not any srobj@reduction), genes to include in centroid calculation}
}
\value{
the input seurat object with pvclust tree in srobj@misc$pvclust
}
\description{
Calculate pvclust() tree (a binary tree of distances between end-clusters) for ARBOL results
tree based on euclidean distance between cluster centroids based on gene medians with complete linkage
}
\examples{
centroids <- getCentroids(srobj = srobj, tree_reduction = tree_reduction, reduction_dims = reduction_dims,
                        centroid_method = centroid_method, centroid_assay = centroid_assay, gene_list = gene_list)
}
