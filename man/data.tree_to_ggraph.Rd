% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ARBOLhelper.R
\name{data.tree_to_ggraph}
\alias{data.tree_to_ggraph}
\title{Convert data.tree to ggraph
Converts a data.tree structure to a ggraph object for easier plotting, by first converting the data.tree structure to a dendrogram via a custom as.dendrogram.Node function, then to a ggraph object, and finally restoring annotations by joining a data frame of data.tree attributes.
This function is used to convert annotated binary phylogeny trees to ggraph objects, offering a more efficient alternative to previous versions that used Newick text conversion.}
\usage{
data.tree_to_ggraph(
  data.tree,
  categories,
  diversities,
  counts,
  totals,
  heightAttribute = "plotHeight"
)
}
\arguments{
\item{data.tree}{A data.tree object representing an annotated binary phylogeny tree.}

\item{categories}{A character vector of category attributes present in the data.tree 
object to be included in the ggraph.}

\item{diversities}{A character vector of diversity attributes present in the data.tree 
object to be included in the ggraph.}

\item{counts}{A character vector of count attributes present in the data.tree object 
to be included in the ggraph.}

\item{totals}{A character vector of total attributes present in the data.tree object 
to be included in the ggraph.}

\item{heightAttribute}{The name of the attribute used to define the height of the nodes 
in the resulting ggraph object. Defaults to 'plotHeight'.}
}
\value{
A tbl_graph object suitable for plotting with ggraph, containing the structure 
        of the data.tree and annotations for nodes and edges.
}
\description{
Convert data.tree to ggraph
Converts a data.tree structure to a ggraph object for easier plotting, by first converting the data.tree structure to a dendrogram via a custom as.dendrogram.Node function, then to a ggraph object, and finally restoring annotations by joining a data frame of data.tree attributes.
This function is used to convert annotated binary phylogeny trees to ggraph objects, offering a more efficient alternative to previous versions that used Newick text conversion.
}
\examples{
# Assuming `myTree` is a data.tree object with necessary annotations:
ggraph_tree <- data.tree_to_ggraph(data.tree = myTree,
                                   categories = c("Category1", "Category2"),
                                   diversities = c("Diversity1", "Diversity2"),
                                   counts = c("Count1", "Count2"),
                                   totals = c("Total1", "Total2"))
# Now `ggraph_tree` can be plotted with ggraph:
ggraph(ggraph_tree) + geom_edge_link() + geom_node_point()
}
