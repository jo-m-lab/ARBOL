#' @name subclusteringTree
#' @title Creates data tree object for ARBOL run
#' 
#' @description Creates data tree object for ARBOL run, adds it to seurat object along with a graph of it
#'   calls \code{prepSubclusteringMetadata()} and \code{propagateTree()}
#' 
#' @param srobj a seurat object with ARBOL 'tierNident', 'CellID', and 'sample' columns. 
#' @param categories columns in metadata for which you want to track categories present per node. Also finds the majority per node
#' @param diversities columns in metadata for which you want to calculate diversity per node 
#' @param diversity_metric one of 'shannon', 'simpson', or 'invsimpson'
#' @param counts columns in metadata for which you want to count each value per node
#' @param totals attributes for which to sum values per node.
#' 
#' @return the input seurat object with tiered clustering tree in srobj@@misc$subclusteringTree, 
#' plot of tree to srobj@@misc$subclusteringViz, and ggraph object to srobj@@misc$cluster_ggraph
#' @examples
#' srobj <- subclusteringTree(srobj)
#' @export
#' 
#' @importFrom tidygraph as_tbl_graph activate
#' @importFrom ape as.phylo
#' 
subclusteringTree <- function(
        srobj,
        categories = 'sample',
        diversities = 'sample',
        diversity_metric = 'simpson',
        counts = 'sample',
        totals = 'nCount_RNA') {
    
    if(any(srobj@meta.data$tierNident %like% '/')) {
        message('tierNident column cannot contain "/"')
        stop()
    }
    
    if (!is.element('sample',diversities)) {
        diversities = c('sample',diversities)
    }
    
    if (!is.element('sample',categories)) {
        categories = c('sample',categories)
    }
    
    if (!is.element('sample',counts)) {
        counts = c('sample',counts)
    }
    
    srobj <- suppressMessages(tierN_diversity(
        srobj,
        diversity_attributes = diversities,
        diversity_metric = diversity_metric))
    
    treemeta <- suppressMessages(prepSubclusteringMetadata(
        srobj,
        categorical_attributes = categories, 
        diversity_attributes = diversities, 
        numerical_attributes = counts,
        total_attributes = totals))
    
    ARBOLtree <- suppressMessages(as.Node(treemeta))
    
    Atree <- suppressMessages(propagateTree(
        ARBOLtree,
        srobj = srobj,
        diversity_attributes = diversities, 
        categorical_attributes = categories,
        numerical_attributes = counts,
        total_attributes = totals))
    
    #obtain counts columns from propagated data tree
    attrs <- Atree$attributesAll
    count_cols <- attrs[grep(sprintf('^(%s)_n',paste0(counts,collapse='|')),attrs)]
    
    ARBOLdf <- do.call(
        allotedTreeToDF,
        c(Atree,'n','pathString', categories, paste0(categories,'_majority'),
          paste0(diversities,'_diversity'),count_cols,totals,'numChildren'))
    
    # simple code call of translation to df disallowing custom diversity_attributes
    # ARBOLdf <- allotedTreeToDF(Atree, 'tier1', 'n', 'pathString', 'sample_diversity')
    ARBOLphylo <- as.phylo(ARBOLtree)
    
    # for original ARBOL clustering tree, set all edge lengths to 1
    ARBOLphylo$edge.length <- rep(1, length(ARBOLphylo$edge.length))
    
    # convert to tbl_graph object to allow easy plotting with ggraph
    x <- as_tbl_graph(ARBOLphylo, directed = TRUE)
    x <- activate(x, nodes)
    nm_keep <- c("levelName", "n", "numChildren",
                 categories, paste0(categories, "_majority"),
                 paste0(diversities, "_diversity"),
                 count_cols, totals)
    
    x <- x %>% left_join(ARBOLdf[, nm_keep]) ## TODO Need to identify by which columns
    
    x <- activate(x, edges)
    x$tier <- sapply(gregexpr("\\.", name), function(x) sum(x >= 0))
    
    bt0 <- ggraph(x, layout = 'tree', circular=T) + 
        geom_edge_diagonal() + geom_node_point(size=0.3) + theme_classic()
    
    srobj@misc$subclusteringTree <- Atree
    srobj@misc$subclusteringViz <- bt0
    srobj@misc$cluster_ggraph <- x
    
    return(srobj)
}
