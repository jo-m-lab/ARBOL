#! /usr/bin/Rscript

#' data.tree to phyloseq object conversion with custom ToNewick function that uses pathString as node names. 
#' @param node A data tree object
#' @return Newick format tree
#' @examples
#' txt <- ToNewickPS(divtree)
#' ARBOLphylo <- ape::read.tree(text=txt)
#' @export
#' Functions for data.tree to phyloseq object conversion with custom ToNewick function that uses pathString as node names. This is useful when node names are redundant throughout a tree, which happens in binary tree calculation with pvclust

ToNewickPS <- function(node, heightAttribute = DefaultPlotHeight, ...) {

  deparse <- function(x) {
    name <- stringi::stri_replace_all_fixed(x$pathString, " ", "_")
    name <- stringi::stri_replace_all_fixed(x$pathString, "/", ".")
    name <- stringi::stri_replace_all_fixed(name, ",", "")
    if(!isRoot(x) && length(heightAttribute) > 0) {
      edge <- GetAttribute(x$parent, heightAttribute, ...) - GetAttribute(x, heightAttribute, ...) 
      me <- paste0(name, ":", edge)
    } else {
      me <- name
    }
    return(me)
  }
  
  Newick <- function(x) {
    if(x$isLeaf) {
      return (deparse(x))
    }
    chNewick <- sapply(x$children, Newick)
    chNewickStr <- paste(chNewick, collapse = ",")
    res <- paste0("(", chNewickStr, ")", deparse(x))
  }
  
  res <- Newick(node)
  res <- paste0(res, ";")
  return (res)
  
}

#Very slow, newick way
as.phylo.NodeNW <- function(x, heightAttribute = DefaultPlotHeight, ...) {
  txt <- ToNewickPS(x, heightAttribute)
  return (ape::read.tree(text = txt))
}

#fast data.tree.Node to ape tree object conversion
as.phylo.NodePS <- function(x, heightAttribute = 'plotHeight') {
  datadend <- as.dendrogram.NodePS(x, heightAttribute = heightAttribute)
  order.dendrogram(datadend) <- seq(1,length(na.omit(get_nodes_attr(datadend, "label"))))
  result = as.phylo(datadend)
  labels(result) <- labels(result) %>% str_replace_all('/','.')
  result$node.label <- datadend %>% get_nodes_attr('label',include_leaves=F) %>% str_replace_all('/','.') %>% na.omit
  return(result)
}

#data.tree.Node to ape tree object conversion with internal nodes named by pathString and leaf nodes named without pathString (only endname)
as.phylo.NodeI <- function(x, heightAttribute = 'plotHeight') {
  result = as.phylo.Node(x)
  result$node.label <- allotedTreeToDF(x,'pathString','isLeaf') %>% filter(!isLeaf) %>% pull(pathString) %>% str_replace_all('/','.')
  return(result)
}


#' data.tree to ggraph object conversion with custom Node -> dend object conversion function that uses pathString as node names. 
#' @param object A data tree object
#' @param heightAttribute plot height attribute. Default is plotHeight
#' @return dendrogram object
#' @examples
#' txt <- ToNewickPS(divtree)
#' ARBOLphylo <- ape::read.tree(text=txt)
#' @export
#' Functions for data.tree to dendrogram object conversion using Node pathString as node names. This is useful when node names are redundant throughout a tree, which happens in taxonomy calculation with pvclust

as.dendrogram.NodePS <- function(object, heightAttribute = 'plotHeight', edgetext = FALSE,
                               namefield = 'pathString', ...) {
  node <- object

  height <- as.vector(GetAttribute(node, heightAttribute))

  if (node$isLeaf) {
    res <- node$value
    if (is.null(res)) res <- 0
    res <- structure(res, 
                     label = node[[namefield]], 
                     members = 1,
                     height = height,
                     leaf = node$isLeaf,
                     class = "dendrogram")
    
  } else {

    res <- unname(lapply(node$children, FUN = function(x) as.dendrogram.NodePS(x, heightAttribute, ...)))
    res <- structure(res, 
                     label = node[[namefield]],
                     members = node$leafCount,
                     midpoint = node$midpoint,
                     height = height,
                     class = "dendrogram")
    
    if (edgetext) attr(res, "edgetext") <- node[[namefield]]
    
  }
  
  return (res)
  
}


#' Load 'endclusters' folder from R ARBOL output as a dataframe
#' @param folder R ARBOL endclusts folder
#' @param maxtiers max_tiers parameter used in ARBOL
#' @return dataframe with tier membership per cell including tierNident, the endcluster for that cell. cell barcode column CellID
#' @examples
#' tiers <- LoadTiersAsDF(folder='./endclusts',maxtiers=10)
#' #metadata from Seurat object
#' metadata <- metadata %>% left_join(tiers, by='CellID')
#' @export
LoadTiersAsDF <- function(folder='./endclusts',maxtiers=10) {
    tierfiles <- list.files(folder,full.names=TRUE,recursive=TRUE)
    sample_strings <- sub('\\..*$','',basename(tierfiles))
    sample_strings <- sub('_T0C0_', '',sample_strings)
    tiers <- map2(tierfiles, sample_strings, ~fread(.x,header=FALSE) %>% mutate(id = .y))
    tiers <- rbindlist(tiers) %>% data.frame
    tiers$id <- gsub('_','.',tiers$id)
    tiers <- tiers %>% separate(id,into=paste0('tier',1:maxtiers),sep='\\.',remove=FALSE)
    tiers <- tiers %>% dplyr::rename(CellID=V1)
    tiers <- tiers %>% data.frame %>% dplyr::rename(tierNident=id)
    return(tiers)
}

#' Calculate diversity index per group in a dataframe (adds '%s_diversity' column),
#' diversity metric options same as those available in vegan diversity() function
#' @param df a dataframe with species and groups columns
#' @param species species column to stratify per group
#' @param group groups column for which to calculate diversity
#' @param diversity_metric one of 'shannon', 'simpson', or 'invsimpson'
#' @return dataframe with species and diversity columns
#' @examples
#' dataframe <- dataframe %>% left_join(diversityPerGroup(dataframe, species=pathString, group='sample')) %>% 
#'                            suppressMessages
#' metadata <- metadata %>% left_join(diversityPerGroup(metadata, species=tierNident, group='sample')) %>% 
#'                          suppressMessages
#' @export
diversityPerGroup <- function(df, species, group, diversity_metric = 'simpson') {
  #enquo parameters to allow dplyr calls
    divcols <- enquo(species)
    #count groups per species
    tierNcount <- df %>% group_by_at(vars(!!divcols)) %>% dplyr::count(.data[[group]])
    #obtain total N per group
    nums <- tierNcount %>% summarize(nums=sum(n))
    #add N to count dataframe
    tierNcount <- tierNcount %>% left_join(nums) %>% suppressMessages
    #calculate fraction of sample per species
    tierNcount <- tierNcount %>% mutate(presence = n/nums)
    #Pivot table to allow vegan::diversity call
    tierNwide <- tierNcount %>% select(-presence,-nums) %>% pivot_wider(names_from=.data[[group]],values_from=n) %>% ungroup %>% data.frame
    #rownames are necessary, here hacking rownames from enquo with alphanumeric gsub
    tierNwide <- tierNwide %>% column_to_rownames(str_replace_all(as.character(vars(!!divcols)),"[^[:alnum:]]", ""))
    #NA to 0 necessary for vegan::diversity
    tierNwide[is.na(tierNwide)] <- 0
    #calculate diversity. Can change simpson to other types. Second use of enquo hack just to be able to name it the diversity name, idk if this is necessary
    simps <- diversity(tierNwide,diversity_metric) %>% data.frame %>% rownames_to_column(str_replace_all(as.character(vars(!!divcols)),"[^[:alnum:]]", ""))
    #Add 2 column names: your species category (i usually use pathString pulled from the node on the tree), and diversity. This joins to the dataframe by species category.
    colnames(simps) <- c(str_replace_all(as.character(vars(!!divcols)),"[^[:alnum:]]", ""),sprintf('%s_diversity',group))
    return(simps)
}

#' Calculate diversity index for a set of IDs in a dataframe
#' @param df a dataframe with group information per cell
#' @param group groups column for which to calculate diversity
#' @param diversity_metric one of 'shannon', 'simpson', or 'invsimpson'
#' @return a list of diversity per group per cell
#' @examples
#' in a data tree
#' node[[sprintf('%s_diversity',y)]] <- SIperIDs(metadata, group=y)
#' directly on a dataframe
#' metadata$condition_diversity <- SIperIDs(metadata, group=condition)
#' @export
SIperIDs <- function(df, group, diversity_metric = 'simpson') {
    #count groups
    tierNcount <- df %>% dplyr::count(.data[[group]])
    #obtain total N per group
    nums <- tierNcount %>% summarize(nums=sum(n))
    #add N to count dataframe
    tierNcount$nums <- nums$nums
    #calculate fraction of each sample
    tierNcount <- tierNcount %>% mutate(presence = n/nums)
    #Pivot table to allow vegan::diversity call
    tierNwide <- tierNcount %>% select(-presence,-nums) %>% pivot_wider(names_from=.data[[group]],values_from=n) %>% ungroup %>% data.frame
    #NA to 0 necessary for vegan::diversity
    tierNwide[is.na(tierNwide)] <- 0
    #calculate diversity. Can change simpson to other types. Second use of enquo hack just to be able to name it the diversity name, idk if this is necessary
    return(diversity(tierNwide,diversity_metric))
}


#' Prepare ARBOL seurat object metadata for tree building
#' @param srobj a seurat object with ARBOL 'tierNident' and 'sample' columns
#' @param maxtiers max_tiers parameter used in ARBOL
#' @param numerical_attributes numerical attributes are propagated up a tree as counts of each unique value per node.
#' creates a new column per unique value per attribute. counts are added as node$attribute_n_value
#' @param categorical_attributes categorical attributes are propagated up a tree as unique values per node. also calculates a majority per node
#' majority added as node$attribute_majority
#' @param diversity_attributes attributes for which to calculate diversity in each node. 
#' Supports simpson's, inverse simpson's, and shannon index as implemented in vegan diversity()
#' @return a metadata dataframe ready for tree building
#' @examples
#' prepSubclusteringMetadata(srobj, maxtiers=10)
#' @export
prepSubclusteringMetadata <- function(srobj,maxtiers=10,categorical_attributes,diversity_attributes,numerical_attributes) {
    meta <- srobj@meta.data

    meta <- spread_tierN(meta,max_tiers=maxtiers)

    jointb <- meta %>% group_by(tierNident) %>% mutate(n=n()) %>% 
          dplyr::select(CellID,sample,tierNident,n,all_of(categorical_attributes),
                        all_of(paste0(diversity_attributes,'_diversity')))

    categorydf <- jointb %>% summarize(across(categorical_attributes, ~ list(paste(unique(.x),collapse=', '))))
    divdf <- jointb %>% summarize_at(paste0(diversity_attributes,'_diversity'),unique)

    jointb <- jointb %>% select(-all_of(categorical_attributes)) %>% 
                summarize(ids = list(CellID),n=unique(n))

    numericaltb <- meta %>% dplyr::select(CellID,sample,tierNident,all_of(numerical_attributes))

    for (x in numerical_attributes) {
      attribute <- enquo(x)
      tierNcount <- numericaltb %>% group_by(tierNident) %>% dplyr::count(.data[[attribute]])
      vals <- unique(tierNcount[[x]])
      countWide <- tierNcount %>% pivot_wider(names_from=all_of(x),values_from=n)
      colnames(countWide) <- c('tierNident',sprintf('%s_n_%s',x,vals))
      numericaltb <- numericaltb %>% left_join(countWide)
    }

    numericaltb[is.na(numericaltb)] <- 0

    numericaltb <- numericaltb %>% select(-CellID,-sample,-all_of(numerical_attributes)) %>% distinct

    treemeta <- meta %>% select(-all_of(categorical_attributes)) %>% left_join(jointb) %>% left_join(categorydf) %>% left_join(divdf) %>% left_join(numericaltb)
    
    return(treemeta)
}

#' Creates data tree object for ARBOL run, adds it to seurat object along with a graph of it
#' calls prepSubclusteringMetadata() and propagateTree()
#' @param srobj a seurat object with ARBOL 'tierNident', 'CellID', and 'sample' columns. 
#' @param categories columns in metadata for which you want to track categories present per node. Also finds the majority per node
#' @param diversities columns in metadata for which you want to calculate diversity per node 
#' @param diversity_metric one of 'shannon', 'simpson', or 'invsimpson'
#' @param counts columns in metadata for which you want to count each value per node
#' @return the input seurat object with tiered clustering tree in srobj@@misc$subclusteringTree, 
#' plot of tree to srobj@@misc$subclusteringViz, and ggraph object to srobj@@misc$cluster_ggraph
#' @examples
#' srobj <- subclusteringTree(srobj)
#' @export
subclusteringTree <- function(srobj, categories = 'sample', diversities = 'sample', diversity_metric = 'simpson', counts = 'sample') {

  if (!is.element('sample',diversities)) {
    diversities = c('sample',diversities)
  }

  if (!is.element('sample',categories)) {
    categories = c('sample',categories)
  }

  if (!is.element('sample',counts)) {
    counts = c('sample',counts)
  }

  srobj <- suppressMessages(tierN_diversity(srobj, diversity_attributes = diversities, diversity_metric = diversity_metric))
  
  treemeta <- suppressMessages(prepSubclusteringMetadata(srobj, categorical_attributes = categories, diversity_attributes = diversities, numerical_attributes = counts))

  ARBOLtree <- suppressMessages(as.Node(treemeta))

  Atree <- suppressMessages(propagateTree(ARBOLtree, srobj=srobj, diversity_attributes = diversities, categorical_attributes = categories, numerical_attributes = counts))

  #obtain counts columns from propagated data tree
  attrs <- Atree$attributesAll
  count_cols <- attrs[grep(sprintf('^(%s)_n',paste0(counts,collapse='|')),attrs)]

  ARBOLdf <- do.call(allotedTreeToDF, c(Atree,'n','pathString', categories, paste0(categories,'_majority'), paste0(diversities,'_diversity'),count_cols))

  #simple code call of translation to df disallowing custom diversity_attributes
  #ARBOLdf <- allotedTreeToDF(Atree, 'tier1', 'n', 'pathString', 'sample_diversity')

  ARBOLphylo <- as.phylo(ARBOLtree)
  #for original ARBOL clustering tree, set all edge lengths to 1
  ARBOLphylo$edge.length <- rep(1,ARBOLphylo$edge.length %>% length)

  #convert to tbl_graph object to allow easy plotting with ggraph
  x <- tidygraph::as_tbl_graph(ARBOLphylo,directed=T) %>% activate(nodes) %>% 
      left_join(ARBOLdf %>% select(name=levelName,n,all_of(categories),all_of(paste0(categories,'_majority')),
                                    all_of(paste0(diversities,'_diversity')),all_of(count_cols)))

  x <- x %>% activate(edges) #%>% left_join(ARBOLdf %>% select(to=i))

  x <- x %>% activate(nodes) %>% mutate(tier = str_count(name, "\\."))

  bt0 <- ggraph(x, layout = 'tree', circular=T) + 
  geom_edge_diagonal() + geom_node_point(size=0.3) + theme_classic()

  srobj@misc$subclusteringTree <- Atree
  srobj@misc$subclusteringViz <- bt0
  srobj@misc$cluster_ggraph <- x

  return(srobj)
}


#' Propagate attributes up a data.tree object to enable manipulation and visualization of all nodes
#' Automatically run when building trees using subclusteringTree or ARBOLcentroidTaxonomy
#' @param srobj a seurat object with ARBOL 'tierNident' and 'sample' columns
#' @param ARBOLtree data.tree object
#' @param numerical_attributes numerical attributes are propagated up a tree as counts of each unique value per node.
#' creates a new column per unique value per attribute. counts are added as node$attribute_n_value
#' @param total_attributes total attributes are numerical attributes that you want to propagate up the tree by summing
#' the value of each node.
#' @param categorical_attributes categorical attributes are propagated up a tree as unique values per node. 
#' also calculates a majority per node - majority added as node$attribute_majority
#' @param diversity_attributes attributes for which to calculate diversity in each node. 
#' Supports simpson's, inverse simpson's, and shannon index as implemented in vegan diversity()
#' @param diversity_metric one of 'shannon', 'simpson', or 'invsimpson'
#' added as node$attribute_diversity
#' @return data.tree object with attributes propagated to all nodes
#' @examples
#' arbolTree <- propagateTree(arbolTree,srobj=srobj, categorical_attributes = categories, diversity_attributes = diversities, numerical_attributes = NA)
#' @export
propagateTree <- function(ARBOLtree, srobj, numerical_attributes = NA, categorical_attributes = NA, total_attributes = NA,
                          diversity_attributes = 'sample', diversity_metric = 'simpson') {

    #calculate number of children per node
    ARBOLtree$Do(function(node) node$numChildren <- node$children %>% length)

    #calculate total n samples for diversity propagation
    totalsamples <- srobj@meta.data$sample %>% unique %>% length

    #add number of cells per node to totals
    if (!is.element('n',total_attributes)) {
      total_attributes = c('n',total_attributes)
    }

    #propagate totals variables up tree by summing the end-nodes
    if(is.character(total_attributes)) {
      for (y in total_attributes){
        ARBOLtree$Do(function(node) node[[y]] <- Aggregate(node, attribute = y, aggFun = sum), traversal = "post-order")
      }
    }

    #propagate numerical variables up tree by counting occurrences of each variant z of category y in each node
    if(is.character(numerical_attributes)) {                     
        for (y in numerical_attributes){
            uniqs <- unique(srobj@meta.data[[y]])
            attrs = sprintf('%s_n_%s',y,uniqs)
            for (attr in attrs) {
                ARBOLtree$Do(function(node) node[[attr]] <- Aggregate(node, attribute = attr, aggFun = sum), traversal = "post-order")
            }            
        }
    }

    #propagate list of cell barcodes per node up the tree
    ARBOLtree$Do(function(node) node[['ids']] <- Aggregate(node, attribute = 'ids', aggFun = c), traversal = "post-order")
    #also samples. introducing this line to the function enforces "sample" column in srobj metadata
    ARBOLtree$Do(function(node) node[['samples']] <- Aggregate(node, attribute = 'samples', aggFun = c), traversal = "post-order")

    #propagate additional categorical variables
    if(is.character(categorical_attributes)) { 
        for (y in categorical_attributes) {
            ARBOLtree$Do(function(node) {
              ids <- node$ids %>% unlist; meta <- srobj@meta.data %>% filter(CellID %in% ids);
                                node[[y]] <- unique(meta[[y]])
            })
            ARBOLtree$Do(function(node) {
              ids <- node$ids %>% unlist; meta <- srobj@meta.data %>% filter(CellID %in% ids);
                                node[[sprintf('%s_majority',y)]] <- names(which.max(table(meta[[y]])))
            })
        }
    }
                         
    #make samples only unique values
    ARBOLtree$Do(function(node) node$samples <- unique(unlist(node$samples)))

    #calculate diversity attributes per node
    for (y in diversity_attributes) {
        ARBOLtree$Do(function(node) {ids <- node$ids %>% unlist; meta <- srobj@meta.data %>% filter(CellID %in% ids);
                                node[[sprintf('%s_diversity',y)]] <- SIperIDs(meta, group=y, diversity_metric = diversity_metric)
        })
    }

    #call internal nodes as certain types depending on most prevalent major type annotation
    ARBOLtree$Do(function(node) {node$tier1 <- names(which.max(table(unlist(node$tier1))))})

    return(ARBOLtree)
}

#' Calculate pvclust() tree (a binary tree of distances between end-clusters) for ARBOL results
#' tree based on euclidean distance between cluster centroids based on gene medians with complete linkage
#' @param srobj a seurat object with ARBOL 'tierNident' column
#' @param tree_reduction either 'centroids', which calculates centroids among all genes, or any reduction slot in srobj
#' @param hclust_method any hierarchical clustering method implemented in pvclust::pvclust(method.hclust), defaults to 'complete'
#' @param distance_method any distance method implemented in pvclust::pvclust(method.dist) - one of "correlation", "abscor", "uncentered", "euclidean" -
#' or cosine (no quotes) as implemented in ARBOL::cosine. you may also write your own function that returns a dist object, as in pvclust::pvclust()
#' @param centroid_method the function used to calculate centroids in the tree_reduction matrix, as implemented in Matrix.utils::aggregate.Matrix(fun)
#' currently, sum, count, mean, and median are supported
#' @param centroid_assay if using cell x gene data (not any srobj@@reduction), the assay within which to calculate centroids
#' @param gene_list if using cell x gene data (not any srobj@@reduction), genes to include in centroid calculation. passes to getCentroids
#' @param reduction_dims the dimensions of the reduction slot to use for centroid calculation. defaults to 1:25
#' @return only pvclust object
#' @export
centroidTaxonomy_noSeurat <- function(srobj, tree_reduction = 'centroids', hclust_method = 'complete',
                                distance_method = 'euclidean', centroid_method = 'mean', 
                               centroid_assay = 'SCT', reduction_dims = 1:25, gene_list = rownames(srobj[["RNA"]]@data), nboot = 1) {

    if (tree_reduction == 'centroids' | tree_reduction %in% names(srobj@reductions)) {
      centroids <- getCentroids(srobj = srobj, tree_reduction = tree_reduction, reduction_dims = reduction_dims,
                         centroid_method = centroid_method, centroid_assay = centroid_assay, gene_list = gene_list)

      result <- pvclust(centroids, method.dist=distance_method, method.hclust=hclust_method, nboot=nboot)

      return(result)
    }

    else {
      message(sprintf('%s isnt an ARBOL implemented reduction for tree building, or does not exist in srobj@reductions',tree_reduction))
    }
  
}


#' Calculate pvclust() tree (a binary tree of distances between end-clusters) for ARBOL results
#' default method creates a tree based on euclidean distance between cluster centroids based on gene medians with complete linkage
#' @param srobj a seurat object with ARBOL 'tierNident' column
#' @param tree_reduction either 'centroids', which calculates centroids among all genes, or any reduction slot in srobj
#' @param hclust_method any hierarchical clustering method implemented in pvclust::pvclust(method.hclust), defaults to 'complete'
#' @param distance_method any distance method implemented in pvclust::pvclust(method.dist) - one of "correlation", "abscor", "uncentered", "euclidean" -
#' or cosine (no quotes) as implemented in ARBOL::cosine. you may also write your own function that returns a dist object, as in pvclust::pvclust()
#' @param centroid_method the function used to calculate centroids in the tree_reduction matrix, as implemented in Matrix.utils::aggregate.Matrix(fun)
#' currently, sum, count, mean, and median are supported
#' @param centroid_assay if using cell x gene data (not any srobj@@reduction), the assay within which to calculate centroids
#' @param gene_list if using cell x gene data (not any srobj@@reduction), genes to include in centroid calculation. passes to getCentroids
#' @param reduction_dims the dimensions of the reduction slot to use for centroid calculation. defaults to 1:25
#' @return the input seurat object with pvclust tree in srobj@@misc$pvclust
#' @examples
#' srobj <- centroidTaxonomy(srobj = srobj, tree_reduction = 'pca', hclust_method = 'complete',
#'                          distance_method = 'euclidean', centroid_method = 'median', 
#'                          centroid_assay = 'SCT', reduction_dims = 1:25)
#' @export
centroidTaxonomy <- function(srobj, tree_reduction = 'centroids', hclust_method = 'complete',
                                distance_method = 'euclidean', centroid_method = 'mean', 
                               centroid_assay = 'SCT', reduction_dims = 1:25, gene_list = rownames(srobj[["RNA"]]@data), nboot = 1) {

    if (tree_reduction == 'centroids' | tree_reduction %in% names(srobj@reductions)) {
      centroids <- getCentroids(srobj = srobj, tree_reduction = tree_reduction, reduction_dims = reduction_dims,
                         centroid_method = centroid_method, centroid_assay = centroid_assay, gene_list = gene_list)

      result <- pvclust(centroids, method.dist=distance_method, method.hclust=hclust_method, nboot=nboot)

      srobj@misc$pvclust <- result
    }

    else {
      message(sprintf('%s isnt an ARBOL implemented reduction for tree building, or does not exist in srobj@reductions',tree_reduction))
    }


    return(srobj)
}

#' Calculate pvclust() tree (a binary tree of distances between end-clusters) for ARBOL results
#' tree based on euclidean distance between cluster centroids based on gene medians with complete linkage
#' @param srobj a seurat object with ARBOL 'tierNident' column
#' @param centroid_method function to calculate gene centroids per ARBOL end cluster, see Matrix.utils::aggregate.Matrix for options
#' @param centroid_assay if using cell x gene data (not any srobj@@reduction), the assay within which to calculate centroids
#' @param gene_list if using cell x gene data (not any srobj@@reduction), genes to include in centroid calculation
#' @param tree_reduction cell x gene reduction space to work from. 'centroids' uses full cell x gene.
#' @param reduction_dims the dimensions of the reduction slot to use for centroid calculation. defaults to 1:25
#' @return the input seurat object with pvclust tree in srobj@@misc$pvclust
#' @examples

#' @export
getCentroids <- function(srobj = srobj, tree_reduction = tree_reduction, reduction_dims = reduction_dims,
                         centroid_method = centroid_method, centroid_assay = centroid_assay, gene_list = rownames(srobj[["RNA"]]@data)) {

    srobj_meta <- srobj@meta.data
    if (tree_reduction %in% names(srobj@reductions)) {
      scaled.data.mtx <- Embeddings(srobj, reduction=tree_reduction)[,reduction_dims]
    }
    else {
      gene_list = match(gene_list,rownames(srobj[[centroid_assay]]@data))
      gene_list <- gene_list[!is.na(gene_list)]
      submtx <- srobj[[centroid_assay]]@data[gene_list,]
      scaled.data.mtx <- Matrix(t(as.matrix(submtx)),sparse=TRUE)
    }
        #pull seurat object's cell ID + tierNident
    t2 <- srobj_meta[,c('CellID','tierNident')]
    #create numerical representation of tierNident per ident 
    t3 <- data.frame(tierNident=t2$tierNident %>% unique,i = seq(1:length(unique(t2$tierNident))))
    #place numerical representation in cell matrix t2
    t2$i <- match(t2$tierNident,t3$tierNident) %>% as.double
    #turn numerical representation matrix of tierNident to a sparse matrix
    t4 <- Matrix(t2 %>% select(i) %>% unlist,sparse=TRUE)
    #add numerical representation of tierNident to a sparse matrix of actual data (to allow aggregation)
    scaled.data.t <- cbind(t4,scaled.data.mtx)

    #remember rownames
    rows = row.names(scaled.data.t)

    #aggregate sparse matrix by cluster, taking median of each gene per cluster
    agg.clst.cntrs<- aggregate.Matrix(scaled.data.t[rows,-1],groupings=scaled.data.t[rows,1,drop=FALSE],fun=centroid_method)
    clst.cntrs <- agg.clst.cntrs

    #transpose cluster centers dataframe
    g <- clst.cntrs %>% as.matrix %>% t %>% data.frame

    #remember actual tierNident names
    colnames(g) <- t3$tierNident
    return(g)

}

#' Creates binary tree object for ARBOL run, adds it to seurat object along with a graph of it
#' calls centroidTaxonomy() in which assay + methods for tree building are called
#' Also adds metrics used in building tree to seurat object in case needed for downstream applications
#' @param srobj a seurat object with ARBOL 'tierNident', 'CellID', 'sample' columns. 
#' @param categories categorical attributes are propagated up a tree as unique values per node. also calculates a majority per node
#' @param diversities attributes for which to calculate diversity in each node. Currently only supports gini-simpson's index.
#' @param diversity_metric one of 'shannon', 'simpson', or 'invsimpson'
#' @param counts attributes for which to count unique values per node.
#' @param totals attributes for which to sum values per node.
#' @param tree_reduction either 'centroids', which calculates centroids among all genes, or any reduction slot in srobj
#' @param hclust_method any hierarchical clustering method implemented in pvclust::pvclust(method.hclust), defaults to 'complete'
#' @param distance_method any distance method implemented in pvclust::pvclust(method.dist) - one of "correlation", "abscor", "uncentered", "euclidean" -
#' or cosine (no quotes) as implemented in ARBOL::cosine. you may also write your own function that returns a dist object, as in pvclust::pvclust()
#' @param centroid_method the function used to calculate centroids in the tree_reduction matrix, as implemented in Matrix.utils::aggregate.Matrix(fun)
#' currently, sum, count, mean, and median are supported
#' @param centroid_assay if using cell x gene data (not any srobj@@reduction), the assay within which to calculate centroids
#' @param gene_list if using cell x gene data (not any srobj@@reduction), list of genes to include in centroid calculation. 
#' defaults to all genes in centroid_assay if no input is given
#' @param reduction_dims the dimensions of the reduction slot to use for centroid calculation. defaults to 1:25
#' @param nboot number of bootstraps for pvclust() function. defaults to 1 for speed
#' @return the input seurat object with binary tree pvclust object in srobj@@misc$pvclust, 
#' @examples
#' srobj <- ARBOLcentroidTaxonomy(srobj, categories = c('celltype','disease'))
#' @export
ARBOLcentroidTaxonomy <- function(srobj, categories = 'sample', diversities = 'sample', 
                                diversity_metric = 'simpson', counts = 'sample',
                                totals = 'nCount_RNA', tree_reduction = 'centroids', 
                                hclust_method = 'complete', distance_method = 'euclidean', 
                                centroid_method = 'mean', centroid_assay = 'SCT', 
                                reduction_dims = 1:25, gene_list = rownames(srobj[["RNA"]]@data), 
                                nboot = 1) {


  if (!is.character(gene_list)) {
    gene_list = rownames(srobj[[centroid_assay]]@data)
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
  
  srobj <- centroidTaxonomy(srobj = srobj, tree_reduction = tree_reduction, hclust_method = hclust_method,
                          distance_method = distance_method, centroid_method = centroid_method, 
                          centroid_assay = centroid_assay, gene_list = gene_list, reduction_dims = reduction_dims, nboot = nboot)

  binarydf <- binaryTreeToDF(pvclust_tree=srobj@misc$pvclust)

  divtree <- treeAllotment(srobj, treedf = binarydf, categories = categories, diversities = diversities, 
                                diversity_metric = diversity_metric, counts = counts, totals = totals)

  divtree <- propagateTree(divtree,srobj = srobj, categorical_attributes = categories, 
    diversity_attributes = diversities, numerical_attributes = counts, total_attributes = totals)

  x <- data.tree_to_ggraph(divtree, categories, diversities, counts, totals)

  x <- x %>% activate(nodes) %>% mutate(tier = str_count(label, "\\."))

  x <- x %>% activate(nodes) %>% mutate(string = label, name = basename(label) %>% str_replace_all('T0C0.',''))

  bt1 <- ggraph(x, layout = 'dendrogram', height=plotHeight*10) +
    geom_edge_elbow() + 
    geom_node_point(size=0) + 
    geom_node_text(aes(filter = leaf, label = name), nudge_y=-0.75,vjust=0.5,hjust=0,angle=270) + 
    geom_node_text(aes(filter = leaf, label = n),color='grey30',nudge_y=-0.2,vjust=0.5,hjust=0,size=3,angle=270) +
    theme_void() +
    geom_node_point(aes(filter = leaf,color=sample_diversity),size=4,shape='square') + 
    scale_color_gradient(low='grey90',high='grey10',limits=c(0,1)) +
    expand_limits(y=-5)

  srobj@misc$taxViz <- bt1

  srobj@misc$taxTree <- divtree
  
  srobj@misc$tax_ggraph <- x

  #add metrics used for tree creation to seurat object
  srobj@misc$tree_metrics[['diversities']] <- diversities
  srobj@misc$tree_metrics[['diversity_metric']] <- diversity_metric
  srobj@misc$tree_metrics[['categories']] <- categories
  srobj@misc$tree_metrics[['counts']] <- counts
  srobj@misc$tree_metrics[['tree_reduction']] <- tree_reduction
  srobj@misc$tree_metrics[['hclust_method']] <- hclust_method
  srobj@misc$tree_metrics[['distance_method']] <- distance_method
  srobj@misc$tree_metrics[['centroid_method']] <- centroid_method
  srobj@misc$tree_metrics[['centroid_assay']] <- centroid_assay
  srobj@misc$tree_metrics[['reduction_dims']] <- reduction_dims
  srobj@misc$tree_metrics[['gene_list']] <- gene_list

  return(srobj)
}


#' custom cosine distance method for pvclust
cosine <- function(x) {
         x <- as.matrix(x)
         y <- t(x) %*% x
         res <- 1 - y / (sqrt(diag(y)) %*% t(sqrt(diag(y))))
         res <- as.dist(res)
         attr(res, "method") <- "cosine"
         return(res)
     }

#' data.tree to ggraph conversion, by converting data.tree structure to ape via newick (doesn't carry annotations)
#' then ape (still doesn't carry annotations) to ggraph, 
#' then join data frame of data.tree to restore annotations.
#' requires 'n' and 'pathString' annotations in data.tree
#' Used to convert annotated binary phylogeny tree to ggraph for easier plotting 
#' 1) write data.tree object to Newick using custom Newick function,
#' 2) read Newick into ape tree object
#' 3) tidygraph::as_tbl_graph to convert from ape to ggraph
#' 4) data.tree to node-level dataframe 
#' 5) join node-level dataframe to tbl_graph nodes
#' The use of ToNewick is a major bottleneck in computation speed right now
data.tree_to_ggraphNW <- function(data.tree, categories, diversities, counts) {
  txt <- ToNewickPS(data.tree)
  apeTree <- ape::read.tree(text=txt)
  attrs <- data.tree$attributesAll
  count_cols <- attrs[grep(sprintf('^(%s)_n',paste0(counts,collapse='|')),attrs)]
  treeDF <- do.call(allotedTreeToDF, c(data.tree, 'n','pathString','plotHeight', 'pval',categories, paste0(categories,'_majority'),paste0(diversities,'_diversity'),count_cols))
  treeDF <- treeDF %>% select(name=pathString,n,i,plotHeight,pval,all_of(categories),all_of(paste0(categories,'_majority')),all_of(paste0(diversities,'_diversity')),all_of(count_cols))
  x <- tidygraph::as_tbl_graph(apeTree,directed=T) %>% activate(nodes) %>% left_join(treeDF)
  x <- x %>% activate(edges) %>% left_join(treeDF %>% select(to=i,height=plotHeight))
  return(x)
}

#' data.tree to ggraph conversion, by converting data.tree structure to dendrogram via custom PS function
#' then ape (still doesn't carry annotations) to ggraph,
#' then join data frame of data.tree to restore annotations.
#' requires 'n' and 'pathString' annotations in data.tree
#' Used to convert annotated binary phylogeny tree to ggraph for easier plotting 
#' 1) write data.tree object to dend object using custom as.dendrogram.Node function,
#' 2) convert dendrogram structure object to ggraph using tidygraph::as_tbl_graph
#' 3) left join heights, categories, diversities, and counts to tbl_graph object
#' previous version data.tree_to_ggraphNW used highly inefficient newick text conversion for tree structure
data.tree_to_ggraph <- function(data.tree, categories, diversities, counts, totals, heightAttribute = 'plotHeight') {
  datadend <- as.dendrogram.NodePS(data.tree, heightAttribute = heightAttribute)
  attrs <- data.tree$attributesAll
  count_cols <- attrs[grep(sprintf('^(%s)_n',paste0(counts,collapse='|')),attrs)]
  treeDF <- do.call(allotedTreeToDF, c(data.tree, 'n','pathString','plotHeight', 'pval',categories, paste0(categories,'_majority'),paste0(diversities,'_diversity'),count_cols,totals))
  treeDF <- treeDF %>% select(label=pathString,n,i,plotHeight,pval,all_of(categories),all_of(paste0(categories,'_majority')),all_of(paste0(diversities,'_diversity')),all_of(count_cols),all_of(totals))
  x <- tidygraph::as_tbl_graph(datadend,directed=T) %>% activate(nodes) %>% left_join(treeDF)
  x <- x %>% activate(edges) %>% left_join(treeDF %>% select(to=i,height=plotHeight))
  return(x)
}

#' Allots seurat object metadata attributes to a data.tree object
#' Once attributes are alloted, propagateTree should be run to propagate attributes up the tree.
#' @param srobj a seurat object with ARBOL 'tierNident', 'CellID', 'sample' columns. 
#' @param categories categorical attributes are propagated up a tree as unique values per node. also calculates a majority per node
#' @param diversities attributes for which to calculate diversity in each node. Currently only supports gini-simpson's index.
#' @param diversity_metric one of 'shannon', 'simpson', or 'invsimpson'
#' @param counts attributes for which to count unique values per node.
#' @param totals attributes for which to sum values per node.
#' @return tree with attributes alloted to internal nodes
#' @export 
treeAllotment <- function(srobj, treedf, categories, diversities, diversity_metric, counts, totals) {
  srobj <- tierN_diversity(srobj, diversity_attributes = diversities, diversity_metric = diversity_metric)

  jointb <- srobj@meta.data %>% group_by(tierNident) %>% mutate(n=n()) %>% 
          dplyr::select(CellID,sample,tierNident,n,all_of(categories),all_of(paste0(diversities,'_diversity')),all_of(totals))

  categorydf <- jointb %>% summarize(across(categories, ~ list(strsplit(paste(unique(.x),collapse=','),','))))

  divdf <- jointb %>% summarize_at(paste0(diversities,'_diversity'),unique)

  jointb <- jointb %>% select(-all_of(categories)) %>% 
            summarize(ids = list(CellID),n=unique(n))

  totalsdf <- srobj@meta.data %>% dplyr::select(CellID,tierNident,all_of(totals))

  for (x in totals) {
    attribute <- enquo(x)
    totalstb <- totalsdf %>% dplyr::select(tierNident,x) %>% group_by(tierNident) %>%
        summarize(tierNident=unique(tierNident),sum(.data[[attribute]]))
    colnames(totalstb) <- c('tierNident',x)
    totalsdf <- totalsdf %>% dplyr::select(-x) %>% left_join(totalstb)
  }

  totalsdf <- totalsdf %>% select(-CellID) %>% distinct

  numericaltb <- srobj@meta.data %>% dplyr::select(CellID,sample,tierNident,all_of(counts))

  for (x in counts) {
  attribute <- enquo(x)
  tierNcount <- numericaltb %>% group_by(tierNident) %>% dplyr::count(.data[[attribute]])
  vals <- unique(tierNcount[[x]])
  countWide <- tierNcount %>% pivot_wider(names_from=all_of(x),values_from=n)
  colnames(countWide) <- c('tierNident',sprintf('%s_n_%s',x,vals))
  numericaltb <- numericaltb %>% left_join(countWide)
  }

  numericaltb[is.na(numericaltb)] <- 0

  numericaltb <- numericaltb %>% select(-CellID,-sample,-all_of(counts)) %>% distinct

  finaltreedf <- treedf %>% left_join(jointb) %>% left_join(categorydf) %>% left_join(divdf) %>% left_join(numericaltb) %>% left_join(totalsdf)

  tree <- as.Node(finaltreedf)

  return(tree)

}

#' Merges tierNidents with their nearest neighbors in a binary tree if 
#' their sample diversity and number of cells do not meet specified thresholds.
#' Directly prunes srobj-attached binary tree. 
#' We suggest re-calculating tree (or just ggraph for viz) from here using centroidTaxonomy() or remake_ggraph(), 
#' so that new endclusters are treated as leaf nodes
#' @param srobj a seurat object with a binarytree calculated in slot srobj@@misc$taxTree, 
#' typically calculated using ARBOLcentroidTaxonomy
#' @param sample_diversity_threshold sample diversity below which to prune nodes from tree
#' @param size_threshold cluster size below which to prune nodes from tree
#' @return the input seurat object with merged tierNidents in a new metadata column, mergedIdent
#' @examples
#' srobj <- mergeEndclusts(srobj, sample_diversity_threshold = 0.1, size_threshold = 10)
#' @export
mergeEndclusts <- function(srobj, sample_diversity_threshold, size_threshold) {

  srobj@misc$rawTaxTree <- Clone(srobj@misc$taxTree)
  #DataTree::Prune chops all nodes that don't meet a threshold
  Prune(srobj@misc$taxTree, pruneFun = function(x) x$sample_diversity > sample_diversity_threshold)
  Prune(srobj@misc$taxTree, pruneFun = function(x) x$n > size_threshold)

  #remove unnecessary nodes that have only 1 child - these are created in binary tree threshold merging
  Prune(srobj@misc$taxTree, pruneFun = function(x) any(x$children %>% length > 1 || x$children %>% length == 0))

  divtestdf <- allotedTreeToDF(srobj@misc$taxTree, 'height', "pathString", 'ids')
  divdf2 <- divtestdf %>% mutate(ids = strsplit(ids, ", ")) %>% unnest(cols = c(ids))
  divdf2 <- divdf2 %>% group_by(ids) %>% slice_min(height) %>% ungroup

  divdf2 <- divdf2 %>% dplyr::rename(CellID=ids,binIdent = pathString)
  divdf2$binIdent <- divdf2$binIdent %>% str_replace_all('\\/','.')

  srobj@meta.data <- srobj@meta.data %>% left_join(divdf2 %>% select(CellID,mergedIdent=binIdent)) %>%
                     rename(tierNident=mergedIdent,rawIdent=tierNident)

  return(srobj)
}

#' Merge based on sample diversity and endclust size, only outputting new end-clusters per cell
#' @param srobj a seurat object with a binarytree calculated in slot srobj@@misc$taxTree, typically calculated using ARBOLcentroidTaxonomy
#' @param sample_diversity_threshold sample diversity below which to prune nodes from tree
#' @param size_threshold cluster size below which to prune nodes from tree
#' @return dataframe with 2 columns, cellid and new endcluster
#' @examples
#' srobj <- mergeEndclustsIdents(srobj, sample_diversity_threshold = 0.1, size_threshold = 10)
#' @export
mergeEndclustsIdents <- function(srobj, sample_diversity_threshold, size_threshold) {

  workingTree <- Clone(srobj@misc$taxTree)
  #DataTree::Prune chops all nodes that don't meet a threshold
  Prune(workingTree, pruneFun = function(x) x$sample_diversity > sample_diversity_threshold)
  Prune(workingTree, pruneFun = function(x) x$n > size_threshold)

  #remove unnecessary nodes that have only 1 child - these are created in binary tree threshold merging
  Prune(workingTree, pruneFun = function(x) any(x$children %>% length > 1 || x$children %>% length == 0))

  divtestdf <- allotedTreeToDF(workingTree, 'height', "pathString", 'ids')
  divdf2 <- divtestdf %>% mutate(ids = strsplit(ids, ", ")) %>% unnest(cols = c(ids))
  divdf2 <- divdf2 %>% group_by(ids) %>% slice_min(height) %>% ungroup

  divdf2 <- divdf2 %>% dplyr::rename(CellID=ids,mIdent = pathString)
  divdf2$mIdent <- divdf2$mIdent %>% str_replace_all('\\/','.')
  divdf3 <- divdf2 %>% dplyr::select(mIdent,CellID)
        
  return(divdf3)
}

#' Merges tierNidents with their nearest neighbors in a srobj-attached binary tree by custom thresholds
#' We suggest re-calculating tree (or just ggraph for viz) from here using centroidTaxonomy() or remake_ggraph(), 
#' so that new endclusters are treated as leaf nodes
#' @param srobj a seurat object with a binarytree calculated in slot srobj@@misc$taxTree, typically calculated using ARBOLcentroidTaxonomy
#' @param threshold_attributes list of srobj metadata columns to threshold on
#' @param thresholds list of threshold values to prune, in same order as threshold_attributes
#' @return the input seurat object with merged tierNidents in a new metadata column, mergedIdent 
#' and a merged data.tree object in srobj@@misc$workingTree
#' @examples
#' srobj <- mergeEndclustsCustom(srobj, threshold_attributes = c('sample_diversity','n'), thresholds = c(0.2,50))
#' @export
mergeEndclustsCustom <- function(srobj, threshold_attributes, thresholds) {

 workingTree <- Clone(srobj@misc$taxTree)
  #DataTree::Prune chops all nodes that don't meet a threshold
  for (z in seq(1,length(threshold_attributes))) {
    Prune(workingTree, pruneFun = function(x) x[[threshold_attributes[z]]] > thresholds[z])
  } 

  #remove unnecessary nodes that have only 1 child - these are created in binary tree threshold merging
  Prune(workingTree, pruneFun = function(x) any(x$children %>% length > 1 || x$children %>% length == 0))

  divtestdf <- allotedTreeToDF(workingTree, 'height', "pathString", 'ids')
  divdf2 <- divtestdf %>% mutate(ids = strsplit(ids, ", ")) %>% unnest(cols = c(ids))
  divdf2 <- divdf2 %>% group_by(ids) %>% slice_min(height) %>% ungroup

  divdf2 <- divdf2 %>% dplyr::rename(CellID=ids,mIdent = pathString)
  divdf2$mIdent <- divdf2$mIdent %>% str_replace_all('\\/','.')
  divdf3 <- divdf2 %>% dplyr::select(mIdent,CellID)

  srobj@meta.data <- srobj@meta.data %>% left_join(divdf3)
  srobj@meta.data$rawIdent <- srobj@meta.data$tierNident
  srobj@meta.data$tierNident <- srobj@meta.data$mIdent
  srobj@meta.data <- srobj@meta.data %>% select(-mIdent)      

  srobj@misc$workingTree <- workingTree

  return(srobj)
}

#' Merges tierNidents with their nearest neighbors in a srobj-attached binary tree by custom thresholds, 
#' outputs dataframe of new idents
#' @param srobj a seurat object with a binarytree calculated in slot srobj@@misc$taxTree, typically calculated using ARBOLcentroidTaxonomy
#' @param threshold_attributes list of srobj metadata columns to threshold on
#' @param thresholds list of threshold values to prune, in same order as threshold_attributes
#' @return dataframe of new idents
#' @examples
#' mIdents <- mergeEndclustsCustomIdents(srobj, threshold_attributes = c('sample_diversity','n'), thresholds = c(0.2,50))
#' @export
mergeEndclustsCustomIdents <- function(srobj, threshold_attributes, thresholds) {

  workingTree <- Clone(srobj@misc$taxTree)
  #DataTree::Prune chops all nodes that don't meet a threshold
  for (z in seq(1,length(threshold_attributes))) {
    Prune(workingTree, pruneFun = function(x) x[[threshold_attributes[z]]] > thresholds[z])
  } 

  #remove unnecessary nodes that have only 1 child - these are created in binary tree threshold merging
  Prune(workingTree, pruneFun = function(x) any(x$children %>% length > 1 || x$children %>% length == 0))

  divtestdf <- allotedTreeToDF(workingTree, 'height', "pathString", 'ids')
  divdf2 <- divtestdf %>% mutate(ids = strsplit(ids, ", ")) %>% unnest(cols = c(ids))
  divdf2 <- divdf2 %>% group_by(ids) %>% slice_min(height) %>% ungroup

  divdf2 <- divdf2 %>% dplyr::rename(CellID=ids,mIdent = pathString)
  divdf2$mIdent <- divdf2$mIdent %>% str_replace_all('\\/','.')
  divdf3 <- divdf2 %>% dplyr::select(mIdent,CellID)

  return(divdf3)
}

#' Converts pvclust tree object to dataframe with row per node and pvclust p-values
#' @param pvclust_tree a pvclust or hclust tree 
#' @return a dataframe with one row per node of the tree
#' @examples
#' binarydf <- binaryTreeToDF(pvclust_tree=srobj@@misc$pvclust)
#' @export
binaryTreeToDF <- function(pvclust_tree) {
  bintree <- as.Node(as.dendrogram(pvclust_tree$hclust))
  binarydf <- data.frame(ToDataFrameTree(bintree, 'pathString','plotHeight','isLeaf'))
  bin2 <- binarydf %>% filter(!isLeaf) %>% mutate(pval = 1-pvclust_tree$edges$au)
  binarydf <- binarydf %>% left_join(bin2)
  colnames(binarydf) <- c('remove','pathString','plotHeight','remove2','pval')
  #following line is heuristic and causes wrapper functions to crash if tierNident has / in it. Can cause problems with curatednames!
  #need a warning line in wrapper functions or something
  binarydf$tierNident <- sub('.*\\/', '', binarydf$pathString)
  binarydf <- binarydf %>% dplyr::select(-remove,-remove2)
  return(binarydf)
}

#' Converts data.tree object to dataframe with row per node, this function aids conversion of highly annotated data.tree objects to dataframes.
#' Used in ARBOL workflow to create a tbl_graph object, a class useful for plotting dendrograms
#' @param tree usually a custom annotated data.tree object
#' @return a dataframe of the tree with one row per node
#' @examples
#' ARBOLdf <- allotedTreeToDF(divtree)
#' txt <- ARBOL::ToNewickPS(divtree)
#' ARBOLphylo <- ape::read.tree(text=txt)
#' 
#' x <- as_tbl_graph(ARBOLphylo, directed=T) %>% activate(nodes) %>% 
#' left_join(ARBOLdf %>% select(name=pathString,sample_diversity,disease_diversity))
#' @export
allotedTreeToDF <- function(tree, ...) {
    ARBOLdf <- ToDataFrameTree(tree, ...)

    #remove graphics from levelName column
    ARBOLdf$levelName <- ARBOLdf$levelName %>% str_replace_all(' ','') %>% str_replace_all('-','') %>%
                            str_replace_all('¦','') %>% str_replace_all('°','')

    #add a row index column for joining with ggraph object                        
    ARBOLdf$i <- as.numeric(row.names(ARBOLdf))

    return(ARBOLdf)
}

#' Calculate Diversity Index for tierNidents for multiple attributes of the data
#' diversity metric options same as those available in vegan's diversity() function
#' @param srobj a seurat object with ARBOL tierNident column
#' @param diversity_attributes the attributes you wish to calculate diversity for (e.g. disease, sample)
#' @param diversity_metric one of 'shannon', 'simpson', or 'invsimpson'
#' @return the seurat object with diversity per attribute per tierNident added to the metadata
#' @examples
#' srobj <- tierN_diversity(srobj, diversity_attributes = c('sample','disease'))
#' @export
tierN_diversity <- function(srobj, diversity_attributes, diversity_metric = 'simpson') {
  #remove existing diversity metrics
  srobj@meta.data <- srobj@meta.data %>% dplyr::select(-any_of(paste0(diversity_attributes, '_diversity')))
    for(z in diversity_attributes) {
        #calculate and join new diversity
        srobj@meta.data <- srobj@meta.data %>% left_join(diversityPerGroup(srobj@meta.data, species=tierNident, group=z, diversity_metric = diversity_metric)) %>% suppressMessages
    }
    return(srobj)
}

#' spread tierNident column in a dataframe into separate columns per tier, with full paths and single tier idents
#' @param df a dataframe with ARBOL tierNident (usually srobj@@meta.data)
#' @param max_tiers max_tiers parameter used in ARBOL
#' @return the dataframe with separate tier columns, useful for producing and plotting dendrograms
#' @examples
#' srobj@@meta.data <- spread_tierN(srobj@@meta.data)
#' @export
spread_tierN <- function(df, max_tiers = 10) {
    df <- df %>% separate(tierNident,into=paste0('tier',seq(1,max_tiers)),remove=F)
    df$tier0 <- 'T0C0'

    #add delimiter to end of tierNident to allow programmatic parsing!
    df$tierNident <- paste0(df$tierNident,'.')

    for(x in seq(1,max_tiers)) {
        df <- df %>% mutate(tierfull = strex::str_before_nth(tierNident, '\\.', n=x))
        # following line was left-over code causing problems in alternative ARBOL tree building. 
        # keeping it commented here for now in case something breaks in subclusteringTree() function, can try adding it back
        #df$tierfull[!grepl(paste0('T',x),df$tierfull)] <- NA
        names(df)[names(df) == 'tierfull'] <- paste0('tier',x,'full')
    }

    #remove delimiter
    df$tierNident = substr(df$tierNident,1,nchar(df$tierNident)-1)

    df2 <- df %>% unite('pathString', tier0:!!paste0('tier',max_tiers,'full'), sep = "/", na.rm = TRUE, remove = FALSE)

    return(df2)
}

#' Calculate curated names for set of clusters as function of their major type
#' @param srobj a seurat object with a celltype metadata column, specified in celltype_col, 
#' CellID (for preserving metadata rownames), and tierNident (i.e. srobj@@meta.data$CellID, srobj@@meta.data$tierNident)
#' @param figdir the figdir from an GenTieredClusters (ARBOL) run 
#' @param max_cells_per_ident maximum cells to use in FindAllMarkers call. defaults to 200
#' @param celltype_col the metadata column corresponding to celltype to assign standard names
#' @param standardname_col metadata column to which to output celltype names. defaults to 'curatedname'
#' @param n_genes number of genes with which to create standard names. defaults to 2
#' @return the seurat object with curatedname column in metadata
#' @examples
#' srobj <- getStandardNames(srobj,figdir='ARBOLoutput/figs')
#' @export
getStandardNames <- function(srobj, figdir, max_cells_per_ident=200, celltype_col = 'celltype', 
                              standardname_col = 'curatedname', n_genes = 2) {

  if (!is.element('CellID',colnames(srobj@meta.data))) {
    message('you are missing the CellID column - function will fail with confusing error message')
  }
  if (!is.element('tierNident',colnames(srobj@meta.data))) {
    message('you are missing the tierNident column - function will fail with confusing error message')
  }

  Idents(srobj) <- srobj@meta.data[[celltype_col]]
  typeobjs <- SplitSrobjOnMeta(srobj, meta=celltype_col,'typeobjects')

  if (!file.exists(sprintf('%s/EndClustersAmongTier1DE.csv',figdir))) {
    cellStateMarkers <- lapply(typeobjs,function(obj) {Idents(obj) <- obj@meta.data$tierNident;
      #if only one cell state, calculate markers against all other cells
          if(length(unique(Idents(obj)))==1) {
            ident = unique(obj@meta.data$tierNident)
            type = unique(obj@meta.data[[celltype_col]])
            tmpobj <- srobj
            Idents(tmpobj) <- ifelse(tmpobj@meta.data$tierNident==ident,ident,'other')
            message('found celltype with only one cluster. Calculating markers by wilcoxon test against all other cells')
            tmp <- FindAllMarkers(tmpobj,only.pos=TRUE,min.pct = 0.25,logfc.threshold = 0.25) %>% 
            #remove markers for all the other cells from the table
              filter(cluster!='other')
            tmp[[celltype_col]] <- type
          }
          else {
            tmp <- FindAllMarkers(obj,only.pos=TRUE,min.pct = 0.25,logfc.threshold = 0.25, max.cells.per.ident = max_cells_per_ident);
            tmp[[celltype_col]] <- unique(obj@meta.data[[celltype_col]])
          }
           return(tmp)})
    write.table(rbindlist(cellStateMarkers), sprintf('%s/EndClustersAmongTier1DE.csv',figdir), sep=",", row.names=F)
  }
  else {
    cellStateMarkers <- fread(sprintf('%s/EndClustersAmongTier1DE.csv',figdir)) %>%
                    split(f=.$cluster)
  }
  
   markersL <- lapply(cellStateMarkers,function(cellStateDF) {
    tryCatch({
           fccol <- grep("FC",colnames(cellStateDF))
            biomarkers <- cellStateDF %>%
          mutate(rnkscr = -log(p_val_adj+1e-310) * sapply(cellStateDF[,!!fccol], as.numeric) * (pct.1 / (pct.2 + 1e-300))) %>%
                    group_by(across(celltype_col),cluster) %>% dplyr::slice_max(n=n_genes, order_by=rnkscr, with_ties=F) %>% 
                    arrange(cluster, desc(rnkscr));

              endname <- biomarkers %>% dplyr::mutate(markers = paste(gene, collapse=".")) %>% 
                    dplyr::select(-p_val,-avg_log2FC,-pct.1,-pct.2,-p_val_adj,-rnkscr)
        
              endname <- endname %>% dplyr::select(-gene) %>% unite({{standardname_col}},-cluster,sep='.',remove=F)        
                
            endname <- endname %>% distinct
        
              return(endname)
                 },
           error = function(e) { message(sprintf('standard name calculation failed for a %s cluster... reason: %s',unique(cellStateDF[[celltype_col]])[1],e)) 
           })
               })

  suppressWarnings( bind_rows(markersL) -> markersAsList)

  message(sprintf('number of tierNident-celltype pairs for which at least %s biomarkers were found: ',n_genes),
       length(markersAsList$cluster))
  message('total number of end clusters: ',length(unique(srobj@meta.data$tierNident)))

  markersAsList <- markersAsList %>% dplyr::rename(tierNident=cluster)

  srobj@meta.data <- left_join(srobj@meta.data,markersAsList,by=c("tierNident",celltype_col))
  
  #the following lines will cause the function to give a standard name only for the majority celltype_col per cellstate
  majority <- srobj@meta.data %>% dplyr::count(tierNident,across(celltype_col),markers) %>% group_by(tierNident) %>% 
                             slice_max(n) %>% dplyr::select(-n)

  colnames(majority) <-  c('tierNident',paste0('majority_',celltype_col),'markers')
  majority <- majority %>% unite({{standardname_col}},-tierNident,sep='.')

  #remove duplicate lines per cellID
  srobj@meta.data <- srobj@meta.data %>% dplyr::select(-markers) %>% distinct

  srobj@meta.data <- srobj@meta.data %>% left_join(majority)

  #Hoping to include an if-statement around the previous block to enable users to 
  #allow or disallow multiple standard names per cellstate

  #replace rownames of seurat object metadata
  row.names(srobj@meta.data) <- srobj@meta.data$CellID

  return(srobj)
}

#' Calculate curated names for set of clusters as function of their major type
#' @param srobj a seurat object with a celltype metadata column, specified in celltype_col, 
#' CellID (for preserving metadata rownames), and tierNident (i.e. srobj@@meta.data$CellID, srobj@@meta.data$tierNident)
#' @param figdir the figdir from an ARBOL run 
#' @param max_cells_per_ident maximum cells to use in FindAllMarkers call. defaults to 200
#' @param celltype_col the metadata column corresponding to celltype to assign standard names
#' @param standardname_col metadata column to which to output celltype names. defaults to 'curatedname'
#' @param n_genes number of genes with which to create standard names. defaults to 2
#' @return a dataframe with standard names per end-cluster
#' @examples
#' srobj <- outputStandardNames(srobj,figdir='ARBOLoutput/figs')
#' @export
outputStandardNames <- function(srobj, figdir, max_cells_per_ident=200, celltype_col = 'celltype', 
                              standardname_col = 'curatedname', n_genes = 2) {

  if (!is.element('CellID',colnames(srobj@meta.data))) {
    message('you are missing the CellID column - function will fail with confusing error message')
  }
  if (!is.element('tierNident',colnames(srobj@meta.data))) {
    message('you are missing the tierNident column - function will fail with confusing error message')
  }

  Idents(srobj) <- srobj@meta.data[[celltype_col]]
  typeobjs <- SplitSrobjOnMeta(srobj, meta=celltype_col,'typeobjects')

  if (!file.exists(sprintf('%s/EndClustersAmongTier1DE.csv',figdir))) {
    cellStateMarkers <- lapply(typeobjs,function(obj) {Idents(obj) <- obj@meta.data$tierNident;
      #if only one cell state, calculate markers against all other cells
          if(length(unique(Idents(obj)))==1) {
            ident = unique(obj@meta.data$tierNident)
            type = unique(obj@meta.data[[celltype_col]])
            tmpobj <- srobj
            Idents(tmpobj) <- ifelse(tmpobj@meta.data$tierNident==ident,ident,'other')
            message('found celltype with only one cluster. Calculating markers by wilcoxon test against all other cells')
            tmp <- FindAllMarkers(tmpobj,only.pos=TRUE,min.pct = 0.25,logfc.threshold = 0.25) %>% 
            #remove markers for all the other cells from the table
              filter(cluster!='other')
            tmp[[celltype_col]] <- type
          }
          else {
            tmp <- FindAllMarkers(obj,only.pos=TRUE,min.pct = 0.25,logfc.threshold = 0.25, max.cells.per.ident = max_cells_per_ident);
            tmp[[celltype_col]] <- unique(obj@meta.data[[celltype_col]])
          }
           return(tmp)})
    write.table(rbindlist(cellStateMarkers), sprintf('%s/EndClustersAmongTier1DE.csv',figdir), sep=",", row.names=F)
  }
  else {
    cellStateMarkers <- fread(sprintf('%s/EndClustersAmongTier1DE.csv',figdir)) %>%
                    split(f=.$cluster)
  }
  
   markersL <- lapply(cellStateMarkers,function(cellStateDF) {
    tryCatch({
           fccol <- grep("FC",colnames(cellStateDF))
            biomarkers <- cellStateDF %>%
          mutate(rnkscr = -log(p_val_adj+1e-310) * sapply(cellStateDF[,!!fccol], as.numeric) * (pct.1 / (pct.2 + 1e-300))) %>%
                    group_by(across(celltype_col),cluster) %>% dplyr::slice_max(n=n_genes, order_by=rnkscr, with_ties=F) %>% 
                    arrange(cluster, desc(rnkscr));

              endname <- biomarkers %>% dplyr::mutate(markers = paste(gene, collapse=".")) %>% 
                    dplyr::select(-p_val,-avg_log2FC,-pct.1,-pct.2,-p_val_adj,-rnkscr)
        
              endname <- endname %>% dplyr::select(-gene) %>% unite({{standardname_col}},-cluster,sep='.',remove=F)        
                
            endname <- endname %>% distinct
        
              return(endname)
                 },
           error = function(e) { message(sprintf('standard name calculation failed for a %s cluster... reason: %s',unique(cellStateDF[[celltype_col]])[1],e)) 
           })
               })

  suppressWarnings( bind_rows(markersL) -> markersAsList)

  message(sprintf('number of tierNident-celltype pairs for which at least %s biomarkers were found: ',n_genes),
       length(markersAsList$cluster))
  message('total number of end clusters: ',length(unique(srobj@meta.data$tierNident)))

  markersAsList <- markersAsList %>% dplyr::rename(tierNident=cluster)

  return(markersAsList)
}

#' remake ggraph object with new categories and diversities, useful to add curatednames
#' @param srobj a seurat object with tierNident, sample, and category columns in metadata i.e. srobj@@meta.data$tierNident
#' @param categories categories as in tree building functions
#' @param diversities diversities as in tree building functions
#' @param diversity_metric one of 'shannon', 'simpson', or 'invsimpson'
#' @param counts attributes for which to count unique values per node.
#' @return the seurat object with ggraph object remade in srobj@@misc$tax_ggraph
#' @examples
#' srobj <- remake_ggraph(srobj, categories = c('curatedname','celltype'), diversities = c('curatedname','celltype'))
#' @export
remake_ggraph <- function(srobj, categories, diversities, counts, diversity_metric = 'simpson') {

  if (!is.element('sample',diversities)) {
    diversities = c('sample',diversities)
  }

  if (!is.element('sample',categories)) {
    categories = c('sample',categories)
  } 

  if (!is.element('sample',counts)) {
    counts = c('sample',counts)
  } 

  binarydf <- binaryTreeToDF(pvclust_tree=srobj@misc$pvclust)

  srobj <- tierN_diversity(srobj, diversity_attributes = diversities, diversity_metric = diversity_metric)

  jointb <- srobj@meta.data %>% group_by(tierNident) %>% mutate(n=n()) %>% 
            dplyr::select(CellID,sample,tierNident,n,all_of(categories),all_of(paste0(diversities,'_diversity')))

  categorydf <- jointb %>% summarize(across(categories, ~ list(strsplit(paste(unique(.x),collapse=','),','))))

  divdf <- jointb %>% summarize_at(paste0(diversities,'_diversity'),unique)

  jointb <- jointb %>% select(-all_of(categories)) %>% 
              summarize(ids = list(CellID),n=unique(n))

  numericaltb <- srobj@meta.data %>% dplyr::select(CellID,sample,tierNident,all_of(counts))

  for (x in counts) {
    attribute <- enquo(x)
    tierNcount <- numericaltb %>% group_by(tierNident) %>% dplyr::count(.data[[attribute]])
    vals <- unique(tierNcount[[x]])
    countWide <- tierNcount %>% pivot_wider(names_from=all_of(x),values_from=n)
    colnames(countWide) <- c('tierNident',sprintf('%s_n_%s',x,vals))
    numericaltb <- numericaltb %>% left_join(countWide)
  }

  numericaltb[is.na(numericaltb)] <- 0   

  numericaltb <- numericaltb %>% select(-CellID,-sample,-all_of(counts)) %>% distinct           

  finaltreedf <- binarydf %>% left_join(jointb) %>% left_join(categorydf) %>% left_join(divdf) %>% left_join(numericaltb)

  divtree <- as.Node(finaltreedf)

  divtree <- propagateTree(divtree,srobj=srobj, categorical_attributes = categories, 
    diversity_attributes = diversities, numerical_attributes = counts)

  srobj@misc$tax_ggraph <- data.tree_to_ggraph(divtree, categories, diversities, counts)

  return(srobj)

}





###Sparse matrix aggregation functions pulled from the (now defunct) Matrix.utils package
#' Data.frame-Like Operations on Sparse and Dense \code{Matrix} Objects
#' 
#' Implements data manipulation methods such as cast, aggregate, and merge/join
#' for \code{\link{Matrix}} and matrix-like objects.
#' 
#' @name Matrix.utils
#' @docType package
#' @import Matrix
#' @import grr
#' @importFrom stats aggregate as.formula terms contrasts na.omit na.pass
#' @importFrom methods is as
NULL

#' Casts or pivots a long \code{data frame} into a wide sparse matrix
#' 
#' Similar in function to \code{\link[reshape2]{dcast}}, but produces a sparse 
#' \code{\link{Matrix}} as an output. Sparse matrices are beneficial for this 
#' application because such outputs are often very wide and sparse. Conceptually
#' similar to a \code{pivot} operation.
#' 
#' Casting formulas are slightly different than those in \code{dcast} and follow
#' the conventions of \code{\link{model.matrix}}. See \code{\link{formula}} for 
#' details.  Briefly, the left hand side of the \code{~} will be used as the 
#' grouping criteria.  This can either be a single variable, or a group of 
#' variables linked using \code{:}.  The right hand side specifies what the 
#' columns will be. Unlike \code{dcast}, using the \code{+} operator will append
#' the values for each variable as additional columns.  This is useful for 
#' things such as one-hot encoding.  Using \code{:} will combine the columns as 
#' interactions.
#' 
#' @param data a data frame
#' @param formula casting \code{\link[stats]{formula}}, see details for specifics.
#' @param fun.aggregate name of aggregation function.  Defaults to 'sum'
#' @param value.var name of column that stores values to be aggregated numerics
#' @param as.factors if TRUE, treat all columns as factors, including
#' @param factor.nas if TRUE, treat factors with NAs as new levels.  Otherwise, 
#'  rows with NAs will receive zeroes in all columns for that factor
#' @param drop.unused.levels should factors have unused levels dropped? Defaults to TRUE, 
#'  in contrast to \code{\link{model.matrix}}
#' @return a sparse \code{Matrix}
#' @seealso \code{\link[reshape]{cast}}
#' @seealso \code{\link[reshape2]{dcast}}
#' @export
#' @examples
#' #Classic air quality example
#' melt<-function(data,idColumns)
#' {
#'   cols<-setdiff(colnames(data),idColumns)
#'   results<-lapply(cols,function (x) cbind(data[,idColumns],variable=x,value=as.numeric(data[,x])))
#'   results<-Reduce(rbind,results)
#' }
#' names(airquality) <- tolower(names(airquality))
#' aqm <- melt(airquality, idColumns=c("month", "day"))
#' dMcast(aqm, month:day ~variable,fun.aggregate = 'mean',value.var='value')
#' dMcast(aqm, month ~ variable, fun.aggregate = 'mean',value.var='value') 
#' 
#' #One hot encoding
#' #Preserving numerics
#' dMcast(warpbreaks,~.)
#' #Pivoting numerics as well
#' dMcast(warpbreaks,~.,as.factors=TRUE)
#' 
#' \dontrun{
#' orders<-data.frame(orderNum=as.factor(sample(1e6, 1e7, TRUE)), 
#'    sku=as.factor(sample(1e3, 1e7, TRUE)), 
#'    customer=as.factor(sample(1e4,1e7,TRUE)), 
#'    state = sample(letters, 1e7, TRUE),
#'    amount=runif(1e7)) 
#' # For simple aggregations resulting in small tables, dcast.data.table (and
#'    reshape2) will be faster
#' system.time(a<-dcast.data.table(as.data.table(orders),sku~state,sum,
#'    value.var = 'amount')) # .5 seconds 
#' system.time(b<-reshape2::dcast(orders,sku~state,sum,
#'    value.var = 'amount')) # 2.61 seconds 
#' system.time(c<-dMcast(orders,sku~state,
#'    value.var = 'amount')) # 8.66 seconds 
#'    
#' # However, this situation changes as the result set becomes larger 
#' system.time(a<-dcast.data.table(as.data.table(orders),customer~sku,sum,
#'    value.var = 'amount')) # 4.4 seconds 
#' system.time(b<-reshape2::dcast(orders,customer~sku,sum,
#'    value.var = 'amount')) # 34.7 seconds 
#'  system.time(c<-dMcast(orders,customer~sku,
#'    value.var = 'amount')) # 14.55 seconds 
#'    
#' # More complicated: 
#' system.time(a<-dcast.data.table(as.data.table(orders),customer~sku+state,sum,
#'    value.var = 'amount')) # 16.96 seconds, object size = 2084 Mb 
#' system.time(b<-reshape2::dcast(orders,customer~sku+state,sum,
#'    value.var = 'amount')) # Does not return 
#' system.time(c<-dMcast(orders,customer~sku:state,
#'    value.var = 'amount')) # 21.53 seconds, object size = 116.1 Mb
#' 
#' system.time(a<-dcast.data.table(as.data.table(orders),orderNum~sku,sum,
#'    value.var = 'amount')) # Does not return 
#' system.time(c<-dMcast(orders,orderNum~sku,
#'    value.var = 'amount')) # 24.83 seconds, object size = 175Mb
#' 
#' system.time(c<-dMcast(orders,sku:state~customer,
#'    value.var = 'amount')) # 17.97 seconds, object size = 175Mb
#'        
#' }
dMcast<-function(data,formula,fun.aggregate='sum',value.var=NULL,as.factors=FALSE,factor.nas=TRUE,drop.unused.levels=TRUE)
{
  values<-1
  if(!is.null(value.var))
    values<-data[,value.var]
  alltms<-terms(formula,data=data)
  response<-rownames(attr(alltms,'factors'))[attr(alltms,'response')]
  tm<-attr(alltms,"term.labels")
  interactionsIndex<-grep(':',tm)
  interactions<-tm[interactionsIndex]
  simple<-setdiff(tm,interactions)
  i2<-strsplit(interactions,':')
  newterms<-unlist(lapply(i2,function (x) paste("paste(",paste(x,collapse=','),",","sep='_'",")")))
  newterms<-c(simple,newterms)
  newformula<-as.formula(paste('~0+',paste(newterms,collapse='+')))
  allvars<-all.vars(alltms)
  data<-data[,c(allvars),drop=FALSE]
  if(as.factors)
    data<-data.frame(lapply(data,as.factor))
  characters<-unlist(lapply(data,is.character))
  data[,characters]<-lapply(data[,characters,drop=FALSE],as.factor)
  factors<-unlist(lapply(data,is.factor))
  #Prevents errors with 1 or fewer distinct levels
  data[,factors]<-lapply(data[,factors,drop=FALSE],function (x) 
  {
    if(factor.nas)
      if(any(is.na(x)))
      {
        levels(x)<-c(levels(x),'NA')
        x[is.na(x)]<-'NA'
      }
    if(drop.unused.levels)
        if(nlevels(x)!=length(na.omit(unique(x))))
          x<-factor(as.character(x))
    y<-contrasts(x,contrasts=FALSE,sparse=TRUE)
    attr(x,'contrasts')<-y
    return(x)
  })
  #Allows NAs to pass
  attr(data,'na.action')<-na.pass
  result<-sparse.model.matrix(newformula,data,drop.unused.levels = FALSE,row.names=FALSE)
  brokenNames<-grep('paste(',colnames(result),fixed = TRUE)
  colnames(result)[brokenNames]<-lapply(colnames(result)[brokenNames],function (x) {
    x<-gsub('paste(',replacement='',x=x,fixed = TRUE) 
    x<-gsub(pattern=', ',replacement='_',x=x,fixed=TRUE) 
    x<-gsub(pattern='_sep = \"_\")',replacement='',x=x,fixed=TRUE)
    return(x)
  })

  result<-result*values
  if(isTRUE(response>0))
  {
    responses=all.vars(terms(as.formula(paste(response,'~0'))))
    result<-aggregate.Matrix(result,data[,responses,drop=FALSE],fun=fun.aggregate)
  }
  return(result)
}

#' Compute summary statistics of a Matrix
#' 
#' Similar to \code{\link[stats]{aggregate}}.  Splits the matrix into groups as 
#' specified by groupings, which can be one or more variables. Aggregation 
#' function will be applied to all columns in data, or as specified in formula. 
#' Warning: groupings will be made dense if it is sparse, though data will not.
#' 
#' \code{aggregate.Matrix} uses its own implementations of functions and should
#' be passed a string in the \code{fun} argument.
#' 
#' @param x a \code{\link{Matrix}} or matrix-like object
#' @param groupings an object coercible to a group of factors defining the
#'   groups
#' @param form \code{\link[stats]{formula}}
#' @param fun character string specifying the name of aggregation function to be
#'   applied to all columns in data.  Currently "sum", "count", and "mean"
#'   are supported.
#' @param ... arguments to be passed to or from methods.  Currently ignored
#' @return A sparse \code{Matrix}.  The rownames correspond to the values of the
#'   groupings or the interactions of groupings joined by a \code{_}.
#'   
#'   There is an attribute \code{crosswalk} that includes the groupings as a
#'   data frame.  This is necessary because it is not possible to include
#'   character or data frame groupings in a sparse Matrix.  If needed, one can 
#'   \code{cbind(attr(x,"crosswalk"),x)} to combine the groupings and the
#'   aggregates.
#'   
#' @seealso \code{\link[dplyr]{summarise}}
#' @seealso \code{\link[plyr]{summarise}}
#' @seealso \code{\link[stats]{aggregate}}
#' @export
#' @export aggregate.Matrix
#' @examples
#' skus<-Matrix(as.matrix(data.frame(
#'    orderNum=sample(1000,10000,TRUE),
#'    sku=sample(1000,10000,TRUE),
#'    amount=runif(10000))),sparse=TRUE)
#'#Calculate sums for each sku
#' a<-aggregate.Matrix(skus[,'amount'],skus[,'sku',drop=FALSE],fun='sum')
#'#Calculate counts for each sku
#' b<-aggregate.Matrix(skus[,'amount'],skus[,'sku',drop=FALSE],fun='count')
#'#Calculate mean for each sku
#' c<-aggregate.Matrix(skus[,'amount'],skus[,'sku',drop=FALSE],fun='mean')
#' 
#' m<-rsparsematrix(1000000,100,.001)
#' labels<-as.factor(sample(1e4,1e6,TRUE))
#' b<-aggregate.Matrix(m,labels)
#' 
#' \dontrun{
#' orders<-data.frame(orderNum=as.factor(sample(1e6, 1e7, TRUE)),
#'    sku=as.factor(sample(1e3, 1e7, TRUE)),
#'    customer=as.factor(sample(1e4,1e7,TRUE)),
#'    state = sample(letters, 1e7, TRUE), amount=runif(1e7))
#' system.time(d<-aggregate.Matrix(orders[,'amount',drop=FALSE],orders$orderNum))
#' system.time(e<-aggregate.Matrix(orders[,'amount',drop=FALSE],orders[,c('customer','state')]))
#' }
aggregate.Matrix<-function(x,groupings=NULL,form=NULL,fun='sum',...)
{
  if(!is(x,'Matrix'))
    x<-Matrix(as.matrix(x),sparse=TRUE)
  if(fun=='count')
    x<-x!=0
  groupings2<-groupings
  if(!is(groupings2,'data.frame'))
    groupings2<-as(groupings2,'data.frame')
  groupings2<-data.frame(lapply(groupings2,as.factor))
  groupings2<-data.frame(interaction(groupings2,sep = '_'))
  colnames(groupings2)<-'A'
  if(is.null(form))
    form<-as.formula('~0+.')
  form<-as.formula(form)
  mapping<-dMcast(groupings2,form)
  colnames(mapping)<-substring(colnames(mapping),2)
  result<-t(mapping) %*% x
  if(fun=='mean')
    result@x<-result@x/(aggregate.Matrix(x,groupings2,fun='count'))@x
  attr(result,'crosswalk')<-grr::extract(groupings,match(rownames(result),groupings2$A))
  return(result)
}

aggregate2.Matrix<-function(x,groupings=NULL,form=NULL,fun=sum,...)
{
  #if(!is(x,'Matrix'))
    x<-as.matrix(x)
  groupings2<-groupings
  if(!is(groupings2,'data.frame'))
    groupings2<-as(groupings2,'data.frame')
  groupings2<-data.frame(lapply(groupings2,as.factor))
  groupings2<-interaction(groupings2,sep = '_')
  index<-grr::order2(groupings2)
  breaks<-which(!duplicated(groupings2[index]))
  results1<-fun(x[1:breaks[1],])
  results<-matrix(results1,ncol=length(results1),nrow=length(breaks))
  for (i in seq_len(length(breaks)-1)) 
  {
    results[i+1]<-fun(x[(breaks[i]+1):breaks[i+1],])
  }
  return(results)
}
  

#'Merges two Matrices or matrix-like objects
#'
#'Implementation of \code{\link{merge}} for \code{\link{Matrix}}.  By explicitly
#'calling \code{merge.Matrix} it will also work for \code{matrix}, for 
#'\code{data.frame}, and \code{vector} objects as a much faster alternative to 
#'the built-in \code{merge}.
#'
#'#' \code{all.x/all.y} correspond to the four types of database joins in the 
#'following way:
#'
#'\describe{ \item{left}{\code{all.x=TRUE}, \code{all.y=FALSE}} 
#'\item{right}{\code{all.x=FALSE}, \code{all.y=TRUE}} 
#'\item{inner}{\code{all.x=FALSE}, \code{all.y=FALSE}} 
#'\item{full}{\code{all.x=TRUE}, \code{all.y=TRUE}} }
#'
#'Note that \code{NA} values will match other \code{NA} values.
#'
#'@param x,y \code{Matrix} or matrix-like object
#'@param by.x vector indicating the names to match from \code{Matrix} x
#'@param by.y vector indicating the names to match from \code{Matrix} y
#'@param all.x logical; if \code{TRUE}, then each value in \code{x} will be 
#'  included even if it has no matching values in \code{y}
#'@param all.y logical; if \code{TRUE}, then each value in \code{y} will be 
#'  included even if it has no matching values in \code{x}
#'@param out.class the class of the output object.  Defaults to the class of x. 
#'  Note that some output classes are not possible due to R coercion
#'  capabilities, such as converting a character matrix to a Matrix.
#'@param fill.x,fill.y the value to put in merged columns where there is no match.
#'  Defaults to 0/FALSE for sparse matrices in order to preserve sparsity, NA for
#'  all other classes
#'@param ... arguments to be passed to or from methods.  Currently ignored
#'@export
#'@export merge.Matrix
#'@examples
#' 
#' orders<-Matrix(as.matrix(data.frame(orderNum=1:1000, 
#'  customer=sample(100,1000,TRUE)))) 
#'  cancelledOrders<-Matrix(as.matrix(data.frame(orderNum=sample(1000,100), 
#'  cancelled=1))) 
#' skus<-Matrix(as.matrix(data.frame(orderNum=sample(1000,10000,TRUE), 
#' sku=sample(1000,10000,TRUE), amount=runif(10000)))) 
#' a<-merge(orders,cancelledOrders,orders[,'orderNum'],cancelledOrders[,'orderNum'])
#' b<-merge(orders,cancelledOrders,orders[,'orderNum'],cancelledOrders[,'orderNum'],all.x=FALSE)
#' c<-merge(orders,skus,orders[,'orderNum'],skus[,'orderNum'])
#' 
#' #The above Matrices could be converted to matrices or data.frames and handled in other methods.  
#' #However, this is not possible in the sparse case, which can be handled by this function:
#' sm<-cbind2(1:200000,rsparsematrix(200000,10000,density=.0001))
#' sm2<-cbind2(sample(1:200000,50000,TRUE),rsparsematrix(200000,10,density=.01))
#' sm3<-merge.Matrix(sm,sm2,by.x=sm[,1],by.y=sm2[,1])
#' 
#'  \dontrun{
#' #merge.Matrix can also handle many other data types, such as data frames, and is generally fast.
#' orders<-data.frame(orderNum=as.character(sample(1e5, 1e6, TRUE)),
#'    sku=sample(1e3, 1e6, TRUE),
#'    customer=sample(1e4,1e6,TRUE),stringsAsFactors=FALSE)
#' cancelledOrders<-data.frame(orderNum=as.character(sample(1e5,1e4)),
#'    cancelled=1,stringsAsFactors=FALSE)
#' system.time(a<-merge.Matrix(orders,cancelledOrders,orders[,'orderNum'],
#'    cancelledOrders[,'orderNum']))
#' system.time(b<-merge.data.frame(orders,cancelledOrders,all.x = TRUE,all.y=TRUE))
#' system.time(c<-dplyr::full_join(orders,cancelledOrders))
#'system.time({require(data.table);
#'d<-merge(data.table(orders),data.table(cancelledOrders),
#'    by='orderNum',all=TRUE,allow.cartesian=TRUE)})
#'
#'orders<-data.frame(orderNum=sample(1e5, 1e6, TRUE), sku=sample(1e3, 1e6,
#'TRUE), customer=sample(1e4,1e6,TRUE),stringsAsFactors=FALSE) 
#'cancelledOrders<-data.frame(orderNum=sample(1e5,1e4),cancelled=1,stringsAsFactors=FALSE)
#'system.time(b<-merge.Matrix(orders,cancelledOrders,orders[,'orderNum'], 
#'cancelledOrders[,'orderNum'])) 
#'system.time(e<-dplyr::full_join(orders,cancelledOrders)) 
#'system.time({require(data.table);
#'  d<-merge(data.table(orders),data.table(cancelledOrders),
#'  by='orderNum',all=TRUE,allow.cartesian=TRUE)})
#'
#'#In certain cases, merge.Matrix can be much faster than alternatives. 
#'one<-as.character(1:1000000) 
#'two<-as.character(sample(1:1000000,1e5,TRUE)) 
#'system.time(b<-merge.Matrix(one,two,one,two)) 
#'system.time(c<-dplyr::full_join(data.frame(key=one),data.frame(key=two))) 
#'system.time({require(data.table);
#'  d<-merge(data.table(data.frame(key=one)),data.table(data.frame(key=two)),
#'  by='key',all=TRUE,allow.cartesian=TRUE)})
#'}
#'
merge.Matrix<-function(x,y,by.x,by.y,all.x=TRUE,all.y=TRUE,out.class=class(x)[1],
                       fill.x=ifelse(is(x,'sparseMatrix'),FALSE,NA),fill.y=fill.x,...)
{
  requireNamespace('grr')
  if(is.null(dim(x)))
    return(grr::matches(by.x,by.y,all.x,all.y,indexes=FALSE))
  indices<-grr::matches(by.x,by.y,all.x,all.y,nomatch = NULL)
  x<-rbind(x,fill.x)
  x<-as(grr::extract(x,indices$x),out.class)
  
  y<-rbind(y,fill.y)
  if(!is.null(colnames(x)) & !is.null(colnames(y)))
    colnames(y)[colnames(y) %in% colnames(x)]<-paste('y',colnames(y)[colnames(y) %in% colnames(x)],sep='.')

  y<-as(grr::extract(y,indices$y),out.class)
  result<-cbind2(x,y)
  return(result)
}

#' @rdname merge.Matrix
#' @export
join.Matrix<-merge.Matrix

#' Combine matrixes by row, fill in missing columns
#' 
#' rbinds a list of Matrix or matrix like objects, filling in missing columns.
#' 
#' Similar to \code{\link[plyr]{rbind.fill.matrix}}, but works for
#' \code{\link{Matrix}} as well as all other R objects.  It is completely
#' agnostic to class, and will produce an object of the class of the first input
#' (or of class \code{matrix} if the first object is one dimensional).
#' 
#' The implementation is recursive, so it can handle an arbitrary number of 
#' inputs, albeit inefficiently for large numbers of inputs.
#' 
#' This method is still experimental, but should work in most cases.  If the
#' data sets consist solely of data frames, \code{\link[plyr]{rbind.fill}} is 
#' preferred.
#' 
#' @param x,... Objects to combine.  If the first argument is a list and 
#'   \code{..} is unpopulated, the objects in that list will be combined.
#' @param fill value with which to fill unmatched columns
#' @param out.class the class of the output object.  Defaults to the class of x. 
#'  Note that some output classes are not possible due to R coercion
#'  capabilities, such as converting a character matrix to a Matrix.
#' @return a single object of the same class as the first input, or of class
#'   \code{matrix} if the first object is one dimensional
#' @seealso \code{\link[plyr]{rbind.fill}}
#' @seealso \code{\link[plyr]{rbind.fill.matrix}}
#' @export
#' @examples
#' df1 = data.frame(x = c(1,2,3), y = c(4,5,6))
#' rownames(df1) = c("a", "b", "c")
#' df2 = data.frame(x = c(7,8), z = c(9,10))
#' rownames(df2) = c("a", "d")
#' rBind.fill(df1,df2,fill=NA)
#' rBind.fill(as(df1,'Matrix'),df2,fill=0)
#' rBind.fill(as.matrix(df1),as(df2,'Matrix'),c(1,2),fill=0)
#' rBind.fill(c(1,2,3),list(4,5,6,7))
#' rBind.fill(df1,c(1,2,3,4))
#' 
#' m<-rsparsematrix(1000000,100,.001)
#' m2<-m
#' colnames(m)<-1:100
#' colnames(m2)<-3:102
#' system.time(b<-rBind.fill(m,m2))
#' 
rBind.fill<-function(x,...,fill=NULL,out.class=class(rbind(x,x))[1])
{
  if (is.list(x) && !is.data.frame(x) && missing(...)) {
    Reduce(function (x,y) rBind.fill.internal(x,y,fill,out.class),x)
  }
  else {
    Reduce(function (x,y) rBind.fill.internal(x,y,fill,out.class),list(x,...))
  }
}

rBind.fill.internal<-function(x,y,fill,out.class)
{
  out.class<-force(out.class)
  fillMissing<-is.null(fill)
  if(fillMissing)
    fill<-if(is(x,'Matrix')) 0 else NA
  if (is.null(nrow(x)))
       x<-matrix(x,nrow=1,dimnames=list(NULL,names(x)))
  if (is.null(nrow(y)))
      y<-matrix(y,nrow=1,dimnames=list(NULL,names(y)))
  
  nullNames<-FALSE
  #Cannot currently handle duplicate column names
  if(is.null(colnames(x)))
    colnames(x)<-make.names(colnames(y)[1:ncol(x)],unique = TRUE)
  if(is.null(colnames(y)))
     colnames(y)<-make.names(colnames(x)[1:ncol(y)],unique = TRUE)
  if(is.null(colnames(x)))
  {
    nullNames<-TRUE
    colnames(x)<-1:ncol(x)
    colnames(y)<-1:ncol(y)
  }
  ymiss<-colnames(x)[which(is.na(match(colnames(x),colnames(y))))]
  ybind<-rsparsematrix(nrow=nrow(y),ncol=length(ymiss),0)
  colnames(ybind)<-ymiss
  if(!fillMissing)
    ybind[seq_along(ybind)]<-fill
  xmiss<-colnames(y)[which(is.na(match(colnames(y),colnames(x))))]
  xbind<-rsparsematrix(nrow=nrow(x),ncol=length(xmiss),0)
  colnames(xbind)<-xmiss
  if(!fillMissing)
    xbind[seq_along(xbind)]<-fill
  if (ncol(xbind>0))
    x<-cbind2(x,xbind)
  if(ncol(ybind)>0)
    y<-cbind2(y,ybind)
  y<-as(y,out.class)
  x<-as(x,out.class)
  result<-rbind2(x,y[,order(match(colnames(y),colnames(x)))])
  if(nullNames)
    colnames(result)<-NULL
  return(result)
}

len<-function (data) 
{
  result <- ifelse(is.null(nrow(data)), length(data), nrow(data))
  return(result)
}

setAs('Matrix','data.frame',function (from) as.data.frame(as.matrix(from)),)

setAs('Matrix','data.frame',function (from) as.data.frame(as.matrix(from)))

setAs('data.frame','dgeMatrix', function (from) as(as.matrix(from),'dgeMatrix'))

setAs('data.frame','dgCMatrix', function (from) as(as.matrix(from),'dgCMatrix'))

setAs('matrix','data.frame',function (from) as.data.frame(from))

setAs('vector','data.frame',function (from) data.frame(from))

setMethod("cbind2",c('data.frame','Matrix'),function (x,y) 
{
  y<-as.matrix(y)
  cbind2(x,y)
}
)

setMethod("cbind2",c('Matrix','data.frame'),function (x,y) 
{
  y<-as.matrix(y)
  cbind2(x,y)
}
)