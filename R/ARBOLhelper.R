#! /usr/bin/Rscript

#require(vegan)
#require(pvclust)
#require(reshape2)
#packages <- c("Seurat","tidyverse","data.table","vegan","Matrix","pvclust",
#              "reshape2","Matrix.utils","tictoc","ggfortify","ggpubr","ggrepel",
#              "tidygraph","ggraph","viridis","data.tree","ape","phyloseq","ggalluvial",'scatterpie','ggnewscale')

#loadpak <- function(pkg){  
#  suppressPackageStartupMessages(sapply(pkg, require, character.only = TRUE))
#}

#loadpak(packages)

#' data.tree to phyloseq object conversion with custom ToNewick function that uses pathString as node names. 
#' @param node A data tree object
#' @return Newick format tree
#' @examples
#' txt <- ToNewickPS(divtree)
#' ARBOLphylo <- ape::read.tree(text=txt)
#' @export
#' Functions for data.tree to phyloseq object conversion with custom ToNewick function that uses pathString as node names. 
#' This is useful when node names are redundant throughout a tree, which happens in binary tree calculation with pvclust

ToNewickPS <- function(node, heightAttribute = DefaultPlotHeight, ...) {

  deparse <- function(x) {
    name <- stringi::stri_replace_all_fixed(x$pathString, " ", "_")
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

#' data.tree to ggraph object conversion with custom Node->dendrogram function that uses pathString as node names. 
#' @param object A data tree object
#' @param heightAttribute plot height attribute. Default is plotHeight
#' @return dendrogram object
#' @examples
#' txt <- ToNewickPS(divtree)
#' ARBOLphylo <- ape::read.tree(text=txt)
#' @export
#' Functions for data.tree to phyloseq object conversion with custom ToNewick function that uses pathString as node names. 
#' This is useful when node names are redundant throughout a tree, which happens in binary tree calculation with pvclust

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
    tsvlist <- list()
    sample_strings <- sub('\\..*$','',basename(tierfiles))
    sample_strings <- sub('_T0C0_', '',sample_strings)
    tiers <- map2(tierfiles, sample_strings, ~fread(.x,header=FALSE) %>% mutate(id = .y))
    tiers <- rbindlist(tiers) %>% data.frame
    tiers$id <- gsub('_','.',tiers$id)
    tiers <- tiers %>% separate(id,into=paste0('tier',1:maxtiers),sep='\\.',remove=FALSE)
    #lapply(3:maxtiers+2,function(x){cols<-colnames(tiers[,3:x]);print(cols)})
    tiers <- tiers %>% rename(CellID=V1)
    tiers <- tiers %>% data.frame %>% rename(tierNident=id)
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
#' dataframe <- dataframe %>% left_join(DiversityPerGroup(dataframe, species=pathString, group='sample')) %>% 
#'                            suppressMessages
#' metadata <- metadata %>% left_join(DiversityPerGroup(metadata, species=tierNident, group='sample')) %>% 
#'                          suppressMessages
#' @export
DiversityPerGroup <- function(df, species, group, diversity_metric = 'simpson') {
  #enquo parameters to allow dplyr calls
    divcols <- enquo(species)
    #count groups per species
    tierNcount <- df %>% group_by_at(vars(!!divcols)) %>% count(.data[[group]])
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
    tierNcount <- df %>% count(.data[[group]])
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
#' prepARBOLmeta_tree(srobj, maxtiers=10)
#' @export
prepARBOLmeta_tree <- function(srobj,maxtiers=10,categorical_attributes,diversity_attributes,numerical_attributes) {
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
      tierNcount <- numericaltb %>% group_by(tierNident) %>% count(.data[[attribute]])
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
#' calls prepARBOLmeta_tree() and prepTree()
#' @param srobj a seurat object with ARBOL 'tierNident', 'CellID', and 'sample' columns. 
#' @param categories columns in metadata for which you want to track categories present per node. Also finds the majority per node
#' @param diversities columns in metadata for which you want to calculate diversity per node 
#' @param diversity_metric one of 'shannon', 'simpson', or 'invsimpson'
#' @param counts columns in metadata for which you want to count each value per node
#' @return the input seurat object with tiered clustering tree in srobj@@misc$ARBOLclustertree, 
#' plot of tree to srobj@@misc$ARBOLclustertreeviz, and ggraph object to srobj@@misc$ARBOLclustertreeggraph
#' @examples
#' srobj <- sr_ARBOLclustertree(srobj)
#' @export
sr_ARBOLclustertree <- function(srobj, categories = 'sample', diversities = 'sample', diversity_metric = 'simpson', counts = 'sample') {

  if (!is.element('sample',diversities)) {
    diversities = c('sample',diversities)
  }

  if (!is.element('sample',categories)) {
    categories = c('sample',categories)
  }

  if (!is.element('sample',counts)) {
    counts = c('sample',counts)
  }

  srobj <- tierN_diversity(srobj, diversity_attributes = diversities, diversity_metric = diversity_metric)
  
  treemeta <- prepARBOLmeta_tree(srobj, categorical_attributes = categories, diversity_attributes = diversities, numerical_attributes = counts)

  ARBOLtree <- as.Node(treemeta) 

  Atree <- prepTree(ARBOLtree, srobj=srobj, diversity_attributes = diversities, categorical_attributes = categories, numerical_attributes = counts)

  #obtain counts columns from propagated data tree
  attrs <- Atree$attributesAll
  count_cols <- attrs[grep(sprintf('^(%s)_n',paste0(counts,collapse='|')),attrs)]

  ARBOLdf <- do.call(preppedTree_toDF, c(Atree,'n','pathString', categories, paste0(categories,'_majority'), paste0(diversities,'_diversity'),count_cols))

  #simple code call of translation to df disallowing custom diversity_attributes
  #ARBOLdf <- preppedTree_toDF(Atree, 'tier1', 'n', 'pathString', 'sample_diversity')

  ARBOLphylo <- as.phylo(ARBOLtree)
  #for original ARBOL clustering tree, set all edge lengths to 1
  ARBOLphylo$edge.length <- rep(1,ARBOLphylo$edge.length %>% length)

  #convert to tbl_graph object to allow easy plotting with ggraph
  x <- as_tbl_graph(ARBOLphylo,directed=T) %>% activate(nodes) %>% 
      left_join(ARBOLdf %>% select(name=levelName,n,all_of(categories),all_of(paste0(categories,'_majority')),all_of(paste0(diversities,'_diversity')),all_of(count_cols)))

  x <- x %>% activate(edges) #%>% left_join(ARBOLdf %>% select(to=i))

  x <- x %>% activate(nodes) %>% mutate(tier = str_count(name, "\\."))

  bt0 <- ggraph(x, layout = 'tree', circular=T) + 
  geom_edge_diagonal() + geom_node_point(size=0.3) + theme_classic()

  srobj@misc$ARBOLclustertree <- Atree
  srobj@misc$ARBOLclustertreeviz <- bt0
  srobj@misc$ARBOLclustertreeggraph <- x

  return(srobj)
}


#' Propagate attributes up a data.tree object to enable manipulation and visualization of all nodes
#' Automatically run when building trees using sr_ARBOLclustertree or sr_ARBOLbinarytree
#' @param srobj a seurat object with ARBOL 'tierNident' and 'sample' columns
#' @param ARBOLtree data.tree object
#' @param numerical_attributes numerical attributes are propagated up a tree as counts of each unique value per node.
#' creates a new column per unique value per attribute. counts are added as node$attribute_n_value
#' @param categorical_attributes categorical attributes are propagated up a tree as unique values per node. also calculates a majority per node
#' majority added as node$attribute_majority
#' @param diversity_attributes attributes for which to calculate diversity in each node. 
#' Supports simpson's, inverse simpson's, and shannon index as implemented in vegan diversity()
#' @param diversity_metric one of 'shannon', 'simpson', or 'invsimpson'
#' added as node$attribute_diversity
#' @return data.tree object with attributes propagated to all nodes
#' @examples
#' arbolTree <- prepTree(arbolTree,srobj=srobj, categorical_attributes = categories, diversity_attributes = diversities, numerical_attributes = NA)
#' @export
prepTree <- function(ARBOLtree, srobj, numerical_attributes = NA, categorical_attributes = NA, diversity_attributes = 'sample', diversity_metric = 'simpson') {

    #calculate number of children per node
    ARBOLtree$Do(function(node) node$numChildren <- node$children %>% length)

    #calculate total n samples for diversity propagation
    totalsamples <- srobj@meta.data$sample %>% unique %>% length

    #propagate node size up tree
    ARBOLtree$Do(function(node) node[['n']] <- Aggregate(node, attribute = 'n', aggFun = sum), traversal = "post-order")

    #propagate numerical variables up tree by counting occurrences of each variant z of category y in each node
    if(!is.na(numerical_attributes)) {                     
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
    if(!is.na(categorical_attributes)) { 
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
#' @param gene_list if using cell x gene data (not any srobj@@reduction), genes to include in centroid calculation. passes to get_Centroids
#' @param reduction_dims the dimensions of the reduction slot to use for centroid calculation. defaults to 1:25
#' @return the input seurat object with pvclust tree in srobj@@misc$pvclust
#' @examples
#' srobj <- sr_binarytree(srobj = srobj, tree_reduction = 'pca', hclust_method = 'complete',
#'                          distance_method = 'euclidean', centroid_method = 'median', 
#'                          centroid_assay = 'SCT', reduction_dims = 1:25)
#' @export
sr_binarytree <- function(srobj, tree_reduction = 'centroids', hclust_method = 'complete',
                                distance_method = 'euclidean', centroid_method = 'mean', 
                               centroid_assay = 'SCT', reduction_dims = 1:25, gene_list = rownames(srobj[["RNA"]]@data)) {

    if (tree_reduction == 'centroids' | tree_reduction %in% names(srobj@reductions)) {
      centroids <- get_Centroids(srobj = srobj, tree_reduction = tree_reduction, reduction_dims = reduction_dims,
                         centroid_method = centroid_method, centroid_assay = centroid_assay, gene_list = gene_list)

      result <- pvclust(centroids, method.dist=distance_method, method.hclust=hclust_method, nboot=1)

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
get_Centroids <- function(srobj = srobj, tree_reduction = tree_reduction, reduction_dims = reduction_dims,
                         centroid_method = centroid_method, centroid_assay = centroid_assay, gene_list = rownames(srobj[["RNA"]]@data)) {

    srobj_meta <- srobj@meta.data
    if (tree_reduction %in% names(srobj@reductions)) {
      scaled.data.mtx <- Embeddings(srobj, reduction=tree_reduction)[,reduction_dims]
    }
    else {
      gene_list = match(gene_list,rownames(srobj[[centroid_assay]]@data))
      sbmtx <- srobj[[centroid_assay]]@data[gene_list,]
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
#' calls sr_binarytree() in which assay + methods for tree building are called
#' @param srobj a seurat object with ARBOL 'tierNident', 'CellID', 'sample' columns. 
#' @param categories categorical attributes are propagated up a tree as unique values per node. also calculates a majority per node
#' @param diversities attributes for which to calculate diversity in each node. Currently only supports gini-simpson's index.
#' @param diversity_metric one of 'shannon', 'simpson', or 'invsimpson'
#' @param counts attributes for which to count unique values per node.
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
#' @return the input seurat object with binary tree pvclust object in srobj@@misc$pvclust, 
#' @examples
#' srobj <- sr_ARBOLbinarytree(srobj, categories = c('celltype','disease'))
#' @export
sr_ARBOLbinarytree <- function(srobj, categories = 'sample', diversities = 'sample', 
                                diversity_metric = 'simpson', counts = 'sample',
                                tree_reduction = 'centroids', hclust_method = 'complete',
                                distance_method = 'euclidean', centroid_method = 'mean', 
                                centroid_assay = 'SCT', reduction_dims = 1:25, gene_list = rownames(srobj[["RNA"]]@data)) {

  if(is.na(gene_list)) {
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
  
  srobj <- sr_binarytree(srobj = srobj, tree_reduction = tree_reduction, hclust_method = hclust_method,
                          distance_method = distance_method, centroid_method = centroid_method, 
                          centroid_assay = centroid_assay, gene_list = gene_list, reduction_dims = reduction_dims)

  binarydf <- bintree_to_df(pvclust_tree=srobj@misc$pvclust)

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
    tierNcount <- numericaltb %>% group_by(tierNident) %>% count(.data[[attribute]])
    vals <- unique(tierNcount[[x]])
    countWide <- tierNcount %>% pivot_wider(names_from=all_of(x),values_from=n)
    colnames(countWide) <- c('tierNident',sprintf('%s_n_%s',x,vals))
    numericaltb <- numericaltb %>% left_join(countWide)
  }

  numericaltb[is.na(numericaltb)] <- 0

  numericaltb <- numericaltb %>% select(-CellID,-sample,-all_of(counts)) %>% distinct

  finaltreedf <- binarydf %>% left_join(jointb) %>% left_join(categorydf) %>% left_join(divdf) %>% left_join(numericaltb)

  divtree <- as.Node(finaltreedf)

  divtree <- prepTree(divtree,srobj=srobj, categorical_attributes = categories, 
    diversity_attributes = diversities, numerical_attributes = counts)

  x <- data.tree_to_ggraph(divtree, categories, diversities, counts)

  x <- x %>% activate(nodes) %>% mutate(tier = str_count(label, "\\."))

  x <- x %>% activate(nodes) %>% mutate(string = label, name = basename(label) %>% str_replace_all('T0C0.',''))

  bt1 <- ggraph(x, layout = 'dendrogram', height=plotHeight*10) +
    geom_edge_elbow() + 
    geom_node_point(size=0) + 
    geom_node_text(aes(filter = leaf, label = name), nudge_y=-0.75,vjust=0.5,hjust=1.01,angle=90) + 
    geom_node_text(aes(filter = leaf, label = n),color='grey30',nudge_y=-0.2,vjust=0.5,hjust=1.01,size=3,angle=90) +
    theme_void() +
    geom_node_point(aes(filter = leaf,color=sample_diversity),size=4,shape='square') + 
    scale_color_gradient(low='grey90',high='grey10',limits=c(0,1)) +
    expand_limits(y=-5)

  srobj@misc$binarytreeviz <- bt1

  srobj@misc$binarytree <- divtree
  
  srobj@misc$binarytreeggraph <- x

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
#' 3) ggraph::as_tbl_graph to convert from ape to ggraph
#' 4) data.tree to node-level dataframe 
#' 5) join node-level dataframe to tbl_graph nodes
#' The use of ToNewick is a major bottleneck in computation speed right now
data.tree_to_ggraphNW <- function(data.tree, categories, diversities, counts) {
  txt <- ToNewickPS(data.tree)
  apeTree <- ape::read.tree(text=txt)
  attrs <- data.tree$attributesAll
  count_cols <- attrs[grep(sprintf('^(%s)_n',paste0(counts,collapse='|')),attrs)]
  treeDF <- do.call(preppedTree_toDF, c(data.tree, 'n','pathString','plotHeight', categories, paste0(categories,'_majority'),paste0(diversities,'_diversity'),count_cols))
  treeDF <- treeDF %>% select(name=pathString,n,i,plotHeight,all_of(categories),all_of(paste0(categories,'_majority')),all_of(paste0(diversities,'_diversity')),all_of(count_cols))
  x <- as_tbl_graph(apeTree,directed=T) %>% activate(nodes) %>% left_join(treeDF)
  x <- x %>% activate(edges) %>% left_join(treeDF %>% select(to=i,height=plotHeight))
  return(x)
}

#' data.tree to ggraph conversion, by converting data.tree structure to dendrogram via custom PS function
#' then ape (still doesn't carry annotations) to ggraph, 
#' then join data frame of data.tree to restore annotations.
#' requires 'n' and 'pathString' annotations in data.tree
#' Used to convert annotated binary phylogeny tree to ggraph for easier plotting 
#' 1) write data.tree object to Newick using custom Newick function,
#' 2) read Newick into ape tree object
#' 3) ggraph::as_tbl_graph to convert from ape to ggraph
#' 4) data.tree to node-level dataframe 
#' 5) join node-level dataframe to tbl_graph nodes
#' The use of ToNewick is a major bottleneck in computation speed right now
data.tree_to_ggraph <- function(data.tree, categories, diversities, counts, heightAttribute = 'plotHeight') {
  datadend <- as.dendrogram.NodePS(data.tree, heightAttribute = heightAttribute)
  x <- as_tbl_graph(datadend)
  attrs <- data.tree$attributesAll
  count_cols <- attrs[grep(sprintf('^(%s)_n',paste0(counts,collapse='|')),attrs)]
  treeDF <- do.call(preppedTree_toDF, c(data.tree, 'n','pathString','plotHeight', categories, paste0(categories,'_majority'),paste0(diversities,'_diversity'),count_cols))
  treeDF <- treeDF %>% select(label=pathString,n,i,plotHeight,all_of(categories),all_of(paste0(categories,'_majority')),all_of(paste0(diversities,'_diversity')),all_of(count_cols))
  x <- as_tbl_graph(datadend,directed=T) %>% activate(nodes) %>% left_join(treeDF)
  x <- x %>% activate(edges) %>% left_join(treeDF %>% select(to=i,height=plotHeight))
  return(x)
}

#' Merges tierNidents with their nearest neighbors in a binary tree if their sample diversity and number of cells do not meet thresholds
#' @param srobj a seurat object with a binarytree calculated in slot srobj@@misc$binarytree, typically calculated using sr_ARBOLbinarytree
#' @return the input seurat object with merged tierNidents in a new metadata column, mergedIdent
#' @examples
#' srobj <- MergeEndclusts(srobj, sample_diversity_threshold = 0.1, size_threshold = 10)
#' @export
MergeEndclusts <- function(srobj, sample_diversity_threshold, size_threshold) {

  srobj@misc$rawBinaryTree <- Clone(srobj@misc$binarytree)
  #DataTree::Prune chops all nodes that don't meet a threshold
  Prune(srobj@misc$binarytree, pruneFun = function(x) x$sample_diversity > sample_diversity_threshold)
  Prune(srobj@misc$binarytree, pruneFun = function(x) x$n > size_threshold)

  #remove unnecessary nodes that have only 1 child - these are created in binary tree threshold merging
  Prune(srobj@misc$binarytree, pruneFun = function(x) any(x$children %>% length > 1 || x$children %>% length == 0))

  divtestdf <- preppedTree_toDF(srobj@misc$binarytree, 'height', "pathString", 'ids')
  divdf2 <- divtestdf %>% mutate(ids = strsplit(ids, ", ")) %>% unnest
  divdf2 <- divdf2 %>% group_by(ids) %>% slice(which.min(height)) %>% ungroup

  divdf2 <- divdf2 %>% rename(CellID=ids,binIdent = pathString)
  divdf2$binIdent <- divdf2$binIdent %>% str_replace_all('\\/','.')

  srobj@meta.data <- srobj@meta.data %>% left_join(divdf2 %>% select(CellID,mergedIdent=binIdent)) %>%
                     rename(tierNident=mergedIdent,rawIdent=tierNident)

  return(srobj)
}

#' Merges tierNidents with their nearest neighbors in a binary tree by custom thresholds
#' @param srobj a seurat object with a binarytree calculated in slot srobj@@misc$binarytree, typically calculated using sr_ARBOLbinarytree
#' @return the input seurat object with merged tierNidents in a new metadata column, mergedIdent
#' @examples
#' srobj <- MergeEndclusts(srobj, sample_diversity_threshold = 0.1, size_threshold = 10)
#' @export
MergeEndclustsCustom <- function(srobj, threshold_attributes, thresholds) {

  srobj@misc$rawBinaryTree <- Clone(srobj@misc$binarytree)
  #DataTree::Prune chops all nodes that don't meet a threshold
  for (z in seq(1,length(threshold_attributes))) {
    Prune(srobj@misc$binarytree, pruneFun = function(x) x[[threshold_attributes[z]]] > thresholds[z])
  } 

  #remove unnecessary nodes that have only 1 child - these are created in binary tree threshold merging
  Prune(srobj@misc$binarytree, pruneFun = function(x) any(x$children %>% length > 1 || x$children %>% length == 0))

  divtestdf <- preppedTree_toDF(srobj@misc$binarytree, 'height', "pathString", 'ids')
  divdf2 <- divtestdf %>% mutate(ids = strsplit(ids, ", ")) %>% unnest
  divdf2 <- divdf2 %>% group_by(ids) %>% slice(which.min(height)) %>% ungroup

  divdf2 <- divdf2 %>% rename(CellID=ids,binIdent = pathString)
  divdf2$binIdent <- divdf2$binIdent %>% str_replace_all('\\/','.')

  srobj@meta.data <- srobj@meta.data %>% left_join(divdf2 %>% select(CellID,mergedIdent=binIdent)) %>%
                     rename(tierNident=mergedIdent,rawIdent=tierNident)

  return(srobj)
}

#' Converts pvclust tree object to dataframe with row per node
#' @param pvclust_tree a pvclust or hclust tree 
#' @return a dataframe with one row per node of the tree
#' @examples
#' binarydf <- bintree_to_df(pvclust_tree=srobj@@misc$pvclust)
#' @export
bintree_to_df <- function(pvclust_tree) {
  bintree <- as.Node(as.dendrogram(pvclust_tree$hclust))
  binarydf <- data.frame(ToDataFrameTree(bintree, 'pathString','plotHeight'))
  colnames(binarydf) <- c('remove','pathString','plotHeight')
  binarydf$tierNident <- sub('.*\\/', '', binarydf$pathString)
  binarydf <- binarydf %>% dplyr::select(-remove)
  return(binarydf)
}

#' Converts data.tree object to dataframe with row per node, this function aids conversion of highly annotated data.tree objects to dataframes.
#' Used in ARBOL workflow to create a tbl_graph object, a class useful for plotting dendrograms
#' @param tree usually a custom annotated data.tree object
#' @return a dataframe of the tree with one row per node
#' @examples
#' ARBOLdf <- preppedTree_toDF(divtree)
#' txt <- ARBOL::ToNewickPS(divtree)
#' ARBOLphylo <- ape::read.tree(text=txt)
#' 
#' x <- as_tbl_graph(ARBOLphylo, directed=T) %>% activate(nodes) %>% 
#' left_join(ARBOLdf %>% select(name=pathString,sample_diversity,disease_diversity))
#' @export
preppedTree_toDF <- function(tree, ...) {
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
        srobj@meta.data <- srobj@meta.data %>% left_join(DiversityPerGroup(srobj@meta.data, species=tierNident, group=z, diversity_metric = diversity_metric)) %>% suppressMessages
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
        df$tierfull[!grepl(paste0('T',x),df$tierfull)] <- NA
        names(df)[names(df) == 'tierfull'] <- paste0('tier',x,'full')
    }

    #remove delimiter
    df$tierNident = substr(df$tierNident,1,nchar(df$tierNident)-1)

    df2 <- df %>% unite('pathString', tier0:!!paste0('tier',max_tiers,'full'), sep = "/", na.rm = TRUE, remove = FALSE)

    return(df2)
}

#' Calculate curated names for set of clusters as function of their major type
#' @param srobj a seurat object with a celltype metadata column, specified in celltype_col, and tierNident i.e. srobj@@meta.data$tierNident
#' @param figdir the figdir from an GenTieredClusters (ARBOL) run 
#' @param max_cells_per_ident maximum cells to use in FindAllMarkers call. defaults to 200
#' @param celltype_col the metadata column corresponding to celltype to assign standard names
#' @param standardname_col metadata column to which to output celltype names. defaults to 'curatedname'
#' @param n_genes number of genes with which to create standard names. defaults to 2
#' @return the seurat object with curatedname column in metadata
#' @examples
#' srobj <- GetStandardNames(srobj,figdir='ARBOLoutput/figs')
#' @export
GetStandardNames <- function(srobj,figdir,max_cells_per_ident=200,celltype_col = 'celltype',standardname_col = 'curatedname',n_genes = 2) {

  Idents(srobj) <- srobj@meta.data[[celltype_col]]
  typeobjs <- SplitSrobjOnMeta(srobj, meta=celltype_col,'typeobjects')

  if (!file.exists(sprintf('%s/EndClustersAmongTier1DE.csv',figdir))) {
    CuratedtypeAmongType1markers <- lapply(typeobjs,function(obj) {Idents(obj) <- obj@meta.data$tierNident;
          tmp <- FindAllMarkers(obj,only.pos=TRUE,min.pct = 0.25,logfc.threshold = 0.25, max.cells.per.ident = max_cells_per_ident);
           return(tmp)})
    write.table(rbindlist(CuratedtypeAmongType1markers), sprintf('%s/EndClustersAmongTier1DE.csv',figdir), sep=",", row.names=F)
  }
  else {
    CuratedtypeAmongType1markers <- fread(sprintf('%s/EndClustersAmongTier1DE.csv',figdir)) %>%
                    split(f=.$cluster)
  }
  
  t <- lapply(CuratedtypeAmongType1markers,function(listmarkers) {
    tryCatch({
           fccol <- grep("FC",colnames(listmarkers))
            biomarkers <- listmarkers %>%
          mutate(rnkscr = -log(p_val_adj+1e-310) * sapply(listmarkers[,!!fccol], as.numeric) * (pct.1 / (pct.2 + 1e-300))) %>%
                    group_by(cluster) %>% top_n(n_genes, wt=rnkscr) %>% 
                    arrange(cluster, desc(rnkscr));

            endname <- biomarkers %>% summarize(markers = paste(gene, collapse="."))

              endname <- endname %>% select(cluster,markers)
              return(endname)
                 },
           error = function(e) {
                       #print(unique(listmarkers$cluster));
                       #endname = data.frame(cluster=unique(listmarkers$cluster))
                       #endname$markers <- NA
                               #return(endname)
                                     })
               })

  suppressWarnings( bind_rows(t,.id='subtype_id') -> markersAsList)

  message('number of end clusters for which at least two biomarkers were found: ',
       length(markersAsList$cluster))
  message('total number of end clusters: ',length(unique(srobj@meta.data$tierNident)))

  markersAsList <- markersAsList %>% rename(tierNident=cluster)

  srobj@meta.data$CellID <- row.names(srobj@meta.data)

  srobj@meta.data <- left_join(srobj@meta.data,markersAsList,by="tierNident")
  srobj@meta.data[[standardname_col]] <- paste(srobj@meta.data[[celltype_col]],srobj@meta.data$markers,sep='.')


  row.names(srobj@meta.data) <- srobj@meta.data$CellID

  return(srobj)
}

#' remake ggraph object with new categories and diversities, useful to add curatednames
#' @param srobj a seurat object with tierNident, sample, and category columns in metadata i.e. srobj@@meta.data$tierNident
#' @param categories categories as in tree building functions
#' @param diversities diversities as in tree building functions
#' @param diversity_metric one of 'shannon', 'simpson', or 'invsimpson'
#' @param counts attributes for which to count unique values per node.
#' @return the seurat object with ggraph object remade in srobj@@misc$binarytreeggraph
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

  binarydf <- bintree_to_df(pvclust_tree=srobj@misc$pvclust)

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
    tierNcount <- numericaltb %>% group_by(tierNident) %>% count(.data[[attribute]])
    vals <- unique(tierNcount[[x]])
    countWide <- tierNcount %>% pivot_wider(names_from=all_of(x),values_from=n)
    colnames(countWide) <- c('tierNident',sprintf('%s_n_%s',x,vals))
    numericaltb <- numericaltb %>% left_join(countWide)
  }

  numericaltb[is.na(numericaltb)] <- 0   

  numericaltb <- numericaltb %>% select(-CellID,-sample,-all_of(counts)) %>% distinct           

  finaltreedf <- binarydf %>% left_join(jointb) %>% left_join(categorydf) %>% left_join(divdf) %>% left_join(numericaltb)

  divtree <- as.Node(finaltreedf)

  divtree <- prepTree(divtree,srobj=srobj, categorical_attributes = categories, 
    diversity_attributes = diversities, numerical_attributes = counts)

  srobj@misc$binarytreeggraph <- data.tree_to_ggraph(divtree, categories, diversities, counts)

  return(srobj)

}




