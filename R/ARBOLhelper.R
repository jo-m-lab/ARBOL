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

#' Calculate Simpson's index per group in a dataframe (only creates 'diversity column', can't do multiple kinds),
#' @param df a dataframe with species and groups columns
#' @param species species column to stratify per group
#' @param group groups column for which to calculate diversity
#' @return dataframe with species and diversity columns
#' @examples
#' dataframe <- dataframe %>% left_join(SIperGroup(dataframe, species=pathString, group='sample')) %>% 
#'                            suppressMessages
#' metadata <- metadata %>% left_join(SIperGroup(metadata, species=tierNident, group='sample')) %>% 
#'                          suppressMessages
#' @export
SIperGroup <- function(df, species, group) {
  #enquo parameters to allow dplyr calls
    divcols <- enquo(species)
    groupcols <- enquo(group)
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
    simps <- diversity(tierNwide,'simpson') %>% data.frame %>% rownames_to_column(str_replace_all(as.character(vars(!!divcols)),"[^[:alnum:]]", ""))
    #Add 2 column names: your species category (i usually use pathString pulled from the node on the tree), and diversity. This joins to the dataframe by species category.
    colnames(simps) <- c(str_replace_all(as.character(vars(!!divcols)),"[^[:alnum:]]", ""),sprintf('%s_diversity',group))
    return(simps)
}

#' Calculate Simpson's index for a set of IDs in a dataframe
#' @param df a dataframe with group information per cell
#' @param group groups column for which to calculate diversity
#' @return a list of diversity per group per cell
#' @examples
#' in a data tree
#' node[[sprintf('%s_diversity',y)]] <- SIperIDs(metadata, group=y)
#' directly on a dataframe
#' metadata$condition_diversity <- SIperIDs(metadata, group=condition)
#' @export
SIperIDs <- function(df, group) {
  #enquo parameters to allow dplyr calls
    groupcols <- enquo(group)
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
    return(diversity(tierNwide,'simpson'))
}


#' Prepare ARBOL seurat object metadata for tree building
#' @param srobj a seurat object with ARBOL 'tierNident' and 'sample' columns
#' @param maxtiers max_tiers parameter used in ARBOL
#' @return a metadata dataframe ready for tree building
#' @examples
#' prepARBOLmeta_tree(srobj, maxtiers=10)
#' @export
prepARBOLmeta_tree <- function(srobj,maxtiers=10,categorical_attributes,diversity_attributes) {
    meta <- srobj@meta.data

    meta <- spread_tierN(meta,max_tiers=maxtiers)

    jointb <- srobj@meta.data %>% group_by(tierNident) %>% mutate(n=n()) %>% 
          dplyr::select(CellID,sample,tierNident,n,all_of(categorical_attributes),all_of(paste0(diversity_attributes,'_diversity')))

    categorydf <- jointb %>% summarize(across(categorical_attributes, ~ list(paste(unique(.x),collapse=', '))))
    divdf <- jointb %>% summarize_at(paste0(diversity_attributes,'_diversity'),unique)

    jointb <- jointb %>% select(-all_of(categorical_attributes)) %>% 
                summarize(ids = list(CellID),n=unique(n))

    treemeta <- meta %>% select(-all_of(categorical_attributes)) %>% left_join(jointb) %>% left_join(categorydf) %>% left_join(divdf)
    
    return(treemeta)
}

#' Prepare ARBOL seurat object metadata for tree building
#' @param srobj a seurat object with ARBOL 'tierNident' and 'sample' columns
#' @param maxtiers max_tiers parameter used in ARBOL
#' @return a metadata dataframe ready for tree building
#' @examples
#' prepARBOLmeta_tree(srobj, maxtiers=10)
#' @export
prepTree <- function(ARBOLtree, srobj, numerical_attributes = NA, categorical_attributes = NA, diversity_attributes = 'sample') {

    #calculate number of children per node
    ARBOLtree$Do(function(node) node$numChildren <- node$children %>% length)

    #calculate total n samples for diversity propagation
    totalsamples <- srobj@meta.data$sample %>% unique %>% length

    #propagate node size up tree
    ARBOLtree$Do(function(node) node[['n']] <- Aggregate(node, attribute = 'n', aggFun = sum), traversal = "post-order")

    #propagate additional numerical variables
    if(!is.na(numerical_attributes)) {                     
        for (y in numerical_attributes){
            ARBOLtree$Do(function(node) node[[y]] <- ifelse(is.na(node[[y]]) & node[['numChildren']]==0,0,node[[y]]))
            ARBOLtree$Do(function(node) node[[y]] <- Aggregate(node, attribute = y, aggFun = sum), traversal = "post-order")
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
                                node[[sprintf('%s_diversity',y)]] <- SIperIDs(meta, group=y)
        })
    }

    #call internal nodes as certain types depending on most prevalent major type annotation
    ARBOLtree$Do(function(node) {node$tier1 <- names(which.max(table(unlist(node$tier1))))})

    return(ARBOLtree)
}


#' Calculate pvclust() tree (a binary tree of distances between end-clusters) for ARBOL results
#' tree based on euclidean distance between cluster centroids based on gene medians with complete linkage
#' @param srobj a seurat object with ARBOL 'tierNident' column
#' @param assay the seurat object assay in which to calculate gene medians. Defaults to 'SCT'
#' @return the input seurat object with pvclust tree in srobj@@misc$pvclust
#' @examples
#' srobj <- sr_binarytree(srobj)
#' @export
sr_binarytree <- function(srobj,assay='SCT') {
    srobj_meta <- srobj@meta.data
    scaled.data.mtx <- Matrix(t(as.matrix(srobj[[assay]]@data)),sparse=TRUE)
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
    agg.clst.cntrs<- aggregate.Matrix(scaled.data.t[rows,-1],groupings=scaled.data.t[rows,1,drop=FALSE],fun='median')
    clst.cntrs <- agg.clst.cntrs

    #transpose cluster centers dataframe
    g <- clst.cntrs %>% as.matrix %>% t %>% data.frame

    #remember actual tierNident names
    colnames(g) <- t3$tierNident

    #hclust call. nboot = 1. bootstrapping will be performed by a different method later. 
    result <- pvclust(g, method.dist="euclidean", method.hclust="complete", nboot=1)

    srobj@misc$pvclust <- result
    return(srobj)
}

#' Creates data tree object for ARBOL run, adds it to seurat object along with a graph of it
#' calls prepARBOLmeta_tree() and prepTree()
#' @param srobj a seurat object with ARBOL 'tierNident', 'CellID', and 'sample' columns. 
#' @param diversity_attributes columns in metadata for which you wish to calculate diversity per node 
#' @return the input seurat object with tiered clustering tree in srobj@@misc$ARBOLclustertree, 
#' plot of tree to srobj@@misc$ARBOLclustertreeviz, and ggraph object to srobj@@misc$ARBOLclustertreeggraph
#' @examples
#' srobj <- sr_ARBOLclustertree(srobj)
#' @export
sr_ARBOLclustertree <- function(srobj, categories = 'sample', diversities = 'sample') {

  if (!is.element('sample',diversities)) {
    diversities = c('sample',diversities)
  }

  if (!is.element('sample',categories)) {
    categories = c('sample',categories)
  }

  srobj <- tierN_SI(srobj, diversity_attributes = diversities)
  
  treemeta <- prepARBOLmeta_tree(srobj, categorical_attributes = categories, diversity_attributes = diversities)

  ARBOLtree <- as.Node(treemeta) 

  Atree <- prepTree(ARBOLtree, srobj=srobj, diversity_attributes = diversities, categorical_attributes = categories)

  ARBOLdf <- do.call(preppedTree_toDF, c(Atree,'n','pathString', categories, paste0(categories,'_majority'), paste0(diversities,'_diversity')))

  #simple code call of translation to df disallowing custom diversity_attributes
  #ARBOLdf <- preppedTree_toDF(Atree, 'tier1', 'n', 'pathString', 'sample_diversity')

  ARBOLphylo <- as.phylo(ARBOLtree)
  #for original ARBOL clustering tree, set all edge lengths to 1
  ARBOLphylo$edge.length <- rep(1,ARBOLphylo$edge.length %>% length)

  #convert to tbl_graph object to allow easy plotting with ggraph
  x <- as_tbl_graph(ARBOLphylo,directed=T) %>% activate(nodes) %>% 
      left_join(ARBOLdf %>% select(name=levelName,n,all_of(categories),all_of(paste0(categories,'_majority')),all_of(paste0(diversities,'_diversity'))))

  x <- x %>% activate(edges) #%>% left_join(ARBOLdf %>% select(to=i))

  x <- x %>% activate(nodes) %>% mutate(tier = str_count(name, "\\."))

  bt0 <- ggraph(x, layout = 'tree', circular=T) + 
  geom_edge_diagonal() + geom_node_point(size=0.3)

  srobj@misc$ARBOLclustertree <- Atree
  srobj@misc$ARBOLclustertreeviz <- bt0
  srobj@misc$ARBOLclustertreeggraph <- x

  return(srobj)
}

#' Creates binary tree object for ARBOL run, adds it to seurat object along with a graph of it
#' calls sr_binarytree() in which assay + methods for tree building are called
#' @param srobj a seurat object with ARBOL 'tierNident', 'CellID', 'sample' columns. 
#' @return the input seurat object with binary tree pvclust object in srobj@@misc$pvclust, 
#' bina
#' @examples
#' srobj <- sr_ARBOLbinarytree(srobj, categories = c('tier1','type','disease'))
#' @export
sr_ARBOLbinarytree <- function(srobj, categories = 'sample', diversities = 'sample') {

  if (!is.element('sample',diversities)) {
    diversities = c('sample',diversities)
  }

  if (!is.element('sample',categories)) {
    categories = c('sample',categories)
  }
  
  srobj <- sr_binarytree(srobj)

  binarydf <- bintree_to_df(pvclust_tree=srobj@misc$pvclust)

  srobj <- tierN_SI(srobj, diversity_attributes = diversities)

  jointb <- srobj@meta.data %>% group_by(tierNident) %>% mutate(n=n()) %>% 
            dplyr::select(CellID,sample,tierNident,n,all_of(categories),all_of(paste0(diversities,'_diversity')))

  categorydf <- jointb %>% summarize(across(categories, ~ list(strsplit(paste(unique(.x),collapse=','),','))))

  divdf <- jointb %>% summarize_at(paste0(diversities,'_diversity'),unique)

  jointb <- jointb %>% select(-all_of(categories)) %>% 
              summarize(ids = list(CellID),n=unique(n))

  finaltreedf <- binarydf %>% left_join(jointb) %>% left_join(categorydf) %>% left_join(divdf)

  divtree <- as.Node(finaltreedf)

  divtree <- prepTree(divtree,srobj=srobj, categorical_attributes = categories, 
    diversity_attributes = diversities)

  x <- data.tree_to_ggraph(divtree, categories, diversities)

  x <- x %>% activate(nodes) %>% mutate(tier = str_count(name, "\\."))

  x <- x %>% activate(nodes) %>% mutate(string = name, name = basename(name) %>% str_replace_all('T0C0.',''))

  bt1 <- ggraph(x, layout = 'dendrogram') +
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

#' data.tree to ggraph conversion, by converting data.tree structure (doesn't carry annotations) to ggraph, 
#' then join data frame of data.tree
#' requires 'n' and 'pathString' annotations in data.tree
#' Used to convert annotated binary phylogeny tree to ggraph for easier plotting 
#' 1) write data.tree object to Newick using custom Newick function,
#' 2) read Newick into ape tree object
#' 3) ggraph::as_tbl_graph to convert from ape to ggraph
#' 4) data.tree to node-level dataframe 
#' 5) join node-level dataframe to tbl_graph nodes
#' 6) join node-level dataframe to tbl_graph edges
data.tree_to_ggraph <- function(data.tree, categories, diversities) {
  txt <- ToNewickPS(data.tree)
  apeTree <- ape::read.tree(text=txt)
  treeDF <- do.call(preppedTree_toDF, c(data.tree, 'n','pathString', categories, paste0(categories,'_majority'),paste0(diversities,'_diversity')))
  treeDF <- treeDF %>% select(name=pathString,n,i,all_of(categories),all_of(paste0(categories,'_majority')),all_of(paste0(diversities,'_diversity')))
  x <- as_tbl_graph(apeTree,directed=T) %>% activate(nodes) %>% left_join(treeDF)
  x <- x %>% activate(edges) %>% left_join(treeDF %>% select(to=i))
  return(x)
}

#' Merges tierNidents with their nearest neighbors in a binary tree if their sample diversity and number of cells do not meet thresholds
#' @param srobj a seurat object with a binarytree calculated in slot srobj@@misc$binarytree, typically calculated using sr_ARBOLbinarytree
#' @return the input seurat object with merged tierNidents in a new metadata column, mergedIdent
#' @examples
#' srobj <- MergeEndclusts(srobj, sample_diversity_threshold = 0.1, size_threshold = 10)
#' @export
MergeEndclusts <- function(srobj, sample_diversity_threshold, size_threshold) {

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



#' Converts pvclust tree object to dataframe with row per node
#' calls prepARBOLmeta_tree() and prepTree()
#' @param pvclust_tree a pvclust or hclust tree 
#' @return a dataframe with one row per node of the tree
#' @examples
#' binarydf <- bintree_to_df(pvclust_tree=srobj@@misc$pvclust)
#' @export
bintree_to_df <- function(pvclust_tree) {
  bintree <- as.Node(as.dendrogram(pvclust_tree$hclust))
  binarydf <- data.frame(ToDataFrameTable(bintree, 'pathString'))
  colnames(binarydf) <- c('pathString')
  binarydf$tierNident <- sub('.*\\/', '', binarydf$pathString)
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

#' Calculate Gini-Simpson's Index for tierNidents for multiple attributes of the data
#' @param srobj a seurat object with ARBOL tierNident column
#' @param diversity_attributes the attributes you wish to calculate diversity for (e.g. disease, sample)
#' @return the seurat object with diversity per attribute per tierNident added to the metadata
#' @examples
#' srobj <- tierN_SI(srobj, diversity_attributes = c('sample','disease'))
#' @export
tierN_SI <- function(srobj, diversity_attributes) {
  #remove existing diversity metrics
  srobj@meta.data <- srobj@meta.data %>% dplyr::select(-any_of(paste0(diversity_attributes, '_diversity')))
    for(z in diversity_attributes) {
        #calculate and join new diversity
        srobj@meta.data <- srobj@meta.data %>% left_join(SIperGroup(srobj@meta.data, species=tierNident, group=z)) %>% suppressMessages
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
#' @param srobj a seurat object with metadata columns 'celltype' and 'tierNident' i.e. srobj@@meta.data$celltype
#' @param figdir the figdir from an GenTieredClusters (ARBOL) run 
#' @param max_cells_per_ident maximum cells to use in FindAllMarkers call. defaults to 200
#' @return the seurat object with curatedname column in metadata
#' @examples
#' srobj <- GetStandardNames(srobj,figdir='ARBOLoutput/figs')
#' @export

GetStandardNames <- function(srobj,figdir,max_cells_per_ident=200) {

  Idents(srobj) <- srobj@meta.data$celltype
  typeobjs <- SplitSrobjOnMeta(srobj, meta="celltype",'typeobjects')

  if(!file.exists(sprintf('%s/curatedtypeAmongtype1markers.rds',figdir))) {
    CuratedtypeAmongType1markers <- lapply(typeobjs,function(obj) {
          if (obj@meta.data$tierNident %>% unique %>% length > 1){
              Idents(obj) <- obj@meta.data$tierNident;
              tmp <- FindAllMarkers(obj,only.pos=TRUE,min.pct = 0.25,logfc.threshold = 0.25, max.cells.per.ident = max_cells_per_ident);
               return(tmp)
              }
          else {
              return(data.frame())
          }})
    saveRDS(CuratedtypeAmongType1markers,sprintf('%s/curatedtypeAmongtype1markers.rds',figdir))
  }
  else(
    CuratedtypeAmongType1markers <- readRDS(sprintf('%s/curatedtypeAmongtype1markers.rds',figdir))
  )

  t <- lapply(CuratedtypeAmongType1markers,function(listmarkers,objs) {
        tryCatch({
                biomarkers <- listmarkers %>%
              mutate(rnkscr = -log(p_val_adj+1e-310) * sapply(avg_logFC, as.numeric) * 
                      (pct.1 / (pct.2 + 1e-300))) %>%
                        group_by(cluster) %>% top_n(2, wt=rnkscr) %>% 
                        arrange(cluster, desc(rnkscr));

                lapply(seq_along(biomarkers$gene), function(x) 
                  paste(biomarkers$gene[x], biomarkers$gene[x+1],sep='.')) %>% unlist -> top2biomarkers
                  top2biomarkers <- top2biomarkers[seq(1,length(top2biomarkers),2)]
                  biomarkers$markers <- rep(top2biomarkers,each=2)
                  endname <- biomarkers %>% select(cluster,markers) %>% summarise_each(funs(max))
                  return(endname)
                     },
               error = function(e) {

               })
             })

  suppressWarnings(bind_rows(t,.id='type_id') -> markersAsList)
  message('number of end clusters for which at least two biomarkers were found: ',
       length(markersAsList$cluster))
  message('total number of end clusters: ',length(unique(srobj@meta.data$tierNident)))
  markersAsList <- markersAsList %>% rename(tierNident=cluster)
  srobj@meta.data <- left_join(srobj@meta.data,markersAsList,by="tierNident") %>% 
  mutate(curatedname = paste(celltype,markers,sep='.'))
  return(srobj)
}












