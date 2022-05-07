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

#' makes data.tree names unique
#' @param l A data tree object
#' @return a data tree object with names forced to be unique. Useful when producing node-level dataframes from raw data tree objects
#' @examples
#' makeDataTreeNamesUnique(ARBOLtree)
#' @export
makeDataTreeNamesUnique <- function(l) {
  l.names <- names(l$children)
  # multiple children types
  tab <- table(l.names)
  t.names <- names(tab)

  # iterate over types
  for(this.type in seq_along(t.names)) {
    # iterate over duplicate names
    # get an index to this type
    idx <- which(l.names == t.names[this.type])
    for(this.element in seq_along(idx)) {
      # make a copy of this chunk of the tree
      l.sub <- l$Children[[idx[this.element]]]
      # if this is a terminal leaf then re-name and continue
      if(is.null(l.sub$Children)) {
        # print('leaf')
        names(l$Children)[idx[this.element]] <- paste0(t.names[this.type], '_', this.element)
      }
      # otherwise re-name and then step into this element and apply this function recursively
      else {
        # print('branch')
        names(l$Children)[idx[this.element]] <- paste0(t.names[this.type], '_', this.element)
        # fix this branch and splice back into tree
        l$Children[[idx[this.element]]] <- makeNamesUnique(l.sub)
      }
    }
  }

  return(l)
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
prepARBOLmeta_tree <- function(srobj,maxtiers=10) {
    meta <- srobj@meta.data
    meta <- meta %>% separate(tierNident,into=paste0('tier',seq(1,maxtiers)),remove=F)
    meta$tier0 <- 'T0C0'

    for(x in seq(1,maxtiers)) {
        meta <- meta %>% mutate(tierfull = strex::str_before_nth(tierNident, '\\.', n=x))
        meta$tierfull[!grepl(paste0('T',x),meta$tierfull)] <- NA
        names(meta)[names(meta) == 'tierfull'] <- paste0('tier',x,'full')
    }

    meta2 <- meta %>% unite('pathString', tier0, tier1full, tier2full,tier3full,
                             tier4full,tier5full,tier6full,tier7full,tier8full,tier9full,tier10full, sep = "/", na.rm = TRUE)

    #make sure sample metadata column exists
    if('sample' %in% colnames(meta)) {message('found sample column')} else {return(message('no sample column'))}

    tierNcount <- meta2 %>%
        dplyr::count(sample,tierNident)

    meta2 <- meta2 %>% left_join(SIperGroup(srobj@meta.data, species=tierNident, group='sample')) %>% suppressMessages

    jointb <- meta2 %>% 
        group_by(tierNident) %>%
        mutate(n = n()) %>% 
        left_join(tierNpresence, by='tierNident') %>%
        dplyr::select(CellID,sample,tierNident,sample_diversity,n) %>% 
        summarize(ids = list(CellID),samples = list(unique(sample)),
                    diversity=unique(sample_diversity),n=unique(n))

    meta3 <- meta2 %>% left_join(jointb)
      
    return(meta3)
}

#' Prepare ARBOL seurat object metadata for tree building
#' @param srobj a seurat object with ARBOL 'tierNident' and 'sample' columns
#' @param maxtiers max_tiers parameter used in ARBOL
#' @return a metadata dataframe ready for tree building
#' @examples
#' prepARBOLmeta_tree(srobj, maxtiers=10)
#' @export
prepTree <- function(ARBOLtree, srobj, numerical_attributes, categorical_attributes = 'samples', diversity_attributes = 'sample') {

    #calculate number of children per node
    ARBOLtree$Do(function(node) node$numChildren <- node$children %>% length)

    #calculate total n samples for diversity propagation
    totalsamples <- srobj@meta.data$sample %>% unique %>% length

    for (y in numerical_attributes){
        ARBOLtree$Do(function(node) node[[y]] <- ifelse(is.na(node[[y]]) & node[['numChildren']]==0,0,node[[y]]))
        ARBOLtree$Do(function(node) node[[y]] <- Aggregate(node, attribute = y, aggFun = sum), traversal = "post-order")
    }

    for (y in categorical_attributes) {
        ARBOLtree$Do(function(node) node[[y]] <- Aggregate(node, attribute = y, aggFun = c), traversal = "post-order")
    }

    #make samples only unique values
    ARBOLtree$Do(function(node) node$samples <- unique(unlist(node$samples)))

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
    scaled.data.mtx <- Matrix(t(srobj[[assay]]@data),sparse=TRUE)
        #pull seurat object's cell ID + tierNident
    t2 <- srobj_meta[,c('CellID','tierNident')] #%>% column_to_rownames('CellID')
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
    g <- clst.cntrs %>% t %>% data.frame

    #remember actual tierNident names
    colnames(g) <- t3$tierNident

    #hclust call. nboot = 1. bootstrapping will be performed by a different method later. 
    result <- pvclust(g, method.dist="euclidean", method.hclust="complete", nboot=1)

    srobj@misc$pvclust <- result
    return(srobj)
}

#' Creates data tree object for ARBOL run, adds it to seurat object along with a graph of it
#' calls prepARBOLmeta_tree() and prepTree()
#' @param srobj a seurat object with ARBOL 'tierNident' and 'sample' columns. 
#' @return the input seurat object with tiered clustering tree in srobj@@misc$pvclust and ggraph of tree to srobj@@misc$ARBOLclustertreegraph
#' @examples
#' srobj <- sr_ARBOLclustertree(srobj)
#' @export
sr_ARBOLclustertree <- function(srobj) {
  treemeta <- prepARBOLmeta_tree(srobj)

  ARBOLtree <- as.Node(treemeta) 

  Atree <- prepTree(ARBOLtree)

  #create dataframe from the resulting annotation-propagated tree, including internal nodes
  ARBOLdf <- ToDataFrameTree(Atree,  "tier1", 'diversity', 'disease_diversity', 'n', 'pathString')

  #remove graphics from levelName column
  ARBOLdf$levelName <- ARBOLdf$levelName %>% str_replace_all(' ','') %>% str_replace_all('-','') %>%
                          str_replace_all('¦','') %>% str_replace_all('°','')

  #add a row index column for joining with ggraph object                        
  ARBOLdf$i <- as.numeric(row.names(ARBOLdf))

  ARBOLphylo <- as.phylo(ARBOLtree)
  #for original ARBOL clustering tree, set all edge lengths to 1
  ARBOLphylo$edge.length <- rep(1,ARBOLphylo$edge.length %>% length)

  #convert to tbl_graph object to allow easy plotting with ggraph
  x <- as_tbl_graph(ARBOLphylo,directed=T) %>% activate(nodes) %>% left_join(ARBOLdf %>% select(name=levelName,tier1,diversity,n))
  x <- x %>% activate(edges) %>% left_join(ARBOLdf %>% select(to=i,tier1))
  x <- x %>% activate(nodes) %>% mutate(tier = str_count(name, "\\."))

  bt0 <- ggraph(x, layout = 'tree', circular=T) + 
  #color branches and nodes by type or tier1 (must propagate this metadata in prepTree)
    #geom_edge_elbow(aes(colour=factor(tier1))) + geom_node_point(aes(colour = tier1))
  geom_edge_elbow() + geom_node_point(size=0.3)

  srobj@misc$ARBOLclustertree <- Atree
  srobj@misc$ARBOLclustertreegraph <- bt0

  return(srobj)
}

#' Converts pvclust tree object to dataframe with row per node
#' calls prepARBOLmeta_tree() and prepTree()
#' @param pvclust_tree a pvclust or hclust tree 
#' @return a dataframe with one row per node of the tree
#' @examples
#' binarydf <- bintree.to.df(pvclust_tree=tsr@@misc$pvclust)
#' @export
bintree.to.df <- function(pvclust_tree) {
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
#' txt <- ToNewick(divtree)
#' ARBOLphylo <- ape::read.tree(text=txt)
#' 
#' x <- as_tbl_graph(ARBOLphylo, directed=T) %>% activate(nodes) %>% 
#' left_join(ARBOLdf %>% select(name=pathString,sample_diversity,disease_diversity))
#' @export
preppedTree_toDF <- function(tree, ...) {
    ARBOLdf <- ToDataFrameTree(divtree, ...)

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
    for(z in diversity_attributes) {
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
    df <- df %>% separate(tierNident,into=paste0('tier',seq(1,maxtiers)),remove=F)
    df$tier0 <- 'T0C0'

    for(x in seq(1,maxtiers)) {
        df <- df %>% mutate(tierfull = strex::str_before_nth(tierNident, '\\.', n=x))
        df$tierfull[!grepl(paste0('T',x),df$tierfull)] <- NA
        names(df)[names(df) == 'tierfull'] <- paste0('tier',x,'full')
    }

    df2 <- df %>% unite('pathString', tier0, tier1full, tier2full,tier3full,
                             tier4full,tier5full,tier6full,tier7full,tier8full,tier9full,tier10full, sep = "/", na.rm = TRUE)
    return(df2)
}

