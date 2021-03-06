<img src="docs/ARBOLsmall.jpg?raw=true" align="right" width=500px>  

Iteratively cluster single cell datasets using a Seurat v4 object as input. This method identifies and utilizes optimum 
cluster resolution parameters at each tier of clustering. It provides outputs of QC plots for each tier and stage.
Furthermore, it can return output directories figs/ srobjs/ (see examples) containing subset 
seurat objects and QC plots. It can also return an endclusts/ folder which contains tsv files named with the 
path to the endcluster (T0C0_T1C1.tsv) and containing the cellnames (one per line) 
that comprise that cluster.

## Install

within R:
```
devtools::install_github('jo-m-lab/ARBOL')
library(Seurat)
library(ARBOL)
```

or clone the repository and source the functions directly from the script
```
git clone https://github.com/jo-m-lab/ARBOL.git

source("path/to/cloned/git/repo/R/ARBOL.R")
```

there is a docker image available with ARBOL and dependencies preinstalled
https://hub.docker.com/r/kkimler/arbol

## Recommended Usage

ARBOL was developed and used in the paper, "A treatment-naïve cellular atlas of pediatric Crohn’s disease predicts disease severity and therapeutic response"
We include here a tutorial where the FGID atlas figure is reproduced: 
https://jo-m-lab.github.io/ARBOL/ARBOLtutorial.html

This package is meant as a starting point for the way that we approached clustering and and is meant to be edited/customized through community feedback through users such as yourself!

We have dedicated effort to choosing reasonable defaults, but there is no certainty that they are
the best defaults for your data.

We recommend installing the repository through devtools. We have tried to organize the script such that
each processing step is contained in a modular function that can be edited and
inserted into the main iterative function.

The main function of ARBOL is GenTieredClusters() - here is an example call

```
srobj <- readRDS("/path/to/full_seurat_object.rds")
tiers <- GenTieredClusters(srobj,
                           saveSROBJdir = "~/tieredoutput/srobjs",
                           figdir = "~/tieredoutput/figdir",
                           SaveEndNamesDir = "~/tieredoutput/endclusts")
```

**Note** This script can take a long time to run. Running on 20K cells could 
take a few hours. Running on 100k+ cells could take over a day. This timing varies
based on the heterogeneity of your data.

**Note** The script requires approximately 1.2 GB RAM per 1k cells, meaning on a local machine with 16GB RAM, one could reasonably run 12k cells. The current RAM/time bottleneck is the SCTransform() call, which is run at each tier to renormalize to the input subset. 

## Params

* *srobj* Seurat v4 object
* *cluster_assay* assay to use for clustering defaults to "SCT"
* *cells* cellnames if tiered clustering should start on subset of object
* *tier* starting level defaults to 0
* *clustN* cluster starting point. default to 0
* *PreProcess_fun* function to use for preproccessing defaults to PreProcess_sctransform.  PreProcess_sctransform_harmony now available
* *min_cluster_size* minimum number of cells to allow further clustering
* *max_tiers* maximum number of tiers to allow further clustering
* *EnoughDiffUp* minimum number of up-regulated genes to call clusters unique. Differential expression is performed when clustering finds 2 clusters
* *EnoughDiffDown* minimum number of down-regulated genes. If either up or down is not met, the 2 clusters are joined, and further clustering is stopped
* *tierAllowedRecomb* minimum tier where differential expression can be called to decide on recombination. defaults to 0. clustering may stop early when clustering finds 2 clusters with high cell numbers, as Wilcoxon effect sizes may be low in their DE.
* *ChooseOptimalClustering_fun* function that returns srobj with clusters in `srobj$Best.Clusters` after choosing optimal clustering resolution
* *harmony_var* variable over which to iterate harmony integration if using PreProcess_sctransform_harmony
* *saveSROBJdir* where to save seurat objects for each tier and cluster, if null does not save
* *figdir* where to save QC figures for each tier and cluster, if null does not save
* *saveEndNamesDir* where to save directory of end clusters, if null does not save
* *SaveEndFileName* prefix for all end cluster files

## Returns

* list of lists with all seurat objects (highly recommend using folder arguments for saving outputs)
