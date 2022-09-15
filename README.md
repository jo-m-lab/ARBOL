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
library(ARBOL)
library(tidyverse)
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

Here is a vignette where ARBOL visualization and analysis is performed: 
https://jo-m-lab.github.io/ARBOL/ARBOLtutorial_22_9_11.html

This package is meant as a starting point for the way that we approached clustering and and is meant to be edited/customized through community feedback through users such as yourself!  We have tried to organize the script such that
each processing step is contained in a modular function that can be edited and
inserted into the main iterative function.

We have dedicated effort to choosing reasonable defaults, but there is no certainty that they are
the best defaults for your data.

We recommend installing the repository through devtools.

The main function of ARBOL is ARBOL() - here is an example call

```
srobj <- readRDS("/path/to/full_seurat_object.rds")
tiers <- ARBOL(srobj,
                           saveSROBJdir = "~/tieredoutput/srobjs",
                           figdir = "~/tieredoutput/figdir",
                           SaveEndNamesDir = "~/tieredoutput/endclusts")
```

**Note** ARBOL() can take a long time to run. Running on 20K cells could 
take a few hours on a low-power machine. Running on 100k+ cells could take over a day. This timing varies
based on the heterogeneity of your data. 

**Note** There is a native python version of the ARBOL() function which acts on a scanpy anndata object 

https://github.com/jo-m-lab/ARBOLpy

**Note** The R script requires approximately 1.2 GB RAM per 1k cells, meaning on a local machine with 16GB RAM, one could reasonably run 12k cells. The current RAM/time bottleneck is the SCTransform() normalization, which is run at each tier to re-normalize to the input subset. 

## Params

* *srobj* v4 seurat object
* *cluster_assay* assay to use for clustering defaults to "SCT"
* *cells* cellnames if tiered clustering should start on subset of object
* *tier* starting level defaults to 0
* *clustN* cluster starting from default to 0
* *PreProcess_fun* function to use for preproccessing defaults to PreProcess_sctransform. PreProcess_sctransform_harmony now available.
* *min_cluster_size* minimum number of cells to allow further clustering. defaults to 100
* *max_tiers* maximum number of tiers to allow further clustering. defaults to 10
* *EnoughDiffUp* minimum number of up-regulated genes to call clusters unique. Differential expression is performed when clustering finds 2 clusters. defaults to 5
* *EnoughDiffDown* minimum number of down-regulated genes. If either up or down is not met, the 2 clusters are joined, and further clustering is stopped. defaults to 5
* *tierAllowedRecomb* minimum tier where differential expression can be called to decide on recombination. defaults to 0. clustering may stop early when clustering finds 2 clusters with high cell numbers, as Wilcoxon effect sizes may be low.
* *ChooseOptimalClustering_fun* function that returns srobj with clusters in `srobj$Best.Clusters` after choosing optimal clustering resolution 
* *res_scan_step* number of resolutions of decreasing score before an early silhouette analysis resolution scan stop
* *res_scan_min* the smallest resolution to scan
* *res_scan_max* the highest resolution to scan (default is 3 due to legacy, we recommend lowering this)
* *res_scan_n* the number of resolutions to scan. A logarithmic sequence of numbers between min and max are scanned
* *saveSROBJdir* where to save seurat objects for each tier and cluster, if null does not save
* *figdir* where to save QC figures for each tier and cluster, if null does not save
* *saveEndNamesDir* where to save directory of end clusters, if null does not save
* *SaveEndFileName* prefix for all end cluster files
* *harmony_var* variable over which to iterate harmony integration

## Returns

* dataframe with tierN cluster membership per cell
