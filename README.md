# SCTieredClustering

Iterively clusters v3 seurat object from single cell datasets, choosing optimum 
resolution parameters at each stage of clustering. Outputs QC plots for each tier and stage.
can return output directories figs/ srobjs/ (see examples) containing subset 
seurat objects and QC plots matching tiered structure of found in dataset
can also return endclusts/ folder which contains tsv files named with the 
path to the endcluster (T0C0_T1C1.tsv) and contains the cellnames (one per line) 
that comprise that cluster

## Install

```
git clone https://github.com/ShalekLab/SCTieredClustering.git
```

## Recommended Usage

This package is meant as a starting point and to be edited/customized by YOU!

There are a lot of decisions that go into this analysis. I have put a lot of 
work into choosing reasonable defaults, but there is no certainty that they are
the best defaults for your data.

I recommend cloning the git repository, and looking directly at the
`./R/SCTieredClustering.R` script. I have tried to organize the script such that
each processing step is contained in a modular function that can be edited and
inserted into the larger clustering steps.

As a starting point, and to use the defaults here is some example code.


```
source("path/to/cloned/git/repo/R/SCTieredClustering.R")

srobj <- readRDS("/path/to/full_seurat_object.rds")
tiers <- GenTieredClusters(srobj,
                           saveSROBJdir = "~/tieredoutput/srobjs",
                           figdir = "~/tieredoutput/figdir",
                           SaveEndNamesDir = "~/tieredoutput/endclusts")
```

**Note** This script can take a long time to run. running on 20K cells could 
take a few hours. Running on 100k+ cells could take over a day. This timing varies
based on the heterogeneity of your data.

**Note** RAM is also a consideration, for running on ~100k cells I needed 256GB of RAM.

A guided tutorial for analysis is available: https://shaleklab.github.io/SCTieredClustering/

## Params

* *srobj* v3 seurat object
* *cluster_assay* assay to use for clustering defaults to "SCT"
* *cells* cellnames if tiered clustering should start on subset of object
* *tier* starting level defaults to 0
* *clustN* cluster starting point. default to 0
* *PreProcess_fun* function to use for preproccessing defaults to PreProcess_sctransform
* *min_cluster_size* minimum number of cells to allow further clustering
* *max_tiers* maximum number of tiers to allow further clustering
* *EnoughDiffUp* minimum number of up-regulated genes to call clusters unique. Differential expression is performed when clustering finds 2 clusters
* *EnoughDiffDown* minimum number of down-regulated genes. If either up or down is not met, the 2 clusters are joined, and further clustering is stopped
* *ChooseOptimalClustering_fun* function that returns srobj with clusters in `srobj$Best.Clusters` after choosing optimal clustering resolution
* *saveSROBJdir* where to save seurat objects for each tier and cluster, if null does not save
* *figdir* where to save QC figures for each tier and cluster, if null does not save
* *saveEndNamesDir* where to save directory of end clusters, if null does not save
* *SaveEndFileName* prefix for all end cluster files

## Returns

* list of lists with all seurat objects (highly recommend using folder arguments for saving outputs)
