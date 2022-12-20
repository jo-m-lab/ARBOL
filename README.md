<img src="docs/ARBOLsmall.jpg?raw=true" align="right" width=500px>  

Iteratively cluster single cell datasets using a Seurat v4 object as input. ARBOL() clustering identifies and utilizes optimum 
cluster resolution parameters at each tier of clustering. It provides outputs of QC plots for each tier and stage.
Furthermore, it can return output directories figs/ srobjs/ (see examples) containing subset 
seurat objects and QC plots. It can also return an endclusts/ folder which contains tsv files named with the 
path to the endcluster (T0C0_T1C1.tsv) and containing the cellnames (one per line) 
that comprise that cluster. 

the ARBOL package also comes with taxonomy building and visualization methods that enable users to graph taxonomies and subclustering trees using ggplot grammar. This is enabled by the ggraph, tidygraph, and ggforce packages. Currently, ggraph and tree objects, as well as tree information built and used by ARBOL visualization functions are added to the Seurat object miscellaneous slot srobj@misc. 

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
https://jo-m-lab.github.io/ARBOL/ARBOLtutorial_22_12_20.html

This package is meant as a starting point for the way that we approached clustering and is meant to be edited/customized through community feedback through users such as yourself!  We have tried to organize the script such that
each processing step is contained in a modular function that can be edited and
inserted into the main iterative function.

We have dedicated effort to choosing reasonable defaults, but there is no certainty that they are
the best defaults for your data.

We recommend installing the repository through devtools.

The subclustering function of ARBOL is ARBOL() - here is an example call

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
