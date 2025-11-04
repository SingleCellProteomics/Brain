This repository contains the code and data necessary to reproduce the key findings and figures presented in the paper "Single Cell Proteomics in the Developing Human Brain."

## Data Resources

*   Essential data files are located in the `data/` directory.

## Analysis Scripts

The `scripts/` directory contains the scripts used for data analysis, clustering, and visualization.

P.S. there will be two step for protein umap: 

step 1. To label the blood and brain cells

``` R
Rscript 1.cluster_all_proteome.R 1.single_cell_protein_diann_log.csv 1.clean_meta.csv
```
you will get the `all_cell.cluster.qs`

`1.single_cell_protein_diann_log.csv` and `1.clean_meta.csv` could be downloaded from the online browser http://xxx 


step 2. re-cluster to get seperated cells cluster  

``` R
Rscript 2.sub_cluster_proteome.R all_cell.cluster.qs
```
then the `brain.cluster.qs` and  the `blood.cluster.qs` will return for visulization.
which could be load in R through `qs::qread`


## Interactive Notebook

A precompiled notebook demonstrating the generation of the main figures is available in the `notebook/` directory. This notebook can be run using Pluto.jl, allowing for interactive exploration of the analysis.

For complete details, please refer to the paper: Single Cell Proteomics in the Developing Human Brain.
