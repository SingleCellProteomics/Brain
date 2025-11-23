This repository contains the code and data necessary to reproduce the key findings and figures presented in the paper "Single-cell proteomic landscape of the developing human brain"

## Data Resources

*   Essential data files are located in the `data/` directory.
*   raw data for protein and RNA could be downloaded from the online browser http://xxx


## Analysis Scripts

The `scripts/` directory contains the scripts used for data analysis, clustering, and visualization.

P.S. there will be two step for the proteome umap: 

step 1. To label the blood and brain cells

``` R
Rscript 1.cluster_all_proteome.R \
        1.single_cell_protein_diann_log.csv \
        1.clean_meta.csv
```
*raw data availability indicated in the manuscript 


This will output seurat obj in qs format `all_cell.cluster.qs`, which could be load in R through `qs::qread`


step 2. re-cluster to get seperated blood and brain clusters

``` R
Rscript 2.sub_cluster_proteome.R all_cell.cluster.qs
```
then the seurat objs for brain and blood cells will sperate into `brain.cluster.qs` and  the `blood.cluster.qs` for visulization.


## Interactive Notebook

A precompiled notebook demonstrating the generation of the main figures is available in the `notebook/` directory. This notebook can be run using Pluto.jl, allowing for interactive exploration of the analysis.

For complete details, please refer to the paper: Single Cell Proteomics in the Developing Human Brain.
