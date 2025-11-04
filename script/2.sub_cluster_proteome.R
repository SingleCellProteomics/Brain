#!/usr/bin/env Rscript

## args:
library(Seurat)
library(dplyr)
library(harmony)
library(qs)
set.seed(42)
args <- commandArgs(trailingOnly = TRUE)

obj_all = qread(file = args[1])

obj <- subset(obj_all,subset = celltype == "Blood")
DefaultAssay(obj) <- "protein"
obj <- RunHarmony(obj,
                  group.by.vars = "Sample",
                  reduction = "pca",
                  assay.use = "protein",
                  reduction.save = "harmony")

graphName = c( "nn","snn")
clusterName = "protein_cluster"
reduction =   "umap"
redu_name = "harmony"
obj <- FindNeighbors(object = obj,
                     dims = 1:50,
                     reduction = redu_name,
                     graph.name = graphName )

obj = FindClusters(object = obj,
                   resolution = 0.5,
                   algorithm = 1,
                   graph.name = graphName[2],
                   cluster.name = clusterName )

obj<- RunUMAP(object = obj,
              reduction = redu_name,
              n.neighbors = 50L,
              min.dist = 0.5, #0.5,
              spread = 1,
              dims = 1:50,
              reduction.name =reduction)

obj$celltype <- "unknown"
renameCelltype <- function(obj,from,to) {
    idx = obj$seurat_clusters %in% from
    obj$celltype[idx] = to
    obj
}

qsave(obj,"blood.cluster.qs")


obj <- subset(obj_all,subset = celltype == "Blood",  invert = T)

#######################
## intergration step ##
#######################
DefaultAssay(obj) <- "protein"
set.seed(42)

obj <- RunHarmony(obj,
                  group.by.vars = "Sample",
                  reduction = "pca",
                  reduction.save = "harmony",
                  assay.use = "protein")


###############
## umap step ##
###############
DefaultAssay(obj) <- "SCT"
graphName = c( "nn","snn")
clusterName = "protein_cluster"
reduction =   "umap"
redu_name = "harmony"
obj <- FindNeighbors(object = obj,
                     dims = 1:50,
                     reduction = redu_name,
                     graph.name = graphName )

obj = FindClusters(object = obj,
                   resolution = 2,
                   algorithm = 1,
                   graph.name = graphName[2],
                   cluster.name = clusterName )

obj<- RunUMAP(object = obj,
              reduction = redu_name,
              n.neighbors = 50L,
              min.dist = 0.5,
              spread = 1,
              dims = 1:50,
              reduction.name =reduction)

qsave(obj,"brain.cluster.qs")

