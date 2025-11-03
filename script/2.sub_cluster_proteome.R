#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(harmony)
library(qs)
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
obj_all = qread(file = args[1]) # the first cluster with all cell.

##############################################
## subset it and re-cluster the Brain cells ##
##############################################

obj <- subset(obj_all,subset = celltype == "Blood",  invert = T) # remove all Blood cell type

obj <- RunHarmony(obj,
                  group.by.vars = "Sample",
                  reduction = "pca",
                  reduction.save = "harmony",
                  assay.use = "SCT")
DefaultAssay(obj) <- "protein"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "SCT"

###############
## umap step ##
###############
graphName = c( "SCT_nn","SCT_snn")
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
              min.dist = 0.8, #0.5,
              spread = 3,
              dims = 1:50,
              reduction.name =reduction)

qsave(obj,"brain.cluster.SCT.qs")


############################################
## subset it and re-cluster the Blood cells  ##
############################################

obj <- subset(obj_all,subset = celltype == "Blood")

obj[["protein"]] <- split(obj[["protein"]], f = obj$orig.ident)
obj = SCTransform(obj, method = "glmGamPoi", assay = "protein")

obj <- RunHarmony(obj,
                  group.by.vars = "Sample",
                  reduction = "pca",
                  assay.use = "SCT",
                  reduction.save = "harmony")

DefaultAssay(obj) <- "SCT"
graphName = c( "SCT_nn","SCT_snn")
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

qsave(obj,"blood.cluster.qs")
