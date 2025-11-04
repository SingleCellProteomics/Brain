#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(harmony)
library(qs)
set.seed(42)


args <- commandArgs(trailingOnly = TRUE)

df = read.csv(file = args[1]) # single cell proteome data
mtt <- read.csv(args[2]) # meta data

obj = CreateSeuratObject(counts = df,
                         assay = "protein",
                         meta.data = mtt,
                         min.cell=1, min.features=1)


obj[["protein"]] <- split(obj[["protein"]], f = obj$orig.ident)
obj = SCTransform(obj, method = "glmGamPoi", assay = "protein")

obj <- RunPCA(obj, npcs = 50, verbose = F)

#######################
## intergration step ##
#######################
DefaultAssay(obj) <- "SCT"
obj <- IntegrateLayers(
    object = obj, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "seurat.harmony",
    normalization.method = "SCT",
    verbose = FALSE)

obj <- RunHarmony(obj,
                  group.by.vars = "batch",
                  reduction = "pca",
                  reduction.save = "harmony",
                  assay.use = "SCT")
DefaultAssay(obj) <- "SCT"
obj <- JoinLayers(obj)

###############
## umap step ##
###############
graphName = c( "SCT_nn","SCT_snn")
clusterName = "protein_cluster"
reduction =   "umap"
redu_name = "harmony"
obj <- FindNeighbors(object = obj,
                     dims = 1:15,
                     reduction = redu_name,
                     graph.name = graphName )

obj = FindClusters(object = obj,
                   resolution = 0.8,
                   algorithm = 1,
                   graph.name = graphName[2],
                   cluster.name = clusterName )

obj<- RunUMAP(object = obj,
              reduction = redu_name,
              n.neighbors = 50L,
              min.dist = 0.8,
              spread = 1,
              dims = 1:50,
              reduction.name =reduction)

obj <- subset(obj, subset = seurat_clusters %in% c("4","11"), invert = T)

obj$celltype <- "Undefined"
renameCelltype <- function(obj,from,to) {
    idx = obj$seurat_clusters %in% from
    obj$celltype[idx] = to
    obj
}

obj <- renameCelltype(obj, c("1","8","10","2"), "Blood"  )
obj <- renameCelltype(obj, c("3","9","0","6","5","7","13","12","14"), "Brain"  )

qsave(obj,"all_cell.cluster.qs")
