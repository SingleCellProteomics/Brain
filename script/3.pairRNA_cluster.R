#!/usr/bin/env Rscript


library(Seurat)
library(future)
library(dplyr)
library(ggplot2)
library(dplyr)
library(harmony)
library(qs)

intergration = "harmony"
random.seed = 42
set.seed(random.seed)
args <- commandArgs(trailingOnly = TRUE)

obj <- qread(args[1]) # RNA obj

Idents(obj) <- "sample"
obj <- subset(obj, subset = sample %in% c("GW19","GW15","GW13"))

npcs = 50
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(object = obj)
obj <- FindVariableFeatures(object = obj)

obj <- CellCycleScoring(
  object = obj,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)

obj <- ScaleData(object = obj)
obj <- RunPCA(obj, npcs = npcs, reduction.name = "nonSCT.pca", verbose = F)
obj <- RunHarmony(obj,
                  group.by.vars = "batch",
                  reduction = "nonSCT.pca",
                  reduction.save = "nonSCT.harmony")

## for SCT
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)
obj = SCTransform(obj, method = "glmGamPoi", assay = "RNA")
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "SCT"
obj <- RunPCA(obj, npcs = npcs, reduction.name = "SCT.pca", verbose = F)

obj <- RunHarmony(obj,
                  group.by.vars = "batch",
                  reduction = "SCT.pca",
                  reduction.save = "SCT.harmony")

npcs = 50
for (redu_name in c("nonSCT.harmony", "SCT.harmony" )){
    graphName = c( paste0(redu_name, ".nn"), paste0(redu_name, ".snn"))
    clusterName = paste0(redu_name, ".recluster")
    reduction =  paste0(redu_name, "reumap")

    obj <- FindNeighbors(object = obj,
                         dims = 1:npcs,
                         reduction = redu_name,
                         graph.name = graphName )

    obj = FindClusters(object = obj,
                       resolution = 0.5,
                       algorithm = 1,
                       graph.name = graphName[2],
                       cluster.name =  clusterName )


    obj<- RunUMAP(object = obj,
                  reduction = redu_name,
                  n.neighbors = 50L,
                  min.dist = 0.5,
                  spread = 1,
                  dims = 1:npcs,
                  reduction.name =reduction)
}

Idents(obj) <- "SCT.harmony.recluster"
obj = FindSubCluster(
    object = obj,
    cluster = "2",
    graph.name = graphName[1],
    subcluster.name = "cluster",
    resolution = 0.3,
    algorithm = 1
)

## rename clusters
obj$celltype <- "unknown"
renameCelltype <- function(obj,from,to) {
    idx = obj$cluster %in% from
    obj$celltype[idx] = to
    obj
}

obj <- renameCelltype(obj, c("13","6"), "RG"  )
obj <- renameCelltype(obj, c("2_1","2_0"), "oRG"  )
obj <- renameCelltype(obj, c("26"), "OPC"  )
obj <- renameCelltype(obj, c("2_3","18"), "Astrocyte" )
obj <- renameCelltype(obj, c("21","7","2_2","4"), "IPC-EN"  )
obj <- renameCelltype(obj, c("0","22"), "EN-newborn"  )
obj <- renameCelltype(obj, c("10","1","3"), "EN"  )
obj <- renameCelltype(obj, c("8"), "IN-CGE"  )
obj <- renameCelltype(obj, c("9"), "IN-MGE"  )
obj <- renameCelltype(obj, c("15"), "Microglia"  )
obj <- renameCelltype(obj, c("14"), "Vascular"  )

qsave(obj,"RNA_pairsample.qs")
