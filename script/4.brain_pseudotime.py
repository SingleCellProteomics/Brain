#!/usr/bin/env ipython

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import scanpy.external as sce
adata = sc.read_h5ad("brain.cluster.h5ad") # convert from brain.cluster.qs to h5ad
sc.pl.embedding(adata, basis = "umap", color = "celltype")
adata = adata[adata.obs['celltype'].isin(['RG', 'IPC-EN', 'EN'])]
sc.pp.neighbors(adata, n_neighbors=30, use_rep = "harmony")

sc.tl.paga(adata, groups="celltype")
sc.pl.paga(adata, color="celltype")

adata.uns["iroot"] = np.flatnonzero(adata.obs["celltype"] == "RG")[0]


sc.tl.diffmap(adata,n_comps=10)
sc.tl.dpt(adata,n_dcs=5)
sc.pl.embedding(adata, basis = "X_diffmap", color = "celltype")
sc.pl.embedding(adata, basis = "X_diffmap", color = "dpt_pseudotime")
sc.pl.embedding(adata, basis = "umap", color = "dpt_pseudotime")

adata.write("brain_pseudotime.h5ad")
