import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import os

import matplotlib
matplotlib.use('Agg')

data_dir = "PATH/data4anadata"
out_dir =  "/PATH/RNA_velocity"

loom_dir = "/PATH/velocyto_out"

os.chdir(out_dir)

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor = 'white', dpi = 300, frameon = False)

# load the data
adata = sc.read_h5ad(data_dir + '/young_seurat_scvelo.h5ad')


# propotions of spliced/unspliced
scv.pl.proportions(adata, groupby='cell_type', save = "_spliced_by_celltype.pdf")

# pre-process to compute RNA velocity
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# compute velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)


# Visualize velocity fields

# get a broad visualization of RNA velocity across all genes and all cels by bisualizing avector field on top of 2D dimensional reduction.
scv.pl.velocity_embedding(adata, basis = "umap", frameon = False, save = "embedding.pdf", palette = ["#CDAF95", "#87CEFA", "#1E90FF", "#8B4513", "#FF69B4", "#8B8682"] )

scv.pl.velocity_embedding_grid(adata, basis = "umap", color = 'cell_type', save = 'embedding_grid.pdf', title = '', scale = 0.25, palette = ["#CDAF95", "#87CEFA", "#1E90FF", "#8B4513", "#FF69B4", "#8B8682"] )

scv.pl.velocity_embedding_stream(adata, basis = 'umap', color = 'cell_type', palette = ["#CDAF95", "#87CEFA", "#1E90FF", "#8B4513", "#FF69B4", "#8B8682"] , save = 'embedding_stream.svg', title= '')

scv.pl.velocity_embedding_stream(adata, basis = 'umap', color = ['cell_type', 'orig.ident'], save='embedding_stream_sample.svg', title = '', palette = ["#CDAF95", "#87CEFA", "#1E90FF", "#8B4513", "#FF69B4", "#8B8682"] )

# visualize the dynamics of favorie genes. Here we show the ratio of unspiced to spliced transcripts for Bric5, Top2a and Mki67, as well as the velocity and expression values as feature plots.
scv.pl.velocity(adata, var_names = ['Birc5', 'Top2a', 'Mki67'], color = 'cell_type', save = 'dynamics_of_Birc5_Top2a_Mki67.pdf')

adata.write(out_dir + '/young_velocity.h5ad')
