import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

import matplotlib
matplotlib.use('Agg')

data_dir = "/PATH/data4anadata"

os.chdir(data_dir)
# load sparse matrix
X = io.mmread(data_dir + "/counts.mtx")

# create anndata object
adata = anndata.AnnData(X = X.transpose().tocsr())

# load cell metadata
cell_meta = pd.read_csv(data_dir + "/metadata.csv")

# load gene names
with open(data_dir + "/gene_names.csv", "r") as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction
pca = pd.read_csv(data_dir + "/pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['umap_1'].to_numpy(), adata.obs['umap_2'].to_numpy())).T

# plot a UMAP colored bysampleID to test
sc.pl.umap(adata, color = ['cell_type'], frameon = False, save = 'umap_scanpy.svg')

# save dataset as anndata format
adata.write(data_dir + '/young.h5ad')
