import scvelo as scv
import scanpy as sc
import os
import matplotlib
matplotlib.use('Agg')


data_dir = "/PATH/data4anadata"

loom_dir = "/PATH/velocyto_out"

os.chdir(data_dir)

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor = 'white', dpi = 100, frameon = False)

adata = sc.read_h5ad(data_dir + '/young.h5ad')

# load loom files for spliced/unspliced matrics for each sample
ldata1 = scv.read(loom_dir + "/1_FKDL210007809-1a.loom")

# rename barcodes in order to merge
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0 : len(bc) -1] + '-1' for bc in barcodes]
ldata1.obs.index = barcodes


# make variable names unique
ldata1.var_names_make_unique()

# concatenate all sample
ldata = ldata1

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

# plot umap to check
sc.pl.umap(adata, color='cell_type', frameon=False, legend_loc='on data', title='', save='_celltypes.pdf', palette = ["#CDAF95", "#87CEFA", "#1E90FF", "#8B4513", "#FF69B4", "#8B8682"] )

# colors
# "seashell4", "chocolate4", "peachpuff3", "lightskyblue", "dodgerblue", "hotpink"
# "#8B8682", "#8B4513", "#CDAF95", "#87CEFA", "#1E90FF", "#FF69B4"


# save dataset as anndata format
adata.write(data_dir + '/young_seurat_scvelo.h5ad')
