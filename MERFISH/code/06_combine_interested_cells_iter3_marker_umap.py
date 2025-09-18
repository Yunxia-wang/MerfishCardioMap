# %%
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq

import os
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'
import subprocess
# Labeling
plt.rcParams.update({'font.size': 20})
import decoupler as dc
import monkeybread as mb
from matplotlib.cm import _colormaps as colormaps 
# from .autonotebook import tqdm as notebook_tqdm
print(sc.__version__)

# %%
import squidpy as sq
print(dc.__version__)

# %%
def deal_str(data):
    data = str(data)+'\t'
    return data

# --- User Configuration (relative paths) ---
# Set the base directory for all data. This path should point to the root of
# your project, which contains the 'cell_seg_DAPI_mouse_rabbit' folder.

base_data_dir = Path('../data')
marker_data_dir = base_data_dir / 'marker_data/'

# -------------------------------------
# Parameters which are needed to change as each iteration
output_path_plot = base_data_dir / '08_combined_cells_cluster_plot_iter3_umap2'
os.makedirs(output_path_plot, exist_ok=True)

last_cell_meta = 'cell_metadata_combine_iter2.csv'
last_cell_meta_with_cluster = 'cell_metadata_combine_with_cluster_iter2.csv'
new_cell_meta = 'cell_metadata_combine_iter3.csv'
new_cell_by_gene = 'cell_by_gene_combine_iter3.csv'
new_cell_meta_with_cluster = 'cell_metadata_combine_with_cluster_iter3.csv'

# -------------------------------------
filenames = ['slide1region5','slide1region7','slide1region9','slide2region2','slide2region5','slide1region4','slide1region6','slide1region8','slide2region0','slide2region4']
output_path = base_data_dir / '04_combined_cells_iter1'

# -------------------------------------
cellmeta_iter1 = pd.read_csv(output_path / last_cell_meta_with_cluster)
print(cellmeta_iter1.shape)
cellmeta_iter2 = cellmeta_iter1[~cellmeta_iter1['clusters'].isin([5])]

cellmeta_combine_iter1 = pd.read_csv(output_path / last_cell_meta)
cellmeta_combine_iter2 = cellmeta_combine_iter1[cellmeta_combine_iter1['EntityID'].isin(cellmeta_iter2['EntityID'].tolist())]
cellmeta_combine_iter2.to_csv(output_path / new_cell_meta, index=False) 
print(cellmeta_combine_iter2.shape)

# %%
# Cell type annotation
gene_panel = marker_data_dir / 'makers_annotated_merfish.csv'
df_ref_panel_ini = pd.read_csv(gene_panel)
df_ref_panel_all = df_ref_panel_ini.iloc[1:, :2]
df_ref_panel_all.columns = ['Markers', 'cell_type']
print(len(df_ref_panel_all))
df_ref_panel = df_ref_panel_all
marker_cell = dict(zip(df_ref_panel_all['Markers'], df_ref_panel_all['cell_type']))

# ----------------------------------------------------------------------------
cell_by_gene_combine_path = base_data_dir / '01_combined_cells/cell_by_gene_combine.csv'
cell_by_gene_combine_data = pd.read_csv(cell_by_gene_combine_path)
print(cell_by_gene_combine_data.shape)

selected_cell_lists = cellmeta_combine_iter2['EntityID'].tolist()

for index, s_cell in enumerate(selected_cell_lists):
    cell_by_gene_singel = cell_by_gene_combine_data[cell_by_gene_combine_data['cell'] == s_cell]

    if index == 0:
        cell_by_gene_com_new = cell_by_gene_singel
    else:
        cell_by_gene_com_new = pd.concat([cell_by_gene_com_new, cell_by_gene_singel], axis=0) 
print(cell_by_gene_com_new)

# ----------------------------------------------------------------------------
list_127 = ['cell'] + df_ref_panel['Markers'].tolist()
cell_by_gene_com_new_127 = cell_by_gene_com_new[list_127]
cell_by_gene_com_new_127

# %%
# Saving the new cell by gene file
cell_by_gene_com_new_127.to_csv(output_path / new_cell_by_gene, index=False)

# Read the data for the next steps
adata = sq.read.vizgen(
    path=str(output_path),
    counts_file=new_cell_by_gene,
    meta_file=new_cell_meta,
    transformation_file="micron_to_mosaic_pixel_transform.csv",
)

TF_ENABLE_ONEDNN_OPTS=0
sc.pp.calculate_qc_metrics(adata, percent_top=(10, 50, 100, 120), inplace=True)
print("before cell filter: {}".format(adata.X.shape))

# filter the cell
sc.pp.filter_cells(adata, min_counts=3, inplace=True)
sc.pp.filter_genes(adata, min_cells=1, inplace=True)
print("after cell filter: {}".format(adata.X.shape))

# %%
# Preprocess the adata including normalizing and scale
adata.layers["counts"] = adata.X.copy()
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=110)
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(
    adata,
    key_added="clusters",
    resolution=0.5,
    n_iterations=2,
    directed=False,
)

len_cluster = len(list(set(adata.obs['clusters'].tolist())))

def generate_distinct_colors(num_colors):
    colors = plt.cm.get_cmap('Set1', num_colors)
    return [colors(i) for i in range(num_colors)]

def rgba_to_hex(rgba):
    return '#{:02x}{:02x}{:02x}'.format(int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255))

num_colors = len_cluster
distinct_colors = generate_distinct_colors(num_colors)
hex_colors = [rgba_to_hex(color) for color in distinct_colors]
print(hex_colors)
# Plot the colors to visualize them
plt.figure(figsize=(10, 2))
for i in range(num_colors):
    plt.plot(i, 0, marker='o', markersize=20, color=hex_colors[i])
plt.xlim(-1, num_colors)
plt.yticks([])
plt.title(str(len_cluster) + ' Distinct Colors from Set1 Colormap')
plt.show()

custom_colors = dict(zip(range(len_cluster), hex_colors))
custom_colors = {0: '#00aa00', 1: '#67009a', 2: '#ffaa00', 3: '#00aaff', 4: '#ff557f'}

# %%
adata.uns['clusters_colors'] = list(custom_colors.values())

fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata, color=['clusters'],
            legend_fontsize=10,
            legend_fontoutline=2,
            title="clustering of pooled cells",
            show=False, ax=ax)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.set_xscale('linear')
ax.set_yscale('linear')

sq.pl.spatial_scatter(
    adata,
    shape=None,
    color=["clusters"],
    wspace=0.4,
    size=8,
    return_ax=True,
)

# %%
# Use the correct path for saving plots
resultpath = output_path_plot

nrow = 5
ncol = 3
ctypes_dict_new = {'0': 'LPM-2', '1': 'CM', '2': 'LPM-1b', '3': 'LPM-3', '4': 'LPM-1a'}
dict_cell_color = {'CM': '#67009a', 'LPM-1a': '#ff557f', 'LPM-1b': '#ffaa00', 'LPM-2': '#00aa00', 'LPM-3': '#00aaff'}
genes = ['fsta','efnb2a','prdm1a','tbx1','nkx2.5','myl7','hand2','gata5','meis3','hoxb3a','hoxb1b','meis2a','irx3a']

fig, axes = plt.subplots(nrow, ncol, figsize=(15, 20))
desired_order = ['LPM-1a','LPM-1b','LPM-2','LPM-3', 'CM']
axes = axes.flatten()

for idx, gene in enumerate(genes):
    ax = axes[idx]
    sc.pl.umap(adata, color=gene, ax=ax, cmap='Blues', show=False, s=130)
    ax.set_title(gene)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel("")
    ax.set_ylabel("")

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_edgecolor("grey")
    
for i in range(len(genes), len(axes)):
    fig.delaxes(axes[i])

plt.tight_layout()
plt.savefig(resultpath / 'second heart field genes_gene expression_umap.png', dpi=600)
plt.show()

# %%
# Use the correct path for saving plots
resultpath = output_path_plot

nrow = 5
ncol = 3
ctypes_dict_new = {'0': 'LPM-2', '1': 'CM', '2': 'LPM-1b', '3': 'LPM-3', '4': 'LPM-1a'}
dict_cell_color = {'CM': '#67009a', 'LPM-1a': '#ff557f', 'LPM-1b': '#ffaa00', 'LPM-2': '#00aa00', 'LPM-3': '#00aaff'}
genes = ['fsta','efnb2a','prdm1a','tbx1','nkx2.5','myl7','hand2','gata5','meis3','hoxb3a','hoxb1b','meis2a','irx3a','nkx2.7']
fig, axes = plt.subplots(nrow, ncol, figsize=(15, 20))
desired_order = ['LPM-1a','LPM-1b','LPM-2','LPM-3', 'CM']
axes = axes.flatten()
groupby = 'cell type'
for idx, gene in enumerate(genes):
    ax = axes[idx]
    sc.pl.umap(adata, color=gene, ax=ax, cmap='Blues', show=False, s=130)
    ax.set_title(gene)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel("")
    ax.set_ylabel("")

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_edgecolor("grey")
    
for i in range(len(genes), len(axes)):
    fig.delaxes(axes[i])

plt.tight_layout()
plt.savefig(resultpath / 'second heart field genes_gene expression_umap_nkx2.7.png', dpi=600)
plt.show()

# %%
# Use the correct path for saving plots
resultpath = output_path_plot

nrow = 1
ncol = 1
ctypes_dict_new = {'0': 'LPM-2', '1': 'CM', '2': 'LPM-1b', '3': 'LPM-3', '4': 'LPM-1a'}
dict_cell_color = {'CM': '#67009a', 'LPM-1a': '#ff557f', 'LPM-1b': '#ffaa00', 'LPM-2': '#00aa00', 'LPM-3': '#00aaff'}

genes = ['nkx2.7']
fig, ax = plt.subplots(nrow, ncol, figsize=(5, 4))
sc.pl.umap(adata, color=genes, ax=ax, cmap='Blues', show=False, s=130)
ax.set_title(genes[0])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel("")
ax.set_ylabel("")

for spine in ax.spines.values():
    spine.set_linewidth(1.5)
    spine.set_edgecolor("grey")
    
plt.tight_layout()
plt.savefig(resultpath / 'gene_expression_umap_nkx2.7.png', dpi=600)
plt.show()

# %%
# Use the correct path for saving plots
resultpath = output_path_plot

nrow = 2
ncol = 1
ctypes_dict_new = {'0': 'LPM-2', '1': 'CM', '2': 'LPM-1b', '3': 'LPM-3', '4': 'LPM-1a'}
dict_cell_color = {'CM': '#67009a', 'LPM-1a': '#ff557f', 'LPM-1b': '#ffaa00', 'LPM-2': '#00aa00', 'LPM-3': '#00aaff'}
genes = ['rbfox1l','tnnt2a','cyp26c1','pltp','rdh10a']
fig, axes = plt.subplots(nrow, ncol, figsize=(15, 8))
axes = axes.flatten()
for idx, gene in enumerate(genes):
    ax = axes[idx]
    sc.pl.umap(adata, color=gene, ax=ax, cmap='Blues', show=False, s=130)
    ax.set_title(gene)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_edgecolor("grey")
    
for i in range(len(genes), len(axes)):
    fig.delaxes(axes[i])

plt.tight_layout()
plt.savefig(resultpath / 'gene_expression_umap_01.png', dpi=600)
plt.show()

# %%
# Use the correct path for saving plots
resultpath = output_path_plot

nrow = 1
ncol = 1
ctypes_dict_new = {'0': 'LPM-2', '1': 'CM', '2': 'LPM-1b', '3': 'LPM-3', '4': 'LPM-1a'}
dict_cell_color = {'CM': '#67009a', 'LPM-1a': '#ff557f', 'LPM-1b': '#ffaa00', 'LPM-2': '#00aa00', 'LPM-3': '#00aaff'}
genes = ['aldh1a2']
fig, ax = plt.subplots(nrow, ncol, figsize=(5, 4))
sc.pl.umap(adata, color=genes, ax=ax, cmap='Blues', show=False, s=130)
ax.set_title(genes[0])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel("")
ax.set_ylabel("")
for spine in ax.spines.values():
    spine.set_linewidth(1.5)
    spine.set_edgecolor("grey")
    
plt.tight_layout()
plt.savefig(resultpath / 'gene_expression_umap_aldh1a2.png', dpi=600)
plt.show()

# %%
print(resultpath)
