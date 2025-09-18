# %%
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq

import os
# Set the environment variable. No need for subprocess import here.
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

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

# --- User Configuration (Set the base directory here) ---
# This path should be the root of your project.

base_data_dir = Path('../data')

# -------------------------------------
# parameters which are needed to change as each iteration
output_path_plot = base_data_dir / '08_combined_cells_cluster_plot_iter3_umap2'
os.makedirs(output_path_plot, exist_ok = True)

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
cellmeta_combine_iter2.to_csv(output_path / new_cell_meta, index = False) 
print(cellmeta_combine_iter2.shape)

# %%
# cell type annotation
gene_panel = base_data_dir / 'marker_data/makers_annotated_merfish.csv'
df_ref_panel_ini = pd.read_csv(gene_panel)
df_ref_panel_all = df_ref_panel_ini.iloc[1:, :2]
df_ref_panel_all.columns = ['Markers', 'cell_type']
print(len(df_ref_panel_all))
df_ref_panel = df_ref_panel_all
marker_cell = dict(zip(df_ref_panel_all['Markers'],df_ref_panel_all['cell_type']))

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
        cell_by_gene_com_new = pd.concat([cell_by_gene_com_new,cell_by_gene_singel], axis = 0) 
print(cell_by_gene_com_new)

# ----------------------------------------------------------------------------
list_127 = ['cell'] + df_ref_panel['Markers'].tolist()
cell_by_gene_com_new_127 = cell_by_gene_com_new[list_127]
cell_by_gene_com_new_127.to_csv(output_path / new_cell_by_gene, index=False)


# %%
adata = sq.read.vizgen(
    path=str(output_path),
    counts_file=new_cell_by_gene,
    meta_file=new_cell_meta,
    transformation_file="micron_to_mosaic_pixel_transform.csv",
)

sc.pp.calculate_qc_metrics(adata, percent_top=(10, 50, 100, 120), inplace=True)
print("before cell filter: {}".format(adata.X.shape))

# filter the cell
sc.pp.filter_cells(adata, min_counts=3,inplace=True)
sc.pp.filter_genes(adata, min_cells=1,inplace=True)
print("after cell filter: {}".format(adata.X.shape))

# %%
n_top_genes=110
# preprocess the adata including normalizing and scale
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
    # flavor="ipgraph",
    directed=False,
)

len_cluster = len(list(set(adata.obs['clusters'].tolist())))

def generate_distinct_colors(num_colors):
    colors = plt.cm.get_cmap('Set1', num_colors)
    return [colors(i) for i in range(num_colors)]

def rgba_to_hex(rgba):
    return '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))

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
plt.title(str(len_cluster) + ' Distinct Colors from tab20 Colormap')
plt.show()

custom_colors = dict(zip(list(range(len_cluster)),hex_colors))
custom_colors = {0: '#00aa00', 1: '#67009a', 2: '#ffaa00', 3: '#00aaff', 4: '#ff557f'}

# %%
adata.uns['clusters_colors'] = list(custom_colors.values())

fig,ax=plt.subplots(figsize=(10,8))
sc.pl.umap(adata, color=['clusters'],            
            legend_fontsize=10,
            legend_fontoutline=2,
            title="clustering of pooled cells",
            show=False,ax=ax
            )
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.set_xscale('linear')
ax.set_yscale('linear')


sq.pl.spatial_scatter(
    adata,
    shape=None,
    color=[
        "clusters",
    ],
    wspace=0.4,
    size = 8,
    return_ax=True,
)
plt.show()
# output_path_plot

# %%
ctypes_dict_new = {'0': 'LPM-2', '1': 'CM', '2': 'LPM-1b', '3': 'LPM-3', '4': 'LPM-1a'}
adata.obs["cell type"] = adata.obs["clusters"].map(ctypes_dict_new).astype("category")

dict_cell_color = {'CM': '#67009a', 'LPM-1a': '#ff557f', 'LPM-1b': '#ffaa00', 'LPM-2': '#00aa00', 'LPM-3': '#00aaff'}
# Add the custom color palette to the AnnData object
adata.uns['cell type_colors'] = list(dict_cell_color.values())

fig,ax=plt.subplots(figsize=(10,8))
sc.pl.umap(
    adata,
    color="cell type", 
    legend_fontsize=10,
    legend_fontoutline=2,
    title="Annotation of cells_" + str(n_top_genes),
    show=False,ax=ax,
    palette = dict_cell_color,   
)
plt.savefig(output_path_plot / ('Annotation_umap_top_gene_'+ str(n_top_genes) +'.png'), dpi=600)
plt.show()


# %%
adata

# %%
savepath = base_data_dir / 'imputation_tangram'
os.makedirs(savepath, exist_ok=True)
# Extract gene names
genes = list(adata.var.index)

# Save as CSV
pd.DataFrame(genes, columns=["gene"]).to_csv(savepath / "gene_list_127.csv", index=False)

# %%
adata.write(savepath / "ST_1032_127genes.h5ad")

# %%
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score, silhouette_samples

# choose embedding: use UMAP if available, otherwise PCA
if "X_umap" in adata.obsm:
    X_emb = adata.obsm["X_umap"]
elif "X_pca" in adata.obsm:
    X_emb = adata.obsm["X_pca"]
else:
    raise ValueError("No 'X_umap' or 'X_pca' found in adata.obsm; compute an embedding first.")

# ensure labels are categorical and have no missing values
labels = adata.obs["cell type"].astype("category")
mask = labels.notna()
if mask.sum() != len(labels):
    print(f"Warning: {len(labels)-mask.sum()} cells have NA labels and will be ignored for silhouette calculations.")

X_emb_valid = X_emb[mask.values]
label_codes = labels[mask].cat.codes.values   # integer labels for sklearn
label_names = list(labels.cat.categories)     # category order

# global silhouette score (mean)
global_sil = silhouette_score(X_emb_valid, label_codes, metric="euclidean")
print("Global Silhouette Score (mean):", global_sil)

# per-sample silhouette values
sil_vals = silhouette_samples(X_emb_valid, label_codes, metric="euclidean")

# add back to adata.obs (NaNs for excluded cells if any)
sil_series = pd.Series(index=adata.obs.index, dtype=float)
sil_series.loc[mask.index[mask]] = sil_vals  # careful indexing
adata.obs["silhouette"] = sil_series


# %%
# compute mean silhouette per cluster (using the category names)
df_sil = adata.obs.loc[mask.index[mask], ["cell type"]].copy()
df_sil["silhouette"] = sil_vals
per_cluster = df_sil.groupby("cell type")["silhouette"].mean().reindex(label_names)
print("Per-cluster mean silhouette:")
print(per_cluster)
# save per-cluster table
per_cluster.to_csv(savepath / "per_cluster_silhouette_mean.csv")

# %%
print(output_path_plot)
