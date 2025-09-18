# ---------------------------------------------------------------------
# Import necessary libraries
# ---------------------------------------------------------------------
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq

import os
import subprocess

# This line ensures a specific TensorFlow option is set
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

import decoupler as dc
import monkeybread as mb
from matplotlib.cm import _colormaps as colormaps

print(f"Scanpy version: {sc.__version__}")
print(f"Squidpy version: {sq.__version__}")
print(f"DeCoupler version: {dc.__version__}")

# ---------------------------------------------------------------------
# Define paths using pathlib for improved portability
# ---------------------------------------------------------------------
# This assumes the script is located in a parent directory of 'Hakan_SP'
pri_path = Path('../data')

# Parameters which are needed to change for each iteration
plot_result = '08_combined_cells_cluster_plot_iter3_umap2'
output_path_plot = pri_path / plot_result

# Filenames for input and output files
last_cell_meta_with_cluster_file = 'cell_metadata_combine_with_cluster_iter2.csv'
new_cell_meta_file = 'cell_metadata_combine_iter3.csv'
new_cell_by_gene_file = 'cell_by_gene_combine_iter3.csv'
new_cell_meta_with_cluster_file = 'cell_metadata_combine_with_cluster_iter3.csv'
gene_panel_file = pri_path / 'Data' / 'marker_data' / 'makers_annotated_merfish.csv'
cell_by_gene_combine_path = pri_path.parent / '01_combined_cells' / 'cell_by_gene_combine.csv'
plots_vizgen_combine_path = pri_path / 'plots_vizgen_combine'

# Paths for the analysis
output_path = pri_path / '04_combined_cells_iter1'

# Create the output directories if they don't exist
os.makedirs(output_path_plot, exist_ok=True)
os.makedirs(plots_vizgen_combine_path, exist_ok=True)

# ---------------------------------------------------------------------
# Prepare metadata and filter cells
# ---------------------------------------------------------------------
cellmeta_iter1 = pd.read_csv(output_path / last_cell_meta_with_cluster_file)
print(f"Shape of initial combined metadata: {cellmeta_iter1.shape}")

# Filter out clusters
cellmeta_iter2 = cellmeta_iter1[~cellmeta_iter1['clusters'].isin([5])]

cellmeta_combine_iter1 = pd.read_csv(output_path / 'cell_metadata_combine_iter2.csv')
cellmeta_combine_iter2 = cellmeta_combine_iter1[cellmeta_combine_iter1['EntityID'].isin(cellmeta_iter2['EntityID'].tolist())]
cellmeta_combine_iter2.to_csv(output_path / new_cell_meta_file, index=False)
print(f"Shape of filtered combined metadata: {cellmeta_combine_iter2.shape}")

# ---------------------------------------------------------------------
# Gene panel and data preparation
# ---------------------------------------------------------------------
df_ref_panel_ini = pd.read_csv(gene_panel_file)
df_ref_panel_all = df_ref_panel_ini.iloc[1:, :2]
df_ref_panel_all.columns = ['Markers', 'cell_type']
print(f"Number of markers in gene panel: {len(df_ref_panel_all)}")

cell_by_gene_combine_data = pd.read_csv(cell_by_gene_combine_path)
print(f"Shape of combined cell-by-gene data: {cell_by_gene_combine_data.shape}")

selected_cell_lists = cellmeta_combine_iter2['EntityID'].tolist()
cell_by_gene_com_new = cell_by_gene_combine_data[cell_by_gene_combine_data['cell'].isin(selected_cell_lists)]

list_127 = ['cell'] + df_ref_panel_all['Markers'].tolist()
cell_by_gene_com_new_127 = cell_by_gene_com_new[list_127]
cell_by_gene_com_new_127.to_csv(output_path / new_cell_by_gene_file, index=False)
print(f"Shape of new cell-by-gene data: {cell_by_gene_com_new_127.shape}")

# ---------------------------------------------------------------------
# Load data into AnnData object and perform QC
# ---------------------------------------------------------------------
adata = sq.read.vizgen(
    path=output_path,
    counts_file=new_cell_by_gene_file,
    meta_file=new_cell_meta_file,
    transformation_file="micron_to_mosaic_pixel_transform.csv",
)

sc.pp.calculate_qc_metrics(adata, percent_top=(10, 50, 100, 120), inplace=True)
print(f"Shape before cell filtering: {adata.X.shape}")

# Filter the cells and genes
sc.pp.filter_cells(adata, min_counts=3, inplace=True)
sc.pp.filter_genes(adata, min_cells=1, inplace=True)
print(f"Shape after cell filtering: {adata.X.shape}")

# Plot QC distributions
fig, axs = plt.subplots(1, 3, figsize=(15, 4))
axs[0].set_title("Total transcripts per cell")
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
axs[1].set_title("Unique transcripts per cell")
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1])
axs[2].set_title("Volume of segmented cells")
sns.histplot(adata.obs["volume"], kde=False, ax=axs[2])
plt.tight_layout()
plt.savefig(plots_vizgen_combine_path / 'the distribution of total transcripts.png', dpi=600)
plt.show()

# ---------------------------------------------------------------------
# Preprocess data and perform clustering
# ---------------------------------------------------------------------
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

# Define custom colors for clusters
len_cluster = len(adata.obs['clusters'].unique())
custom_colors = {0: '#00aa00', 1: '#67009a', 2: '#ffaa00', 3: '#00aaff', 4: '#ff557f'}

# Plot the UMAP with clusters
adata.uns['clusters_colors'] = list(custom_colors.values())
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(
    adata,
    color=['clusters'],
    legend_fontsize=10,
    legend_fontoutline=2,
    title="Clustering of pooled cells",
    show=False,
    ax=ax
)
plt.savefig(output_path_plot / 'cluster_umap_pooled_cells.png', dpi=600)
plt.show()

# ---------------------------------------------------------------------
# Machine learning annotation with DeCoupler
# ---------------------------------------------------------------------
df_ref_panel = df_ref_panel_all
dc.run_ora(
    mat=adata,
    net=df_ref_panel,
    source='cell_type',
    target='Markers',
    min_n=3,
    verbose=True,
    use_raw=False
)

acts = dc.get_acts(adata, obsm_key='ora_estimate')
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

df = dc.rank_sources_groups(acts, groupby='clusters', reference='rest', method='t-test_overestim_var')
n_ctypes = 3
ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
ctypes_dict_new = {k: v[0] for k, v in ctypes_dict.items()}
ctypes_dict_new = {'0': 'LPM-2', '1': 'CM', '2': 'LPM-1b', '3': 'LPM-3', '4': 'LPM-1a'}

dict_cell_color = {
    'CM': '#67009a',
    'LPM-1a': '#ff557f',
    'LPM-1b': '#ffaa00',
    'LPM-2': '#00aa00',
    'LPM-3': '#00aaff'
}

acts.obs["cell type"] = acts.obs["clusters"].map(ctypes_dict_new).astype("category")
acts.uns['cell type_colors'] = list(dict_cell_color.values())

fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(
    acts,
    color="cell type",
    legend_fontsize=10,
    legend_fontoutline=2,
    title="Annotation of cells",
    show=False,
    ax=ax,
    palette=dict_cell_color,
)
plt.savefig(output_path_plot / 'Annotation_umap.png', dpi=600)
plt.show()

# ---------------------------------------------------------------------
# Violin plots for gene expression
# ---------------------------------------------------------------------
adata.obs["clusters"] = acts.obs["clusters"]
adata.obs["cell type"] = adata.obs["clusters"].map(ctypes_dict_new).astype("category")

nrow = 4
ncol = 4
genes = ['fsta', 'efnb2a', 'prdm1a', 'tbx1', 'nkx2.5', 'myl7', 'hand2', 'gata5', 'meis3', 'hoxb3a', 'hoxb1b', 'meis2a', 'irx3a']
fig, axes = plt.subplots(nrow, ncol, figsize=(15, 10))
desired_order = ['LPM-1a', 'LPM-1b', 'LPM-2', 'LPM-3', 'CM']
axes = axes.flatten()
groupby = 'cell type'

for idx, gene in enumerate(genes):
    ax = axes[idx]
    if gene in adata.var_names:
        sc.pl.violin(adata, gene, groupby=groupby, ax=ax, order=desired_order, palette=dict_cell_color, bw=0.2, show=False)
        ax.set_title(gene)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlabel("")
        ax.set_ylabel("Expression level")
        
        group_means = adata.to_df()[gene].groupby(adata.obs[groupby]).mean()
        group_means = group_means.loc[desired_order]
        for i, (group, mean_val) in enumerate(group_means.items()):
            ax.hlines(y=mean_val, xmin=i - 0.1, xmax=i + 0.1, color='black', linewidth=2, label="Mean" if i == 0 else "")
        ax.set_ylim(-0.2, 8)
    else:
        ax.set_title(f"Gene '{gene}' not found")
        ax.axis('off')

for i in range(len(genes), len(axes)):
    fig.delaxes(axes[i])

plt.tight_layout()
plt.savefig(plots_vizgen_combine_path / 'second_heart_field_genes_gene_expression.png', dpi=600)
plt.show()
