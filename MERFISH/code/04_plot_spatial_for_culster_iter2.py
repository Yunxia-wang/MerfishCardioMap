# %%
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq

import os

import decoupler as dc
import monkeybread as mb

# --- User Configuration ---
# Set the base directory for all data. This path should point to the root of
# your project, which contains the 'cell_seg_DAPI_mouse_rabbit' folder.
base_data_dir = Path('../data')

# Define the paths for the cluster data and plotting results relative to the base path
cluster_datapath = base_data_dir / '04_combined_cells_iter1'
plot_result_dir = base_data_dir / '07_combined_cells_cluster_plot_iter2_umap2'

# %%
# Get the cell meta with cluster using the relative path
annotated_cell = pd.read_csv(cluster_datapath / 'cell_metadata_combine_with_cluster_iter2.csv')
dict_cell_clu = annotated_cell.set_index('EntityID')['clusters'].to_dict()

len_cluster = len(list(set(annotated_cell['clusters'].tolist())))
cell_list = [int(i) for i in annotated_cell['EntityID']]

def generate_distinct_colors(num_colors):
    """Generates a list of distinct colors from the 'tab20' colormap."""
    colors = plt.cm.get_cmap('tab20', num_colors)
    return [colors(i) for i in range(num_colors)]

def rgba_to_hex(rgba):
    """Converts an RGBA tuple to a hex color string."""
    return '#{:02x}{:02x}{:02x}'.format(int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255))

num_colors = len_cluster
distinct_colors = generate_distinct_colors(num_colors)
hex_colors = [rgba_to_hex(color) for color in distinct_colors]
print(hex_colors)

custom_colors = dict(zip(range(len_cluster), hex_colors))
custom_colors[num_colors] = '#cccccc' # Color for un-clustered cells
print(custom_colors)

# Plot the colors to visualize them
plt.figure(figsize=(12, 2))
for i in range(len(custom_colors)):
    plt.plot(i, 0, marker='o', markersize=20, color=custom_colors[i])
plt.xlim(-1, len(custom_colors))
plt.yticks([])
plt.title(f'{len_cluster} Distinct Colors from tab20 Colormap')
plt.show()

# %%
filenames = ['slide1region5','slide1region7','slide1region9','slide2region2','slide2region5','slide1region4','slide1region6','slide1region8','slide2region0','slide2region4']

for filename in filenames:
    print('===========================================')
    print(filename)

    vizgen_dir = base_data_dir / filename
    resultpath = plot_result_dir / filename
    os.makedirs(resultpath, exist_ok=True)
    
    cellpath = "cell_by_gene.csv"
    cellmetapath = "cell_metadata_pooled_iter1.csv"

    # Generate new cell_metadata file
    datapath_meta = vizgen_dir / 'cell_metadata.csv'
    datameta = pd.read_csv(datapath_meta)
    cell_clusters = []
    
    count_cell = 0
    for i in datameta['EntityID'].tolist():
        if i in cell_list:
            count_cell += 1
            cell_clusters.append(dict_cell_clu[i])
        else:
            cell_clusters.append(num_colors)
    
    datameta.loc[:, ('cluster')] = cell_clusters
    datameta.to_csv(vizgen_dir / 'cell_metadata_pooled_iter1.csv', index=False)
    
    print(f'count_cell:{count_cell}')
    print(f'datameta.shape:{datameta.shape}')
    print(datameta.head())

    # Load AnnData object with relative paths
    adata = sq.read.vizgen(
        path=str(vizgen_dir),
        counts_file=cellpath,
        meta_file=cellmetapath,
        transformation_file="micron_to_mosaic_pixel_transform.csv",
    )

    # QC and filtering
    sc.pp.calculate_qc_metrics(adata, percent_top=(10, 50, 100, 120), inplace=True)
    fdr = adata.obsm["blank_genes"].to_numpy().sum() / adata.var["total_counts"].sum() * 100
    print(f"fdr:{fdr}")

    sc.pp.filter_cells(adata, min_counts=3, inplace=True)
    sc.pp.filter_genes(adata, min_cells=1, inplace=True)
    print(f"after cell filter: {adata.X.shape}")
    
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    
    # Set cell type and custom colors
    cell_types = [str(i) for i in adata.obs['cluster'].tolist()]
    adata.obs['cell_type'] = pd.Categorical(cell_types)
    
    custom_colors_region = {str(cell_type): custom_colors[int(cell_type)] for cell_type in set(cell_types)}
    
    adata.uns['cell_type_colors'] = list(custom_colors_region.values())    

    # Plot spatial figure for all clusters
    fig, axs = plt.subplots(1, 1, figsize=(16, 20))
    sc.pl.embedding(
        adata,
        "spatial",
        color='cell_type',
        groups=cell_types,
        wspace=0.4,
        size=7,
        ax=axs,
        palette=custom_colors_region,
        show=False
    )
    axs.yaxis.set_inverted(True)
    plt.savefig(resultpath / f'{filename}_all_cluster_kernel_density_spatial.png', dpi=600)
    
    # Plot spatial figure for each cluster
    plot_cluster = set(cell_types)
    if str(num_colors) in plot_cluster:
        plot_cluster.remove(str(num_colors))
    
    for i in list(plot_cluster):
        fig, axs = plt.subplots(1, 1, figsize=(16, 20))
        sc.pl.embedding(
            adata,
            "spatial",
            color='cell_type',
            groups=str(i),
            wspace=0.4,
            ax=axs,
            size=7,
            palette=custom_colors_region,
            show=False,
            na_color='#f1f2f1'
        )
        axs.yaxis.set_inverted(True)
        plt.savefig(resultpath / f'{filename}_cluster{i}_kernel_density_spatial.png', dpi=600)     
    
# %%
# This part of the code was separate, but it still works with the new Path objects
print(datameta['EntityID'].tolist())
print(dict_cell_clu)
